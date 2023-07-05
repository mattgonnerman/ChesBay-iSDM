lapply(c("dplyr", "sf", "auk", "lubridate", "tidyr", "ggplot2"), require, character.only = T)

### Shared Information
year.list <- 2016:2021
nyear <- length(year.list)

'%notin%' <- Negate('%in%')

### Create grid for estimating underlying intensity surface
#Use state outlines for Maryland, Virginia, Deleware
states <- st_read("./State SHP/cb_2018_us_state_5m.shp") %>% 
  filter(STUSPS %in% c("MD", "VA", "DE")) %>%
  st_transform(32619)
#Create Hexagonal grid across study area. Estimate occupancy/abundance intensity surface within each.
grid.surface <- st_as_sf(st_make_grid(states, what = "polygons", cellsize = 10000, square = F)) %>%
  mutate(GridID = row_number())
grid.centers <- st_as_sf(st_make_grid(states, what = "centers", cellsize = 10000, square = F)) %>%
  mutate(GridID = row_number())

n.cells <- max(grid.surface$GridID)

grid.neighbors <- st_intersection(grid.surface, grid.surface) %>% st_drop_geometry() %>%
  rename(GridID = GridID.1, Touching = GridID) %>%
  filter(GridID != Touching) %>%
  arrange(GridID, Touching) %>%
  group_by(GridID) %>%
  mutate(row = row_number())

adj <- grid.neighbors$Touching
nadj <- length(adj)
grid.n.neighbors.df <- grid.neighbors %>%
  summarize(Num = max(row))
num <- grid.n.neighbors.df$Num
weights <- rep(1, length(adj))

st_write(grid.surface, dsn = "./SamplingGrid.shp", delete_layer = T)

###List of species to track
bird.codes.all <- read.csv("./BirdCodes.csv")
birds.use <- sort(c("MALL", "ABDK", "ABDU", "CAGO", "BWTE", "AMWI", "NOPI", 
                    "GWTE", "AGWT","GADW", "NSHO", "WODU", "GSGO", "GWFG",
                    "BLSC", "SUSC", "LTDU", "LESC"))
birds.codes <- bird.codes.all %>% filter(Alpha %in% birds.use)


### EBird
#Data accessed directly through the EBird download portal
#https://ebird.org/data/download
setwd("./EBird/")
dir.create("data", showWarnings = FALSE)

#Data only comes in txt files, need to use the auk package
#first provide file names for data and where you want filtered outputs to be stored
f_in <- c("ebd_US-MD_relNov-2022.txt","ebd_US-DE_relNov-2022.txt","ebd_US-VA_relNov-2022.txt")
f_sample <- "ebd_sampling_relNov-2022.txt"
f_out <- c("EBird_MD_Waterfowl.txt", "EBird_DE_Waterfowl.txt", "EBird_VA_Waterfowl.txt")
f_sample_out <- c("EBird_MD_Waterfowl.txt", "EBird_DE_Waterfowl.txt", "EBird_VA_Waterfowl.txt")
ebird.raw.list <- ebird.sampling.list <-  ebird.obs.list <- vector("list", 3)

# Use the auk package to parse the raw EBird data and pull out the data of interest
# Only need to run once, afterwards just load object
for(i in 1:3){
  ebird.raw.list[[i]] <- auk_ebd(f_in[i], file_sampling = f_sample) %>%
    # 2. define filters
    auk_date(date = c("2016-01-01", "2022-12-31")) %>%
    auk_complete() %>%
    auk_species(unique(birds.codes$Scientific)) %>%
    # 3. run filtering
    auk_filter(file = f_out[i], overwrite = T,
               file_sampling = f_sample_out[i]) %>%
    # 4. read text file into r data frame
    auk_zerofill(unique = T, collapse = T)
}

#Provides lat/longs to filter any unusual points (incorrect lat/long entry?). only a few
ebird.bbox <- st_bbox(grid.surface %>% st_transform(4326))

#Clean the raw data 
ebird.full.df <- do.call("rbind", ebird.raw.list) %>%
  select(ID = checklist_id, Scientific = scientific_name.x, Count = observation_count.x,
         latitude, longitude, observation_date, Minutes = duration_minutes, 
         DistKM = effort_distance_km, NObs = number_observers) %>%
  filter(latitude >= ebird.bbox[2] & latitude <= ebird.bbox[4] &
         longitude >= ebird.bbox[1] & longitude <= ebird.bbox[3]) %>%
  distinct() %>%
  group_by(ID, Scientific) %>%
  slice(1L) %>% #Some surveys (~6000) have two locations for a single checklist because it spanned two states. Just taking the first location for now
  ungroup()

#Corresponding grid ID for each survey location
ebird.grid.df <- ebird.full.df %>% select(ID, latitude, longitude) %>%
  distinct() %>%
  group_by(ID) %>%
  slice(1L) %>%
  ungroup() %>%  
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(32619) %>%
  st_intersection(., grid.surface) %>%
  arrange(ID) %>%
  st_drop_geometry()
ebird.grid <- ebird.grid.df$GridID

#Bird count data
ebird.counts.df <- ebird.full.df %>%
  select(ID, Scientific, Count) %>%
  filter(ID %in% ebird.grid.df$ID) %>%
  mutate(Count = as.numeric(ifelse(Count == "X", "0", Count))) %>%
  merge(., birds.codes %>% select(Alpha, Scientific), by = "Scientific", all.x = T) %>%
  select(-Scientific) %>%
  arrange(Alpha) %>%
  pivot_wider(values_from = "Count", names_from = "Alpha", values_fill = 0) %>%
  arrange(ID)
ebird.counts <- ebird.counts.df %>% select(-ID) %>% as.matrix()

#Effort Covariates
ebird.covs.df <- ebird.full.df %>% select(ID, Minutes, DistKM, NObs)  %>%
  filter(ID %in% ebird.grid.df$ID) %>%
  distinct() %>%
  mutate(Minutes = ifelse(is.na(Minutes), 1, Minutes),
         DistKM = ifelse(is.na(DistKM), 0, DistKM), 
         NObs = ifelse(is.na(NObs), 1, NObs)) 
ebird.covs <- ebird.covs.df%>%
  select(-ID) %>% as.matrix()

#Check that the rows are correctly aligned.
which(ebird.counts.df$ID != ebird.covs.df$ID)



#Indexing constants for NIMBLE
ebird.nsurvey <- nrow(ebird.counts)
ebird.nspecies <- ncol(ebird.counts)
ebird.ncovs <- ncol(ebird.covs)


### BBS
# downloaded directly from sciencebase
# https://www.sciencebase.gov/catalog/item/625f151ed34e85fa62b7f926
setwd("./../BBS/50-StopData")

#Process raw data files
bbs.files <- list.files(".")
bbs.raw50 <- vector("list", length(bbs.files))
for(i in 1:length(bbs.files)){
  bbs.raw50[[i]] <- read.csv(bbs.files[[i]]) %>%
    filter(Year %in% year.list,
           StateNum %in% c(21, 46, 88)
    )
}

#Vector of column names for subsetting df
na.to.0 <- vector("list", 50)
na.to.0[1:50] <- 0
names(na.to.0) <- paste0("Stop", 1:50)

#Clean the raw data
bbs.raw50.df <- do.call("rbind", bbs.raw50) %>% 
  mutate(RouteID = paste(CountryNum, StateNum, Route, sep = "_")) %>%
  select(-RouteDataID, -CountryNum, -StateNum, -Route, -RPID) %>%
  arrange(AOU, RouteID, Year) %>%
  complete(AOU, RouteID, Year) %>%
  mutate(SurveyID = paste(RouteID, Year, sep = "_")) %>%
  replace_na(na.to.0) %>%
  arrange(SurveyID)

# Count by route
bbs.counts.df <- bbs.raw50.df %>%
  merge(., birds.codes %>% select(Species = Alpha, AOU = Number), by = "AOU", all.x = T) %>%
  select(-AOU) %>%
  filter(!is.na(Species)) %>%
  mutate(SumCount = rowSums(.[,which(colnames(.) %in% paste0("Stop", 1:50))])) %>%
  select(SurveyID, RouteID, Species, SumCount) %>%
  pivot_wider(names_from = Species, values_from = SumCount) %>%
  arrange(SurveyID)
bbs.counts <- bbs.counts.df %>% select(-SurveyID, -RouteID) %>% as.matrix()

#Route Covariates - number of stops with positive id for a species in a given survey
bbs.nposid.df1 <- bbs.raw50.df %>%
  merge(., birds.codes %>% select(Species = Alpha, AOU = Number), by = "AOU", all.x = T) %>%
  select(-AOU) %>%
  filter(!is.na(Species)) %>%
  select(SurveyID, Species, paste0("Stop", 1:50))
bbs.nposid.df2 <- bbs.nposid.df1 %>% select(SurveyID, Species)
bbs.nposid.df <- bbs.nposid.df1 %>% select(paste0("Stop", 1:50)) %>%
  mutate_all(function(x) ifelse(x> 0,1,0)) %>%
  mutate(NPosID = rowSums(.)) %>%
  select(NPosID) %>% cbind(bbs.nposid.df2, .) %>%
  pivot_wider(names_from = Species, values_from = NPosID) %>%
  arrange(SurveyID)
bbs.nposid <- bbs.nposid.df %>% select(-SurveyID) %>% as.matrix()

#Route Covariates - number of species observed during entire survey
bbs.nseen.df <- bbs.raw50.df %>%
  select(SurveyID, AOU, paste0("Stop", 1:50)) %>%
  mutate(SumCount = rowSums(.[,which(colnames(.) %in% paste0("Stop", 1:50))])) %>%
  mutate(Obs = ifelse(SumCount > 0, 1, 0)) %>%
  select(SurveyID, AOU, Obs) %>%
  group_by(SurveyID) %>%
  summarize(NSeen = sum(Obs)) %>%
  arrange(SurveyID)
bbs.nseen <- bbs.nseen.df$NSeen

#Grid Cell Index Object
bbs.surveys <- read.csv("./../routes.csv") %>% 
  mutate(RouteID = paste(CountryNum, StateNum, Route, sep = "_")) %>%
  filter(RouteID %in% bbs.raw50.df$RouteID) %>%
  select(RouteID, Latitude, Longitude) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4269) %>%
  st_transform(32619) %>%
  st_intersection(., grid.surface) %>%
  st_drop_geometry() %>%
  merge(bbs.counts.df %>% select(SurveyID, RouteID), ., by = "RouteID", all.x = T) %>%
  arrange(SurveyID)
bbs.grid <- bbs.surveys$GridID  

#Indexing constants for NIMBLE
bbs.nspecies <- ncol(bbs.counts)
bbs.nsurveys <- nrow(bbs.counts)


### BBL
# downloaded directly from sciencebase
# https://www.sciencebase.gov/catalog/item/632b2d7bd34e71c6d67bc161
setwd("./../../BBL/")

#Process the raw files
bbl.files <- list.files(".")
bbl.list.data <- list()
for(i in 1:length(bbl.files)){
  bbl.list.data[[i]] <- read.csv(bbl.files[i]) %>%
    select(Date = EVENT_DATE, State = ISO_SUBDIVISION, Lat = LAT_DD, Lon = LON_DD, Number = SPECIES_ID, HOW_OBTAINED) %>%
    mutate(Date = as.Date(Date, format = "%m/%d/%Y", origin = "1970-01-01")) %>%
    filter(year(Date) %in% year.list, 
           State %in% c("US-MD", "US-VA", "US-DE"))
}

#Clean the raw data
bbl.count.raw <- do.call("rbind", bbl.list.data) %>%
  # filter(Number %in% birds.codes$Number) %>%
  group_by(Date, Lat, Lon, Number) %>%
  filter(!is.na(Lat) & Lat != 0) %>%
  summarize(Count = n()) %>%
  merge(., bird.codes.all %>% select(Number, Species = Alpha), by = "Number", all.x = T) %>%
  select(Species, Date, Lat, Lon, Count) %>%
  arrange(Species) %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4269) %>%
  st_transform(32619) %>%
  st_intersection(., grid.surface) %>%
  st_drop_geometry() %>%
  group_by(Species, GridID) %>%
  summarize(TotalCap = sum(Count)) %>%
  pivot_wider(names_from = "Species", values_from = "TotalCap", values_fill = 0) %>%
  merge(., grid.surface %>% st_drop_geometry() %>% select(GridID), by = "GridID", all = T) %>%
  arrange(GridID)
bbl.count.raw[is.na(bbl.count.raw)]<- 0

#Site occupancy status
bbl.count <- bbl.count.raw %>% 
  select(which(colnames(bbl.count.raw) %in% birds.use)) %>%
  select(which(colSums(.) > 0)) %>%
  as.matrix()

#Effort Covariates
#Number of captures conducted within a grid cell
bbl.covs1 <- do.call("rbind", bbl.list.data) %>%
  select(Date, Lat, Lon) %>%
  distinct() %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4269) %>%
  st_transform(32619) %>%
  st_intersection(., grid.surface) %>%
  st_drop_geometry() %>%
  group_by(GridID) %>%
  summarize(NCapEvent = n()) %>%
  merge(., grid.surface %>% st_drop_geometry() %>% select(GridID), by = "GridID", all = T) %>%
  arrange(GridID) %>%
  mutate(NCapEvent = ifelse(is.na(NCapEvent), 0, NCapEvent))

#Number of birds captured
bbl.covs.df <- do.call("rbind", bbl.list.data) %>%
  select(Date, Lat, Lon) %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4269) %>%
  st_transform(32619) %>%
  st_intersection(., grid.surface) %>%
  st_drop_geometry() %>%
  group_by(GridID) %>%
  summarize(NCaptured = n()) %>%
  merge(., grid.surface %>% st_drop_geometry() %>% select(GridID), by = "GridID", all = T)%>%
  mutate(NCaptured = ifelse(is.na(NCaptured), 0, NCaptured)) %>%
  merge(bbl.covs1, ., by = "GridID", all = T) %>%
  arrange(GridID) 
bbl.covs <- bbl.covs.df %>% select(-GridID) %>% as.matrix()


#Nimble Indexing constants
bbl.nspecies <- ncol(bbl.count)
bbl.nsurveys <- nrow(bbl.count)
bbl.ncovs <- ncol(bbl.covs)


### MWS
# Directly aquired (through Jeff)
setwd("./../MWS/")

#Clean and organize original data
mws.raw <- read.csv("MWS_ChesBay.csv") %>%
  filter(Include == "Yes", 
         SurveyYear != "AVG") 
mws.surveyinfo <- mws.raw %>%
  select(Name, SurveyYear, Total_Minutes, )

#Generate Raw Count Data
mws.countonly <- mws.raw %>%
  select(Name, SurveyYear, 27:ncol(.))
mws.countonly[is.na(mws.countonly)] <- 0
mws.counts.df <- mws.countonly %>%
  select(Name, SurveyYear, 3:ncol(mws.countonly)) %>%
  pivot_longer(names_to = "Alpha", values_to = "Count", cols = 3:ncol(.)) %>%
  arrange(Alpha, Name, SurveyYear) %>%
  filter(Alpha %in% birds.codes$Alpha) %>%
  pivot_wider(names_from = "SurveyYear", names_prefix = "Year", values_from = "Count")
mws.counts.list <- mws.counts.df %>%
  select(-Name) %>%
  group_split(Alpha, .keep = F)

#3 dimensional data (survey, species, year) requires array formatting 
mws.counts <- array(unlist(mws.counts.list), dim = c(nrow(mws.counts.list[[1]]),
                                                     ncol(mws.counts.list[[1]]),
                                                     length(mws.counts.list)))

#Aerial surveys cover much larger area than other surveys
#Instead of assessing a single grid cell, assess all overlapping grid cells simultaneously
#Create Grid Cell Index vector/matrix combo for working with list
mws.survey.sf <- st_read("ChesapeakeBay_Combined_WithBuffer.shp") %>%
  select(Name, SurveyLength = Shape_Leng ) %>%
  st_transform(32619) %>%
  arrange(Name)

mws.grid.df <- st_intersection(mws.survey.sf, grid.surface) %>%
  mutate(intersect_area = st_area(.)) %>%
  group_by(Name) %>%
  filter(max(intersect_area) == intersect_area) %>%
  st_drop_geometry() %>% 
  select(Name, GridID) %>%
  arrange(Name) %>%
  ungroup() %>%
  merge(mws.counts.df %>% filter(Alpha == "AMWI") %>% select(Name), ., by = "Name", all.x = T)

mws.gridid <- mws.grid.df$GridID

# #This is code I was trying to use to aggregate grid cells within the model. 
# #Decided to use greatest overlap for now.
# mws.ngrids <- mws.grid.df %>% summarize(NGrid = n()) %>% .$NGrid
# 
# mws.gridid <- mws.grid.df %>%
#   mutate(row = row_number()) %>%
#   pivot_wider(values_from = GridID, names_from = row, names_prefix = "Grid") %>%
#   as.data.frame() %>%
#   select(-Name) %>% as.matrix()

#Effort Covariates
mws.covs.df <- mws.countonly %>% rowwise() %>%
  mutate(NObs = sum(c_across(3:ncol(.)))) %>%
  select(Name, SurveyYear, NObs) %>%
  merge(., mws.survey.sf %>% st_drop_geometry(), by = "Name", all.x = T)

mws.covs <- array(NA, dim = c(length(unique(mws.covs.df$Name)),
                              length(unique(mws.covs.df$SurveyYear)),
                              2))
#Number of observers
mws.covs[,,1] <- mws.covs.df %>% select(Name, SurveyYear, NObs) %>%
  pivot_wider(names_from = "SurveyYear", values_from = "NObs", values_fill = 0) %>%
  select(-Name) %>% as.matrix()
#Survey distance
mws.covs[,,2] <- mws.covs.df %>% select(Name, SurveyYear, SurveyLength) %>%
  pivot_wider(names_from = "SurveyYear", values_from = "SurveyLength") %>%
  mutate(`2018` = coalesce(`2016`, `2018`)) %>%
  select(-Name) %>% as.matrix()

#Indexing objects for Nimble
mws.nsurveys <- dim(mws.covs)[1]
mws.nyears <- dim(mws.covs)[2]
mws.nspecies <- dim(mws.counts)[3]
mws.ncovs <- dim(mws.covs)[3]


### Christmas Bird Count Data Prep
#Data request submitted through https://survey123.arcgis.com/share/7dc33b4fff77468a8bba855291f86527
setwd("./../CBC/")

#Load raw data
cbc.count.raw <- read.csv("20230110_Matthew_Gonnerman_CB_Circle_Species_Report.csv")
cbc.effort <- read.csv("20230110_Matthew_Gonnerman_CBC_Effort_Many_Types_Report.csv")
cbc.effortsummary <- read.csv("20230110_Matthew_Gonnerman_CBC_Effort_Summary_Report.csv")

#Create count matrix
cbc.count.df <- cbc.count.raw %>%
  mutate(ID = paste(abbrev, count_yr, sep = "_")) %>%
  select(ID, Scientific = sci_name, 
         Count = how_many, NSpecies = totalspecies) %>%
  distinct() %>%
  complete(ID, Scientific) %>%
  group_by(ID) %>%
  fill(NSpecies, .direction = "downup") %>%
  filter(Scientific %in% birds.codes$Scientific) %>%
  mutate(Count = ifelse(is.na(Count), 0, Count)) %>%
  merge(., birds.codes %>% select(Alpha, Scientific), by = "Scientific", all.x = T) %>%
  select(-Scientific)%>%
  pivot_wider(values_from = "Count", names_from = "Alpha") %>%
  arrange(ID) 
cbc.count <- cbc.count.df %>% select(-ID, -NSpecies) %>% as.matrix()

#Effort covariates
#Number of species observed during a survey
cbc.obsspecies <- cbc.count.df %>% select(ID, NSpecies)

# Total distance covered during survey
cbc.covs.df <- cbc.effort %>%
  select(Site = abbrev, Year = count_yr, distance) %>%
  mutate(ID = paste(Site, Year, sep = "_")) %>%
  group_by(ID) %>%
  summarize(TDist = sum(distance)) %>%
  ungroup() %>% arrange(ID) %>%
  merge(., cbc.obsspecies, by = "ID", all = T)
cbc.covs <- cbc.covs.df %>% select(-ID) %>% as.matrix()
cbc.ncovs <- ncol(cbc.covs)

#GridID indexing object
cbc.grid.df <- cbc.count.raw %>%
  mutate(ID = paste(abbrev, count_yr, sep = "_")) %>%
  select(ID, latitude, longitude) %>%
  distinct() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(32619) %>%
  st_intersection(., grid.surface) %>%
  st_drop_geometry() %>%
  merge(cbc.count.df %>% select(ID), ., by = "ID", all.x = T) %>%
  arrange(ID)
cbc.gridid <- cbc.grid.df$GridID

#Nimble indexing objects
cbc.nspecies <- ncol(cbc.count)
cbc.nsurveys <- nrow(cbc.count)


### Final Species List
#full list of unique species across surveys
final.species.list <- sort(unique(c(colnames(ebird.counts), colnames(bbs.counts), colnames(bbl.count),
         mws.counts.df$Alpha, colnames(cbc.count))))
ntotspecies <- length(final.species.list)
# match(colnames(ebird.counts), final.species.list) #1 to 1, so don't need to adjust indexing 
bbs.sp <- match(colnames(bbs.counts), final.species.list)
bbl.sp <- match(colnames(bbl.count), final.species.list)
mws.sp <- match(unique(mws.counts.df$Alpha),final.species.list)
cbc.sp <- match(colnames(cbc.count), final.species.list)


### Grid Cell Covariates
#Used Zonal Statistics as Table in arcgis to generate these values
setwd("./../data/")
ag.grid <- read.csv("Ag_grid.csv") %>% select(GridID, Ag = SUM) 
dev.grid <- read.csv("Dev_grid.csv") %>% select(GridID, Dev = SUM)  
body.grid <- read.csv("Body_grid.csv") %>% select(GridID, Body = SUM)  
flow.grid <- read.csv("Flow_grid.csv") %>% select(GridID, Flow = SUM)  
wetland.grid <- read.csv("Wetland_grid.csv") %>% select(GridID, Wetland = SUM)  
  
grid.covs.df <- merge(ag.grid, dev.grid, by = "GridID", all.x = T) %>%
  merge(., body.grid, by = "GridID", all.x = T) %>%
  merge(., flow.grid, by = "GridID", all.x = T) %>%
  merge(., wetland.grid, by = "GridID", all.x = T) %>%
  arrange(GridID) %>% 
  mutate(Ag = ifelse(is.na(Ag), 0, Ag)) %>% mutate(Ag = scale(Ag)[,1]) %>% 
  mutate(Dev = ifelse(is.na(Dev), 0, Dev)) %>% mutate(Dev = scale(Dev)[,1])%>% 
  mutate(Body = ifelse(is.na(Body), 0, Body)) %>% mutate(Body = scale(Body)[,1])%>% 
  mutate(Flow = ifelse(is.na(Flow), 0, Flow)) %>% mutate(Flow = scale(Flow)[,1])%>% 
  mutate(Wetland = ifelse(is.na(Wetland), 0, Wetland)) %>% mutate(Wetland = scale(Wetland)[,1])
grid.covs <- grid.covs.df %>% select(-GridID) %>% as.matrix()
ncovs.grid <-  ncol(grid.covs)

### Save relevant objects used in NIMBLE Model
setwd("./../")
save(ebird.counts, ebird.covs, ebird.grid, ebird.nsurvey, ebird.nspecies, ebird.ncovs, #EBird
     bbs.counts, bbs.nspecies, bbs.nsurveys, bbs.nposid, bbs.nseen, bbs.grid, bbs.sp, #BBS
     bbl.count, bbl.covs, bbl.nsurveys, bbl.nspecies, bbl.ncovs, bbl.sp, #BBL
     mws.counts, mws.covs, mws.gridid, mws.nsurveys, mws.nyears, mws.nspecies, mws.ncovs, mws.sp,#MWS
     cbc.count, cbc.covs, cbc.nsurveys, cbc.nspecies, cbc.ncovs, cbc.gridid, cbc.sp, #CBC
     grid.covs, ncovs.grid, ntotspecies, n.cells, adj, num, weights, nadj,
     file = "./iSDM_Nimble_Objects.R")

