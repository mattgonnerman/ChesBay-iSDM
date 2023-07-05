lapply(c("dplyr", "sf", "auk", "lubridate", "tidyr", "ggplot2"), require, character.only = T)

#### TO DO LIST #####
### Double check all surveys for miss-specified names
### Distance to coast metric would be good
### Filter EBIRD by best practices

### Set working directory depending on which computer you are using
# setwd("E:/GitHub/ChesBay iSDM/") #Home
setwd("C:/Users/mgonnerman/Box/Prosser_USGS/Matt Projects/ChesBay iSDM/") #Laptop

### Years of Interest
year.list <- 2016:2021

### Specify Study Area Boundaries
studyarea <- st_read("./Data/StudyArea.shp") %>%
  st_transform(32619)

###List of species to track
bird.codes.all <- read.csv("./BirdCodes.csv") %>% arrange(Common)
birds.use <- sort(c("MALL", "ABDK", "ABDU", "CAGO", "BWTE", "AMWI", "NOPI", 
                    "GWTE", "AGWT","GADW", "NSHO", "WODU", "GSGO", "GWFG",
                    "BLSC", "SUSC", "LTDU", "LESC", "RUDU"))
birds.codes <- bird.codes.all %>% filter(Alpha %in% birds.use)

### EBird
### EBIRD DATA BEST PRACTICES
### Strimas-Mackey et al. 2020
###  (1) only completed checklists (i.e., all the birds observed are reported); 
###  (2) surveys conducted between 5:00 a.m. and 9:00 p.m.; 
###  (3) covering <5 km; 
###  (4) collected over a period of no longer than 5 h; and 
###  (5) with no more than 10 observers.

#Data accessed directly through the EBird download portal https://ebird.org/data/download
setwd("E:/ChesBay iSDM/EBird/")
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
# for(i in 1:3){
#   ebird.raw.list[[i]] <- auk_ebd(f_in[i], file_sampling = f_sample) %>%
#     # 2. define filters
#     auk_date(date = c("2016-01-01", "2022-12-31")) %>%
#     auk_complete() %>%
#     auk_species(unique(birds.codes$Scientific)) %>%
#     # 3. run filtering
#     auk_filter(file = f_out[i], overwrite = T,
#                file_sampling = f_sample_out[i]) %>%
#     # 4. read text file into r data frame
#     auk_zerofill(unique = T, collapse = T)
# }

ebird.raw.list <- lapply(f_out, FUN = read_ebd)

#Provides lat/longs to filter any unusual points (incorrect lat/long entry?). only a few
ebird.bbox <- st_bbox(studyarea %>% st_transform(4326))

#Clean the raw data 
ebird.full.df <- do.call("rbind", ebird.raw.list) %>%
  #Filter to all species reported so you can treat as count, keep only approved
  filter(approved == T, all_species_reported == T, effort_distance_km <= 5) %>%
  #Filter by best practices (see above) %>%
  filter(duration_minutes < 60*5, number_observers < 11) %>%
  filter(!is.na(time_observations_started)) %>%
  filter(hms(time_observations_started) > hms("05:00:00") & 
           hms(time_observations_started) < hms("17:00:00")) %>%
  select(ID = checklist_id, Scientific = scientific_name, Count = observation_count,
         latitude, longitude, observation_date, Minutes = duration_minutes, 
         DistKM = effort_distance_km, NObs = number_observers) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_intersection(., studyarea %>% st_transform(4326)) %>%
  mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  distinct() %>%
  group_by(ID, Scientific) %>%
  arrange(ID, Scientific) %>%
  ungroup %>% distinct(ID, Scientific, Count, .keep_all = T) %>% #Some surveys (~6000) have two locations for a single checklist because it spanned two states. Just taking the first location for now
  #If a person didn't record a count but was present, its a "X" for present (and keep track for detection?). 
  mutate(RealCount = ifelse(Count == "X", F, T), 
         Count = as.numeric(ifelse(Count == "X", "1", Count))) %>%
  merge(., birds.codes %>% select(AOU = Alpha, Scientific), by = "Scientific", all.x = T) %>%
  select(-Scientific, ID, AOU, X, Y, Date = observation_date, Minutes, DistKM, NObs) %>%
  # pivot_wider(values_from = "Count", names_from = "Alpha", values_fill = 0) %>%
  # arrange(ID) %>%
  mutate(Minutes = ifelse(is.na(Minutes), 1, Minutes),
         DistKM = ifelse(is.na(DistKM), 0, DistKM), 
         NObs = ifelse(is.na(NObs), 1, NObs)) %>%
  mutate(Occ = ifelse(Count > 0, 1, 0),
         Survey = "EBird") %>%
  select(Survey, AOU, X, Y, ChecklistID = ID, Date, Count, Occ, DistKM, NObs, Minutes)

ebird.surv.info <- ebird.full.df %>% 
  select(ChecklistID, Survey, X, Y, Date, DistKM, NObs, Minutes) %>%
  distinct()

# Not all species observed in BBS, add 0s to represent no observations
ebird.obs.spp <- ebird.full.df %>% select(ChecklistID, AOU, Count, Occ) %>%
  complete(ChecklistID, AOU, fill = list(Count = 0, Occ = 0))
ebird.data.full <- merge(ebird.obs.spp, ebird.surv.info, by = c("ChecklistID"), all.x = T) %>%
  select(Survey, AOU, X, Y, ChecklistID, Date, Count, Occ, DistKM, NObs, Minutes)

write.csv(ebird.full.df, 
          "C:/Users/mgonnerman/Box/Prosser_USGS/Matt Projects/ChesBay iSDM/Data/EBird_Counts.csv",
          row.names = F)


### BBS
# downloaded directly from sciencebase
# https://www.sciencebase.gov/catalog/item/625f151ed34e85fa62b7f926
setwd("E:/ChesBay iSDM/BBS/50-StopData")

bbs.codes.all <- bbs.codes <- read.csv("./../../BBL/BirdCodes.csv") %>%
  select(Number = Species.Number, AOU = Alpha.Code, Common = Common.Name, Scientific = Scientific.Name)
bbs.codes <- bbs.codes.all %>%
  filter(AOU %in% birds.use)

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
  rename(Number = AOU) %>%
  mutate(RouteID = paste(CountryNum, StateNum, Route, sep = "_")) %>%
  select(-RouteDataID, -CountryNum, -StateNum, -Route, -RPID) %>%
  arrange(Number, RouteID, Year) %>%
  complete(Number, RouteID, Year) %>%
  mutate(SurveyID = paste(RouteID, Year, sep = "_")) %>%
  replace_na(na.to.0) %>%
  arrange(SurveyID) 

#Consolidate information across stops within a given survey route
bbs.data.obsspp <- read.csv("./../routes.csv") %>% 
  mutate(RouteID = paste(CountryNum, StateNum, Route, sep = "_")) %>%
  select(RouteID, Latitude, Longitude) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4269) %>%
  st_intersection(., studyarea %>% st_transform(4269)) %>%
  st_transform(4326) %>%
  mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  arrange(RouteID) %>%
  select(RouteID, X, Y) %>%
  merge(., bbs.raw50.df, by = "RouteID", all.x = T) %>%
  filter(Number %in% bbs.codes$Number) %>%
  merge(., bbs.codes %>% select(Number, AOU), by = "Number", all.x = T)%>%
  mutate(Count = rowSums(.[,which(colnames(.) %in% paste0("Stop", 1:50))]))  %>%
  relocate(SurveyID, AOU) %>%
  select(-paste0("Stop", 1:50), -Number, -SurveyID) %>%
  arrange(RouteID, AOU, Year)

# Not all species observed in BBS, add 0s to represent no observtions 
bbs.data.absspp <- expand.grid(AOU = bbs.codes$AOU,
                               RouteID = unique(bbs.data.obsspp$RouteID),
                               Year = unique(bbs.data.obsspp$Year)) %>%
  filter(!AOU %in% bbs.data.obsspp$AOU) %>%
  merge(., bbs.data.obsspp %>% select(RouteID, X, Y), by = "RouteID", all.x = T) %>%
  mutate(Count = 0)

#Route Covariates - number of stops with positive id for a species in a given survey
bbs.nposid <- bbs.raw50.df %>%
  select(RouteID, Year, Number, paste0("Stop", 1:50)) %>%
  mutate_at(paste0("Stop", 1:50), function(x) ifelse(x> 0,1,0)) %>%
  mutate(NPosID = rowSums(across(paste0("Stop", 1:50)))) %>%
  select(RouteID, Year, Number, NPosID)

bbs.nspp <- bbs.nposid %>% 
  mutate(NPosID = ifelse(NPosID > 0, 1, 0)) %>%
  group_by(RouteID, Year) %>%
  summarize(NSppObs = sum(NPosID))

bbs.covs <- merge(bbs.nposid, bbs.nspp, by = c("RouteID", "Year"), all.x = T) %>%
  merge(., bbs.codes.all %>% select(Number, AOU), by = "Number", all.x = T) %>%
  select(RouteID, Year, AOU, NStopPos = NPosID, NSppObs) %>%
  filter(AOU %in% bbs.codes$AOU)

# Extract Date from the Weather file
bbs.date <- read.csv("./../weather.csv") %>%
  mutate(RouteID = paste(CountryNum, StateNum, Route, sep = "_"),
         Date = paste(Year, formatC(Month, digits = 1, format = "d", flag = "0"), 
                      formatC(Day, digits = 1, format = "d", flag = "0"), sep = "-")) %>%
  filter(as.numeric(Year) %in% year.list) %>% select(RouteID, Date, Year)

# Combine observed and absent df's
bbs.data.final <- rbind(bbs.data.obsspp, bbs.data.absspp) %>%
  merge(., bbs.date, by = c("RouteID", "Year")) %>%
  merge(., bbs.covs, by = c("RouteID", "Year", "AOU")) %>%
  mutate(RowID = row_number(), Survey = "BBS",
         Occ = ifelse(Count > 0, 1, 0)) %>%
  select(Survey, RowID, AOU, X, Y, RouteID, Date, Count, Occ, NStopPos, NSppObs)
write.csv(bbs.data.final, 
          "C:/Users/mgonnerman/Box/Prosser_USGS/Matt Projects/ChesBay iSDM/Data/BBS_Counts.csv",
          row.names = F)

### BBL
# downloaded directly from sciencebase
# https://www.sciencebase.gov/catalog/item/632b2d7bd34e71c6d67bc161
setwd("E:/ChesBay iSDM/BBL/")

bbl.codes <- read.csv("BirdCodes.csv") %>%
  select(Number = Species.Number, AOU = Alpha.Code, Common = Common.Name, Scientific = Scientific.Name) %>%
  filter(AOU %in% birds.use)

#Process the raw files
bbl.files <- list.files(".")[!grepl("BirdCodes", list.files("."))]
bbl.list.data <- list()
for(i in 1:length(bbl.files)){
  bbl.list.data[[i]] <- read.csv(bbl.files[i]) %>%
    select(Date = EVENT_DATE, State = ISO_SUBDIVISION, Lat = LAT_DD, Lon = LON_DD, Number = SPECIES_ID, HOW_OBTAINED) %>%
    mutate(Date = as.Date(Date, format = "%m/%d/%Y", origin = "1970-01-01")) %>%
    filter(year(Date) %in% year.list, 
           State %in% c("US-MD", "US-VA", "US-DE"))
}

#Clean the raw data
bbl.count <- do.call("rbind", bbl.list.data) %>%
  # filter(Number %in% birds.codes$Number) %>%
  group_by(Date, Lat, Lon, Number) %>%
  filter(!is.na(Lat) & Lat != 0) %>%
  summarize(Count = n()) %>%
  merge(., bbl.codes %>% select(Number, AOU), by = "Number", all.x = T) %>%
  filter(!is.na(AOU)) %>%
  select(AOU, Date, Lat, Lon, Count) %>%
  group_by(AOU, Date, Lat, Lon) %>%
  summarize(Count = sum(Count, na.rm = T)) %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4269) %>%
  st_intersection(., studyarea %>% st_transform(4269)) %>%
  st_transform(4326) %>%
  mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>% select(AOU, Date, X, Y, Count) %>%
  mutate(RowID = row_number(), Survey = "BBL")

bbl.final <- do.call("rbind", bbl.list.data) %>%
  merge(., bbl.codes %>% select(Number, AOU), by = "Number", all.x = T) %>%
  select(X = Lon, Y = Lat, Date, AOU) %>%
  group_by(Y, X, Date) %>%
  summarize(NCap = n()) %>%
  merge(bbl.count, ., by = c("X", "Y", "Date"), all.x = T) %>%
  mutate(Occ = ifelse(Count > 0, 1, 0)) %>%
  select(Survey, RowID, AOU, X, Y, Date, Count, Occ, NCap)

write.csv(bbl.final, 
          "C:/Users/mgonnerman/Box/Prosser_USGS/Matt Projects/ChesBay iSDM/Data/BBL_Counts.csv",
          row.names = F)

# ### MWS
# # Directly aquired (through Jeff)
# setwd("E:/ChesBay iSDM/MWS/")
# 
# #Clean and organize original data
# mws.counts <- read.csv("MWS_ChesBay.csv") %>%
#   filter(Include == "Yes", 
#          SurveyYear != "AVG") %>%
#   select(Name, Year = SurveyYear, Date = SurveyDate, 27:ncol(.)) %>%
#   pivot_longer(names_to = "AOU", values_to = "Count", cols = 4:ncol(.)) %>%
#   arrange(AOU, Name, Year) %>%
#   filter(AOU %in% birds.codes$Alpha, Year %in% year.list) %>%
#   mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
#   group_by(Year) %>%
#   mutate(Year = as.numeric(Year),
#          Date = as.Date(ifelse(is.na(Date), mean.Date(Date, na.rm = T), Date), origin = "1970-01-01")) %>%
#   filter(!is.na(Count))
# 
# 
# #Aerial surveys cover much larger area than other surveys
# #Instead of assessing a single grid cell, assess all overlapping grid cells simultaneously
# #Create Grid Cell Index vector/matrix combo for working with list
# mws.xy <- st_read("ChesapeakeBay_Combined_WithBuffer.shp") %>%
#   select(Name)  %>%
#   mutate(SurvArea = st_area(.)/1000000) %>%
#   st_point_on_surface() %>% st_transform(4326) %>%
#   mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>%
#   st_drop_geometry()
# 
# mws.data.final <- merge(mws.counts, mws.xy, by = "Name", all.x = T) %>%
#   arrange(AOU, Name, Year) %>%
#   mutate(RowID = row_number(), Survey = "MWS",
#          Occ = ifelse(Count > 0, 1, 0),
#          SurvArea = as.numeric(SurvArea)) %>%
#   select(Survey, RowID, AOU, X, Y, SurveyID = Name, Date, Count, Occ, SurvArea)
# 
# write.csv(mws.data.final, 
#           "C:/Users/mgonnerman/Box/Prosser_USGS/Matt Projects/ChesBay iSDM/Data/MWS_Counts.csv",
#           row.names = F)

### Christmas Bird Count Data Prep
#Data request submitted through https://survey123.arcgis.com/share/7dc33b4fff77468a8bba855291f86527
setwd("E:/ChesBay iSDM/CBC/")

#CBC has slightly different names, so adjust accordingly. 
#If something is unexpectedly missing, check the names
cbc.species.filter <- data.frame(Scientific = c(birds.codes$Scientific, "Anas strepera", "Anas americana", "Anas discors", "Anas clypeata", "Chen caerulescens"),
                                 AOU = c(birds.codes$Alpha, "GADW", "AMWI", "BWTE", "NSHO", "GSGO"))

#Load raw data
cbc.count.raw <- read.csv("20230110_Matthew_Gonnerman_CB_Circle_Species_Report.csv") %>%
  mutate(Date = as.POSIXct(cnt_dt, format = "%m/%d/%Y %H:%M:%S")) %>% 
  mutate(Year = year(Date)) %>%
  filter(Year %in% year.list)

cbc.effort.raw <- read.csv("20230110_Matthew_Gonnerman_CBC_Effort_Many_Types_Report.csv")

cbc.surveys <- unique(cbc.effort.raw$abbrev)

cbc.count <- cbc.count.raw %>%
  select(SiteID = abbrev, X = longitude, Y = latitude, Year, Date, Scientific = sci_name, Count = how_many) %>%
  merge(., cbc.species.filter, by = "Scientific", all.x = T) %>% select(-Scientific) %>%
  complete(AOU, SiteID, Year, fill = list(Count = 0)) %>%
  arrange(AOU, SiteID, Year) %>%
  group_by(SiteID) %>%
  mutate(X = unique(X[!is.na(X)]),
         Y = unique(Y[!is.na(Y)])) %>%
  group_by(SiteID, Date) %>%
  filter(!all(is.na(Date))) %>%
  mutate(Date = unique(Date[!is.na(Date)])) %>%
  filter(!is.na(AOU)) 

#Number of species observed during a survey
cbc.obsspecies <- cbc.count.raw %>%
  select(SiteID = abbrev, Date, NSppObs = totalspecies) %>%
  distinct()

# Total distance covered during survey
cbc.distance <- cbc.effort.raw %>%
  select(SiteID = abbrev, Year = count_yr, distance) %>%
  mutate(Year = Year + 1900) %>%
  group_by(SiteID, Year) %>%
  summarize(TDist = sum(distance)) 

  
cbc.data.final <- merge(cbc.count, cbc.obsspecies, by = c("SiteID", "Date"), all.x=T)  %>%
  merge(., cbc.distance, by = c("SiteID", "Year"), all.x=T)  %>%
  st_as_sf(., coords = c("X", "Y"), crs = 4326) %>%
  st_intersection(., studyarea %>% st_transform(4326)) %>%
  mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  mutate(RowID = row_number(), Survey = "CBC",
         Occ = ifelse(Count > 0, 1, 0)) %>%
  select(Survey, RowID, AOU, X, Y, SiteID, Date, Count, Occ, TDist, NSppObs)

write.csv(cbc.data.final, 
          "C:/Users/mgonnerman/Box/Prosser_USGS/Matt Projects/ChesBay iSDM/Data/CBC_Counts.csv",
          row.names = F)


### Create object with all locations to extract covariates
## Ensure you aren't extracting for duplicates
setwd("C:/Users/mgonnerman/Box/Prosser_USGS/Matt Projects/ChesBay iSDM/Data/")
ebird.sf <- read.csv("EBird_Counts.csv") %>% select(X, Y) %>% distinct()
bbs.sf <- read.csv("BBS_Counts.csv") %>% select(X, Y) %>% distinct()
bbl.sf <- read.csv("BBL_Counts.csv") %>% select(X, Y) %>% distinct()
# mws.sf <- read.csv("MWS_Counts.csv") %>% select(X, Y) %>% distinct()
cbc.sf <- read.csv("CBC_Counts.csv") %>% select(X, Y) %>% distinct()

isdm.sf <- rbind(ebird.sf, bbs.sf, bbl.sf, cbc.sf) %>%
  distinct() %>%
  st_as_sf(., coords = c("X", "Y"), crs = 4326)

st_write(isdm.sf, dsn = "E:/ChesBay iSDM/SurveyLocationsAll.shp", delete_layer = T)
#Use extract tool in ArcGIS to get values of interest. 
