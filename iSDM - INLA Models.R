#Load relevant packages
lapply(c("dplyr", "sf", "lubridate", "INLA", "ggplot2", "fields", "splancs"), require, character.only =T)

#Set Working Directory
setwd("C:/Users/mgonnerman/Box/Prosser_USGS/Matt Projects/ChesBay iSDM/") 

#Clear Environment
rm(list = ls())

#Study Area
studyarea <- st_read("Data/StudyArea.shp") %>%
  st_transform(4326)

birds.use <- sort(c("MALL", "ABDK", "ABDU", "CAGO", "BWTE", "AMWI", "NOPI", 
                    "GWTE", "AGWT","GADW", "NSHO", "WODU", "GSGO", "GWFG",
                    "BLSC", "SUSC", "LTDU", "LESC", "RUDU"))

#Count Data
ebird.data <- read.csv("Data/EBird_Counts.csv")
bbs.data <- read.csv("Data/BBS_Counts.csv")
bbl.data <- read.csv("Data/BBL_Counts.csv")
# mws.data <- read.csv("Data/MWS_Counts.csv")
cbc.data <- read.csv("Data/CBC_Counts.csv")

#Environmental Covariates
isdm.sf <- st_read("E:/ChesBay iSDM/SurveyLocationsAll.shp") %>%
  mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  relocate(X, Y) %>%
  mutate_all(~ ifelse(. == -9999, 0, .)) %>% select(-FID_1)

#For within belore function. Provide temporal periods
temp.periods <- function(df){df %>% select(RowID, Date) %>%
    mutate(Date = as.Date(Date)) %>%
    mutate(Week2 = ceiling(yday(Date)/14), Month = month(Date)) %>%
    mutate(Season = with(., case_when((Month >= 12 | Month <= 2)~1,
                                      (Month > 2 & Month <= 5) ~ 2,
                                      (Month >= 6 & Month <= 8) ~ 3,
                                      (Month >= 9 | Month <= 11) ~ 4,
                                      is.na(Month)~NA, TRUE~NA)))}

#Prepare SPDE Mesh
### Generate Mesh for SPDE
studyarea.sp <- spTransform(as(studyarea, "Spatial"), CRS("+proj=longlat +datum=WGS84"))
mesh <- inla.mesh.2d(boundary = inla.sp2segment(studyarea.sp),
                     cutoff = 0.01,
                     max.edge = c(0.1, 0.2))

#Create 1d temporal mesh for month, to = # of weeks in analysis
mesh.t <- inla.mesh.1d(seq(from = 1, to = 4, by = 1))
k <- mesh.t$n

# Generate spde components for each Season
spde <- inla.spde2.pcmatern(mesh = mesh,
                            constr = T,
                            prior.range = c(1, 0.5),
                            prior.sigma = c(1, 0.01))

#links temporal and spatial
iset <- inla.spde.make.index('i',
                             n.spde = spde$n.spde,
                             n.group = k)

#Save the mesh locations as a shapefile so you can get covariate values
#Only resave when you update the mesh parameters
mesh.coords <- st_as_sf(data.frame(X = mesh$loc[,1], Y = mesh$loc[,2]), coords = c("X", "Y"), crs = 4326)
# st_write(mesh.coords, dsn = "E:/ChesBay iSDM/IntegrationLocs.shp", delete_layer = T)
# We will still need to get environmental covariates at the mesh coords (integration points),
# So we will export locations as a shapefile and use extract by multipoint tool again
#Then read the file back in
integrat.covs <- st_read("E:/ChesBay iSDM/IntegrationLocs.shp") %>%
  mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  relocate(X, Y) %>%
  mutate_all(~ ifelse(. == -9999, 0, .)) %>%
  select(-FID_1)


### Code adapted from Ahmad Suhaimi et al 2021
### Create integration point objects
min_x <- min(isdm.sf$X)
min_y <- min(isdm.sf$Y)
max_x <- max(isdm.sf$X)
max_y <- max(isdm.sf$Y)
loc.d <- t(matrix(c(min_x,min_y,max_x,min_y,max_x,max_y,min_x,max_y,min_x,min_y), 2))
#make domain into spatial polygon
domainSP <- SpatialPolygons(list(Polygons(list(Polygon(loc.d)), '0')))
#intersection between domain and dual mesh
poly.gpc <- as(domainSP@polygons[[1]]@Polygons[[1]]@coords, "gpc.poly")
#make dual mesh
dd <- deldir::deldir(mesh$loc[, 1], mesh$loc[, 2])
tiles <- deldir::tile.list(dd)
w <- sapply(tiles, function(p) rgeos::area.poly(rgeos::intersect(as(cbind(p$x, p$y), "gpc.poly"), poly.gpc)))
nv <- mesh$n
#diagonal matrix for integration point A matrix
imat <- Diagonal(nv, rep(1, nv))

#Function to save all objects within a funciton to the global environment
allglobal <- function() list2env(mget(ls(name = parent.frame()),
                                      envir = parent.frame()),
                                 envir = .GlobalEnv)

### Filter and format data for a single species
spp_subset <- function(spp){
  ebird.data <- ebird.data %>% filter(AOU == spp)
  bbs.data <- bbs.data %>% filter(AOU == spp)
  bbl.data <- bbl.data %>% filter(AOU == spp)
  cbc.data <- cbc.data %>% filter(AOU == spp)
  # mws.data <- mws.data %>% filter(AOU == spp)
  
  ebird.det <- ebird.data %>% select(RowID, Minutes, DistKM, NObs)
  bbs.det <- bbs.data %>% select(RowID, NStopPos, NSO.BBS = NSppObs) 
  bbl.det <- bbl.data %>% select(RowID, NCap)
  cbc.det <- cbc.data %>% select(RowID, TDist, NSO.CBC = NSppObs)
  # mws.det <- mws.data %>% select(RowID, SurvArea)
  
  #Reduce to Points
  ebird.xy <- ebird.data %>% select(RowID, X, Y)
  bbs.xy <- bbs.data %>% select(RowID, X, Y)
  bbl.xy <- bbl.data %>% select(RowID, X, Y)
  cbc.xy <- cbc.data %>% select(RowID, X, Y)
  # mws.xy <- mws.data %>% select(RowID, X, Y)
  
  #Merge spatial covariates and counts
  ebird.env <- merge(ebird.xy, isdm.sf, by = c("X", "Y"), all.x = T) %>%
    select(-RowID)
  bbs.env <- merge(bbs.xy, isdm.sf, by = c("X", "Y"), all.x = T) %>%
    select(-RowID)
  bbl.env <- merge(bbl.xy, isdm.sf, by = c("X", "Y"), all.x = T) %>%
    select(-RowID)
  # mws.env <- merge(mws.xy, isdm.sf, by = c("X", "Y"), all.x = T) %>%
  #   select(-RowID)
  cbc.env <- merge(cbc.xy, isdm.sf, by = c("X", "Y"), all.x = T) %>%
    select(-RowID)
  
  #Save cov names for function
  env.covs.names <- colnames(isdm.sf)[3:(ncol(isdm.sf))]
  det.covs.names <- c(colnames(bbs.det)[-1],
                      colnames(bbl.det)[-1],
                      colnames(ebird.det)[-1],
                      colnames(cbc.det)[-1])
  
  #Define temporal periods using dates and homebrew function above
  ebird.temp <- temp.periods(ebird.data)
  bbs.temp <- temp.periods(bbs.data)
  bbl.temp <- temp.periods(bbl.data)
  cbc.temp <- temp.periods(cbc.data)
  # mws.temp <- temp.periods(mws.data)
  
  #A Matrix for just the spatial mesh 
  A.sb.ebird <- inla.spde.make.A(mesh, loc = as.matrix(ebird.xy[, c("X", "Y")]))
  A.sb.bbs <- inla.spde.make.A(mesh, loc = as.matrix(bbs.xy[, c("X", "Y")]))
  A.sb.bbl <- inla.spde.make.A(mesh, loc = as.matrix(bbl.xy[, c("X", "Y")]))
  A.sb.cbc <- inla.spde.make.A(mesh, loc = as.matrix(cbc.xy[, c("X", "Y")]))
  # A.sb.mws <- inla.spde.make.A(mesh, loc = as.matrix(mws.xy[, c("X", "Y")]))
  
  #Information for integration of PO data
  # n.ebird <- nrow(ebird.data)
  n.bbl <- nrow(bbl.data)
  
  #change data to include 0s for nodes and 1s for presences
  # y.ebird.pp <- rep(0:1, c(nv*4, n.ebird))
  y.bbl.pp <- rep(0:1, c(nv*4, n.bbl))
  
  #add expectation vector (area for integration points/nodes and 0 for presences)
  # e.ebird.pp <- c(w, w, w, w, rep(0, n.ebird))
  e.bbl.pp <- c(w, w, w, w, rep(0, n.bbl))

  
  #Create all relevant A Matrix
  A.ebird <- inla.spde.make.A(mesh, loc = rbind(as.matrix(ebird.xy[, c("X", "Y")])),
                              group = ebird.temp$Season, group.mesh = mesh.t)
  A.bbs <- inla.spde.make.A(mesh, loc = as.matrix(bbs.xy[, c("X", "Y")]), 
                            group = bbs.temp$Season, group.mesh = mesh.t)
  A.bbl <- inla.spde.make.A(mesh, loc = rbind(mesh$loc[,1:2], mesh$loc[,1:2], mesh$loc[,1:2], 
                                              mesh$loc[,1:2], as.matrix(bbl.xy[, c("X", "Y")])),
                            group = c(rep(1:4, each = nrow(mesh$loc)), bbl.temp$Season), 
                            group.mesh = mesh.t)
  A.cbc <- inla.spde.make.A(mesh, loc = as.matrix(cbc.xy[, c("X", "Y")]),
                            group = cbc.temp$Season, group.mesh = mesh.t)
  # A.mws <- inla.spde.make.A(mesh, loc = as.matrix(mws.xy[, c("X", "Y")]),
  #                           group = mws.temp$Season, 
  #                           group.mesh = mesh.t)
  
  # A.ebird.pp <- rbind(imat, imat, imat, imat, A.sb.ebird)
  A.bbl.pp <- rbind(imat, imat, imat, imat, A.sb.bbl)
  
  ### Create Stacks
  ### PRESENCE-ONLY DATA
  #BBL
  bbl.stk <- inla.stack(data = list(y = cbind(y.bbl.pp, NA, NA, NA),
                                      e = e.bbl.pp),
                          effects = list(list(cbind(data.frame(int.BBL = rep(1, nv*4 + n.bbl),
                                                               Season = c(rep(1:4, each = nv), bbl.temp$Season),
                                                               Week2 = c(rep(c(1,3,7,10)*4, each = nv), bbl.temp$Week2),
                                                               Month = c(rep(c(1,3,7,10), each = nv), bbl.temp$Month)),
                                                    rbind(integrat.covs, integrat.covs, integrat.covs,
                                                          integrat.covs, bbl.env),
                                                    NCap = c(rep(1, nv*4), bbl.det[,2]))),
                                         BBL_field = 1:spde$n.spde,
                                         iset),
                          A=list(1, A.bbl.pp, A.bbl),
                          tag="bbl_data")
  
  ### COUNT DATA
  
  ebird.stk <- inla.stack(data = list(y = cbind(NA, ebird.data$Count, NA, NA)),
                        effects=list(list(cbind(data.frame(int.EBird = rep(1, length(ebird.data$Occ))),
                                                ebird.temp[, 3:ncol(ebird.temp)],
                                                ebird.env[,3:ncol(ebird.env)],
                                                ebird.det[,2:ncol(ebird.det)])),
                                     EBird_field = 1:spde$n.spde,
                                     iset),
                        A=list(1, A.sb.ebird, A.ebird),
                        tag="ebird_data")
  
  bbs.stk <- inla.stack(data = list(y = cbind(NA, NA, bbs.data$Count, NA)),
                        effects=list(list(cbind(data.frame(int.BBS = rep(1, length(bbs.data$Occ))),
                                                bbs.temp[, 3:ncol(bbs.temp)],
                                                bbs.env[,3:ncol(bbs.env)],
                                                bbs.det[,2:ncol(bbs.det)])),
                                     BBS_field = 1:spde$n.spde,
                                     iset),
                        A=list(1, A.sb.bbs, A.bbs),
                        tag="bbs_data")
  
  
  cbc.stk <- inla.stack(data = list(y = cbind(NA, NA, NA, cbc.data$Occ)),
                        effects=list(list(cbind(data.frame(int.CBC = rep(1, length(cbc.data$Occ))),
                                                cbc.temp[, 3:ncol(cbc.temp)],
                                                cbc.env[,3:ncol(cbc.env)],
                                                cbc.det[,2:ncol(cbc.det)])),
                                     CBC_field = 1:spde$n.spde,
                                     iset),
                        A=list(1, A.sb.cbc, A.cbc),
                        tag="cbc_data")
  
  all.stk <- inla.stack(bbl.stk, cbc.stk, bbs.stk, ebird.stk)
  
  all.stk$effects$data <- all.stk$effects$data %>%
    mutate_at(c(env.covs.names, det.covs.names), function(x){scale(x)[,1]})
    
  
  allglobal()
}

### A function to call a specific 
isdm_fun <- function(formfixed){
  #Formula
  isdm.formula = as.formula(paste0("y ~ -1 + int.BBL + int.EBird + int.BBS + int.CBC + ",
    formfixed,
    " + f(i, model = spde, group = i.group, control.group = list(model = 'iid')) +
  f(Month,  model = 'ar1', hyper = prior.week2, cyclic = T) +
  f(BBL_field, model = spde) +
  f(EBird_field, copy = 'BBL_field', fixed = TRUE) +
  f(BBS_field, copy = 'BBL_field', fixed = TRUE) +
  f(CBC_field, copy = 'BBL_field', fixed = TRUE)"))
  
  #INLA Call
  result <- inla(isdm.formula, 
                 family = c("poisson", "poisson", "poisson", "poisson"),
                 data = inla.stack.data(all.stk),
                 control.predictor = list(A = inla.stack.A(all.stk), 
                                          compute=TRUE),
                 control.family = list(list(link = "log"), 
                                       list(link = "log"), 
                                       list(link = "log"), 
                                       list(link = "log")),
                 E = inla.stack.data(all.stk)$e,
                 control.compute = list(dic = FALSE, cpo = FALSE,  waic = TRUE),
                 verbose=T, inla.mode = "experimental", 
                 control.inla = list(int.strategy='eb'),
                 control.fixed = list(expand.factor.strategy = 'inla')
  )
  
  return(result)
}

#Prior objects for INLA (spatiotemporal SPDE related; temporal AR1)
h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))
prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))
prior.week2 <- list(theta = list(prior = 'pc.cor1', param = c(0, 0.9)))


### Run Models
#Pick a species
spp_subset("ABDU")

### Check that a null model runs
null.form <- y ~ -1 + int.BBL + int.EBird + int.BBS + int.CBC + 
  f(i, model = spde, group = i.group, control.group = list(model = 'iid')) +
  f(Month,  model = 'ar1', hyper = prior.week2, cyclic = T) +
  f(BBL_field, model = spde) +
  f(EBird_field, copy = 'BBL_field', fixed = TRUE) +
  f(BBS_field, copy = 'BBL_field', fixed = TRUE) +
  f(CBC_field, copy = 'BBL_field', fixed = TRUE)

#1.23 hours for null model with Month as ar1, Spatiotemp SPDE (ar1 for season), BBL Field, Null
#1.18 hours for null

start <- Sys.time()
null.test <- inla(null.form, 
                  family = c("poisson", "poisson", "poisson", "poisson"),
                  data = inla.stack.data(all.stk),
                  control.predictor = list(A = inla.stack.A(all.stk), 
                                           compute=TRUE),
                  control.family = list(list(link = "log"), 
                                        list(link = "log"), 
                                        list(link = "log"), 
                                        list(link = "log")),
                  E = inla.stack.data(all.stk)$e,
                  control.compute = list(dic = FALSE, cpo = FALSE,  waic = TRUE),
                  verbose=T, inla.mode = "experimental", 
                  control.inla = list(int.strategy='eb'),
                  control.fixed = list(expand.factor.strategy = 'inla')
)
Sys.time() -start

save(null.test, file = "C:/Users/mgonnerman/Desktop/Projects/Test.R")

### EXAMINE MODEL OUTPUTS
# All code is from the SPDE INLA book:https://www.google.com/url?q=https%3A%2F%2Fbecarioprecario.bitbucket.io%2Fspde-gitbook%2F&sa=D&sntz=1&usg=AOvVaw20bjmXuRn21WiCqeN0Lzd7
# link to code download here: https://www.google.com/url?q=https%3A%2F%2Furldefense.com%2Fv3%2F__https%3A%2F%2Finla.r-inla-download.org%2Fshared%2Fspde-book%2Fspde-book-files.zip__%3B!!Nmw4Hv0!xD5QKGbpWOwPGQ-_L-GG4_JlRoOQv4tN1ISruG5DW_h0TC-OrstKdv13y_kJIr_qDjc_LDDalbvWVRKw6uy3Z2Rh%24&sa=D&sntz=1&usg=AOvVaw2zF92fXhPU1h-IwubJxWz4
load("C:/Users/mgonnerman/Desktop/Projects/Test.R")

stepsize <- 2 * 1 / 111
nxy <- round(c(max_x-min_x, max_y-min_y) / stepsize)
projgrid <- inla.mesh.projector(mesh, xlim = c(min_x, max_x), ylim = c(min_y, max_y), dims = nxy)
st_field_mean <- list()
for (j in 1:4){
  st_field_mean[[j]] <- inla.mesh.project(
    projgrid, null.test$summary.random$i$mean[iset$i.group == j])
}

#BBL_Field
bbl_field_mean <- inla.mesh.project(projgrid, 
                            null.test$summary.random$BBL_field$mean)

xy.in <- inout(projgrid$lattice$loc, 
               st_coordinates(studyarea)[,1:2])

#Code in script within R folder of downloaded code zip 
image.plot(x = projgrid$x, y = projgrid$y, z = bbl_field_mean, 
           xlab = "Longitude", ylab = "Latitude")


#Plot ST SPDE
par(c(2,2))
for (j in 1:4) {
  st_field_mean[[j]][!xy.in] <- NA
  image.plot(x = projgrid$x, y = projgrid$y, z = st_field_mean[[j]], 
             xlab = "Longitude", ylab = "Latitude", col = viridis(n = 210))
}
# Look at distribution of raw observations from datasets
# ebird.check <- ebird.data %>% select(Survey, Count, Date, X, Y)
# bbs.check <- bbs.data %>% select(Survey, Count, Date, X, Y)
# bbl.check <- bbl.data %>% select(Survey, Count, Date, X, Y)
# cbc.check <- cbc.data %>% select(Survey, Count, Date, X, Y)
# 
# all.check <- rbind(ebird.check, bbs.check, bbl.check, cbc.check) %>%
#   filter(Count > 0) %>%
#   mutate(Date = as.Date(Date)) %>%
#   mutate(Week2 = ceiling(yday(Date)/14), Month = month(Date)) %>%
#   mutate(Season = with(., case_when((Month >= 12 | Month <= 2)~1,
#                                     (Month > 2 & Month <= 5) ~ 2,
#                                     (Month >= 6 & Month <= 8) ~ 3,
#                                     (Month >= 9 | Month <= 11) ~ 4,
#                                     is.na(Month)~NA, TRUE~NA))) %>%
#   st_as_sf(., coords = c("X", "Y"), crs = 4326)
# 
# ggplot() +
#   geom_sf(data = studyarea, fill = NA) +
#   geom_sf(data = all.check, aes(color = Count)) +
#   theme_void() + 
#   facet_wrap(~Season)

##########################################
###########################################
##########################################

#Create Model list
cov.forms <- c(det.covs, paste(env.covs, det.covs, sep = " + "))

# 2 covs - 3.3 hours

start <- Sys.time()
test <- isdm_fun("D2Coast + Wet_10")
Sys.time()-start

save(test, file = "C:/Users/mgonnerman/Desktop/Projects/Test_2covs.R")

### EXAMINE MODEL OUTPUTS
# All code is from the SPDE INLA book:https://www.google.com/url?q=https%3A%2F%2Fbecarioprecario.bitbucket.io%2Fspde-gitbook%2F&sa=D&sntz=1&usg=AOvVaw20bjmXuRn21WiCqeN0Lzd7
# link to code download here: https://www.google.com/url?q=https%3A%2F%2Furldefense.com%2Fv3%2F__https%3A%2F%2Finla.r-inla-download.org%2Fshared%2Fspde-book%2Fspde-book-files.zip__%3B!!Nmw4Hv0!xD5QKGbpWOwPGQ-_L-GG4_JlRoOQv4tN1ISruG5DW_h0TC-OrstKdv13y_kJIr_qDjc_LDDalbvWVRKw6uy3Z2Rh%24&sa=D&sntz=1&usg=AOvVaw2zF92fXhPU1h-IwubJxWz4
load("C:/Users/mgonnerman/Desktop/Projects/Test_2covs.R")

stepsize <- 2 * 1 / 111
nxy <- round(c(max_x-min_x, max_y-min_y) / stepsize)
projgrid <- inla.mesh.projector(mesh, xlim = c(min_x, max_x), ylim = c(min_y, max_y), dims = nxy)
st_field_mean <- list()
for (j in 1:4){
  st_field_mean[[j]] <- inla.mesh.project(
    projgrid, test$summary.random$i$mean[iset$i.group == j])
}

#BBL_Field
bbl_field_mean <- inla.mesh.project(projgrid, 
                                    test$summary.random$BBL_field$mean)

xy.in <- inout(projgrid$lattice$loc, 
               st_coordinates(studyarea)[,1:2])

#Code in script within R folder of downloaded code zip 
image.plot(x = projgrid$x, y = projgrid$y, z = bbl_field_mean, 
           xlab = "Longitude", ylab = "Latitude")


#Plot ST SPDE
par(c(2,2))
for (j in 1:4) {
  st_field_mean[[j]][!xy.in] <- NA
  image.plot(x = projgrid$x, y = projgrid$y, z = st_field_mean[[j]], 
             xlab = "Longitude", ylab = "Latitude", col = viridis(n = 210))
}

test$summary.fixed
