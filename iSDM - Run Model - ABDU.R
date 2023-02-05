lapply(c("nimble", "dplyr", "parallel", "coda", "MCMCvis", "matrixStats"), require, character.only = T)


### List of loaded data objects
### ebird.counts, ebird.covs, ebird.grid, ebird.nsurvey, ebird.nspecies, ebird.ncovs, #EBird
### bbs.counts, bbs.nspecies, bbs.nsurveys, bbs.nposid, bbs.nseen, bbs.grid, bbs.sp, #BBS
### bbl.count, bbl.covs, bbl.nsurveys, bbl.nspecies, bbl.ncovs, bbl.sp, #BBL
### mws.counts, mws.covs, mws.gridid, mws.nsurveys, mws.nyears, mws.nspecies, mws.ncovs, mws.sp,#MWS
### cbc.count, cbc.covs, cbc.nsurveys, cbc.nspecies, cbc.ncovs, cbc.gridid, cbc.sp, #CBC
### grid.covs, ncovs.grid, ntotspecies, n.cells, adj, num, weights, nadj
load(file = "iSDM_Nimble_Objects.R")

###Load Bird List
birds.use <- birds.use <- sort(c("MALL", "ABDK", "ABDU", "CAGO", "BWTE", "AMWI", "NOPI", 
                                 "GWTE", "AGWT","GADW", "NSHO", "WODU", "GSGO", "GWFG",
                                 "BLSC", "SUSC", "LTDU", "LESC"))

### Simplify to SPECIES OF INTEREST only (can adjust to whichever species you want to use)
#EBird
colnum <- which(colnames(ebird.counts) == "ABDU")
ebird.counts <- ebird.counts[,colnum]

#BBS
colnum <- which(colnames(bbs.counts) == "ABDU")
bbs.counts <- bbs.counts[,colnum]
bbs.nposid <- bbs.nposid[,colnum]

# #BBL
# colnum <- which(colnames(bbl.count) == "ABDU")
# bbl.count <- bbl.count[,colnum]

#CBC
colnum <- which(colnames(cbc.count) == "ABDU")
cbc.count <- cbc.count[,colnum]

#MWS
mws.counts <- read.csv("./MWS/MWS_ChesBay.csv") %>%
  filter(Include == "Yes", 
         SurveyYear != "AVG") %>%
  select(Name, SurveyYear, ABDU) %>% #CHANGE HERE IF DIFF SPECIES
  arrange(Name, SurveyYear) %>%
  tidyr::pivot_wider(names_from = "SurveyYear", 
                     names_prefix = "Year", 
                     values_from = "ABDU") %>% #CHANGE THIS TOO
  select(-Name) %>% as.matrix()
mws.counts[is.na(mws.counts)] <- 0


# ########################################################
# ### SUBSET EBIRD FOR TESTS
# ### CAN REMOVE FOR FINAL RUN
# ebl <- length(ebird.counts)
# Per <- .10 #percent to cut off of ebird
# eb.sub <- sort(sample(1:ebl, floor(ebl*Per)))
# eb.sub1 <- max(eb.sub)
# ebird.counts <- ebird.counts[eb.sub]
# ebird.covs <- ebird.covs[eb.sub,]
# ebird.grid <- ebird.grid[eb.sub]
# ebird.nsurvey <- length(eb.sub)
# ########################################################

### Package Data for NIMBLE
data <- list(

  #BBS
  y.bbs = bbs.counts,
  bbs.Nposid = bbs.nposid,
  bbs.Nseen = bbs.nseen,
  
  #EBird
  y.ebird = ebird.counts,
  ebird.covs = ebird.covs,

  #MWS
  y.mws = mws.counts,
  mws.covs = mws.covs,

  #CBC
  y.cbc = cbc.count,
  cbc.covs = cbc.covs,

  # #BBL
  # y.bbl = bbl.count,
  # bbl.covs = bbl.covs,
  
  #Occupancy/Abundance
  grid.covs = grid.covs
)

constants <- list(
  #BBS
  bbs.nsurvey = bbs.nsurveys,
  grid.bbs = bbs.grid,

  #EBird
  grid.ebird = ebird.grid,
  ebird.nsurvey = ebird.nsurvey,
  ebird.ncovs = ebird.ncovs,

  #MWS
  mws.gridid = mws.gridid,
  mws.nsurveys = mws.nsurveys,
  mws.nyears = mws.nyears,
  mws.ncovs = mws.ncovs,

  #CBC
  cbc.nsurveys = cbc.nsurveys,
  cbc.ncovs = cbc.ncovs,
  grid.cbc = cbc.gridid,

  # #BBL
  # bbl.nsurveys = bbl.nsurveys,
  # bbl.ncovs = bbl.ncovs,

  #dCAR
  adj = adj,
  weights = weights,
  num = num,
  nadj = nadj,
  
  #Occupancy/Abundance
  ncovs.grid = ncovs.grid,
  n.cells = n.cells
)

### Initial Values
### Initial Values
N.max <- rep(1, n.cells)

for(i in 1:length(ebird.grid)){
  if(ebird.counts[i] >= N.max[ebird.grid[i]]){
    N.max[ebird.grid[i]] <- ebird.counts[i] + 1
  }
}

for(i in 1:length(bbs.grid)){
  if(bbs.counts[i] >= N.max[bbs.grid[i]]){
    N.max[bbs.grid[i]] <- bbs.counts[i] + 1
  }
}

for(i in 1:length(cbc.count)){
  if(cbc.count[i] >= N.max[cbc.gridid[i]]){
    N.max[cbc.gridid[i]] <- cbc.count[i] + 1
  }
}

for(i in 1:nrow(mws.counts)){
  if(rowMaxs(mws.counts)[i] >= N.max[mws.gridid[i]]){
    N.max[mws.gridid[i]] <- rowMaxs(mws.counts)[i] + 1
  }
}


gammainit <- function(covX){
  log(.98)/(covX*log(.5))
}

bbs.covmax <- c(max(bbs.nposid), max(bbs.nseen))
ebird.covmax <- colMaxs(ebird.covs)
cbc.covmax <- colMaxs(cbc.covs)
mws.covmax <- c(max(mws.covs[,,1]),max(mws.covs[,,2]))

inits <- list(
  beta0 = log(mean(N.max)),
  beta = rep(0, ncovs.grid),
  beta.car = rep(0, n.cells),
  tau.beta = 1,
  alpha0 = logit(.5),
  alpha = rep(0, ncovs.grid),
  alpha.car = rep(0, n.cells),
  tau.alpha = 1,
  gamma.bbs = gammainit(bbs.covmax),
  gamma.ebird = gammainit(ebird.covmax),
  gamma.cbc = gammainit(cbc.covmax),
  gamma.mws = gammainit(mws.covmax),
  N = N.max
)

### Load Model Code
source("iSDM - Model Code - ABDU.R")
source("iSDM - ZIP Function.R")

## DONT NEED TO RUN ONCE ALL VARIABLES ARE INITIALIZED
## AND CALCULATE RETURNS NON -INF/NA VALUE
## Check model before fully running
## Looking for a non-NA/-Inf value returned from calculate()
# model_test <- nimbleModel( code = code,
#                            constants = constants,
#                            data =  data,
#                            inits = inits,
#                            calculate = F)
# model_test$simulate(c("E.ebird", "p.ebird", "E.bbs", "p.bbs",
#                       "E.cbc", "p.cbc", "E.mws", "p.mws",
#                       "psi", "lambda", "lifted_d1_minus_psi_oBi_cB_L4",
#                       "lifted_CAR_calcNumIslands_oPadj_oB1to21578_cB_comma_num_oB1to3680_cB_cP")) 
# ### When the above throws an error, you probably changed the grid size and need to rename this
# 
# 
# model_test$initializeInfo()
# model_test$calculate()



# ### Assess issues
# wrong.bbs <- which(is.infinite(model_test$logProb_y.bbs))
# wrong.ebird <- which(is.infinite(model_test$logProb_y.ebird))
# wrong.cbc <- which(is.infinite(model_test$logProb_y.cbc))
# wrong.mws <- which(is.infinite(model_test$logProb_y.mws))
# 
# model_test$y.ebird[wrong.ebird]
# model_test$p.ebird[wrong.ebird]
# model_test$ebird.covs[wrong.ebird,]
# ###


#Set the paramaters to monitor
monitor.coeff <- c("alpha0", "alpha", "alpha.car",
                   "beta0", "beta", "beta.car",
                   "gamma.ebird", "gamma.bbs", "gamma.cbc", "gamma.mws")
monitor.est <- c("psi", "lambda")

# Parallel Processing Setup
rm(out.full.predict)
start_time <- Sys.time() # To track runtime
start_time
nc <- detectCores()/2    # number of chains
cl<-makeCluster(nc,timeout=5184000) #Start 3 parallel processing clusters

clusterExport(cl, c("code", "inits", "data", "constants", "monitor.coeff", "monitor.est", "rZIP", "dZIP")) #identify what is to be exported to each cluster

out.full.predict <- clusterEvalQ(cl, {
  require(nimble)
  require(coda)
  
  registerDistributions(list(
    dZIP = list(
      BUGSdist = "dZIP(lambda, zeroProb)",
      discrete = TRUE,
      range = c(0, Inf),
      types = c('value = integer()', 'lambda = double()', 'zeroProb = double()')
    )))
  
  model_test <- nimbleModel( code = code,
                             constants = constants,
                             data =  data,
                             inits = inits )
  model_test$simulate(c("E.ebird", "p.ebird", "E.bbs", "p.bbs",
                        "E.cbc", "p.cbc", "E.mws", "p.mws",
                        "psi", "lambda", "lifted_d1_minus_psi_oBi_cB_L4"))
  model_test$initializeInfo()
  model_test$calculate()
  
  mcmcConf <-  configureMCMC( model_test,   monitors =  monitor.est, monitors2 = monitor.coeff)
  mcmc     <-  buildMCMC( mcmcConf)
  Cmodel   <- compileNimble(model_test)
  Cmcmc    <- compileNimble(mcmc)
  samplesList <- runMCMC(Cmcmc,nburnin = 5000, niter = 10000, thin = 5, thin2 = 5)
  return(samplesList)
})
#Stop parallel cluster
stopCluster(cl)

#Find model runtime
end_time <- Sys.time()
end_time - start_time

samples.param <- list(chain1 =  out.full.predict[[1]]$samples,
                 chain2 =  out.full.predict[[2]]$samples,
                 chain3 =  out.full.predict[[3]]$samples)

mcmcList.param <- as.mcmc.list(lapply(samples.param, mcmc))

samples.coef <- list(chain1 =  out.full.predict[[1]]$samples2,
                 chain2 =  out.full.predict[[2]]$samples2,
                 chain3 =  out.full.predict[[3]]$samples2)

mcmcList.coef <- as.mcmc.list(lapply(samples.coef, mcmc))

#Save Outputs as file
abdu.files <- list(mcmcList.param, mcmcList.coef, code, data)
save(abdu.files, file = "./Output/iSDM_ABDU.rdata")


MCMCtrace(mcmcList.param, filename = "./Output/TraceOut - Estimates.pdf")
MCMCtrace(mcmcList.coef, filename = "./Output/TraceOut - Covariates.pdf")
