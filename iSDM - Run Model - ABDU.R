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

### Simplify to Mallard only (can adjust to whichever species you want to use)
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


### Package Data for NIMBLE
data <- list(
  #EBird
  y.ebird = ebird.counts,
  ebird.covs = ebird.covs,

  #BBS
  y.bbs = bbs.counts,
  bbs.Nposid = bbs.nposid,
  bbs.Nseen = bbs.nseen,

  # #BBL
  # y.bbl = bbl.count,
  # bbl.covs = bbl.covs,

  #MWS
  y.mws = mws.counts,
  mws.covs = mws.covs,

  #CBC
  y.cbc = cbc.count,
  cbc.covs = cbc.covs,
  
  #Occupancy/Abundance
  grid.covs = grid.covs
)

constants <- list(
  #EBird
  grid.ebird = ebird.grid,
  ebird.nsurvey = ebird.nsurvey,
  ebird.ncovs = ebird.ncovs,

  #BBS
  bbs.nsurvey = bbs.nsurveys,
  grid.bbs = bbs.grid,

  # #BBL
  # bbl.nsurveys = bbl.nsurveys,
  # bbl.ncovs = bbl.ncovs,

  #MWS
  mws.gridid = mws.gridid,
  mws.nsurveys = mws.nsurveys,
  mws.nyears = mws.nyears,
  mws.ncovs = mws.ncovs,

  #CBC
  cbc.nsurveys = cbc.nsurveys,
  cbc.ncovs = cbc.ncovs,
  grid.cbc = cbc.gridid,

  # #dCAR
  # adj = adj,
  # weights = weights,
  # num = num,
  # nadj = nadj,
  
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
  beta = rep(0, ncovs.grid),
  beta0 = log(mean(N.max)),
  alpha = rep(0, ncovs.grid),
  alpha0 = logit(.5),
  gamma.bbs = gammainit(bbs.covmax),
  gamma.ebird = gammainit(ebird.covmax),
  gamma.cbc = gammainit(cbc.covmax),
  gamma.mws = gammainit(mws.covmax),
  N = N.max
)

### Load Model Code
source("iSDM - Nimble Model - Mallard.R")
source("iSDM - ZIP Function.R")


### Check model before fully running
#Looking for a non-NA value returned from calculate()
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits,
                           calculate = F)
model_test$simulate(c("E.ebird", "p.ebird", "E.bbs", "p.bbs",
                      "E.cbc", "p.cbc", "E.mws", "p.mws",
                      "psi", "lambda"))
model_test$initializeInfo()
model_test$calculate()


monitor.coeff <- c("beta0.psi", "beta.psi", "beta0.lambda", #"beta.lambda",
                   "alpha.bbs", "alpha.bbl")
monitor.est <- c("psi", "lambda")


mcmcConf <-  configureMCMC( model_test,   monitors =  monitor.est, monitors2 = monitor.coeff)
mcmc     <-  buildMCMC( mcmcConf)
Cmodel   <- compileNimble(model_test)
Cmcmc    <- compileNimble(mcmc)
samplesList <- runMCMC(Cmcmc,nburnin = 5000, niter = 10000, thin = 5, thin2 = 5)

samples1 <- list(chain1 =  samplesList$samples)
mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
samples2 <- list(chain1 =  samplesList$samples2)
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))

MCMCsummary(mcmcList2, 'alpha.bbs')



# Parallel Processing Setup
rm(out.full.predict)
start_time <- Sys.time() # To track runtime
start_time
nc <- detectCores()/2    # number of chains
cl<-makeCluster(nc,timeout=5184000) #Start 3 parallel processing clusters

clusterExport(cl, c("code", "inits", "data", "constants", "monitor.coeff", "monitor.est")) #identify what is to be exported to each cluster

out.full.predict <- clusterEvalQ(cl, {
  require(nimble)
  require(coda)
  
  model_test <- nimbleModel( code = code,
                             constants = constants,
                             data =  data,
                             inits = inits )
  model_test$simulate(c(
    # "E.ebird", "p.ebird", "alpha.ebird", "sig.alpha.ebird", "mu.alpha.ebird",
    # "E.bbs", "p.bbs", "alpha.bbs", "sig.alpha.bbs", "mu.alpha.bbs",
    "E.bbl", "p.bbl", "alpha.bbl", "sig.alpha.bbl", "mu.alpha.bbl",
    # "E.mws", "p.mws", "alpha.mws", "sig.alpha.mws", "mu.alpha.mws",
    # "E.cbc", "p.cbc", "alpha.cbc", "sig.alpha.cbc", "mu.alpha.cbc"
    "psi", "z", "lambda", "N", "r.nb", "p.nb",
    "gamma.psi", "spacetau.psi", #"gamma.lambda", "spacetau.lambda",
    "beta0.psi", "beta.psi", "mu.beta.psi", "sig.beta.psi",
    "beta0.lambda", "beta.lambda", "mu.beta.lambda", "sig.beta.lambda"
  ))
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