lapply(c("nimble", "dplyr", "parallel", "coda", "MCMCvis"), require, character.only = T)

load(file = "iSDM_Nimble_Objects.R")


#########################################################################
### SUBSET FOR TESTS
### List of loaded data objects
### ebird.counts, ebird.covs, ebird.grid, ebird.nsurvey, ebird.nspecies, ebird.ncovs, #EBird
### bbs.counts, bbs.nspecies, bbs.nsurveys, bbs.nposid, bbs.nseen, bbs.grid, bbs.sp, #BBS
### bbl.occ, bbl.covs, bbl.nsurveys, bbl.nspecies, bbl.ncovs, bbl.sp, #BBL
### mws.counts, mws.covs, mws.gridid, mws.nsurveys, mws.nyears, mws.nspecies, mws.ncovs, mws.sp,#MWS
### cbc.count, cbc.covs, cbc.nsurveys, cbc.nspecies, cbc.ncovs, cbc.gridid, cbc.sp, #CBC
### grid.covs, ncovs.grid, ntotspecies, n.cells, adj, num, weights, nadj

set.seed(1234)

#EBird
ebird.nsurvey.og <- ebird.nsurvey
ebird.nsurvey <- round(ebird.nsurvey*.01)
ebird.subset.index <- sort(sample(1:ebird.nsurvey.og, ebird.nsurvey))
ebird.counts <- ebird.counts[ebird.subset.index,]
ebird.covs <- ebird.covs[ebird.subset.index,]
ebird.grid <- ebird.grid[ebird.subset.index]

#BBS
bbs.nsurvey.og <- bbs.nsurveys
bbs.nsurveys <- round(bbs.nsurveys*.25)
bbs.subset.index <- sort(sample(1:bbs.nsurvey.og, bbs.nsurveys))
bbs.counts <- bbs.counts[bbs.subset.index,]
bbs.nposid <- bbs.nposid[bbs.subset.index,]
bbs.nseen <- bbs.nseen[bbs.subset.index]
bbs.grid <- bbs.grid[bbs.subset.index]

# #BBL
# bbl.nsurvey.og <- bbl.nsurveys
# bbl.nsurveys <- round(bbl.nsurveys*.05)
# bbl.subset.index <- sort(sample(1:bbl.nsurvey.og, bbl.nsurveys))
# bbl.occ <- bbl.occ[bbl.subset.index,]
# bbl.covs <- bbl.covs[bbl.subset.index,]

# #MWS
# mws.nsurvey.og <- mws.nsurveys
# mws.nsurveys <- round(mws.nsurveys*.05)
# mws.subset.index <- sort(sample(1:mws.nsurvey.og, mws.nsurveys))
# mws.counts <- mws.counts[mws.subset.index,,]
# mws.covs <- mws.covs[mws.subset.index,,]
# mws.grid <- mws.grid[mws.subset.index]

#CBC
cbc.nsurvey.og <- cbc.nsurveys
cbc.nsurveys <- round(cbc.nsurveys*.15)
cbc.subset.index <- sort(sample(1:cbc.nsurvey.og, cbc.nsurveys))
cbc.count <- cbc.count[cbc.subset.index,]
cbc.covs <- cbc.covs[cbc.subset.index,]
cbc.gridid <- cbc.gridid[cbc.subset.index]

set.seed(NULL)
#########################################################################


data <- list(
  # #EBird
  # y.ebird = ebird.counts,
  # ebird.covs = ebird.covs,
  # 
  # #BBS
  # y.bbs = bbs.counts,
  # bbs.Nposid = bbs.nposid,
  # bbs.Nseen = bbs.nseen,

  #BBL
  y.bbl = bbl.occ,
  bbl.covs = bbl.covs,

  # #MWS
  # y.mws = mws.counts,
  # mws.covs = mws.covs,
  # 
  # #CBC
  # y.cbc = cbc.count,
  # cbc.covs = cbc.covs 
  
  #Occupancy/Abundance
  grid.covs = grid.covs
  

)

constants <- list(
  # #EBird
  # grid.ebird = ebird.grid,
  # ebird.nsurvey = ebird.nsurvey,
  # ebird.nspecies = ebird.nspecies,
  # ebird.ncovs = ebird.ncovs,
  
  # #BBS
  # bbs.nspecies = bbs.nspecies,
  # bbs.nsurvey = bbs.nsurveys,
  # grid.bbs = bbs.grid,
  # bbs.sp = bbs.sp,

  #BBL
  bbl.nspecies = bbl.nspecies,
  bbl.nsurveys = bbl.nsurveys,
  bbl.sp = bbl.sp,
  bbl.ncovs = bbl.ncovs,

  # #MWS
  # mws.gridid = mws.gridid,
  # mws.nsurveys = mws.nsurveys,
  # mws.nyears = mws.nyears,
  # mws.nspecies = mws.nspecies,
  # mws.ncovs = mws.ncovs,
  # mws.sp = mws.sp,
  # 
  # #CBC
  # cbc.nsurveys = cbc.nsurveys,
  # cbc.nspecies = cbc.nspecies,
  # cbc.ncovs = cbc.ncovs,
  # grid.cbc = cbc.gridid,
  # cbc.sp = cbc.spntotspecies = ntotspecies,
  
  #Occupancy/Abundance
  ncovs.psi = ncovs.grid,
  ncovs.lambda = ncovs.grid,
  n.cells = n.cells,
  adj = adj, 
  weights = weights,
  num = num,
  nadj = nadj
)

### Initial Values
alpha.ebird <- matrix(0.1, ebird.ncovs, ebird.nspecies)
p.ebird <- E.ebird <- matrix(NA, ebird.nsurvey, ebird.nspecies)
for(s in 1:ebird.nspecies){
  for(j in 1:ebird.nsurvey){
    E.ebird[j,s] <- inprod(alpha.ebird[,s], ebird.covs[j,])
    p.ebird[j,s] <- 1 - pow(.5, E.ebird[j,s])
  }
}

inits <- list(
  # mu.alpha.ebird = rep(0.1, ebird.ncovs),
  # sig.alpha.ebird = rep(1, ebird.ncovs),
  # alpha.ebird = matrix(0.1, ebird.ncovs, ebird.nspecies),  beta0.psi = rep(0, ntotspecies),
  
  #Occupancy Process
  mu.beta.psi = rep(0, ncovs.grid),
  sig.beta.psi = rep(1, ncovs.grid),
  beta.psi = matrix(0, ntotspecies, ncovs.grid),
  # spacesigma.psi = rep(1, ntotspecies),
  # spacetau.psi = rep(1, ntotspecies),
  # gamma.psi = matrix(0, n.cells, ntotspecies),
  z = matrix(1, n.cells, ntotspecies),
  
  #Abundance Process
  beta0.lambda = log(colMeans(ebird.counts)),
  # spacesigma.lambda = rep(1, ntotspecies),
  # gamma.lambda = matrix(0, n.cells, ntotspecies),
  r.nb = rep(1,ntotspecies),
  lambda = matrix(rep(colMeans(ebird.counts), each = n.cells), ncol = ntotspecies),
  mu.beta.lambda = rep(0, ncovs.grid),
  sig.beta.lambda = rep(1, ncovs.grid),
  beta.lambda = matrix(0, ntotspecies, ncovs.grid)

)

### Load Model Code
source("iSDM - Nimble Model.R")

### Check model before fully running
#Looking for a non-NA value returned from calculate()
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits,
                           calculate = F)
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

monitor.coeff <- c("beta0.psi", "beta.psi", "mu.beta.psi", "sig.beta.psi", "beta0.lambda")
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

MCMCsummary(mcmcList2, 'beta0.lambda')



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