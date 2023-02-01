lapply(c("nimble", "dplyr", "parallel", "coda", "MCMCvis"), require, character.only = T)


### List of loaded data objects
### ebird.counts, ebird.covs, ebird.grid, ebird.nsurvey, ebird.nspecies, ebird.ncovs, #EBird
### bbs.counts, bbs.nspecies, bbs.nsurveys, bbs.nposid, bbs.nseen, bbs.grid, bbs.sp, #BBS
### bbl.occ, bbl.covs, bbl.nsurveys, bbl.nspecies, bbl.ncovs, bbl.sp, #BBL
### mws.counts, mws.covs, mws.gridid, mws.nsurveys, mws.nyears, mws.nspecies, mws.ncovs, mws.sp,#MWS
### cbc.count, cbc.covs, cbc.nsurveys, cbc.nspecies, cbc.ncovs, cbc.gridid, cbc.sp, #CBC
### grid.covs, ncovs.grid, ntotspecies, n.cells, adj, num, weights, nadj
load(file = "iSDM_Nimble_Objects.R")

###Load Bird List
bird.codes.all <- read.csv("./BirdCodes.csv")
birds.use <- c("MALL", "ABDK", "CAGO", "BWTE", "AMWI", "NOPI", "GWTE", "GADW", "NSHO", "WODU", "GSGO", "GWFG")
birds.codes <- bird.codes.all %>% filter(Alpha %in% birds.use)

### Simplify to Mallard only (can adjust to whichever species you want to use)
#EBird
colnum <- which(colnames(ebird.counts) == "MALL")
ebird.counts <- ebird.counts[,colnum]

#BBS
colnum <- which(colnames(bbs.counts) == "MALL")
bbs.counts <- bbs.counts[,colnum]
bbs.nposid <- bbs.nposid[,colnum]

#BBL
colnum <- which(colnames(bbl.occ) == "MALL")
bbl.occ <- bbl.occ[,colnum]

#CBC
colnum <- which(colnames(cbc.count) == "MALL")
cbc.count <- cbc.count[,colnum]

#MWS
mws.counts <- read.csv("./MWS/MWS_ChesBay.csv") %>%
  filter(Include == "Yes", 
         SurveyYear != "AVG") %>%
  select(Name, SurveyYear, MALL) %>% #CHANGE HERE IF DIFF SPECIES
  arrange(Name, SurveyYear) %>%
  tidyr::pivot_wider(names_from = "SurveyYear", 
                     names_prefix = "Year", 
                     values_from = "MALL") %>% #CHANGE THIS TOO
  select(-Name) %>% as.matrix()
mws.counts[is.na(mws.counts)] <- 0


### Package Data for NIMBLE
data <- list(
  # #EBird
  # y.ebird = ebird.counts,
  # ebird.covs = ebird.covs,

  #BBS
  y.bbs = bbs.counts,
  bbs.Nposid = bbs.nposid,
  bbs.Nseen = bbs.nseen,

  #BBL
  y.bbl = bbl.occ,
  bbl.covs = bbl.covs,
# 
#   #MWS
#   y.mws = mws.counts,
#   mws.covs = mws.covs,
# 
#   #CBC
#   y.cbc = cbc.count,
#   cbc.covs = cbc.covs,
  
  #Occupancy/Abundance
  grid.covs = grid.covs
)

constants <- list(
  # #EBird
  # grid.ebird = ebird.grid,
  # ebird.nsurvey = ebird.nsurvey,
  # ebird.ncovs = ebird.ncovs,

  #BBS
  bbs.nsurvey = bbs.nsurveys,
  grid.bbs = bbs.grid,

  #BBL
  bbl.nsurveys = bbl.nsurveys,
  bbl.ncovs = bbl.ncovs,

  # #MWS
  # mws.gridid = mws.gridid,
  # mws.nsurveys = mws.nsurveys,
  # mws.nyears = mws.nyears,
  # mws.ncovs = mws.ncovs,
  # 
  # #CBC
  # cbc.nsurveys = cbc.nsurveys,
  # cbc.ncovs = cbc.ncovs,
  # grid.cbc = cbc.gridid,
  # 
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
z.init <- rep(0, n.cells)
for(i in 1:length(bbs.grid)){ 
  if(z.init[bbs.grid[i]] == 0){
    z.init[bbs.grid[i]] <- ifelse(bbs.counts[i] > 0, 1, 0)
  }
}
for(i in 1:n.cells){ 
  if(z.init[i] == 0){
    z.init[i] <- ifelse(bbl.occ[i] > 0, 1, 0)
  }
}
# 
# for(i in 1:length(ebird.grid)){ 
#   if(z.init[ebird.grid[i]] == 0){
#     z.init[ebird.grid[i]] <- ifelse(ebird.counts[i] > 0, 1, 0)
#   }
# }
# 
# 
# for(i in 1:length(cbc.gridid)){ 
#   if(z.init[cbc.gridid[i]] == 0){
#     z.init[cbc.gridid[i]] <- ifelse(cbc.count[i] > 0, 1, 0)
#   }
# }
# 
# for(i in 1:length(mws.gridid)){ 
#   if(z.init[mws.gridid[i]] == 0){
#     z.init[mws.gridid[i]] <- ifelse(rowSums(mws.counts)[i] > 0, 1, 0)
#   }
# }

inits <- list(
  #Effort
  # alpha.ebird = rep(.01, ebird.ncovs),
  # alpha.bbs = rep(.01, 2),
  # alpha.bbl = rep(.01, bbl.ncovs),
  # alpha.cbc = rep(.01, cbc.ncovs),
  # alpha.mws = rep(.01, mws.ncovs),
  # #dCAR
  # spacetau.psi = 1,
  # gamma.psi = rep(0, n.cells),
  # spacetau.lambda = 1,
  # gamma.lambda = matrix(0, n.cells, ntotspecies),
  #Abundance
  # r.nb = 1,
  # lambda.z = z.init*rep(ceiling(mean(ebird.counts)), n.cells) + 1e-10*(1-z.init),
  # lambda = rep(exp(1), n.cells),
  beta.lambda = rep(0, ncovs.grid),
  beta0.lambda = 1,
  #Occupancy
  # z = z.init,
  # psi = rep(.5, n.cells),
  beta.psi = rep(0, ncovs.grid),
  beta0.psi = log(-log(1-.5))
)

### Load Model Code
source("iSDM - Nimble Model - Mallard.R")

### Check model before fully running
#Looking for a non-NA value returned from calculate()
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           # inits = inits,
                           calculate = F)
model_test$simulate(c(
                      # "E.ebird", "p.ebird", "alpha.ebird",
                      "E.bbs", "p.bbs", "alpha.bbs",
                      "E.bbl", "p.bbl", "alpha.bbl",
                      # "E.mws", "p.mws", "alpha.mws",
                      # "E.cbc", "p.cbc", "alpha.cbc",
                      # "r.nb", "p.nb",
                      # "gamma.psi", "spacetau.psi", "gamma.lambda", "spacetau.lambda",
                      "beta0.psi", "beta.psi", "beta0.lambda", "beta.lambda",
                      "psi", "N", "lambda", "lambda.z","z" 
                      ))
model_test$initializeInfo()
model_test$calculate()
model_test$calculate("lambda")


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