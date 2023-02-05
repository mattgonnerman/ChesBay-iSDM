lapply(c("nimble", "dplyr", "parallel", "coda", "MCMCvis"), require, character.only = T)

load(file = "./iSDM_Nimble_Objects_share.R")

### Package Data for NIMBLE
data <- list(
  #BBS
  y.bbs = as.integer(bbs.counts),
  bbs.Nposid = bbs.nposid,
  bbs.Nseen = bbs.nseen,

  #BBL
  y.bbl = as.integer(bbl.count),
  bbl.covs = bbl.covs,

  #Occupancy/Abundance
  grid.covs = grid.covs
)

constants <- list(
  #BBS
  bbs.nsurvey = bbs.nsurveys,
  grid.bbs = bbs.grid,

  #BBL
  bbl.nsurveys = bbl.nsurveys,
  bbl.ncovs = bbl.ncovs,

  #Occupancy/Abundance
  ncovs.grid = ncovs.grid,
  n.cells = n.cells
)

### Initial Values
N.max <- rep(1, n.cells)
for(i in 1:length(bbs.grid)){
  if(bbs.counts[i] >= N.max[bbs.grid[i]]){
    N.max[bbs.grid[i]] <- bbs.counts[i] + 1
  }
}

for(i in 1:length(bbl.count)){
  if(bbl.count[i] >= N.max[i]){
    N.max[i] <- bbl.count[i] + 1
  }
}

inits <- list(
  N = N.max,
  alpha = rep(0, ncovs.grid),
  alpha0 = logit(.8),
  beta = rep(0, ncovs.grid),
  beta0 = log(mean(N.max)),
  gamma.bbs = rep(.01,2),
  gamma.bbl = rep(.01, bbl.ncovs)
)

source("iSDM - Nimble Model - toshare.R")
source("iSDM - ZIP Function.R")


### Check model before fully running
#Looking for a non-NA value returned from calculate()
model_test <- nimbleModel( code = code,
                           constants = constants,
                           data =  data,
                           inits = inits,
                           calculate = F)
model_test$simulate(c("E.bbl", "p.bbl", "E.bbs", "p.bbs",
                      "psi", "lambda"))
model_test$initializeInfo()
model_test$calculate()


### Assess issues
wrong.bbl <- which(is.infinite(model_test$logProb_y.bbl))
wrong.bbs <- which(is.infinite(model_test$logProb_y.bbs))

bbs.counts[wrong.bbs]

model_test$N[bbs.grid[wrong.bbs]] > model_test$y.bbs[wrong.bbs]
model_test$N[wrong.bbl] > model_test$y.bbl[wrong.bbl]


bbs.grid.wrong <- sort(unique(bbs.grid[wrong.bbs]))

bbs.bbl <- sort(unique(bbs.grid.wrong[which(bbs.grid.wrong %in% wrong.bbl)]))
bbs.only <- sort(unique(bbs.grid.wrong[which(!(bbs.grid.wrong %in% wrong.bbl))]))
bbl.only <- sort(unique(wrong.bbl[which(!(wrong.bbl %in% bbs.grid.wrong))]))

summary(model_test$N[which(1:n.cells %in% c(bbs.only,bbl.only, bbs.bbl))])
summary(model_test$N[which(!(1:n.cells %in% c(bbs.only,bbl.only, bbs.bbl)))])

dZIP(2000, lambda = model_test$lambda[1], zeroProb = model_test$psi[1])

model_test$logProb_alpha0

which(!(model_test$N > model_test$y.bbl))
which(!(model_test$N[bbs.grid] > model_test$y.bbs))

exp(model_test$alpha0)/(1+exp(model_test$alpha0))


model_test$p.bbl[wrong.bbl]
which(model_test$p.bbl != 1) == wrong.bbl

model_test$bbl.covs[wrong.bbl,]


#########################################
### Run the Model
#Set monitors
monitor.coeff <- c("alpha0", "alpha", "beta0", "beta",
                   "gamma.bbs", "gamma.bbl")
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

MCMCsummary(mcmcList2, 'alpha')
MCMCsummary(mcmcList2, 'gamma.bbl')
MCMCsummary(mcmcList2, 'beta')
