set.seed(1234)

N <- 3680 # Number of grid cells 
alpha0 <- exp(.2)/(1+exp(.2)) #Occ intercept
alpha1 <- .5 #Occ beta coefficient
beta0 <- log(5) #Abun intercept
beta1 <- 1 #Abun beta coefficient
x <- rnorm(n=N) # Occ/Abun covariate values

#Real data has multiple sites with orders of magnitude greater abundance
x[sample(1:1000, 70)] <- runif(30, 6,8)


mu.binary <- alpha0*1 + alpha1*x #Psi calculation
pi <- exp(mu.binary)/(1+exp(mu.binary)) #link function
zero <- rbinom(n=N, size=1, prob=pi) #True occ state for each cell
mu.count <- beta0*1 + beta1*x #Lambda regression
lambda <- exp(mu.count) #link function
y <- rpois(n=N, lambda=lambda*(zero)) #Zero-inflated Abundance
  
surv.id1 <- sort(sample(1:1000, 500, replace = F))
z1 <- runif(n=500, 0, 2) #detection covariate
gamma1 <- .5 #Detection coefficient
effort1 <- gamma1*z1 #Detection calculation
p1 <- 1- (.5)^effort1
obs.y1 <- rbinom(n = 500, size = y[surv.id1], prob = p1)

z2 <- runif(n=N, 0, 4) #detection covariate
gamma2 <- .2 #Detection coefficient
effort2 <- gamma2*z2 #Detection calculation
p2 <- 1- (.5)^effort2
obs.y2 <- rbinom(n = N, size = y, prob = p2)


#https://rushinglab.github.io/WILD6900/articles/n-mixture2.html
#https://georgederpa.github.io/teaching/countModels.html
#https://r-nimble.org/nimbleExamples/zero_inflated_poisson.html

#Load zero inflated poisson distribution function
source("iSDM - ZIP Function.R")

ZIPcode <- nimbleCode({  
  for (i in 1:1000){
    cloglog(psi[i]) <- alpha0 + alpha1*cov1[i]
    log(lambda[i]) <- beta0 + beta1*cov1[i]
    
    N[i] ~ dZIP(lambda[i], zeroProb = 1-psi[i])
  }
  alpha0 ~ dnorm(0, sd = 10)
  beta0 ~ dnorm(0, sd = 10)
  alpha1 ~ dnorm(0, sd = 10) 
  beta1 ~ dnorm(0, sd = 10)  
  
  ###Survey 1
  for(i in 1:500){
    E.1[i] <- gamma1*cov2[i]
    p.1[i] <- 1 - (.5)^E.1[i]
    
    y.1[i] ~ dbinom(size = N[surv.id1[i]], prob = p.1[i])
  }
  gamma1 ~ T(dnorm(0,0.0001), 0,1000)  
  
  ###Survey2
  for(i in 1:1000){
    E.2[i] <- gamma2*cov3[i]
    p.2[i] <- 1 - (.5)^E.2[i]
    
    y.2[i] ~ dbinom(size = N[i], prob = p.2[i])  
  }
  gamma2 ~ T(dnorm(0,0.0001), 0,10000)
}) 

data2 <- list(
  y.1 = obs.y1,
  y.2 = obs.y2,
  cov1 = x,
  cov2 = z1,
  cov3 = z2
)

constants2 <- list(
  surv.id1 = surv.id1
)

N.max <- rep(1, 1000)
for(i in 1:500){
  if(obs.y1[i] > N.max[surv.id1[i]]){
    N.max[surv.id1[i]] <- obs.y1[i] + 1
  }
}

for(i in 1:1000){
  if(obs.y2[i] > N.max[i]){
    N.max[i] <- obs.y2[i] + 1
  }
}

inits2 <- list(
  N = N.max,
  alpha0 = logit(.5),
  alpha1 = 0,
  beta0 = log(4),
  beta1 = 0,
  gamma1 = 1,
  gamma2 = 1
)

model_test <- nimbleModel( code = ZIPcode,
                           data =  data2,
                           constants = constants2,
                           inits = inits2)
# model_test$simulate()
model_test$initializeInfo()
model_test$calculate()

monitor.coeff <- c("alpha0", "alpha1", "beta0", "beta1",
                   "gamma1", "gamma2")
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

#Output Estimates
MCMCsummary(mcmcList1, 'psi')
MCMCsummary(mcmcList1, 'lambda')

MCMCsummary(mcmcList2, 'alpha0')
MCMCsummary(mcmcList2, 'alpha1')
MCMCsummary(mcmcList2, 'beta0')
MCMCsummary(mcmcList2, 'beta1')
MCMCsummary(mcmcList2, 'gamma1')
MCMCsummary(mcmcList2, 'gamma2')