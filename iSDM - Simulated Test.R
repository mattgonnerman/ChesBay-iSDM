set.seed(1234)

N <- 1000 # Number of grid cells 
alpha0 <- -.5 #Occ intercept
alpha1 <- -.5 #Occ beta coefficient
beta0 <- 1 #Abun intercept
beta1 <- 1 #Abun beta coefficient
x <- rnorm(n=N) # Occ/Abun covariate values
mu.binary <- alpha0*1 + alpha1*x #Psi calculation
pi <- exp(mu.binary)/(1+exp(mu.binary)) #link function
zero <- rbinom(n=N, size=1, prob=pi) #True occ state for each cell
mu.count <- beta0*1 + beta1*x #Lambda regression
lambda <- exp(mu.count) #link function
N.expected <- rpois(n=N, lambda=lambda) #expected Abundance
y <- N.expected*(zero) #Zero-inflated Abundance
  
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

dat <- cbind(x,z1,z2,obs.y1, obs.y2)

#https://rushinglab.github.io/WILD6900/articles/n-mixture2.html
#https://r-nimble.org/nimbleExamples/zero_inflated_poisson.html

#######################################################################
### Attempt to use nimbleFunction to see if this helps
dZIP <- nimbleFunction(
  run = function(x = integer(), lambda = double(), 
                 zeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(dpois(x, lambda, log = TRUE) + log(1 - zeroProb))
      ## or the probability if log = FALSE
      else return((1 - zeroProb) * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1 - zeroProb) * dpois(0, lambda, log = FALSE)
    if (log) return(log(totalProbZero))
    return(totalProbZero)
  })

rZIP <- nimbleFunction(
  run = function(n = integer(), lambda = double(), zeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if (isStructuralZero) return(0)
    return(rpois(1, lambda))
  })

registerDistributions(list(
  dZIP = list(
    BUGSdist = "dZIP(lambda, zeroProb)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = integer()', 'lambda = double()', 'zeroProb = double()')
  )))

ZIPcode <- nimbleCode({  
  for (i in 1:1000){
    cloglog(psi[i]) <- alpha0 + alpha1*cov1[i]
    lambda[i] <- exp(beta0 + beta1*cov1[i])
    
    N[i] ~ dZIP(lambda[i], zeroProb = 1-psi[i])
    
    E.2[i] <- gamma2*cov3[i]
    p.2[i] <- 1 - (.5)^E.2[i]
    
    y.2[i] ~ dbinom(size = N[i], prob = p.2[i])  
    
  }
  
  for(i in 1:500){
    E.1[i] <- gamma1*cov2[i]
    p.1[i] <- 1 - (.5)^E.1[i]
    
    y.1[i] ~ dbinom(size = N[surv.id1[i]], prob = p.1[i])
    
  }
  
  alpha0 ~ dnorm(0, sd = 10)
  beta0 ~ dnorm(0, sd = 10)
  alpha1 ~ dnorm(0, sd = 10) 
  beta1 ~ dnorm(0, sd = 10)
  gamma1 ~ T(dnorm(0,0.0001), 0,1000)
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
