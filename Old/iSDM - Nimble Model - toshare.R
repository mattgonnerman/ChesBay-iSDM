code <- nimbleCode({
  ### Grid Cell Occupancy/Abundance Intensity surface
  for(i in 1:n.cells){
    #Zero Inflated Poisson
    cloglog(psi[i]) <- alpha0 + inprod(alpha[1:ncovs.grid], grid.covs[i, 1:ncovs.grid])
    log(lambda[i]) <- beta0 + inprod(beta[1:ncovs.grid], grid.covs[i, 1:ncovs.grid])

    N[i] ~ dZIP(lambda[i], zeroProb = 1-psi[i])
  }
  
  alpha0 ~ dnorm(0, sd = 100)
  beta0 ~ dnorm(0, sd = 100)
  for(i in 1:ncovs.grid){
    alpha[i] ~ dnorm(0, sd = 10)
    beta[i] ~ dnorm(0, sd = 10)
  }

  ### Breeding Bird Survey
  for(j in 1:bbs.nsurvey){
    E.bbs[j] <- gamma.bbs[1]*bbs.Nposid[j] + gamma.bbs[2]*bbs.Nseen[j]
    p.bbs[j] <- 1 - (.5)^E.bbs[j]

    y.bbs[j] ~ dbinom(size = N[grid.bbs[j]], prob = p.bbs[j])
  }
  for(i in 1:2){gamma.bbs[i] ~ T(dnorm(0,0.0001), 0,10000)}
  
  ### Bird Banding Lab
  for(j in 1:bbl.nsurveys){
    E.bbl[j] <- inprod(gamma.bbl[1:bbl.ncovs], bbl.covs[j,1:bbl.ncovs])
    p.bbl[j] <- 1 - (.5)^E.bbl[j]

    y.bbl[j] ~ dbinom(size = N[j], prob = p.bbl[j])
  }
  for(i in 1:bbl.ncovs){gamma.bbl[i] ~ T(dnorm(0,0.0001), 0,10000)}
})