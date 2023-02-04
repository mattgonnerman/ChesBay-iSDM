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
  
  # #Conditional Auto-Regressive
  # #https://r-nimble.org/html_manual/cha-spatial.html
  # #https://ecosystems.psu.edu/research/labs/walter-lab/manual/chapter-9-spatial-epidemiology-in-winbugs/link-to-pdf
  # spacetau.psi ~ dgamma(0.001, 0.001)
  # gamma.psi[1:n.cells] ~ dcar_normal(adj = adj[1:nadj], weights = weights[1:nadj], num = num[1:n.cells], tau = spacetau.psi)
  # 
  # spacetau.lambda ~ dgamma(0.001, 0.001)
  # gamma.lambda[1:n.cells] ~ dcar_normal(adj[1:nadj],weights[1:nadj],num[1:n.cells],spacetau.lambda)

  
  ### BBS
  for(j in 1:bbs.nsurvey){
    E.bbs[j] <- gamma.bbs[1]*bbs.Nposid[j] + gamma.bbs[2]*bbs.Nseen[j]
    p.bbs[j] <- 1 - (.5)^E.bbs[j]
    
    y.bbs[j] ~ dbinom(size = N[grid.bbs[j]], prob = p.bbs[j])
  }
  for(i in 1:2){gamma.bbs[i] ~ T(dnorm(0,0.0001), 0,10000)}

  #### EBird
  for(j in 1:ebird.nsurvey){
    E.ebird[j] <- inprod(gamma.ebird[1:ebird.ncovs], ebird.covs[j,1:ebird.ncovs])
    p.ebird[j] <- 1 - (.5)^E.ebird[j]

    y.ebird[j] ~ dbinom(size = N[grid.ebird[j]], prob = p.ebird[j])
  }
  for(i in 1:ebird.ncovs){gamma.ebird[i] ~  T(dnorm(0,0.0001), 0,10000)}

  ### MWS
  for(t in 1:mws.nyears){
    for(j in 1:mws.nsurveys){
      E.mws[j,t] <- inprod(gamma.mws[1:mws.ncovs], mws.covs[j,t,1:mws.ncovs])
      p.mws[j,t] <- 1 - (.5)^E.mws[j,t]

      y.mws[j,t] ~ dbinom(size = N[mws.gridid[j]], prob = p.mws[j,t])
    }
  }
  for(i in 1:mws.ncovs){gamma.mws[i] ~ T(dnorm(0,0.0001), 0,10000)}
  
  ### CBC
  for(j in 1:cbc.nsurveys){
    E.cbc[j] <- inprod(gamma.cbc[1:cbc.ncovs], cbc.covs[j,1:cbc.ncovs])
    p.cbc[j] <- 1 - (.5)^E.cbc[j]

    y.cbc[j] ~ dbinom(size = N[grid.cbc[j]], prob = p.cbc[j])
  }
  for(i in 1:cbc.ncovs){gamma.cbc[i] ~ T(dnorm(0,0.0001), 0,10000)}


  # ### BBL
  # for(j in 1:bbl.nsurveys){
  #   E.bbl[j] <- inprod(alpha.bbl[1:bbl.ncovs], bbl.covs[j,1:bbl.ncovs])
  #   p.bbl[j] <- 1 - (.5)^E.bbl[j]
  #   
  #   y.bbl[j] ~ dbinom(size = N[j], prob = p.bbl[j])
  # }
  # for(i in 1:bbl.ncovs){alpha.bbl[i] ~ T(dnorm(0,0.0001), 0,10000)}
})