


code <- nimbleCode({
  ### Grid Cell Occupancy/Abundance Intensity surface
  for(i in 1:n.cells){
    #Zero Inflated Poisson
    cloglog(psi[i]) <- beta0.psi #+ inprod(beta.psi[1:5], grid.covs[i, 1:5]) #+ gamma.psi[i,s] 
    log(lambda[i]) <- beta0.lambda # + inprod(beta.lambda[1:5], grid.covs[i, 1:5])) # + gamma.lambda[i,s]
    psi.ZIP[i] <- 1-psi[i]
    N[i] ~ dZIP(lambda[i], zeroProb = psi.ZIP[i])
  }
  
  beta0.psi ~ dnorm(0, sd = 10)
  beta0.lambda ~ dnorm(0, sd = 10)
  for(i in 1:ncovs.grid){
    beta.psi[i] ~ dnorm(0, sd = 10) 
    beta.lambda[i] ~ dnorm(0, sd = 10)
  }
  
  # #Conditional Auto-Regressive
  # #https://r-nimble.org/html_manual/cha-spatial.html
  # #https://ecosystems.psu.edu/research/labs/walter-lab/manual/chapter-9-spatial-epidemiology-in-winbugs/link-to-pdf
  # spacetau.psi ~ dgamma(0.001, 0.001)
  # gamma.psi[1:n.cells] ~ dcar_normal(adj = adj[1:nadj], weights = weights[1:nadj], num = num[1:n.cells], tau = spacetau.psi)
  # 
  # spacetau.lambda ~ dgamma(0.001, 0.001)
  # gamma.lambda[1:n.cells] ~ dcar_normal(adj[1:nadj],weights[1:nadj],num[1:n.cells],spacetau.lambda)
  # 
  
  # ### BBS
  # for(j in 1:bbs.nsurvey){
  #   E.bbs[j] <- alpha.bbs[1]*bbs.Nposid[j] + alpha.bbs[2]*bbs.Nseen[j]
  #   p.bbs[j] <- 1 - (.5)^E.bbs[j]
  # 
  #   y.bbs[j] ~ dbinom(size = N[grid.bbs[j]], prob = p.bbs[j])
  # }
  # for(i in 1:2){alpha.bbs[i] ~ T(dnorm(0,0.0001), 0,10000)}

  
  ### BBL
  for(j in 1:bbl.nsurveys){
    E.bbl[j] <- inprod(alpha.bbl[1:bbl.ncovs], bbl.covs[j,1:bbl.ncovs])
    p.bbl[j] <- 1 - (.5)^E.bbl[j]
    
    y.bbl[j] ~ dbinom(size = N[j], prob = p.bbl[j])
  }
  for(i in 1:bbl.ncovs){alpha.bbl[i] ~ T(dnorm(0,0.0001), 0,10000)}
  # 
  # ### EBird
  # for(j in 1:ebird.nsurvey){
  #   E.ebird[j] <- inprod(alpha.ebird[1:ebird.ncovs], ebird.covs[j,1:ebird.ncovs])
  #   p.ebird[j] <- 1 - pow(.5,E.ebird[j])
  #   
  #   y.ebird[j] ~ dbinom(size = N[grid.ebird[j]], prob = p.ebird[j])
  # }
  # for(i in 1:ebird.ncovs){alpha.ebird[i] ~  T(dnorm(0,0.0001), 0,10000)}
  # 
  
  
  # 
  # ### CBC
  # for(j in 1:cbc.nsurveys){
  #   E.cbc[j] <- inprod(alpha.cbc[1:cbc.ncovs], cbc.covs[j,1:cbc.ncovs])
  #   p.cbc[j] <- 1 - (.5)^E.cbc[j]
  #   
  #   y.cbc[j] ~ dbinom(size = N[grid.cbc[j]], prob = p.cbc[j])
  # }
  # for(i in 1:cbc.ncovs){alpha.cbc[i] ~ T(dnorm(0,0.0001), 0,10000)}
  # 
  # 
  # ### MWS
  # for(t in 1:mws.nyears){
  #   for(j in 1:mws.nsurveys){
  #     E.mws[j,t] <- inprod(alpha.mws[1:mws.ncovs], mws.covs[j,t,1:mws.ncovs])
  #     p.mws[j,t] <- 1 - (.5)^E.mws[j,t]
  #     
  #     y.mws[j,t] ~ dbinom(size = N[mws.gridid[j]], prob = p.mws[j,t])
  #   }
  # }
  # for(i in 1:mws.ncovs){alpha.mws[i] ~ T(dnorm(0,0.0001), 0,10000)}
  
})