
code <- nimbleCode({
  ### Grid Cell Occupancy/Abundance Intensity surface
  for(s in 1:ntotspecies){
    for(i in 1:n.cells){
      #Occupancy/Zero-inflation term
      cloglog(psi[i,s]) <- beta0.psi[s] + inprod(beta.psi[s, 1:5], grid.covs[i, 1:5]) #+ gamma.psi[i,s] 
      z[i,s] ~ dbern(psi[i,s])

      #Zero Inflated Negative Binomial
      log(lambda[i,s]) <- beta0.lambda[s] #+ inprod(beta.lambda[s, 1:5], grid.covs[i, 1:5]) # + gamma.lambda[i,s]
      p.nb[i,s] <- r.nb[s]/(r.nb[s] + (z[i,s]*lambda[i,s])) - 1e-10*(1-z[i,s])
      N[i,s] ~ dnegbin(prob = p.nb[i,s], size = r.nb[s])
    }
    
    beta0.psi[s] ~ dnorm(0, sd = 100)
    beta0.lambda[s] ~ dnorm(0, sd = 100)
    r.nb[s] ~ dunif(0,50)
  }
  
  #Occupancy Coefficients
  for(i in 1:ncovs.psi){
    mu.beta.psi[i] ~ dnorm(0, sd = 100)
    sig.beta.psi[i] ~ T(dt(0, pow(2.5, -2), 1), 0, )
    for(s in 1:ntotspecies){
      beta.psi[s, i] ~ dnorm(mu.beta.psi[i], sd = sig.beta.psi[i])
    }
  }
  #Abundance Coefficients
  for(i in 1:ncovs.lambda){
    mu.beta.lambda[i] ~ dnorm(0, sd = 100)
    sig.beta.lambda[i] ~ T(dt(0, pow(2.5, -2), 1), 0, )
    for(s in 1:ntotspecies){
      beta.lambda[s, i] ~ dnorm(mu.beta.lambda[i], sd = sig.beta.lambda[i])
    }
  }
  
  #Conditional Auto-Regressive
  #https://r-nimble.org/html_manual/cha-spatial.html
  #https://ecosystems.psu.edu/research/labs/walter-lab/manual/chapter-9-spatial-epidemiology-in-winbugs/link-to-pdf
  # for(s in 1:ntotspecies){
  #   spacetau.psi[s] ~ dgamma(0.001, 0.001)
  #   gamma.psi[1:n.cells, s] ~ dcar_normal(adj = adj[1:nadj], weights = weights[1:nadj], num = num[1:n.cells], tau = spacetau.psi[s])
  # 
  #   # spacesigma.lambda[s] ~ dunif(0,5)
  #   # spacetau.lambda[s] <- 1/(spacesigma.lambda[s]*spacesigma.lambda[s])
  #   # gamma.lambda[1:n.cells, s] ~ dcar_normal(adj[1:nadj],weights[1:nadj],num[1:n.cells],spacetau.lambda[s])
  # }
  
  ### EBird
  # for(s in 1:ebird.nspecies){
  #   for(j in 1:ebird.nsurvey){
  #     E.ebird[j,s] <- inprod(alpha.ebird[1:ebird.ncovs,s], ebird.covs[j,1:ebird.ncovs])
  #     p.ebird[j,s] <- 1 - pow(.5,E.ebird[j,s])
  #     
  #     y.ebird[j,s] ~ dbinom(N[grid.ebird[j],s], z[grid.ebird[j],s]*p.ebird[j,s])
  #   }
  # }
  
  
  # ### BBS
  # for(s in 1:bbs.nspecies){
  #   for(j in 1:bbs.nsurvey){
  #     E.bbs[j,s] <- alpha.bbs[1,s]*bbs.Nposid[j,s] + alpha.bbs[2,s]*bbs.Nseen[j]
  #     p.bbs[j,s] <- 1 - (.5)^E.bbs[j,s]
  # 
  #     y.bbs[j,s] ~ dbinom(N[grid.bbs[j], bbs.sp[s]], z[grid.bbs[j],bbs.sp[s]]*p.bbs[j,s])
  #   }
  # }


  ### BBL
  for(s in 1:bbl.nspecies){
    for(j in 1:bbl.nsurveys){
      E.bbl[j,s] <- inprod(alpha.bbl[1:bbl.ncovs,s], bbl.covs[j,1:bbl.ncovs])
      p.bbl[j,s] <- 1 - (.5)^E.bbl[j,s]

      y.bbl[j,s] ~ dbinom(p.bbl[j,s], N[j,bbl.sp[s]])
    }
  }


  # ### MWS
  # for(s in 1:mws.nspecies){
  #   for(t in 1:mws.nyears){
  #     for(j in 1:mws.nsurveys){
  #       E.mws[j,t,s] <- inprod(alpha.mws[1:mws.ncovs,s], mws.covs[j,t,1:mws.ncovs])
  #       p.mws[j,t,s] <- 1 - (.5)^E.mws[j,t,s]
  # 
  #       y.mws[j,t,s] ~ dbinom(N[mws.gridid[j], mws.sp[s]], z[mws.gridid[j],mws.sp[s]]*p.mws[j,t,s])
  # 
  #       # #Failed attempt to aggregate grids which overlap a survey.
  #       # for(x in 2:mws.ngrids[j]){
  #       #   z.mws.hold[x-1,j,t,s] <- z[mws.gridid[j,x], mws.sp[s]] + z[mws.gridid[j,x-1], mws.sp[s]]
  #       #   N.mws.hold[x-1,j,t,s] <- N[mws.gridid[j,x], mws.sp[s]]*z[mws.gridid[j,x], mws.sp[s]] + N[mws.gridid[j,x-1], mws.sp[s]]*z[mws.gridid[j,x-1], mws.sp[s]]
  #       # }
  #       #
  #       # #If any of the cells were occupied by a species, then the survey as a whole is a 1
  #       # z.mws.combo[j,t,s] <- ifelse(z.mws.hold[mws.ngrids[j]-1,j,t,s] > 0,1,0)
  #       # #Multiply abundance by true occupancy state to get expected abundance for a survey
  #       # N.mws.combo[j,t,s] <- N.mws.hold[mws.ngrids[j]-1,j,t,s]
  #       # y.mws[j,t,s] ~ dbinom(N.mws.combo[j,t,s], z.mws.combo[j,t,s]*p.mws[j,s])
  #     }
  #   }
  # }
  # 
  # ### CBC
  # for(s in 1:cbc.nspecies){
  #   for(j in 1:cbc.nsurveys){
  #       E.cbc[j,s] <- inprod(alpha.cbc[1:cbc.ncovs,s], cbc.covs[j,1:cbc.ncovs])
  #       p.cbc[j,s] <- 1 - (.5)^E.cbc[j,s]
  # 
  #       y.cbc[j,s] ~ dbinom(N[grid.cbc[j], cbc.sp[s]], z[grid.cbc[j],cbc.sp[s]]*p.cbc[j,s])
  #   }
  # }
  
  #########################################
  ### Priors
  # #EBird
  # for(i in 1:ebird.ncovs){
  #   mu.alpha.ebird[i] ~  T(dnorm(0,0.0001), 0,10000)
  #   sig.alpha.ebird[i] ~ T(dt(0, pow(2.5, -2), 1), 0, )
  #   for(s in 1:ebird.nspecies){
  #     alpha.ebird[i,s] ~ T(dnorm(mu.alpha.ebird[i], sd = sig.alpha.ebird[i]), 0,10000)
  #   }
  # }
  # 
  # #BBS
  # for(i in 1:2){
  #   mu.alpha.bbs[i] ~ T(dnorm(0,0.0001), 0,10000)
  #   sig.alpha.bbs[i] ~ T(dt(0, pow(2.5, -2), 1), 0, )
  #   for(s in 1:bbs.nspecies){
  #     alpha.bbs[i,s] ~ T(dnorm(mu.alpha.bbs[i], sd = sig.alpha.bbs[i]), 0,10000)
  #   }
  # }

  #BBL
  for(i in 1:bbl.ncovs){
    mu.alpha.bbl[i] ~ T(dnorm(0,0.0001), 0,10000)
    sig.alpha.bbl[i] ~ T(dt(0, pow(2.5, -2), 1), 0, )
    for(s in 1:bbl.nspecies){
      alpha.bbl[i,s] ~ T(dnorm(mu.alpha.bbl[i], sd = sig.alpha.bbl[i]), 0,10000)
    }
  }

  # #MWS
  # for(i in 1:mws.ncovs){
  #   mu.alpha.mws[i] ~ T(dnorm(0,0.0001), 0,10000)
  #   sig.alpha.mws[i] ~ T(dt(0, pow(2.5, -2), 1), 0, )
  #   for(s in 1:mws.nspecies){
  #     alpha.mws[i,s] ~ T(dnorm(mu.alpha.mws[i], sd = sig.alpha.mws[i]), 0,10000)
  #   }
  # }
  # 
  # #CBC
  # for(i in 1:cbc.ncovs){
  #   mu.alpha.cbc[i] ~ T(dnorm(0,0.0001), 0,10000)
  #   sig.alpha.cbc[i] ~ T(dt(0, pow(2.5, -2), 1), 0, )
  #   for(s in 1:cbc.nspecies){
  #     alpha.cbc[i,s] ~ T(dnorm(mu.alpha.cbc[i], sd = sig.alpha.cbc[i]), 0,10000)
  #   }
  # }
})