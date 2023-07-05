### Parameter/Index List
# i = surface grid cell
# s = survey
# t = year
# cell = dynamic indexing to track overlap of point/cell0
# Z = true occurrence
# psi = probability of occurence


### PROCESS/SAMPLING MODELS ###
# Occupancy Surface
for(i in 1:n.cells){
  Z[i] ~ dbern(psi[i])
  logit(psi[i]) = beta0 + beta*Covs[i] + SPDE[i]
}

### DATA MODELS ###
# MWS
for(s in 1:mws.n.survey){
  for(t in 1:mws.years){
    #Effort
    E.mws[s,t] = alpha*Covs[s,t]
    p.mws[s,t] = 1 - (1 - p)^E.mws[j,t]
    
    #Data Model
    log(lambda.mws[s]) = gamma*Covs[s]
    N.mws[s,t] ~ dpois(lambda.mws[s])
    y.mws[s,t] ~ dbinom(N.mws[s,t], Z[max(mws.grid[s]])*p.mws[s,t])
  }
}


# BBS
for(s in 1:bbs.n.survey){
  for(t in 1:bbs.years){
    #Effort
    E.bbs[s,t] = alpha*Cov[s,t]
    p.bbs[s,t] = 1 - (1 - p)^E.bbs[j,t]
    
    #Data Model
    log(lambda.mws[s]) = gamma*Covs[s]
    N.mws[s,t] ~ dpois(lambda.mws[s])
    y.bbs[s,t] ~ dbinom(N[s,t], Z[bbs.grid[s]]*p.bbs[j,t])
  }
}


# eBIRD
for(s in 1:ebird.n.survey){
  for(t in 1:ebird.years){
    #Effort
    E.ebird[j,t] = alpha*Cov[s,t]
    p.ebird[j,t] = 1 - (1 - p)^E.ebird[j,t]
    
    #Data Model
    y.ebird[s,t] ~ dpois(E*(Z[cell.ebird[s,t]]*lambda[s]+p.ebird[t]))
  }
}


### PRIORS ###
p = .5