lapply(c("nimble", "dplyr", "parallel", "coda", "MCMCvis", "matrixStats", "ggplot2", "sf"), require, character.only = T)

#Load original file,
#abdu.files <- list(mcmcList.param, mcmcList.coef, code, data)
load("./Output/iSDM_ABDU.rdata")
str(abdu.files)

grid <- st_read(dsn = "./SamplingGrid.shp")

est.psi.raw <- MCMCsummary(mcmcList.param, 'psi')
summary(est.psi.raw)

est.psi <- est.psi.raw %>% 
  select(Psi = mean, Psi.LCL = '2.5%', Psi.UCL = '97.5%', psi.Rhat = Rhat) %>%
  mutate(GridID = 1:nrow(.))
  
est.lambda.raw <- MCMCsummary(mcmcList.param, 'lambda')
summary(est.lambda.raw)

est.lambda <- est.lambda.raw %>% 
  select(Lambda = mean, Lambda.LCL = '2.5%', Lambda.UCL = '97.5%', lambda.Rhat = Rhat) %>%
  mutate(GridID = 1:nrow(.))

grid.par <- grid %>%
  merge(., est.psi, by = "GridID", all = T) %>%
  merge(., est.lambda, by = "GridID", all = T)

ggplot(data = grid.par) +
  geom_sf(aes(fill = Psi))

ggplot(data = grid.par) +
  geom_sf(aes(fill = log(Lambda)))
