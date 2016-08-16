## modified from Bolker et al wildflower model.




## sessionInfo()

## rawData <- read.csv("../DATA/wildflower_pods_complete.csv",
##                     header=TRUE, row.names = 1)
## year <- rawData$Year - min(rawData$Year) + 1  ## convert to an index
## plant <- as.integer(rawData$PlantID)
## stage <- as.integer(factor(rawData$State,
##                            levels = c("F","L","M","S","D"))) ## match ADMB/R ordering
## ## levels = c("D","S","M","L","F"))) ## Intuitive ordering of states
## Nyear <- length(unique(year))
## Nplant <- length(unique(plant))
## Ndata <- nrow(rawData)
## Nstage <- length(unique(stage))
## data <- list(year = year, plant = plant, stage = stage, Pods = rawData$Pods, Nstage = Nstage, Nyear = Nyear, Nplant = Nplant, Ndata = Ndata, toF = rawData$toF)
## saveRDS(data, file='data.RDS')

library(R2jags)
data <- readRDS(file='data.RDS')
pars <-
    c("yearInterceptSD", "plantInterceptSD", "plantSlopeSD", "intercept",
      "slope", "yearInterceptEffect", "plantSlopeEffect", "plantInterceptEffect")
inits <-
    list(list(yearInterceptSD = 1,
              plantInterceptSD = 1,
              plantSlopeSD = 1,
              intercept = rep(0,dataForJags$Nstage), slope = 0,
              yearInterceptEffect = rep(0, data$Nyear),
              plantInterceptEffect = rep(0, data$Nplant),
              plantSlopeEffect = rep(0, data$Nplant)))
jags.fit <- jags(model="wildflower_BUGS.txt",
                 data = data, parameters.to.save = pars,
                 n.chains = 1, inits = inits)

set.seed(81571)
