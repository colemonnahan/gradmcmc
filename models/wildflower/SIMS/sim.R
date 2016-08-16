#library(bbmle)
library(lme4)

fruitdat1 <- read.csv("../../../DATA/wildflower_pods_complete.csv")
fruitdat1$State = factor(fruitdat1$State, levels = c("F", "L", "M", "S", "D"))

if (FALSE) {
  parms <- read.table("../../true.dat",header=TRUE)
  rownames(parms) <- parms$parameter

  nplant <- length(levels(fruitdat1$PlantID))
  nyear <- length(unique(fruitdat1$Year))

  library(Matrix)

  X <- model.matrix(~State+Pods,data=fruitdat1)
  Z.plantid <- sparse.model.matrix(~PlantID-1,data=fruitdat1)
  Z.plantid_pods <- sparse.model.matrix(~(PlantID-1):Pods,data=fruitdat1)
  Z.year <- sparse.model.matrix(~Year-1,data=fruitdat1)
  
  beta_fix <- parms$value[1:6]
  u_plantid <- rnorm(nplant,sd=sqrt(parms$value["plant.int.var"]))
  u_plantid_pods <- rnorm(nplant,sd=sqrt(parms$value["plant.pods.var"]))
  u_year <- rnorm(nyear,sd=sqrt(parms$value["year.int.var"]))
}

load("../../wildflower.fit.RData")
d <- as.data.frame(fruitdat1)

simdata <- function(seed) {
  if (!missing(seed)) set.seed(seed)
  transform(d,toF=unlist(simulate(Rfit)))
}

dat<-simdata()

# Write data in the same filename and format as the real data 

write.csv(dat,file="wildflower_pods_complete.csv",row.names=FALSE)
