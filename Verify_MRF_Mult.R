
rm(list = ls())

source("mixmrf_mult.R")
source("Utility.R")

#Theta
source("thetaInit.R")
source("theta.optim.R")

##Phi
source("loglikPhi.R")
source("logPhi0.R")
source("logPhi1.R")
source("grrLoglikPhi.R")
source("grrLogPhi0.R")
source("grrLogPhi1.R")
source("X0Init.R")

#################################################################
#read gaussian data
all <- read.csv("../Data/gs.data.csv", string = F)
ground.truth <- all$Key

all <- all[, -c(2, 3,5, 9,19)]



all <- cbind(price = all[,1], Intercept = rep(1, nrow(all)), all[,-1])

pred <- as.matrix(all[,-1])
resp <- all[,1]



#read multinomial data
mult <- read.csv("../Data/mult.data.csv", string = F)
mult <- as.matrix(mult)

mult[mult!=0] <- 1


#read neighboring data
load("../Data/nb.Rdata")
nb <- final.mrf.mat

#################################################################
#browser()
mrf.mult.start <- proc.time()
rslt <- mxmrf.mult(resp, pred, mult, nb)
mrf.mult.time <- proc.time() - mrf.mult.start
browser()
save(rslt, file = "gamma.mrf.mult.RData")

gamma <- cbind(ground.truth, rslt$postProb0, resp, pred, mult)
write.csv(gamma, "rslt.mrf.mult.csv", row.names = F, quote = T)



