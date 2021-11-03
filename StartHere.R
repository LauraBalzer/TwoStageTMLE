##############
#  Laura B. Balzer
# Wrapper R code for simulation studies - 
#"Two-Stage TMLE to Reduce Bias and Improve Efficiency in Cluster Randomized Trials"

rm(list=ls())

library('nbpMatching')
library('MASS')
library('SuperLearner')
library('ltmle')

# load functions to generate the data
source('MainFunction.R')
source('GenData.R')

# load functions to do estimation 
source('Estimators.R')

# the following code is adapted from or taken directly from 
# https://github.com/LauraBalzer/SEARCH_Analysis_Adults
source('Stage2_Functions_sim.R')
source('Adapt_Functions.R')

set.seed(1)
# whether there is an intervention or not
effect <- T
# do.complex = T for main simulation study
# do.complex = F for supplementary simulation study
do.complex <- T
# dropM = T explore performance of GEEs, mixed models, & CARE when not adjusting for M
# (only applicable to the complex DGP)
dropM <- T

# number of clusters
J <- 30
# number of participants per cluster
n <- c(100, 150, 200)
# set number of simulation repetitions
nReps <- 500
# SuperLearner library for individual-level TMLE to control for differential missingness
SL.library <- c('SL.mean', 'SL.glm', 'SL.gam')

# a bunch of data frames to save the results 
truth<- data.frame(matrix(NA, nrow=nReps, ncol=7) )
aug.RR<-  data.frame(matrix(NA, nrow=nReps, ncol=9))
colnames(aug.RR) <- c("psi", "psi.hat", "se", "tstat", "CI.lo", "CI.hi", "pval", "cover", "reject")
ttest.b <- ttest.p <- gee.b <- gee.p <- care.b <- care.p <- mixed.b <- mixed.p <- aug.p <- aug.RR
tmle.RR.b <- tmle.RR.p <- tmle.RD.b <- tmle.RD.p <-  data.frame(matrix(NA, nrow=nReps, ncol=19))

# calculate the true values of the population-means and effects
pop.truth <- get.truth.pop(do.complex=do.complex, effect=effect, n=n, J=5000)

# for nReps 
for(i in 1:nReps){
  
  out<- getTrialData(do.complex=do.complex, effect=effect, n=n, J=J, pop.truth=pop.truth,
                    SL.library=SL.library, dropM=dropM, verbose=F)
  truth[i,]<- out$truth
  
  ttest.b[i,] <- out$ttest.b
  ttest.p[i,] <- out$ttest.p
  
  gee.b[i,] <- out$gee.b
  gee.p[i,] <- out$gee.p
  
  care.b[i,] <- out$care.b
  care.p[i,] <- out$care.p
  
  mixed.b[i,] <- out$mixed.b
  mixed.p[i,] <- out$mixed.p

  aug.RR[i,] <- out$aug.RR
  
  tmle.RR.b[i,] <- out$tmle.RR.b
  tmle.RR.p[i,] <- out$tmle.RR.p
  tmle.RD.b[i,] <- out$tmle.RD.b
  tmle.RD.p[i,] <- out$tmle.RD.p

  print(i)
  
}

colnames(truth)<- colnames(out$truth)
colnames(tmle.RR.b) <- colnames(tmle.RR.p) <- colnames(tmle.RD.b) <- colnames(tmle.RD.p) <- colnames(out$tmle.RD.b)

# quick function to summarize the simulation results 
make.pretty <- function(ttest, care, mixed, gee, augRR,  tmleRR, tmleRD){
  make.pretty.mini <- function(est, RD=T){
      these <- c("psi", "psi.hat", "se", "tstat", "CI.lo", "CI.hi", "pval", "cover", "reject")
      if(RD){
        bias<- mean(est$psi.hat - est$psi)
        se.mc <-  sqrt(var(est$psi.hat))
      }else{
        bias<- mean(est$psi.hat/est$psi)
        se.mc <-  sqrt(var(log(est$psi.hat)) )
      }
      c(colMeans(est[,these]), bias=bias, se.mc=se.mc )
  }
  
  rbind(ttest.RD = make.pretty.mini(ttest),
        CARE.RD = make.pretty.mini(care),
        tmle.RD = make.pretty.mini(tmleRD),
        mixed.RR = make.pretty.mini(mixed, RD=F), 
        gee.RR = make.pretty.mini(gee, RD=F), 
        aug.RR = make.pretty.mini(augRR, RD=F), 
        tmle.RR = make.pretty.mini(tmleRR, RD=F)
  )
}


# Results breaking the matches used for randomization 
B <- make.pretty(ttest.b, care.b, mixed.b, gee.b, aug.RR, tmle.RR.b, tmle.RD.b)
round(B, 2)

# Results keeping the matches used for randomization 
P <- make.pretty(ttest.p, care.p, mixed.p, gee.p, aug.p,  tmle.RR.p, tmle.RD.p)
round(P,2)

colMeans(truth, na.rm=T)


file.name <- paste('MAIN', do.complex, 'J', J, 'effect', effect,
                   ifelse(!do.complex | dropM, 'noM', 'wM'),
                   'reps', nReps, format(Sys.time(),"%d%b%Y"), 'Rdata', sep='.' )
save(truth, J, n, ttest.b, care.b, mixed.b, gee.b, aug.RR,tmle.RD.b, tmle.RR.b,
     ttest.p, care.p, mixed.p, gee.p, aug.p, tmle.RD.p, tmle.RR.p, SL.library,
     file=file.name)
