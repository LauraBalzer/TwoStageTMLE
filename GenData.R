##############
# Laura B. Balzer
# Code to generate the data for the simulation studies in
#"Two-Stage TMLE to Reduce Bias and Improve Efficiency in Cluster Randomized Trials"


# MAIN SIMULATION STUDY - BASELINE & POST-BASLINE CAUSES OF MISSINGNESS
gen.data.complex <- function(effect, n, J, verbose=F){
  
  # number of individuals per cluster J
  n <- sample(rep(n, times=J), J) 
  # total number of people
  N <- sum(n)
  # get baseline covariates and cluster-level randomization
  df <- get.baseline.cov.and.treatment(n=n, J=J, N=N)
  
  X1 <- df$X1
  X2 <- df$X2
  X1.c <- df$X1.c
  X2.c <- df$X2.c
  U.c <- df$U3.c
  
  A.ind <- list(1,0, df$A.c)
  
  #-------------------
  # POST-INTERVENTION COVARIATE  
  M <- list()
  U.M <- runif(N)
  # generate under both treatment arms & observed treatment
  for (i in seq_along(A.ind)) {
    p <-plogis(-1 + 2*A.ind[[i]]  + 1*(X1+X2) +
                 0.2*(1-A.ind[[i]])*(X1.c + X2.c) + 0.25*U.c)
    if(verbose){ print(summary(p)) }
    M[[i]] <- as.numeric(U.M <   p)
  }
  
  #-------------------
  # UNDERLYING OUTCOME 
  # strength of impact for treatment & mediator depend on if there is an effect or not
  beta.A <- ifelse(effect, -2.5, 0) 
  beta.M <- ifelse(effect, 4 , 0)
  # generate under both treatment arms & observed treatment
  Y <- list() 
  U.Y <- runif(N)
  for (i in seq_along(A.ind)) {
    p <- plogis(1 + beta.A*A.ind[[i]]  + beta.M*M[[i]] +
                  .5*(X1 + X2) + .2*(X1.c + X2.c) +.25*U.c)
    if(verbose){ print(summary(p))}
    Y[[i]] <- as.numeric(U.Y <  p)
  }
  
  #------------------
  # MEASUREMENT 
  # depends on the cluster-level treatment + baseline & time-dependent covariates
  U.Delta <- runif(N)
  txt <- plogis(3- 3*M[[3]] -0.5*(X1+X2)  )
  con <- plogis(-2 + 3*M[[3]] + 0.5*(X1+X2) ) 
  pscore <- A.ind[[3]]*txt + (1-A.ind[[3]])*con
  if(verbose){ 
    print( summary(pscore[A.ind[[3]]==1]) ) 
    print( summary(pscore[A.ind[[3]]==0]) ) 
  }
  # measurement indicator
  Delta <- as.numeric(U.Delta < pscore)
  
  # return the observed data 
  O<- data.frame(id=df$id, pair=df$pair, U=1, X1=X1, X2=X2,
                 A = A.ind[[3]], M= M[[3]],
                 Delta=Delta, 
                 Y = Y[[3]], Y.1=Y[[1]], Y.0=Y[[2]])
  
  # do some aggregation to get cluster-level variables
  O <- do.aggregation(O)
  
  O
}

# get baseline covariates and cluster-level randomization
# in the MAIN simulation study 
get.baseline.cov.and.treatment <- function(n, J, N){
  
  # cluster indicator
  id <-  rep(1:J, times = n)
  
  # CLUSTER-LEVEL ERROR
  U1.cc <- runif(J, -1, 1)
  U2.cc <- runif(J, -1, 1)
  U3.cc <- rnorm(J, 0, 1)
  #-----------------------------------------------
  # CONVERT CLUSTER COVARIATES TO LONG-FORM INDV
  U1.c <- rep(U1.cc, times=n)
  U2.c <- rep(U2.cc, times=n)
  U3.c <- rep(U3.cc, times=n)
  
  #------------------------------
  # INDV COVARIATES
  X1 <- rnorm(N, U1.c, .5)
  X2 <- rnorm(N, U2.c, .5)
  
  # now make observed clust-cov the mean of the ind
  X1.c <- rep(aggregate(X1, by=list(id), mean)$x, times=n)
  X2.c <- rep(aggregate(X2, by=list(id), mean)$x, times=n)
  
  #-------------------------
  # RANDOMIZE INTERVENTION within matched pairs of clusters
  temp <- match.me(U3.cc)
  A.c <- temp$A.c
  pair <- temp$pairs
  A.c <- rep(A.c, times=n)
  pair <- rep(pair, times=n)
  
  data.frame(cbind(id=id, U1.c, U2.c, U3.c, X1.c, X2.c, X1, X2, A.c, pair))
}


# ---------------------------------------
# SIMPLIFIED SIMULATION STUDY - BASELINE (only) CAUSES OF MISSINGNESS
gen.data.simple <- function(effect, n, J, verbose=T){
  
  # number of individuals per cluster J
  n <- sample(rep(n, times=J), J) 
  # total number of people
  N <- sum(n)
  # cluster indicator
  id <-  rep(1:J, times = n)
  
  # cluster-level factors 
  U1.c <- runif(J, 1.75, 2.25)
  U2.c <- rnorm(J, 0, 1)
  U3.c <- rnorm(J, 0, 1)
  
  # RANDOMIZE TREATMENT within matched pairs of clusters 
  temp <- match.me(U3.c)
  A.c <- temp$A.c
  pair <- temp$pairs
  pair <- rep(pair, times=n)
  
  A.ind <- list(1,0, rep(A.c, times = n))
  #-----------------------------------------------
  # at indv-level
  U1.c <- rep(U1.c, times=n)
  U2.c <- rep(U2.c, times=n)
  U3.c <- rep(U3.c, times=n)
  
  # indv covariate
  X1 <- rnorm(N, U1.c, 1)
  X2 <- rnorm(N, U2.c, 1)
  
  # now make observed clust-cov the mean of the ind
  X1.c <- rep(aggregate(X1, by=list(id), mean)$x, times=n)
  X2.c <- rep(aggregate(X2, by=list(id), mean)$x, times=n)
  
  # UNDERLYING OUTCOME 
  # impact of treatment depends on whether or not there is an effect 
  if(effect){
    beta.A <- 0.15
    beta.A.X <- 0.15
  } else{
    beta.A <- beta.A.X <- 0 
  }
  U.Y <- runif(N)
  Y <- list() 
  for (i in seq_along(A.ind)) {
    p <- plogis(-4 + beta.A*A.ind[[i]] + 0.3*X1.c + 0.3*X2.c +
                  0.4*X1 +  0.2*X2 + 0.5*X1.c*X1+ beta.A.X*X1*A.ind[[i]] +
                  .3*U3.c )
    Y[[i]] <- as.numeric(U.Y <  p)
  }
  
  # MEASUREMENT
  U.Delta <- runif(N)
  #  observation depends on cluster-level txt & baseline covariates
  pscore <-  plogis(4 -0.25*A.ind[[3]] -.5*X1.c -.1*X2.c -.1*X2 - 0.75*X1 -.75*X1*A.ind[[3]])
  Delta <- as.numeric(U.Delta < pscore)
  
  
  O<- data.frame(id=id, pair, U=1, X1=X1, X2=X2,
                 A = A.ind[[3]],
                 Delta=Delta,
                 Y = Y[[3]], Y.1=Y[[1]], Y.0=Y[[2]])
  
  # do some aggregation to get cluster-level variables
  O <- do.aggregation(O)
  
  O
}

#################################################
# OTHER HELPER FUNCTIONS TO GENERATE THE DATA 
# function to get cluster-level variables as aggregates of indv-level variables
do.aggregation <- function(O){
  
  # cluster-level covariates
  X.c <- aggregate(O[, c('X1','X2')], by=list(O$id), mean)
  # number of participants
  n.C <- aggregate(O[,'U'], by=list(O$id), sum)$x
  # making long form (i.e., at individual-level)
  X1.c <- rep(X.c[,'X1'], times=n.C)
  X2.c <- rep(X.c[,'X2'], times=n.C)
  nIndv <- rep(n.C, times=n.C)
  new <- data.frame(cbind(O, nIndv, X1.c, X2.c))
  colnames(new) <- c(colnames(O), 'nIndv','X1.c', 'X2.c')
  new
}

# function to finalize the observed data 
get.obs <- function(O, exclude){
  
  # set observed outcome Y as (measurement indicator)x(underlying outcome)
  O$Y <- O$Y*O$Delta  
  if(exclude){
    # if subsetting the data to those who are measured
    O <- O[O$Delta==1, ]
    O <- subset(O, select=-c(Y.1,Y.0, X1.c, X2.c))
    O <- do.aggregation(O=O)
    
  }else{
    # if not subsetting on complete observations, just drop counterfactual outcomes
    O <- subset(O, select=-c(Y.1,Y.0))
    
  }
  O 
}


#########################################
# Pair match clusters, and then randomize within pairs 
match.me <- function(matchData){
  
  dist<- distancematrix(gendistance(data.frame(matchData)))
  matches<- nonbimatch(dist)
  # matches contains ids for the pair match as well as the distance measure
  grpA<- as.numeric(matches$halves[,'Group1.Row'])
  grpB<- as.numeric(matches$halves[,'Group2.Row'])	
  
  J <- length(matchData)
  n.pairs <- J/2
  pairs <- A.c <- rep(NA, J)
  pairs[grpA] <- 1:n.pairs
  pairs[grpB] <- 1:n.pairs
  
  A <- rbinom(n.pairs, 1, .5)
  A2 <- ifelse(A,0,1)
  A.c[grpA] <- A
  A.c[grpB] <- A2
  data.frame(pairs, A.c)
}



########################################################
# CALCULATING THE TRUTH 

# get the counterfactual outcomes for the population 
get.truth.pop<- function(do.complex, effect, n, J){
  if(do.complex){
    # main simulation study: post-baseline causes of measurement
    X.all <- gen.data.complex(effect=effect, n=n, J=J, verbose=F)
  }  else{
    # supplementary study: only baseline causes of measurement
    X.all <- gen.data.simple(effect=effect, n=n, J=J, verbose=F)
  }
  truth <- get.truth(X.all=X.all)
  truth
}


# get the true values for the treatment-specific means & effect measures
get.truth <- function(X.all){
  
  # aggregate to the data to the cluster-level 
  truth.C <- aggregate(X.all[,c('U','Y.1', 'Y.0')], by=list(X.all$id), sum)
  # cluster-level counterfactual outcomes
  Y1.c <- truth.C$Y.1/truth.C$U
  Y0.c <- truth.C$Y.0/truth.C$U
  # treatment-specific means
  Risk1.C <- mean(Y1.c)
  Risk0.C <-  mean(Y0.c)
  # effect on the absolute scale and relative scale
  RD <- Risk1.C - Risk0.C
  aRR <-  Risk1.C / Risk0.C
  
  # measures of variability by arm + with/without pairing
  km<- getMeasuresVariability(X0=Y0.c, X1=Y1.c,
                              pairs= aggregate(X.all$pair, by=list(X.all$id), mean)[,2]  )
  
  data.frame(Risk1.C= Risk1.C, Risk0.C, RD, aRR, km)
  
}


# simple measures of w/in cluster variability by arm + with/without pairing
getMeasuresVariability<- function(X1, X0, pairs){	
  
  k.con<- sd(X0)/mean(X0)
  k.txt<- sd(X1)/mean(X1)
  
  temp<- unique(pairs)
  n.pairs<- length(temp)
  varm<- pi <- rep(NA, n.pairs)
  for(i in 1:n.pairs){		
    varm[i]<- var(X0[ pairs==temp[i]])
    pi[i] <- mean(X0[ pairs==temp[i]])
  }
  km <- mean(sqrt(varm)/pi)
  data.frame(k.con, k.txt, km)
}
