##############
# Laura B. Balzer
# Main R code to generate the data and implement the estimators
# in the simulation studies for 
#"Two-Stage TMLE to Reduce Bias and Improve Efficiency in Cluster Randomized Trials"


getTrialData<- function(do.complex, effect, n, J, pop.truth, SL.library='glm', 
                        dropM=F, verbose=F){
  
  # GENERATE THE DATA   
  if(do.complex){
    # main simulation study 
    # post-baseline causes of missingness
    X.all <- gen.data.complex(effect=effect, n=n, J=J, verbose=verbose)
    # indv-level covariates
    ind.cov <- c('X1', 'X2','M')
  }  else{
    # supplementary simulation study 
    # only baseline causes of missingness
    X.all <- gen.data.simple(effect=effect, n=n, J=J, verbose=verbose)
    # indv-level covariates
    ind.cov <- c('X1', 'X2')
  }
  # cluster-level covariates
  clust.cov <- c('X1.c', 'X2.c')

  
  #TRUTH
  if(is.na(pop.truth[1]) ){
    # can also calculate the sample-specific means & effects 
    truth <- get.truth(X.all=X.all)
  }else{
    truth <- pop.truth
  }
  
  #OBSERVED DATA
  O <- get.obs(O=X.all, exclude=F)
  
  # STAGE 1 ESTIMATION if the cluster-level endpoint Y^c in each 
  # cluster separately 
  # (1) the mean among measured (unadjusted approach)
  # (2) TMLE to control for differential measurement (adjusted appproach)
  data.C <- suppressWarnings(getYc(X.all=X.all[,c('A','Y.1','Y.0')], O=O, 
                                   ind.cov=ind.cov, clust.cov=clust.cov,
                                   SL.library = SL.library ))
  # colMeans(data.C[data.C$A==1,])
  # colMeans(data.C[data.C$A==0,])
  
  # T-TEST on the cluster-level means (among those measured)
  # breaking the match
  ttest.b <- output.ttest(psi=truth$RD, 
                          g=t.test(x=data.C[data.C$A==1, 'unadj'],
                                   y=data.C[data.C$A==0, 'unadj'], var.equal=T), 
                          gRR=F, paired=F)
  # keeping the match
  temp<- unique(data.C$pair)
  n.pairs <- length(temp)
  Y.0p <- Y.1p <- rep(NA, n.pairs)
  for(i in 1:n.pairs){		
    Y.0p[i]<- data.C[ data.C$pair==temp[i] & data.C$A==0, 'unadj']
    Y.1p[i]<- data.C[ data.C$pair==temp[i] & data.C$A==1, 'unadj']
  }
  ttest.p <- output.ttest(psi=truth$RD,
                          g=t.test(y=Y.0p, x=Y.1p, paired=T),
                          gRR=F, paired=T)
  
  #---------------------------- 
  # STAGE 2 TMLEs using Adaptive Prespecification for intervention effects
  # Risk difference 
  TT.RD <- do.tmles(goal='RD', psi=truth$RD, data.C=data.C, clust.cov=clust.cov)
  # Risk ratio
  TT.RR <- do.tmles(goal='aRR', psi=truth$aRR, data.C=data.C, clust.cov=clust.cov)
  #----------------------------
  
  # if dropping M from the adjustment set in the complex setting
  if(do.complex & dropM){
    ind.cov <- c('X1','X2')
  }    
  # all the covariates
  cov <- c(ind.cov, clust.cov)
  
  # DR-GEE for risk ratio
  aug.RR <- do.gee.aug(train=O, psi=truth$aRR, link="poisson", do.complex=do.complex, dropM=dropM)

  # COMPLETE CASE ANALYSES
  ## subset on fully observed observations 
  train <- get.obs(O=X.all, exclude=T)

  # COVARIATE ADJUSTED RESIDUALS ESTIMATORS (CARE) for risk difference
  care.b <- do.care( train=train,  ind.cov=cov, psi=truth$RD, pair=F)
  care.p <- do.care( train=train,  ind.cov=cov, psi=truth$RD, pair=T)

  # GEE for risk ratio 
  gee.b <- do.gee(train=train, ind.cov=cov, psi=truth$aRR, paired=F, link='poisson')
  gee.p <- do.gee(train=train, ind.cov=cov, psi=truth$aRR, paired=T, link='poisson')

  # MIXED MODELS for risk ratio
  mixed.b <- do.mixed.models(train=train, psi=truth$aRR, paired=F, link='poisson',
                             do.complex = do.complex, dropM=dropM)
  mixed.p <- do.mixed.models(train=train, psi=truth$aRR, paired=T, link='poisson',
                             do.complex = do.complex, dropM=dropM)

  # RETURN for all estimators
  RETURN <- list(ttest.b=ttest.b,
                 ttest.p=ttest.p,
                 gee.b=gee.b, gee.p=gee.p,
                 care.b=care.b,care.p=care.p,
                 mixed.b=mixed.b,mixed.p=mixed.p,
                 aug.RR=aug.RR,
                 tmle.RD.b = TT.RD['ignore',],
                 tmle.RD.p = TT.RD['pair',],
                 tmle.RR.b = TT.RR['ignore',],
                 tmle.RR.p = TT.RR['pair',],
                 truth=truth)
  RETURN
  
}	

