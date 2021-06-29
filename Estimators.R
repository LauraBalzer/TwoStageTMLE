##############
# Laura B. Balzer
# Main R code to implement the estimators in the simulation studies for 
#"Two-Stage TMLE to Reduce Bias and Improve Efficiency in Cluster Randomized Trials"


# STAGE 1: estimation of cluster-specific endpoints
getYc<- function(X.all, O, ind.cov, clust.cov, SL.library='glm'){
  
  clusters <- unique(O$id)
  nClust <- length(clusters)
  
  clust.var <- c('id', 'pair', 'nIndv','U', 'A', clust.cov)
  data.clust <- data.frame(matrix(NA, nrow=nClust, ncol=(length(clust.var)+2)))
  
  for(j in 1:nClust){ # for each cluster separately
    
    # subset the data on the cluster of interest 
    these.units <- clusters[j]==O$id
    OC <- O[these.units, ]
    XC <- colMeans( X.all[these.units, ] )

    # unadjusted estimator: mean among those measured 
    unadj <- mean(OC[OC$Delta==1, 'Y'])
    
    # individual-level TMLE to adjust for differential measurement 
    OC.sub <- OC[,c(ind.cov, 'Delta','Y')]
    suppressMessages(
      Yc <-  ltmle(data=OC.sub, Anodes='Delta', Ynodes='Y', abar=1,
                   SL.library=SL.library,
                   estimate.time=F, stratify=T)
    )
    tmle <- Yc$estimates['tmle']
    
    data.clust[j,] <-  c(OC[1, clust.var], unadj, tmle) 
    
  }  
  colnames(data.clust) <- c(clust.var, 'unadj', 'tmle')
  data.clust
}


# STAGE 2 estimation of intervention effect with TMLE using
# adaptive pre-specificaiotn 
do.tmles <- function(goal, psi, data.C, clust.cov){
  
  # wrapper function to put data in proper format for Stage2 functions
  reformat.stage2<- function(data, outcome, indv){
    colnames(data)[grep(outcome, colnames(data) )] <- 'Y'
    colnames(data)[grep(indv, colnames(data))] <- 'nIndv_Y'
    data
  }
  
  data.temp <- reformat.stage2(data.C, outcome='tmle', indv='nIndv') 
  
  # ignoring (i.e., breaking) the matches used for randomization 
  tmle.A.A <- suppressWarnings(Stage2(goal=goal, weighting='clust', data.input=data.temp,
                                      outcome='Y', clust.adj=clust.cov, 
                                      do.data.adapt =T, break.match=T, verbose=F, psi=psi))
  
  # keeping the matches used for randomization 
  tmle.A.A.P <- suppressWarnings(Stage2(goal=goal, weighting='clust', data.input=data.temp,
                                      outcome='Y', clust.adj=clust.cov, 
                                      do.data.adapt =T, break.match=F, verbose=F, psi=psi))
  
  data.frame(rbind(ignore=tmle.A.A, pair=tmle.A.A.P))
}




# COVARIATE ADJUSTED RESIDUALS ESTIMATOR (CARE)
do.care <- function(train, ind.cov, psi, pair){
  
  clust <- aggregate(train, by=list(train$id), mean )[,c('id', 'pair', 'A', 'Y')]
  
  if(pair){
    n.pairs <- sum(!duplicated(clust$pair)) 
    temp <- data.frame(  train[, c('pair',ind.cov, 'Y') ])
    temp$pair <- as.factor(temp$pair)
  }else{
    temp <- train[, c(ind.cov, 'Y') ]
  }
  # pooled indv regression of outcome on covariates but not txt 
  ind.reg <- glm(Y~  ., data=temp, family='binomial')
  ind.risk <- predict(ind.reg, type='response')
  pred <- aggregate(ind.risk, by=list(train$id), mean)$x
  resid <- clust$Y- pred
  clust <- cbind(clust, pred, resid)
  txt <- clust[clust$A==1, ]
  con <- clust[clust$A==0,]
  
  if(pair){
    red1 <- red0 <- rep(NA, n.pairs)
    temp<- unique(clust$pair)
    for(i in 1:n.pairs){		
      red1[i]<- txt[ txt$pair==temp[i], 'resid'] 
      red0[i]<- con[ con$pair==temp[i], 'resid'] 
    }
    out <- output.ttest(psi=psi, g= t.test(red1, red0, paired=T, var.equal=T),
                        gRR=F, paired=T)
  }else{
    red1 <- txt$resid
    red0 <- con$resid
    out <- output.ttest(psi=psi, g= t.test(red1, red0, paired=F, var.equal=T),
                        gRR=F, paired=F)
  }
  out

}

# T-TEST on cluster-level means 
output.ttest<- function(psi, g, gRR=F, paired=F){
	
  if(paired){
    psi.hat<- g$est
  } else{
    psi.hat <- g$est[1]- g$est[2]
  }
	pval <- g$p.value
	reject<- pval< 0.05
	CI.lo <- g$conf.int[1]
	CI.hi <- g$conf.int[2]
	if(gRR){
		psi.hat <- g$est
		CI.lo <-  exp(CI.lo)
		CI.hi <- exp(CI.hi)
	} 
	
	cover<- ( CI.lo <= psi & psi <= CI.hi )
	
  out<- data.frame(psi, psi.hat, se=g$stderr,  tstat=g$stat,  CI.lo, CI.hi, pval,
                   cover, reject)
	out
}


# GEE: Generalized estimating equations 
library(geepack)
do.gee <- function(train, ind.cov, psi, paired=F, link){
  
  N <- sum(!duplicated(train$id))
  id <- train$id
  pair <- train$pair
  if(paired){
    train$pair <- as.factor(train$pair)
    train <- train[, c('pair', ind.cov, 'A', 'Y')]
  }else{
    train <- train[, c(ind.cov, 'A', 'Y')]
  }
  m <- geeglm(Y~ ., data=train, family=link, id=id)
  psi.hat<- coef(m)['A']
  se <- summary(m)$coefficients['A', 'Std.err']
  get.inference(goal='aRR', psi=psi, psi.hat=psi.hat, se=se, df=NA)
}

# MIXED MODELS (a.k.a., random effect models)
library('lme4')
do.mixed.models <- function(train, paired=F, psi, link, do.complex){
  
  N <- sum(!duplicated(train$id))
  if(do.complex){
    if(paired){
      m <-  glmer(Y~ A+ X1 +X2 + M + X1.c + X2.c +(1 | pair),
                  data=train, family=link ) 
    } else{
      m <-  glmer(Y~ A+ X1 +X2 + M + X1.c + X2.c + (1 | id),
                  data=train, family=link ) 
    }
  } else{
    if(paired){
      m <-  glmer(Y~ A+ X1 +X2 + X1.c + X2.c +(1 | pair),
                  data=train, family=link ) 
    } else{
      m <-  glmer(Y~ A+ X1 +X2 + X1.c + X2.c + (1 | id),
                  data=train, family=link ) 
    }
  }
  
  psi.hat <- summary(m)$coef['A','Estimate']
  se <- summary(m)$coef['A','Std. Error']
  get.inference(goal='aRR', psi=psi, psi.hat=psi.hat, se=se, df=NA)
}

# DOUBLE-ROBUST GEE (DR-GEE)
library('CRTgeeDR')
do.gee.aug <- function(train, psi, link , do.complex=T){
  N <- sum(!duplicated(train$id))
  train[train$Delta==0,'Y']<- NA
  if(do.complex){
    out <- geeDREstimation(formula=Y~A, nameY='Y', id='id', nameTRT='A', nameMISS='Delta',
                           data=train, family=link, 
                           model.augmentation.ctrl = Y~ X1 + X2 + M + X1.c + X2.c,
                           model.augmentation.trt = Y ~ X1 + X2 + M + X1.c + X2.c,
                           model.weights = Delta ~ X1 + X2+ M +  X1.c + X2.c + A
    )
  } else{
    out <- geeDREstimation(formula=Y~A, nameY='Y', id='id', nameTRT='A', nameMISS='Delta',
                           data=train, family=link, 
                           model.augmentation.ctrl = Y~ X1 + X2 + X1.c + X2.c,
                           model.augmentation.trt = Y ~ X1 + X2 + X1.c + X2.c,
                           model.weights = Delta ~ X1 + X2+  X1.c + X2.c + A
    )
  }

  psi.hat <- summary(out)$beta[2]
  se <- summary(out)$se.robust[2]    
  get.inference(goal='aRR', psi=psi, psi.hat=psi.hat, se=se, df=NA ) 
  
 }
