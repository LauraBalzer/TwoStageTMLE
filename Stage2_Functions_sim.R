##############
# Stage2_Functions.R 
# R code to implement all Stage2 analyses to compare intervention effects between arms
#
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Lead Statistician for SEARCH
#
# Input: Stage1 estimated cluster-level outcome & cluster-level covariates/exposure 
# Outputs: an estimate of the intervention effect
# Goal option specifies if arithmetic risk ratio (aRR) or risk difference (RD)
# Weighting option specifies to weight clusters or individuals equally when computing
#		the intervention effect. 
##---------------------------
#################



#-----------------------------------------------------#-----------------------------------------------------
# Stage2: function to estimate the intervention effect
#
# input: 
#	  goal;  aRR= arithmetic risk ratio; otherwise risk difference (RD)
#	  weighting; "clust" weights clusters equally; otherwise weights indv equally,
#	  observed data (data.input),
#	  name of column corresponding to the outcome (outcome),
#	  survival: indicator prob survival inputted but outcome is prob of badness
#		  currently applicable to TB/mortality
#	  set of candidate adjustment variables; must be at cluster-level (clust.adj)
#	  prespecified adjustment variables for estimation of the conditional mean outcome with main terms (e.g. Y~A+W1+W2) (QAdj), 
#	  prespecified adjustment variables for estimation of the propensity score with main terms (e.g. A~W1+W2) (gAdj),
#	  indicator to do adaptive prespecification (do.data.adapt),
#	  indicator to break the pair-match (break.match)
#	  indicator to print updates (verbose)
# output: point estimate and inference
#
# note the data.input must have
#		1) a column 'id' indicating cluster membership 
#		2) a column 'U' as dummy adjustment variable (=1 for all)
#		both of these requirements are taken care with preprocess function
#
# note2: if you are running adaptive prespecification to select the candidate C-TMLE 
#   minimizing the variance (which nearly all Stage 2 analyses do), 
#   please set do.data.adapt=T and specify the candidate adjustment variables in clust.adj... 
#   Please do NOT use (QAdj, gAdj) which are reserved for running non-adaptive glms only. 
#
# Note3: when specifying clust.adj, please include the dummy variable U 
#   to make sure that the unadjusted estimator is considered as a candidate. 
#   e.g. clust.adj <- c(‘U’, ‘hiv_prev_0’, ‘mc_cover_0’)
#
# note4:to run the unadjusted estimator than set clust.adj=NULL and do.data.adapt=F.
#
# EXAMPLES
## Stage2 estimation of the arithmetic risk ratio (ARR), weighting clusters equally, & preserving the matches
#
## First with TMLE with adaptive pre-specification
# clust.adj <- c('U', 'chc_prev', 'bmi')
# Yc1.adj <- Stage2(goal=‘aRR’, weighting='clust’, data.input=pri, 
#                   outcome='Yc', clust.adj=clust.adj,  do.data.adapt=T,  break.match=F) 
#                   
## Then with the unadjusted estimator
#  Yc1.unadj <- Stage2(goal=‘aRR’, weighting=‘cluster-level’, data.input=pri, 
#               outcome= 'YcU', clust.adj=NULL, do.data.adapt=F,  break.match=F)
#-----------------------------------------------------#-----------------------------------------------------

Stage2 <- function(goal='aRR', weighting='clust', data.input,
	outcome=NULL,  survival=F,
	clust.adj, QAdj=NULL, gAdj=NULL, 
	do.data.adapt =F, break.match=F, verbose=F, psi=NULL){	

	outcome.string <- make.full.word(outcome)
	names(data.input)[grep(outcome.string, names(data.input) ) ] <- 'Y'

	nIndv.string <- make.full.word( paste('nIndv', outcome, sep='_') )
	names(data.input)[grep(nIndv.string, names(data.input) ) ] <- 'nIndv'

	# if input survival probabilities, but the outcome is probability of badness (e.g. TB)
	if(survival){
		data.input$Y <- 1 - data.input$Y
		print('Flipping outcome to be probability of badness')
	}
	
	# if haven't already, get weights 
	if( sum(grep('alpha', colnames(data.input)))==0){
		data.input <-  get.weights(data.input, weighting=weighting) 
	}


	if(do.data.adapt){
		# implement adaptive prespecificaton
		select<- do.adaptive.prespec(goal=goal, weighting=weighting, break.match = break.match,
			Ldata= data.input, clust.adj=clust.adj, QAdj=QAdj, gAdj=gAdj)
	
		QAdj = select$QAdj
		gAdj = select$gAdj	
	}
	
	# Run full TMLE algorithm with prespecified adjustment set
	est<- do.TMLE(goal=goal, weighting=weighting, train=data.input,  QAdj=QAdj,  
			 gAdj=gAdj, verbose=verbose)
		
	# Get point estimates and inference
	R1=est$R1
	R0=est$R0
	

	
	if( goal=='aRR' ){
	  psi.hat <- log(R1/R0)
	} else if (goal=='RD'){
	  psi.hat <- R1- R0
	} else if (goal=='OR'){
	  psi.hat <- log( R1/(1-R1)*(1-R0)/R0)
	}
	
	n.clust <- length(unique(data.input$id)) 
	n.pair <- length(unique(data.input$pair))
	
	if(break.match){
		# if breaking the match, set df to (#clusters -2)
		df <-n.clust - 2
		se <- sqrt(est$var.break)
		
	} else{
		# if preserving the match, set df to (#pairs-1)
		df <- n.pair-1 
		se <- sqrt(est$var.pair)
		
	}
	inference <- get.inference(goal=goal, psi=psi, psi.hat=psi.hat, se=se, df=df )
	
	
	Txt<- get.CI(psi.hat=R1, var.IC=est$var.R1, df=(n.clust-2) )
	Con <- get.CI(psi.hat=R0, var.IC=est$var.R0, df=(n.clust-2) )

	RETURN<-  data.frame(Txt=Txt, Con=Con, inference,  QAdj=est$QAdj, gAdj=est$gAdj)

	  
	RETURN
}

#-----------------------------------------------------#-----------------------------------------------------
# do.TMLE: master function to run TMLE 
#
# input: 
#		goal - aRR: arithmetic risk ratio; O/W risk difference,
#		training data (train),
#		weighting; "clust" weights clusters equally; O/W weights indv equally,
#		prespecified adjustment variables for the conditional mean outcome (QAdj), 
#		prespecified adjustment variables for the propensity score (gAdj),
#		initial estimator of the conditional mean outcome (Q.out),
#		estimator of the propensity score (p.out),
#		indicator to print updates (verbose)
#
# output: list 
#		training data augmented with estimates,
#		prespecified adjustment variables for the conditional mean outcome (QAdj), 
#		prespecified adjustment variables for the propensity score (gAdj),
#		initial estimator of the conditional mean outcome (Q.out),
#		estimator of the propensity score (p.out),
#		estimted fluctutation coef (epsilon),
#		estimated risk under txt and control (R1, R0),
#		estimated variance preserving the match (var.pair),
# 		estimated variance breaking the match (var.break)
#-----------------------------------------------------#-----------------------------------------------------

do.TMLE <- function(goal, train, weighting='clust', 
	QAdj, gAdj=NULL, Q.out=NULL, p.out=NULL,
	verbose=F){	
	
	J <- length(unique(train$id))
	
	# if haven't already, aggregate to the cluster-level & add weights
	if( nrow(train) >  J){  		
		train <- aggregate(train, by=list(train$id), mean)[,2: (ncol(train)+1)]
	}
	
	# if haven't already, add a column for weights
	if( sum(grep('alpha', colnames(train)))==0){
		train <- get.weights(train, weighting=weighting)
	}
	
	
	#=====================================================
	# Step1 - initial estimation of E(Y|A,W)= Qbar(A,W)
	#=====================================================
	
	# run glm on the adjustment set
	Q<- do.Init.Qbar(train=train, QAdj=QAdj, glm.out=Q.out , verbose=verbose)
	train <- Q$train
	
	#==========================================================
	# Step2: Calculate the clever covariate
	#==========================================================	

	G <- get.clever.cov(train=train, gAdj=gAdj, p.out=p.out,  verbose=verbose)
	train <- G$train

	#==========================================================
	# Step3: Targeting
	#==========================================================
	
	eps <- get.epsilon(train=train, goal=goal, verbose=verbose)

	train <- do.targeting(train=train, eps=eps, goal=goal)
			
	#==========================================================
	# Step4: Parameter estimation
	#==========================================================
	# aggregated data: the weights are 1 for cluster-effect & n_j*J/nTot for individual-level effect 
	# note weights are normalized to sum to J		
	R1 <- mean( train$alpha*train$Qbar1W.star )
	R0 <- mean( train$alpha*train$Qbar0W.star ) 
	
	variance.out <- get.IC.variance(goal=goal, Vdata=train, R1=R1, R0=R0)

	RETURN<- list(train=train, 	
		QAdj=Q$QAdj, Q.out=Q$glm.out,
		gAdj=G$gAdj, p.out=G$p.out, 
		eps=eps, 	R1=R1, R0=R0, 
		var.R1=variance.out$var.R1, 
		var.R0=variance.out$var.R0,
		var.pair=variance.out$var.pair, 
		var.break=variance.out$var.break)	
	RETURN
}

#-----------------------------------------------------#-----------------------------------------------------
# do.Init.Qbar - function to do initial estimation of E[Y|A,W] = Qbar(A,W)
# 	input: 
#		data set, adjustment variable(s), outcome regression fit on training set, verbose
# 	output: 
#		adjustment variable(s)
#		outcome regression fit on training set
#		data set augmented w/ initial predictions: Qbar(A,W), Qbar(1,W) and Qbar(0,W)
#-----------------------------------------------------#-----------------------------------------------------
do.Init.Qbar<- function(train, QAdj, glm.out=NULL, verbose=F){
	
	if( is.null(QAdj) ){
		QAdj<- 'U'
	}

	train.temp <- train[, c(QAdj, 'A', 'Y')]
	
	X1<- X0<- train.temp
	X1$A<-1; X0$A<- 0	
		
	if( is.null(glm.out) ){
		# fit using the training data
		glm.out<- suppressWarnings( glm( Y~. , family='binomial', data=train.temp, weights=train$alpha ) )	
		if(verbose) print(glm.out)
	}	
	
	# get initial predictions
	QbarAW <- predict(glm.out, newdata=train.temp, type='response')
	Qbar1W <- predict(glm.out, newdata=X1, type='response')
	Qbar0W <- predict(glm.out, newdata=X0, type='response')
	
	list(QAdj=QAdj, glm.out=glm.out, train=data.frame(train, QbarAW , Qbar1W , Qbar0W ))
}

#-----------------------------------------------------#-----------------------------------------------------
# get.clever.cov - function calculate the clever covariate
# 	input: 
#		data set, adjustment variable(s), pscore regression fit on training set, verbose
# 	output: 
#		adjustment variable(s), pscore regression fit on training set, 
#		data set augmented with pscore & clever covariate (H.AW, H.1W, H.0W)
#-----------------------------------------------------#-----------------------------------------------------
get.clever.cov<- function(train, gAdj, p.out=NULL, verbose=F){

	if( is.null(gAdj) ){
		gAdj <- 'U'
	}

	train.temp <- train[, c(gAdj, 'A')]  

	if( is.null(p.out) ){
		
		# fit pscore on training set 	
		p.out<-   suppressWarnings( glm( A~. , family='binomial', data= train.temp, weights=train$alpha) )
		if(verbose){ print(p.out)}
	}

	# now use p.out to get estimated pscores
	pscore <- predict(p.out, newdata= train.temp,  type="response")
	# Note if gAdj=U & train$alpha!=1, then this will differ from 0.5

	# bound g - should not apply for a randomized trial
	pscore [pscore < 0.025] <- 0.025
	pscore [pscore > 0.975] <- 0.975

	A.train <- train$A
	# Clever covariate is two-dimensional; 
	H.1W <- A.train/pscore 
	H.0W <- (1-A.train)/(1-pscore )
	# via delta method
	H.AW <- H.1W - H.0W


	list(gAdj=gAdj, p.out=p.out,  train=data.frame(train, pscore, H.1W , H.0W , H.AW) ) 
}	
	
	
#-----------------------------------------------------#-----------------------------------------------------
# get.epsilon - function calculate the fluctuation coefficient
# 	input: 
#		data set, goal with 'aRR'=arithmetic RR, verbose
# 	output: 
#		estimated fluctuation coefficient (eps)
#-----------------------------------------------------#-----------------------------------------------------
get.epsilon <- function(train, goal, verbose=F){
	
	A.train<- train$A
	Y.train<- train$Y
		
	# Skip fitting if outcome=0 for all observations in either txt or control group
	Skip.update <-  mean(Y.train[A.train==1])==0 | mean(Y.train[A.train==0])==0 |  
		mean(Y.train[A.train==1])==1 | mean(Y.train[A.train==0])==1 
			
	if(goal=='RD'){
		
		# if going after RD, then use a 1-dim clever covariate
		if(!Skip.update){
			logitUpdate<- suppressWarnings( 
				glm(Y.train ~ -1 +offset(qlogis(train$QbarAW )) + train$H.AW, family="binomial",  weights=train$alpha))
			eps<-logitUpdate$coef
		} else{
			eps<- 0
		}
		names(eps) <- 'H.AW'
	
	}else{
		# targeting the risk or odds ratio requires a two-dimensional clever covariate
		
		if( !Skip.update  ){
			logitUpdate<- suppressWarnings(
				glm(Y.train ~ -1 +offset(qlogis(train$QbarAW )) + train$H.0W + train$H.1W, family="binomial", weights=train$alpha))
			eps<-logitUpdate$coef
		} else{
			eps <- c(0,0)
		}
		names(eps)<- c('H.0W', 'H.1W')	
	}
	if(verbose) print(eps)

	eps
}

#-----------------------------------------------------#-----------------------------------------------------
# do.targeting - function to update initial estimators of QbarAW
# 	input: 
#		data set (train), fluctuation coefficient (eps), goal (aRR= arithmetic risk ratio; otherwise RD)
# 	output: 
#		data.frame w/ targeted predictions: Qbar*(A,W), Qbar*(1,W), Qbar*(0,W)
#-----------------------------------------------------#-----------------------------------------------------

do.targeting <- function(train, eps, goal){
	
	g1W<- train$pscore
	g0W<- (1 - g1W)
	
	if(goal=='RD'){
		
		# updated QbarAW estimates for training set. 
		QbarAW.star <- plogis( qlogis(train$QbarAW ) + eps*train$H.AW)	
		Qbar1W.star <- plogis( qlogis(train$Qbar1W ) + eps/g1W )
		Qbar0W.star <- plogis( qlogis(train$Qbar0W ) - eps/g0W )
		
	}else{
		# updated QbarAW estimates for training set. 
		QbarAW.star <- plogis( qlogis(train$QbarAW) + eps['H.0W']*train$H.0W + eps['H.1W']*train$H.1W)	
		Qbar0W.star <- plogis( qlogis(train$Qbar0W) + eps['H.0W']/g0W )
		Qbar1W.star <- plogis( qlogis(train$Qbar1W) + eps['H.1W']/g1W )
	}
	train <- data.frame(train, QbarAW.star, Qbar1W.star, Qbar0W.star)		
	train
}
	

#-----------------------------------------------------#-----------------------------------------------------
# get.IC.variance - function to do influence curve-based variance estimate 
# 	input: 
#		goal (aRR= arithmetic risk ratio; otherwise RD)
#		dataset
#		risk estimates under txt and control R1 & R0
# 	output: on log scale for if goal='aRR'
#		estimated IC & variance - preserving/breaking the match
#-----------------------------------------------------#-----------------------------------------------------
get.IC.variance<- function(goal, Vdata, R1, R0){
	
	# aggregated data, wt=1 for cluster-level effect & wt=n_j*J/nTot for indv
	# again noting that weights are normalized to sum to J
	
	DY1 <- Vdata$alpha*Vdata$H.1W*(Vdata$Y - Vdata$Qbar1W.star)
	DY0 <- Vdata$alpha*Vdata$H.0W*(Vdata$Y - Vdata$Qbar0W.star) 	

	
	if(goal=='RD'){
		# going after RD, easy IC
		DY <-  DY1 - DY0
	} else if (goal=='aRR'){ 
		
		# going after aRR, then get IC estimate on log scale
		#	i.e. Delta method for log(aRR) = log(R1) - log(R0)
		DY <- 1/R1*DY1 - 1/R0*DY0
	} else if(goal=='OR'){
	  # Delta method for log(OR)
	  DY <- 1/R1*DY1 + 1/(1-R1)*DY1 - 1/(1-R0)*DY0 - 1/R0*DY0
	}
	
	# estimated variance for txt specific means or if break the match	
	J<- length( unique(Vdata$id) )
	var.R1 <- var(DY1) /J
	var.R0 <- var(DY0) / J
	var.break <- var(DY) /J


	# estimated variance if preserve the match
	pairs <- unique(Vdata$pair)
	n.pairs <- length(pairs)
	DY.paired <-  rep(NA, n.pairs)
	for(i in 1:n.pairs){		
			these<- Vdata$pair== pairs[i] 
			DY.paired[i]<- 0.5*sum(DY[ these] )			
	}
	
	var.pair <- var(DY.paired) / n.pairs

	
	list(var.R1=var.R1, var.R0=var.R0, DY=DY, var.break=var.break, DY.paired=DY.paired, var.pair=var.pair)
}


#-----------------------------------------------------#-----------------------------------------------------
# get.CI: mini function to calculate two-sided (1-sig.level)% confidence intervals 
#	input: 
# 		pt estimate
# 		var.IC (variance estimate)
#		df (degrees of freedom for Student's t-dist ) 
#		sig.level (significance level)
# output: 
#		data.frame with point est, lower/higher 95% CI
#-----------------------------------------------------#-----------------------------------------------------	


get.CI <- function(psi.hat,  var.IC, df, sig.level=0.05){
	
	# cutoff based on t-dist for testing and CI	
  	cutoff <- qt(sig.level/2, df=df, lower.tail=F)
  	
  	# standard error (square root of the variance)
  	se<- sqrt(var.IC)
  	# 95% confidence interval 
  	CI.lo <- (psi.hat - cutoff*se)
  	CI.hi <- (psi.hat + cutoff*se)
	  out<- data.frame(est=psi.hat, CI.lo, CI.hi, var.hat=var.IC)
  	out
}

#-----------------------------------------------------#-----------------------------------------------------
# get.inference: function to calculate (1-sig.level)% confidence intervals & test the null hypothesis
#	input: 
#		goal (aRR= arithmetic risk ratio; otherwise RD)
#   psi (true value)
#   psi.hat (estimate)
#   se (standard error)
#		df (degrees of freedom if using a Student's t-dist ) 
#		sig.level (significance level)
# output: 
#		variance, test statistic, confidence intervals, pval, indicator reject null
# 		note: if goal=aRR, variance & test stat are on log-scale
#-----------------------------------------------------#-----------------------------------------------------	

get.inference <- function(goal, psi=NA, psi.hat, se, df=NA, sig.level=0.05){
  
  # test statistic (on the transformed scale)
  tstat <- psi.hat/se
  
  if(is.na(df)){
    # ass normal
    cutoff <- qnorm(sig.level/2, lower.tail=F)
    pval<- 2*pnorm(abs(tstat), lower.tail=F) 
  }else{
    # cutoff based on t-dist for testing and CI	
    cutoff <- qt(sig.level/2, df=df, lower.tail=F)
    pval<- 2*pt(abs(tstat), df=df, lower.tail=F)
  }

  # reject the null
  reject <- pval < sig.level 
  
  # 95% confidence interval 
  CI.lo <- (psi.hat - cutoff*se)
  CI.hi <- (psi.hat + cutoff*se)
  
  
  if(goal!='RD'){
    psi.hat<- exp(psi.hat)
    CI.lo <- exp(CI.lo)
    CI.hi <- exp(CI.hi)
  }  
  
  # confidence interval coverage
  cover<- ( CI.lo <= psi & psi <= CI.hi )
  #bias <- psi.hat - psi
  data.frame(psi, psi.hat,  se=se, tstat, CI.lo, CI.hi, pval, #bias, 
             cover, reject)
  
}




#----------------------------------------------
# get.weights: function to get appropriate weights for the analysis
# input: data & weighting scheme (clusters equal or indiv equal)
# output: data with a column alpha for weights
#----------------------------------------------

################
get.weights<- function(data, weighting){
	
	# reset any previously inputted weights
	data$alpha <- -99
	
	data.clust <- aggregate(data[,c('id', 'nIndv')], by=list(data$id), mean)
	nTot <- sum(data.clust$nIndv)
	J <- nrow(data.clust) 

	if( nrow(data) >  J ){   
		# if not aggregated to the cluster-level
		if(weighting=='clust'){
			# corresponds to the empirical mean in each cluster
			data$alpha  <- 1/data$nIndv
		} else{
			data$alpha  <- J/nTot
		}
	}else{ 
		# aggregated to cluster-level
		if(weighting=='clust'){
			data$alpha  <- 1
		} else{
			data$alpha  <- data$nIndv*J/nTot
		}
	} 		
	data
	
}

############
make.full.word <- function(string){
	paste("\\<", string, "\\>", sep="")
}




get.worst.pair <- function(data.input, variable.name='hiv_prev_0') {
	temp <- rep(NA, 16)
	 for(j in 1:16){
 		these<- data.input[data.input$pair==j, variable.name]
 		temp[j] <- abs( these[1]-these[2])
 		# print(c(j, temp[j]))
 	}
 	temp
}



