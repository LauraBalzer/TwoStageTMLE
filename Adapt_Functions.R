##############
# Adapt_Functions.R 
# R code to implement adaptive prespecification
#
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Lead Statistician for SEARCH
#---------------------------


#-----------------------------------------------------#-----------------------------------------------------
# do.adaptive.prespec: function to implement adaptive prespecification as described in 
#		Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
# input: goal of analysis ('aRR' or 'RD'), 
#		weighting (clust/indv), 
#		indicator to break the match (break.match)
#		dataset (Ldata)
#		cluster-level adjustment variables (clust.adj),
#		selected adjustment variable for the outcome regression (QAdj)
#		selected adjustment variable for the pscore regression (gAdj)
#	output: selection for candidate TMLE
#-----------------------------------------------------#-----------------------------------------------------
do.adaptive.prespec<- function(goal, weighting, break.match, Ldata, 
	clust.adj, QAdj=NULL, gAdj=NULL){
	
	if( is.null(QAdj) ){
		
		# do adaptive prespecification with cluster-level candidate covariates
		select.Q <- suppressWarnings( CV.selector(goal=goal,  weighting=weighting, break.match=break.match,
				 Ldata=Ldata, clust.adj=clust.adj, forQ=T) )
		QAdj <- select.Q$Adj
		
		# if select unadjusted estimator for QbarAW=E(Y|A,W), then stop
		if( QAdj == 'U' ){ 
			gAdj <- 'U'
			var.CV <- select.Q$var.CV
		} else{
			
			# if did not select the unadjusted for initial estimation of Qbar(A,W), 
			# then need to remove this variable from the candidate for pscore estimation
			clust.adj<- clust.adj[ -which(clust.adj==select.Q$Adj) ]		
		}
	}
	
	if( is.null(gAdj) ){ 		
		select.G<- suppressWarnings( CV.selector(goal=goal, weighting=weighting, break.match=break.match, 
			Ldata=Ldata, clust.adj=clust.adj, forQ=F, QAdj= QAdj) )
		gAdj <- select.G$Adj
		var.CV <- select.G$var.CV		
	}		
	
	list(QAdj=QAdj, gAdj=gAdj, var.CV=var.CV )
}



#-----------------------------------------------------#-----------------------------------------------------
# CV.selector: function to estimate the cross-validated risk
#		Loss function is the squared-IC; Risk is then the variance of the TMLE
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: goal of analysis ('aRR' or 'RD), 
#		weighting (clust/indv), 
#		indicator to break the match (break.match)
#		dataset (Ldata)
#		cluster-level adjustment variables
#		indicator if for the conditional mean outcome (forQ)
#		selected adjustment variable for the outcome regression (QAdj)
#	output: selection for adjustment variable (corresponding to a TMLE)
#-----------------------------------------------------#-----------------------------------------------------
CV.selector <- function(goal, weighting, break.match, Ldata, clust.adj, forQ, QAdj){

	# covariates correspond to different TMLEs  
	num.tmles <- length(clust.adj)
	CV.risk <-  var.CV <-  rep(NA, num.tmles)
	
	for(k in 1: num.tmles){	

		if(forQ){
			# if selecting the adjustment set for the outcome regression
			IC.temp<- get.IC.CV(goal=goal, weighting=weighting, break.match=break.match, 
				 Ldata=Ldata,  QAdj=clust.adj[k], gAdj=NULL)
		} else{
			# if collaboratively selecting the adjustment set for the pscore
			IC.temp<- get.IC.CV(goal=goal, weighting=weighting, break.match=break.match,
				Ldata=Ldata,  QAdj=QAdj,  gAdj= clust.adj[k] )
		}
		# estimating the CV risk for each candidate
		CV.risk[k]<- IC.temp$CV.risk
		# estimating the CV variance for that TMLE
		var.CV[k] <- IC.temp$var.CV

	}
	# select the candidate estimator resulting in the smallest CV-risk
	this<- which.min(CV.risk)
 	list(Adj= clust.adj[this ], var.CV=var.CV[this])
}

#-----------------------------------------------------#-----------------------------------------------------
# getIC.CV: function to obtain a cross-validated estimate of the influence curve
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: goal of analysis ('aRR' or 'RD), 
#		weighting (clust/indv), 
#		indicator to break the match (break.match)
#		dataset (Ldata)
#		selected adjustment variable for the outcome regression (QAdj),
#		selected adjustment variable for the pscore (gAdj),
#	output: cross-validated estimate of the IC for pair
#-----------------------------------------------------#-----------------------------------------------------

get.IC.CV<- function(goal, weighting, break.match, Ldata, QAdj, gAdj){

	# This is implemented for leave-one-out 
	if( !break.match ){
		indpt.unit <- Ldata$pair
	}else{
		indpt.unit <- Ldata$id
	}
	folds <- unique(indpt.unit)
	nFolds <- length(folds)
	DY.CV <- rep(NA, nFolds)
	
	# doing a cross-validated estimate
	for(i in 1: nFolds) {
		
		these.units <- folds[i]==indpt.unit
		valid <- Ldata[these.units, ]
		train <- Ldata[!these.units,]
	
		# run full TMLE algorithm on the training set
		train.out <- do.TMLE(goal=goal, train=train, weighting=weighting, 
			QAdj=QAdj, gAdj=gAdj, verbose=F)	
	
		# get the relevant components of the IC for the validiation set, using fits based on the training set
		valid.out <- do.TMLE.validset(goal=goal, valid=valid, weighting=weighting,
			train.out=train.out)
		
		if(break.match){
			DY.CV[i]<- valid.out$DY
		}else{
			DY.CV[i] <- valid.out$DY.paired
		}
	}
	# estimating the CV risk for each candidate
	CV.risk <- mean( DY.CV^2 )
	# estimating the CV variance for that TMLE
	var.CV <- var(DY.CV)/length(DY.CV)

	list(CV.risk=CV.risk, var.CV=var.CV)

}
	

#-----------------------------------------------------#-----------------------------------------------------
# do.TMLE.for.valid: function to obtain a cross-validated estimate of the influence curve
#	for observations in the validation set
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: goal of analysis ('aRR' or 'RD'),
#		validation dataset ('valid') 
#		weighting (clust/indv), 
#		TMLE-fits from training set (train.out)
#	output: cross-validated estimate of the IC,
#		cross-validated risk estimate (loss=IC^2)
#-----------------------------------------------------#-----------------------------------------------------

do.TMLE.validset <- function(goal, valid, weighting, train.out) {		
	
	J <- length(unique(valid$id) )
	
	# if not already aggregated
	if(nrow(valid) >  J){  		
		valid <- aggregate(valid, by=list(valid$id), mean)[,2: (ncol(valid)+1)]
	}
	# if don't have weights - should NOT be true
	if( sum(grep('alpha', colnames(valid)))==0){
		valid <- get.weights(valid, weighting=weighting)
	}
	

	#=============================================
	# Step1 - initial estimation of E(Y|A,W)= Qbar(A,W)
	#=============================================

	valid<- do.Init.Qbar(train=valid, QAdj=train.out$QAdj, glm.out=train.out$Q.out)$train
	
	#=============================================
	# Step2: Calculate the clever covariate
	#=============================================
	
	valid <- get.clever.cov(train=valid, gAdj=train.out$gAdj, p.out=train.out$p.out)$train

	#=============================================
	# Step3: Targeting - 			
	#=============================================
	
	valid <- do.targeting(train=valid, eps=train.out$eps, goal=goal)
	
	get.IC.variance(goal=goal, Vdata=valid, R1=train.out$R1, R0=train.out$R0)
	
}





