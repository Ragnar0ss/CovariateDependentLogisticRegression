analyse=function(y,x,z,R=200,show.progress=T){
	# analysis of covariate-dependent misclassification
	# calls the function 'estim09'
	
	# y is a 0-1 variable, 1 denotes the event of interest
	# Se = P(y_obs=1 | y_true=1)
	# Sp = P(y_obs=0 | y_true=0)
	
	# model: logit(y_true) = beta0 + beta1*x, logit(Se) = gamma0 + gamma1*z, Sp = const.
	# R is the number of jackknife replications (for the SE of the estimates)
	#    if R=0, no SE is calculated
	
	data=data.frame(y, x, z)
	data.complete=data[complete.cases(data),]
	
	if(show.progress) {cat("\nFitting the model...\n"); flush.console()}
	
	# fit the models
	result=with(data.complete,estim09(y,x,z))
	
	# calculate Se and/or Sp from logitSe and/or logitSp 
	result[8]=exp(result[8])/(1+exp(result[8]))
	result[19]=exp(result[19])/(1+exp(result[19]))
	result[20]=exp(result[20])/(1+exp(result[20]))
	result[26]=exp(result[26])/(1+exp(result[26]))
		
	names(result)=c(
		"p(M5.vs.Sp=1)",
		"p(M5.vs.gamma1=0)",
		"p(M5.vs.beta1=0)",
		"beta0.M5",
		"beta1.M5",
		"gamma0.M5",
		"gamma1.M5",
		"Sp.M5",
		"deviance.M5",
		"converg.M5",
		"beta0.M4",
		"beta1.M4",
		"gamma0.M4",
		"gamma1.M4",
		"deviance.M4",
		"converg.M4",
		"beta0.LZ",
		"beta1.LZ",
		"Se.LZ",
		"Sp.LZ",
		"deviance.LZ",
		"converg.LZ",
		"beta0.M0",
		"gamma0.M0",
		"gamma1.M0",
		"Sp.M0",
		"deviance.M0",
		"converg.M0")
	
	# calculate jackknife SEs and CIs
	N=nrow(data.complete)
	ranjack=sample(1:N,R)  
	jackparam=matrix(nr=R,nc=28)

	jackse=NA
	progress=floor(.05*R); plus=max(3*progress,8)
	if(R>0){
		if(show.progress) {cat("Calculating SEs..."); flush.console()}
		for (r in 1:R){
			if(show.progress) if(r==progress) {
				cat(paste("  ",floor(100*r/R),"%",sep=""))
				progress=progress+plus
				flush.console()
			}
			da=data.complete[-ranjack[r],]
			jackparam[r,]=estim09(da[,1],da[,2],da[,3])
			jackparam[r,8]=exp(jackparam[r,8])/(1+exp(jackparam[r,8]))
			jackparam[r,19]=exp(jackparam[r,19])/(1+exp(jackparam[r,19]))
			jackparam[r,20]=exp(jackparam[r,20])/(1+exp(jackparam[r,20]))
			jackparam[r,26]=exp(jackparam[r,26])/(1+exp(jackparam[r,26]))

		}
		cat("\n")

		jackmean=apply(jackparam,2,mean)
		jackbias=(N-1)*(jackmean-result[4:10])
		jackse=(apply(jackparam,2,var)*(R-1)/R*(N-1))^.5

		names(jackse)=c(
			" ",
			" ",
			" ",
			"SE.beta0.M5",
			"SE.beta1.M5",
			"SE.gamma0.M5",
			"SE.gamma1.M5",
			"SE.Sp.M5",
			" ",
			" ",
			"SE.beta0.M4",
			"SE.beta1.M4",
			"SE.gamma0.M4",
			"SE.gamma1.M4",
			" ",
			" ",
			"SE.beta0.LZ",
			"SE.beta1.LZ",
			"SE.Se.LZ",
			"SE.Sp.LZ",
			" ",
			" ",
			"SE.beta0.M0",
			"SE.gamma0.M0",
			"SE.gamma1.M0",
			"SE.Sp.M0",
			" ",
			" ")
	}
	return(list(pvalues=result[1:3],estimates=result[4:28],jack.SEs=jackse[4:28]))
}

