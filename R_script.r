estim09=function(y, x, z){
	# fits the model M5 and some of its submodels by ML
	
	# y is a 0/1 variable, 1 denotes the event of interest
	# the model is: y ~ x,  sensitivity ~ z,  specificity = const.
	
	# function value:
	# 	1:3 - p-values comparing M4 vs M5; comparing LZ vs M5; testing dependence of y on x
	# 	4:8 - estimated parameters of M5 
	# 	9 - deviance for M5
	#	10 - convergence for M5 (0 means converged)
	# 	11:14 - estimated parameters of M4 
	# 	15 - deviance for M4
	#	16 - convergence for M4 (0 means converged)
	# 	17:20 - estimated parameters of LZ 
	# 	21 - deviance for LZ
	#	22 - convergence for LZ (0 means converged)
	# 	23:26 - estimated parameters of the model beta1=0 (y does not depend on x) 
	# 	27 - deviance for the model beta1=0
	#	28 - convergence for the model beta1=0 (0 means converged)
	
	meanx=mean(x); x=x-meanx
	meanz=mean(z); z=z-meanz
	d=data.frame(y,x,z)
	
	#log-likehood functions:
			loglik5 = function(param,d){
				# M5
				# d[,1:3] - y (1 the event occurs, 0 not), x, z  
				# param[1:5] - beta0, beta1, gamma0, gamma1, logit(Sp)
				szaml.p = exp(param[1]+param[2]*d[,2])
				p=szaml.p/(szaml.p+1)
				szaml.se = exp(param[3]+param[4]*d[,3])
				Se=szaml.se/(szaml.se+1)
				szaml.sp=exp(param[5])
				Sp=szaml.sp/(szaml.sp+1)
				pobs=p*Se+(1-p)*(1-Sp)
				llik=sum(log(pobs)*d[,1] + log(1-pobs)*(1-d[,1]))
				return(llik)
			}

			loglik4 = function(param,d){
				# M4 (specificity = 1)
				# d[,1:3] - y (1 the event occurs, 0 not), x, z  
				# param[1:4] - beta0, beta1, gamma0, gamma1
				szaml.p = exp(param[1]+param[2]*d[,2])
				p=szaml.p/(szaml.p+1)
				szaml.se = exp(param[3]+param[4]*d[,3])
				Se=szaml.se/(szaml.se+1)
				pobs=p*Se
				llik=sum(log(pobs)*d[,1] + log(1-pobs)*(1-d[,1]))
				return(llik)
			}

			loglik3 = function(param,d){
				# M3 (sensitivity = const. & specificity = 1)
				# d[,1:3] - y (1 the event occurs, 0 not), x, z  
				# param[1:3] - beta0, beta1, gamma0
				szaml.p = exp(param[1]+param[2]*d[,2])
				p=szaml.p/(szaml.p+1)
				szaml.se = exp(param[3])
				Se=szaml.se/(szaml.se+1)
				pobs=p*Se
				llik=sum(log(pobs)*d[,1] + log(1-pobs)*(1-d[,1]))
				return(llik)
			}

			loglik2 = function(param,d){
				# M2 (sensitivity = specificity = 1, that is, no misclassification)
				# d[,1:3] - y (1 the event occurs, 0 not), x, z  
				# param[1:2] - beta0, beta1
				szaml.p = exp(param[1]+param[2]*d[,2])
				p=szaml.p/(szaml.p+1)
				llik=sum(log(p)*d[,1] + log(1-p)*(1-d[,1]))
				return(llik)
			}
			
			loglikSeSp = function(param,d){
				# Liu-Zhang model (both sensitivity and specificity are constant)
				# d[,1:3] - y (1 the event occurs, 0 not), x, z  
				# param[1:4] - beta0, beta1, gamma0, logit(Sp)
				szaml.p = exp(param[1]+param[2]*d[,2])
				p=szaml.p/(szaml.p+1)
				szaml.se = exp(param[3])
				Se=szaml.se/(szaml.se+1)
				szaml.sp = exp(param[4])
				Sp=szaml.sp/(szaml.sp+1)
				pobs=p*Se+(1-p)*(1-Sp)
				llik=sum(log(pobs)*d[,1] + log(1-pobs)*(1-d[,1]))
				return(llik)
			}
				
			loglikbeta1 = function(param,d){
				# y does not depend on x
				# d[,1:3] - y (1 the event occurs, 0 not), x, z  
				# param[1:4] - beta0, gamma0, gamma1, logit(Sp)
				szaml.p = exp(param[1])
				p=szaml.p/(szaml.p+1)
				szaml.se = exp(param[2]+param[3]*d[,3])
				Se=szaml.se/(szaml.se+1)
				szaml.sp=exp(param[4])
				Sp=szaml.sp/(szaml.sp+1)
				pobs=p*Se+(1-p)*(1-Sp)
				llik=sum(log(pobs)*d[,1] + log(1-pobs)*(1-d[,1]))
				return(llik)

				
			}


	# M5
	result5 = optim(c(0,0,0,0,4.5), loglik5, method="Nel", control=list(fnscale=-1, maxit=9000, reltol=1e-15), d=d)
	devi5 = -2*result5$value
	
	result5a = optim(c(0,0,0,0,3), loglik5, method="Nel", control=list(fnscale=-1, maxit=9000, reltol=1e-15), d=d)
	devi5a = -2*result5a$value
	if(devi5a<devi5){result5=result5a; devi5=devi5a}

	result5a = optim(c(0,0,0,0,2), loglik5, method="Nel", control=list(fnscale=-1, maxit=9000, reltol=1e-15), d=d)
	devi5a = -2*result5a$value
	if(devi5a<devi5){result5=result5a; devi5=devi5a}

	result5$par[1]= result5$par[1]-result5$par[2]*meanx
	result5$par[3]= result5$par[3]-result5$par[4]*meanz


	# M4 (specificity = 1)
	result4 = optim(c(0,0,0,0), loglik4, method="Nel", control=list(fnscale=-1, maxit=9000, reltol=1e-15), d=d)
	devi4 = -2*result4$value
	pvalue5=1-pchisq(devi4-devi5,df=1)
	result4$par[1]= result4$par[1]-result4$par[2]*meanx
	result4$par[3]= result4$par[3]-result4$par[4]*meanz
	
	
	# Liu-Zhang (both sensitivity and specificity are constant)
	resultSeSp = optim(c(0,0,3,3), loglikSeSp, method="Nel", control=list(fnscale=-1, maxit=9000, reltol=1e-15), d=d)
	deviSeSp = -2*resultSeSp$value
	
	resultSeSp2 = optim(c(0,0,2,2), loglikSeSp, method="Nel", control=list(fnscale=-1, maxit=9000, reltol=1e-15), d=d)
	deviSeSp2 = -2*resultSeSp2$value
	if(deviSeSp2<deviSeSp){resultSeSp=resultSeSp2; deviSeSp=deviSeSp2}
	
	resultSeSp2 = optim(c(0,0,1,1), loglikSeSp, method="Nel", control=list(fnscale=-1, maxit=9000, reltol=1e-15), d=d)
	deviSeSp2 = -2*resultSeSp2$value
	if(deviSeSp2<deviSeSp){resultSeSp=resultSeSp2; deviSeSp=deviSeSp2}
	
	pvalueSeSp=1-pchisq(deviSeSp-devi5,df=1)
	resultSeSp$par[1]= resultSeSp$par[1]-resultSeSp$par[2]*meanx


	# y does not depend on x
	resultbeta1 = optim(c(0,0,0,3), loglikbeta1, method="Nel", control=list(fnscale=-1, maxit=9000, reltol=1e-15), d=d)
	devibeta1 = -2*resultbeta1$value
	pvaluebeta1=1-pchisq(devibeta1-devi5,df=1)

	
	return(c(pvalue5,pvalueSeSp,pvaluebeta1,result5$par,devi5,result5$convergence,
		result4$par,devi4,result4$convergence,
		resultSeSp$par,deviSeSp,resultSeSp$convergence,
		resultbeta1$par,devibeta1,resultbeta1$convergence))
}




analyse=function(y,x,z,swap=F){
	# calls the function 'estim09' (estim09 does the model fitting, see above)
	# and prints the results and makes plots
	# and returns all results of estim09
	
	# y is the outcome, a 0-1 variable, 1 denotes the event of interest
	#  the model is 
	#	- if swap=F:  y~x, sensitivity~z, specificity=const.
	# 	- if swap=T:  y~x, specificity~z, sensitivity=const.
	
	if(swap) y=1-y
	data=data.frame(y, x, z)
	data.complete=data[complete.cases(data),]

	result=with(data.complete,estim09(y,x,z))
	
	cat("\nLogistic regression allowing misclassification of outcome \n") 
	
	if(!swap)
		cat("(sensitivity is covariate-dependent, specificity is constant) \n")
	else 
		cat("(specificity is covariate-dependent, sensitivity is constant) \n")
	
	#cat("\n"); print(round(result,4))
	cat("\nParameters of model M5:\n beta0:",round(result[4],4)," beta1:",round(result[5],4)," gamma0:",round(result[6],4)," gamma1:",round(result[7],4)," ")

	if(!swap) 
		cat(" Sp =",round(exp(result[8])/(exp(result[8])+1),3),"\n\nM5 vs M4 (Sp=1): p-value =",round(result[1],4),"\n") 
	else 
		cat(" Se =",round(exp(result[8])/(exp(result[8])+1),3),"\n\nM5 vs M4 (Se=1): p-value =",round(result[1],4),"\n")
	
	cat("\nM5 vs Y doesn't depend on X: p-value =",round(result[3],4),"\n")
	
	if(!swap) 
		cat("\nLiu-Zhang: Se =",round(exp(result[19])/(exp(result[19])+1),3),
		"Sp =",round(exp(result[20])/(exp(result[20])+1),3),"   M5 vs Liu-Zhang: p-value =",round(result[2],4),"\n")
	else
		cat("\nLiu-Zhang: Se =",round(exp(result[20])/(exp(result[20])+1),3),
		"Sp =",round(exp(result[19])/(exp(result[19])+1),3),"   M5 vs Liu-Zhang: p-value =",round(result[2],4),"\n")
		
	beta0=result[4]; beta1=result[5]; gamma0=result[6]; gamma1=result[7]
	
	par(mfrow=c(3,1),mar=c(5,6,2,2))
	
	#plot y~x
	xx=seq(min(data.complete$x),max(data.complete$x),len=30)
	linpred1=beta0+beta1*xx; y=exp(linpred1)/(exp(linpred1)+1)
	if(!swap) plot(xx,y,xlab="X",ylab="P(Y_true = 1)",ty="l") 
	else plot(xx,1-y, xlab="X", ylab="P(Y_true = 1)",ty="l")
	
	#plot Se~z
	xz=seq(min(data.complete$z),max(data.complete$z),len=30)
	linpred2=gamma0+gamma1*xz; Se=exp(linpred2)/(exp(linpred2)+1)
	if(!swap) plot(xz,Se,xlab="Z",ylab='Sensitivity \n P(Y_obs = 1 | Y_true = 1)',ty="l") 
	else plot(xz,Se,xlab="Z",ylab='Specificity \n P(Y_obs = 0 | Y_true = 0)',ty="l")

	
	#plot y~x from Liu-Zhang model
	beta0LZ=result[17];beta1LZ=result[18]
	linpred3=beta0LZ+beta1LZ*xx
	if(!swap) {
		y=exp(linpred3)/(exp(linpred3)+1)
	} else {
		y=1-exp(linpred3)/(exp(linpred3)+1)	
	}
	plot(xx,y,xlab="X",ylab="P(Y_true = 1)", main="Liu-Zhang model",ty="l") 
	
	cat("\nmodel: Y~X+Z \n")
	glm01<- glm(y ~ x + z, family=binomial(logit), data=data.complete)
	print(round(summary(glm01)$coefficients,4))
	cat("\nAIC(Y~X+Z) - AIC(M5) =",deviance(glm01)-result[9]-4,"\n")
	
	# returns all results of estim09
	return(invisible(result))
}

#with(sample_data, analyse(y,x,z,swap=F))

