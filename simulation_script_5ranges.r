# simulation with 5 ranges

K.datasets.norm.x=function(K,N,x1range,x2range,Sp){

	proby=function(x) exp(x)/(exp(x)+1)
	
	sens=function(x) exp(x)/(exp(x)+1)
	
	results=matrix(nr=K,nc=28)
	k=1; trials=1
	repeat{
		if(k > K) break
		x1=rnorm(N,me=mean(x1range),sd=(x1range[2]-x1range[1])/4)
		x2=rnorm(N,me=mean(x2range),sd=(x2range[2]-x2range[1])/4)
		# about 95% of the X values are in the given range
		
		p=proby(x1)
		Se=sens(x2)
		
		outcome=rbinom(N,si=1,p=p*Se+(1-p)*(1-Sp))
		
		Dat=data.frame(y=outcome, x1=x1, x2=x2)
		param.estim=estim06power(Dat[,1],Dat[,2],Dat[,3])
		if(param.estim[10] == 0) {
			results[k,]=param.estim
			#cat("pvalues",round(param.estim[c(1,2,3,9,10,15,16,21,22,27,28)],4),"\n")
			k=k+1
		}
		trials=trials+1
	}
	trials=trials-1
	cat(trials," ",K," ",N," ",x1range[1]," ",x2range[1]," ",Sp," ",mean(results[,1]<=.05),mean(results[,2]<=.05),mean(results[,3]<=.05),"\n")
	return(invisible(1))
}

# values of predictors are normally distributed

# in the results, first comes the power for Sp<1, 
#   then the power for dependence of Se on x2, 
#   finally the power for dependence of y on x1

# "range1" (-4.6, -.85)    0.01 - 0.3
# "range2" (-1.39, 0)       0.2 - 0.5 
# "range3" (-.62, .62)     0.35 - 0.65
# "range4" (0, 1.39)        0.5 - 0.8
# "range5" (.85, 4.6)     0.7 - 0.99

for (N in c(500,1000,2000,5000,10000))
for (x1 in c("range1","range2","range3","range4","range5"))
for (x2 in c("range1","range2","range3","range4","range5"))
for (Sp in c(1, .999, .9, .8, .7)){
	if(x1=="range1") x1range=c(-4.6, -.85)
	if(x1=="range2") x1range=c(-1.39, 0)  
	if(x1=="range3") x1range=c(-.62, .62) 
	if(x1=="range4") x1range=c(0, 1.39)   
	if(x1=="range5") x1range=c(.85, 4.6)  
	if(x2=="range1") x2range=c(-4.6, -.85)
	if(x2=="range2") x2range=c(-1.39, 0)  
	if(x2=="range3") x2range=c(-.62, .62) 
	if(x2=="range4") x2range=c(0, 1.39)   
	if(x2=="range5") x2range=c(.85, 4.6)  
	K.datasets.norm.x(500,N,x1range,x2range,Sp)
}

