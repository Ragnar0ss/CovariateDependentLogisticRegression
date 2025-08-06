# load estim09.r

R.datasets.norm.x=function(R,N,xrange,zrange,Sp){

	proby=function(x) exp(x)/(exp(x)+1)
	
	sens=function(x) exp(x)/(exp(x)+1)
	
	results=matrix(nr=R,nc=28)
	k=1; trials=1
	repeat{
		if(k > R) break
		x=rnorm(N,me=mean(xrange),sd=(xrange[2]-xrange[1])/4)
		z=rnorm(N,me=mean(zrange),sd=(zrange[2]-zrange[1])/4)
		# about 95% of the X values are in the given range
		
		p=proby(x)
		Se=sens(z)
		
		outcome=rbinom(N,si=1,p=p*Se+(1-p)*(1-Sp))
		
		Dat=data.frame(y=outcome, x=x, z=z)
		param.estim=estim09(Dat[,1],Dat[,2],Dat[,3])
		if(param.estim[10] == 0) {
			results[k,]=param.estim
			#cat("pvalues",round(param.estim[c(1,2,3,9,10,15,16,21,22,27,28)],4),"\n")
			k=k+1
		}
		trials=trials+1
	}
	trials=trials-1
	cat(trials," ",R," ",N," ",xrange[1]," ",zrange[1]," ",Sp," ",mean(results[,1]<=.05),mean(results[,2]<=.05),mean(results[,3]<=.05),"\n")
	
	return(data.frame("1q",trials,R,N,xrange[1],zrange[1],Sp,mean(results[,1]<=.05),mean(results[,2]<=.05),mean(results[,3]<=.05)))
}

# values of predictors are normally distributed

# first comes the power for Sp<1, 
#   then the power for dependence of Se on z, 
#   finally the power for dependence of y on x

# "range1" (-4.61, 4.61) 0.01 - 0.99 
# "range2" (-4.6, .69) 0.01 - 0.66
# "range3" (-.69, 4.6) 0.33 - 0.99



# R-number of replications, N-sample size
R=500; N=500

simresult=as.data.frame(matrix(nr=45, nc=10))
ii=0
for (x in c("range1","range2","range3"))
for (z in c("range1","range2","range3"))
for (Sp in c(1, .999, .9, .8, .7)){
	if(x=="range1") xrange=c(-4.61, 4.61)
	if(x=="range2") xrange=c(-4.6, .69)
	if(x=="range3") xrange=c(-.69, 4.6)
									   
									   
	if(z=="range1") zrange=c(-4.61, 4.61)
	if(z=="range2") zrange=c(-4.6, .69)
	if(z=="range3") zrange=c(-.69, 4.6)
	ii=ii+1
	simresult[ii,]=R.datasets.norm.x(R,N,xrange,zrange,Sp)
}
colnames(simresult)=c("type", "trials", "R", "N", "xrange", "zrange", "Sp", "pow.Sp", "pow.Se", "pow.y")

simresult <- within(simresult, {
  Se.range <- factor(zrange, labels=c('whole','low','high'))
  Y.range <- factor(xrange, labels=c('whole','low','high'))
})
simresult <- within(simresult, {
  f.Sp <- as.factor(Sp)
})
simresult <- within(simresult,{converg.prop=round((R-1)/(trials-1),3)})

save(simresult,file=paste("results_power_",N,".rdata",sep=""))

library(car)
round(Tapply(converg.prop ~ Y.range + Se.range + f.Sp, mean, na.action=na.omit, data=simresult),3)
round(Tapply(pow.Sp ~ Y.range + Se.range + f.Sp, mean, na.action=na.omit, data=simresult),3)
round(Tapply(pow.Se ~ Y.range + Se.range + f.Sp, mean, na.action=na.omit, data=simresult),3)
round(Tapply(pow.y ~ Y.range + Se.range + f.Sp, mean, na.action=na.omit, data=simresult),3)
