# load estim09_M5_only.r

R.datasets.norm.x=function(R,N,xrange,zrange,Sp){

	proby=function(x) exp(x)/(exp(x)+1)
	
	sens=function(x) exp(x)/(exp(x)+1)
	
	results=matrix(nr=R,nc=6)
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
		param.estim=estim09.M5(Dat[,1],Dat[,2],Dat[,3])
		if(param.estim[7] == 0) {
			results[k,]=param.estim[1:6]
			results[k,5]=exp(results[k,5])/(exp(results[k,5])+1)
			k=k+1
		}
		trials=trials+1
	}
	trials=trials-1
	#print(results)
	#cat("-----------\n")
	means=apply(results,2,mean)
	sds=apply(results,2,sd)
	cat(trials," ",R," ",N," ",xrange[1]," ",zrange[1]," ",Sp," ",
		round(means,2)," ",round(sds,2),"\n")
	
	return(data.frame("1q",trials,R,N,xrange[1],zrange[1],Sp,means[1],means[2],means[3],means[4],means[5],means[6],sds[1],sds[2],sds[3],sds[4],sds[5],sds[6]))
}

# values of predictors are normally distributed

# "range1" (-4.61, 4.61) 0.01 - 0.99 
# "range2" (-4.6, .69) 0.01 - 0.66
# "range3" (-.69, 4.6) 0.33 - 0.99



# R-number of replications, N-sample size
R=500; N=500

simresult=as.data.frame(matrix(nr=20, nc=19))
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
colnames(simresult)=c("type", "trials", "R", "N", "xrange", "zrange", "Sp","beta0","beta1","gamma0","gamma1","estSp","dev","SEb0","SEb1","SEg0","SEg1","SESp","SEdev")

round(simresult[,2:19],2)

simresult <- within(simresult, {
  Se.range <- factor(zrange, labels=c('whole','low','high'))
  Y.range <- factor(xrange, labels=c('whole','low','high'))
})
simresult <- within(simresult, {
  f.Sp <- as.factor(Sp)
})
simresult <- within(simresult,{converg.prop=round((R-1)/(trials-1),3)})

rmse = function(meanx, sdx, truex) ((meanx-truex)^2 + sdx^2)^.5

simresult$RMSEb0=with(simresult,rmse(beta0 ,SEb0,0))
simresult$RMSEb1=with(simresult,rmse(beta1 ,SEb1,1))
simresult$RMSEg0=with(simresult,rmse(gamma0,SEg0,0))
simresult$RMSEg1=with(simresult,rmse(gamma1,SEg1,1))
simresult$RMSESp=with(simresult,rmse(estSp ,SESp,Sp))

save(simresult,file=paste("results_bias_",N,".rdata",sep=""))

cbind(round(simresult[,2:19],2),simresult[,20:22],round(simresult[,23:28],2))

library(car)
round(Tapply(converg.prop ~ Y.range + Se.range + f.Sp, mean, na.action=na.omit, data=simresult),2)
round(Tapply(RMSEb0 ~ Y.range + Se.range + f.Sp, mean, na.action=na.omit, data=simresult),2)
round(Tapply(RMSEb1 ~ Y.range + Se.range + f.Sp, mean, na.action=na.omit, data=simresult),2)
round(Tapply(RMSEg0 ~ Y.range + Se.range + f.Sp, mean, na.action=na.omit, data=simresult),2)
round(Tapply(RMSEg1 ~ Y.range + Se.range + f.Sp, mean, na.action=na.omit, data=simresult),2)
round(Tapply(RMSESp ~ Y.range + Se.range + f.Sp, mean, na.action=na.omit, data=simresult),2)


