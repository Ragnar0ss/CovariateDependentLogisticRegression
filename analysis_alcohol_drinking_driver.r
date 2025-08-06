# load d3_age_recoded.rdata, estim09.r, and analyse.r 

# During the past 30 days, what is the largest number of alcoholic drinks you had in a row, that is, within a couple of hours?
# A. I did not drink alcohol during the past 30 days
# B. 1 or 2 drinks,  C. 3 drinks ...

d3$Age=d3$age1

d3$alc.qty.43 = d3[,43]-1; d3$alc.qty.43[d3$alc.qty.43>1]=1

d3$Alcohol.consumption = d3$alc.qty.43 

dat=data.frame(d3["Alcohol.consumption"],d3["Age"],Rode.w.drinking.driver=d3[,9])
data.complete=dat[complete.cases(dat),]

allresults=analyse(data.complete[,1],data.complete[,2],data.complete[,3],200)

print(round(data.frame(allresults[1]),4))
print(round(data.frame(allresults[2],allresults[3]),3))

