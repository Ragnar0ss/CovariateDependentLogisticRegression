# load CzechRep.rdata, estim09.r, and analyse.r 

table(dd$v172)
dd$vote=ifelse(as.numeric(dd$v172)==1,0,1) 
# v172: always, never, usually, thus vote: 0-always  1-not always
# we expect that this is negatively correlated with v142 (How important is it for you
#   to live in a country that is governed democratically? 1-not at all ... 10-absolutely)
# the more important is democracy, the more probable is going always

# covariate dependent: Se = P(says not always | not always)
# we expect that this is positively correlated with v152 (Can it be justified
# that someone accepts a bribe? 1-never ... 10-always)
# the more negative opinion about bribery, the less honest reporting of "not always"
# the less v152, the smaller Se

# constant: Sp = P(says always | always votes)


table(dd$vote)
table(dd$v142)
table(dd$v152)

data=data.frame(y=dd$vote, x=dd$v142, z=dd$v152)
data.complete=data[complete.cases(data),]

allresults=with(data.complete,analyse(y, x, z, R=200, show.progress=T))

print(round(data.frame(allresults[1]),4))

