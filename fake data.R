rm(list=ls(all=TRUE))
set.seed(1)

n=500
p=30
xmat=matrix(rnorm(p*n),n,p)
colnames(xmat)=paste('cov',1:p,sep='')

beta=rnorm(p,mean=0,sd=0.3)
ind=sample(1:p,size=p*0.9)
beta[ind]=0
beta.true=beta

med=xmat%*%beta
true.z=z=rnorm(n,mean=med,sd=1)

nclass=20
true.breaks=breaks=sort(rnorm(nclass-1,mean=0,sd=1))

y=rep(NA,n)
cond=z<breaks[1]; y[cond]=1
for (i in 2:(nclass-1)){
  cond=z<breaks[i] & z>breaks[i-1]
  y[cond]=i
}
cond=z>breaks[nclass-1]
y[cond]=nclass

dat=cbind(y,xmat)

setwd('U:\\modeling abundance\\git_MN_modelsel')
write.csv(dat,'fake data.csv',row.names=F)