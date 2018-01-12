tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  #qnorm can give some imprecise results
  cond=z<lo;    z[cond] = lo[cond]
  cond=z==-Inf; z[cond] = lo[cond]
  cond=z>hi;    z[cond] = hi[cond]
  cond=z==Inf;  z[cond] = hi[cond]
  z
}
#----------------------------------------------------------------------------------------------
sample.z=function(y,cov,betas,b,class1,nobs,nclass){
  z=rep(NA,nobs)
  media=cov%*%betas
  
  #------------------------------
  #get z
  cond=y==class1[1]
  z[cond]=tnorm(sum(cond),lo=-100,hi=b[1],mu=media[cond],sig=1)
  
  for (i in 2:(nclass-1)){
    cond=y==class1[i]
    z[cond]=tnorm(sum(cond),lo=b[i-1],hi=b[i],mu=media[cond],sig=1)
  }

  cond=y==class1[nclass]
  z[cond]=tnorm(sum(cond),lo=b[nclass-1],hi=100,mu=media[cond],sig=1)
  z
}
#----------------------------------------------------------------------------------------------
sample.betas=function(z,cov,prior.var){
  tmp=dim(cov)
  if (is.null(tmp)){
    p=1
    prec=matrix(t(cov)%*%cov+(1/prior.var),1,1)
  }
  if (!is.null(tmp)){
    p=ncol(cov)
    prec=t(cov)%*%cov+diag(1/prior.var,p)
  }
  var1=solve(prec)
  pmedia=var1%*%t(cov)%*%z
  t(rmvnorm(1,pmedia,var1))
}
#----------------------------------------------------------------------------------------------
sample.b=function(z,y,nclass,class1){
  b=rep(NA,nclass-1)
  for (i in 2:nclass){
    cond1=y==class1[i-1]
    cond2=y==class1[i]
    lo=max(z[cond1])
    hi=min(z[cond2]);
    # if (lo==hi) hi=hi+0.001
    b[i-1]=runif(1,lo,hi)
  }
  b
}
#------------------------------
log.marg.likel=function(cov,z,prior.var){
  tmp=dim(cov)
  if (is.null(tmp)){
    p=1
    prec=matrix(t(cov)%*%cov+(1/prior.var),1,1)
  }
  if (!is.null(tmp)){
    p=ncol(cov)
    prec=t(cov)%*%cov+diag(1/prior.var,p)
  }
  var1=solve(prec)
  mu=var1%*%t(cov)%*%z
  
  -(1/2)*(-t(mu)%*%prec%*%mu+t(z)%*%z+p*log(prior.var)-determinant(var1)$modulus[1])
}
#------------------------------------
samp.move=function(indin,indout,maxp,z,xmat.orig,prior.var){
  indin.old=indin
  p=length(indin.old)
  rand1=runif(1)	
  p0=1
  if (p == 1) {
    if (rand1 < 1/2) indin.new=swap(indin,indout)
    if (rand1 > 1/2) indin.new=birth(indin,indout)
    p0=2/3 #death prob 2 -> 1 is (1/3) and birth prob 1 -> 2 (or swap prob 1 -> 1) is 1/2. 
  }
  if (p == maxp) {
      indin.new=death(indin)
      p0=1/3 #birth prob T-1 -> T is (1/3) and death prob T -> T-1 is 1
  }
  if (1 < p & p < maxp) {
    if (rand1 < 1/3) {
      indin.new=birth(indin,indout)
      if (p==maxp-1) p0=3 #death prob from T -> T-1 is 1 and birth prob from T-1 -> T is (1/3)
    }
    if (1/3 < rand1 & rand1 < 2/3) {
      indin.new=death(indin)
      if (p==2) p0=3/2 #birth prob from 1 -> 2 is 1/2 and death prob from 2 -> 1 is 1/3
    }
    if (2/3 < rand1) indin.new=swap(indin,indout)
  }
  pold=log.marg.likel(cov=xmat.orig[,indin.old],z=z,prior.var=prior.var)
  pnew=log.marg.likel(cov=xmat.orig[,indin.new],z=z,prior.var=prior.var)+log(p0)
  prob=exp(pnew-pold)
  rand2=runif(1)
  
  seq1=1:maxp
  k=which(!seq1%in%indin.new)
  indout.new=seq1[k]
  if (rand2<prob) return(list(cov=xmat.orig[,indin.new],indin=indin.new,indout=indout.new))
  return(list(cov=xmat.orig[,indin.old],indin=indin.old,indout=indout))
}
#------------------------------------------
death=function(indinz){
  k=sample(1:length(indinz),size=1) 
  indinz[-k]
}
#---------------------------------------------------------------------------------------------------
swap=function(indinz,indoutz){
  n=length(indinz)
  if (n==1) tmp=numeric()
  if (n >1) {
    k=sample(1:n,size=1)  
    tmp=indinz[-k]
  }
  
  n=length(indoutz)
  if (n==1) include=indoutz
  if (n> 1) include=sample(indoutz,size=1)
  sort(c(tmp,include))
}
#---------------------------------------------------------------------------------------------------
birth=function(indinz,indoutz){
  n=length(indoutz)
  if (n==1) k=indoutz
  if (n>1 ) k=sample(indoutz,size=1)
  sort(c(indinz,k))
}
