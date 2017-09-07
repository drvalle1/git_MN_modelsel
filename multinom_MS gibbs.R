multinom.gibbs=function(dat,ngibbs,covs,burnin,prior.var){
  nobs=nrow(dat)
  
  tmp=sort(unique(dat$y))
  nclass=length(tmp)
  class1=1:nclass
  uni=data.frame(y=tmp,y1=class1)
  
  dat1=merge(dat,uni,all=T); dim(dat); dim(dat1)
  dat2=dat1[,c('y1',covs)]
  
  #get initial values
  z=(dat2$y1-mean(dat2$y1))/sd(dat2$y1)
  b=tmp[1:(nclass-1)]+0.1
  b=(b-mean(dat2$y1))/sd(dat2$y1)
  xmat.orig=cov=data.matrix(dat2[,covs])
  maxp=ncol(cov)
  
  betas=rep(0,maxp)
  indin=1:4
  indout=5:maxp
  betas=betas[indin]
  cov=cov[,indin]

  #gibbs sampler
  store.b=matrix(NA,ngibbs,nclass-1)
  store.betas=matrix(NA,ngibbs,maxp)
  options(warn=2)
  for (i in 1:ngibbs){
    print(c(i,indin))
    
    z=sample.z(y=dat2$y1,cov=cov,betas=betas,b=b,class1=class1,nobs=nobs,nclass=nclass)
    b=sample.b(z=z,y=dat2$y1,nclass=nclass,class1=class1)

    tmp=samp.move(indin,indout,maxp,z,xmat.orig,prior.var)
    indin=tmp$indin
    indout=tmp$indout
    cov=tmp$cov
    
    betas=sample.betas(z=z,cov=cov,prior.var=prior.var)
    
    #store results
    store.b[i,]=b
    
    tmp=rep(0,maxp)
    tmp[indin]=betas
    store.betas[i,]=tmp
  }
  
  after.burn=burnin:ngibbs
  list(b=store.b[after.burn,],
       betas=store.betas[after.burn,])
}

