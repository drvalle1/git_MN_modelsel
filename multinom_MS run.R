rm(list=ls(all=TRUE))
library('mvtnorm')
set.seed(1)

setwd('U:\\modeling abundance\\git_MN_modelsel')
source('multinom_MS functions.R')
source('multinom_MS gibbs.R')
dat=read.csv('fake data.csv',as.is=T)

ngibbs=10000
covs=colnames(dat)[-1]
prior.var=1

res=multinom.gibbs(dat=dat,ngibbs=ngibbs,covs=covs,burnin=5000,prior.var=prior.var)

ebreaks=apply(res$b,2,mean)
rango=range(c(true.breaks,ebreaks))
plot(true.breaks,ebreaks,xlim=rango,ylim=rango)
lines(rango,rango)

betas.estim=colMeans(res$betas)
rango=range(c(beta.true,betas.estim))
plot(beta.true,betas.estim,ylim=rango,xlim=rango)
lines(rango,rango)