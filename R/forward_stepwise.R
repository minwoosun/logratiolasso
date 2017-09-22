custom_fs <- function(expanded_set, resp, k_max, selected_vars, p) {
  #wrapper of myfs function for internal use

  if(k_max <= 0) {
    return(0)
  }
  subs <- myfs(expanded_set, resp, nsteps = k_max, center = FALSE)

  list(theta_vals = subs$beta, theta_ind = subs$ind, subset = which(selected_vars))
}


myfs <- function(x,y,nsteps=min(nrow(x),ncol(x)),center=T, family="gaussian"){
  # fs by minimizing scaled ip;
  # for use in logRatioLasso; assumes predictors correspond to all possible pairs
  # of some raw variables
  # first center x and y
  # note at this point, it always uses Gaussian model for the stepwise calculation.
  #  but with family="binomial", it returns logistic reg coefs
  #returns:
  #  pred,s =predictors and their signs, in order entered 
  #  scor, bhat: scaled ip and ls coef for each predictor entered
  # sigmahat:  est of error variance;
  # rss of each model
  #pss- %var unexplained by each model
  # ind- indices of variable pairs in order entered

  p=ncol(x)
  n=length(y)
  y.orig=y
  if(center){
    x=scale(x,T,F)
    y=y-mean(y)
  }

  nv=.5*(1+sqrt(1+4*2*p));
  
  #construct matrix indices
  thmat=matrix(NA,nv,nv)
  suppressWarnings(thmat[row(thmat)>col(thmat)] <- 1:p)
  thmat=t(thmat)
  yhat=matrix(NA,n,nsteps)
  pred=s=scor=bhat=sigmahat=rss=rep(NA,nsteps)
  ip=t(x)%*%y/sqrt(diag(t(x)%*%x))
  pred[1]=which.max(abs(ip))
  s[1]=sign(sum(x[,pred[1]]*y))
  scor[1]=ip[pred[1]]
  bhat[1]=ip[pred[1]]/sqrt(sum(x[,pred[1]]^2))

  r=lsfit(x[,pred[1]],y)$res
  rss[1]=sum(r^2)
  sigmahat[1]= sqrt(sum(r^2)/(n-1))
  ind=which(thmat ==pred[1], arr.ind = T)

  if(nsteps>1){
    for(j in 2:nsteps){
      mod=pred[1:(j-1)]
      r= lsfit(x[,mod],r,int=center)$res
      sigmahat[j]= sqrt(sum(r^2)/n)
      xr= lsfit(x[,mod],x,int=center)$res
      ip=t(xr)%*%r/sqrt(diag(t(xr)%*%xr))
      ip[mod]=0
      pred[j]=which.max(abs(ip))
      scor[j]=ip[pred[j]]
      s[j]=sign(sum(xr[,pred[j]]*r))
      bhat[j]=ip[pred[j]]/sqrt(sum(xr[,pred[j]]^2))
      rss[j]=sum( lsfit(x[,pred[1:j]],y)$res^2)

      new=which(thmat ==pred[j], arr.ind = T)
      ind=rbind(ind,new)
    }
  }

  #compute ls coefs for all models
  beta=vector("list",nsteps)
  b0=rep(NA,nsteps)
  for(j in 1:nsteps){
    if(family=="gaussian"){
      junk=lsfit(x[,pred[1:j],drop=F],y)
      beta[[j]]=junk$coef[-1]
      b0[j]=junk$coef[j]
    }
    
    #  if(family=="binomial") junk=glm(y.orig~x[,pred[1:j],drop=F],family="binomial")
    if(family=="binomial") {
      junk=glmnet(x[,pred[1:j],drop=F],y.orig,alpha=.05,standardize=FALSE,family="binomial")
      nlam=ncol(junk$beta)
      beta[[j]]=junk$beta[,nlam]
      if(sum(is.na(beta[[j]]))>0) {browser()}
      b0[j]=junk$a0[nlam]
    }
      
  }
  prss=rss/sum( (y-mean(y))^2)
  return(list(pred=pred,ind=ind,s=s,scor=scor,b0=b0,beta=beta,sigmahat=sigmahat,rss=rss,prss=prss))
}

myfs.new=function(x,y,nsteps=min(nrow(x),ncol(x)),center=T, family="gaussian",verbose=FALSE){
p=ncol(x)
n=length(y)
# fs by minimizing scaled ip;
# for use in logRatioLasso; assumes predictors correspond to all possible pairs
# of some raw variables
# first center x and y
# note at this point, it always uses Gaussian model for the stepwise calculation.
#  but with family="binomial", it returns logistic reg coefs
#returns:
#  pred,s =predictors and their signs, in order entered 
#  scor, bhat: scaled ip and ls coef for each predictor entered
# sigmahat:  est of error variance;
# rss of each model
#pss- %var unexplained by each model
# ind- indices of variable pairs in order entered

y.orig=y
if(center){
x=scale(x,T,F)
y=y-mean(y)
}

out=fs(x, y, maxsteps =nsteps, intercept = FALSE, 
    normalize = FALSE, verbose = verbose)



pred=out$act

nv=.5*(1+sqrt(1+4*2*p));
#construct matrix indices
thmat=matrix(NA,nv,nv)
suppressWarnings(thmat[row(thmat)>col(thmat)] <- 1:p)
thmat=t(thmat)
ind=matrix(NA,length(pred),2)
for(j in 1:nrow(ind)){
  ind[j,]=which(thmat ==pred[j], arr.ind = T)

}


#compute ls coefs for all models
nsteps.act=min(nsteps,length(pred))
beta=vector("list",nsteps.act)
b0=rep(NA,nsteps.act)
for(j in 1:nsteps.act){

if(family=="gaussian"){
        junk=lsfit(x[,pred[1:j],drop=F],y)
        beta[[j]]=junk$coef[-1]
          b0[j]=junk$coef[j]
}
   #  if(family=="binomial") junk=glm(y.orig~x[,pred[1:j],drop=F],family="binomial")
    if(family=="binomial") {
        if(j==1) {junk=glm(y.orig~x[,pred[1:j],drop=F],family="binomial")
                   beta[[j]]=junk$coef[-1]
                      b0[j]=junk$coef[j]
              }
         if(j>1){
              junk=glmnet(x[,pred[1:j],drop=F],y.orig,alpha=.01,standardize=FALSE,family="binomial")
        nlam=ncol(junk$beta)
       beta[[j]]=junk$beta[,nlam]
        # if(sum(is.na(beta[[j]]))>0) {browser()}
          b0[j]=junk$a0[nlam]
 }}}

return(list(pred=pred,ind=ind,b0=b0,beta=beta))
}