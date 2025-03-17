###################################################
# Generalized Tweedie SBM with mixed covariates
# Jul 2024
# input:
# Y is a list with length T
# X is a list with the first p1 matrices for the time-invariant
# coefficients, and the later p2 matrices for the time-varying 
# coefficients.
#
# Parameters/output:
# beta is a list with the first v1 scales as the coefficients 
# for the time-invariant covariate, and the later v2 vectors of 
# length T as the coefficients for the time-varying covariates.
#
###################################################

SBM_Tweedie <- function(Y,X,K,p1,
                        initialCommunity="random",
                        tau,
                        beta0=matrix(0,K,K),
                        lambda=0.5,
                        rho=1.5,
                        phi=NULL,
                        saddle=FALSE,
                        max.iter=50,eps=1e-6,
                        optim.iter=200,
                        fp.iter=50,fp.eps=1e-5,
                        step1.max.iter=20,step1.eps=1e-5)
{
  start.time=Sys.time()
  
  # initial community labels
  n=ncol(Y[[1]])
  if(initialCommunity=="one"){
    tau[,1]=1
  }
  
  ########### Step 1: estimate beta and beta(t) via backfitting ############
  
  # initial values of beta: the first p1 entries are scalars with the value 0.5; 
  #                         the next p2 entries are zero vectors of length T
  newbeta=c(as.list(rep(0.5, p1)), replicate(length(X)-p1, vector("numeric", length(Y)), simplify = FALSE))
  
  # B-spline Basis Function Values
  B=bsplineS(seq(0,1,length=(length(Y)+2))[2:(length(Y)+1)], breaks=seq(0,1,length=(length(Y)+2)), norder=4, nderiv=0, returnMatrix=FALSE)
  # B-Spline Penalty Matrix
  Omega=bsplinepen(create.bspline.basis(c(0,1),(length(Y)+4)))
  
  i=0;convg<-F;
  # loop to estimate beta
  while (i < step1.max.iter && !convg)
  {
    beta=newbeta
    # fixed-effect covariate
    newbeta[1:p1]=as.list(fixed.beta.solver.2steps(Y,X,newbeta,calculateZZ(tau),K,rho,p1,max.ite=200))
    
    # time-varying covariate
    for (j in (p1+1):length(X)) {
      theta=tv.beta.solver.2step(Y,X,calculateZZ(tau),rho,lambda,B,Omega,K,newbeta,j,p1,max.ite=optim.iter)
      newbeta[[j]]=B%*%matrix(theta,,1)
    }
    
    convg <- (sum(mapply(function(x, y) if (is.numeric(x) && length(x) == 1) (x - y)^2 else sum((x - y)^2), beta, newbeta))
              <step1.eps*length(Y)*length(X))
    i <- i + 1
    
  }
  
  ########### Step 2: estimate other parameters ############
  i=0;convg<-F;newbeta0=beta0;newtau=tau;newPi=Pi=colMeans(tau);
  beta=newbeta
  
  est.phi=!is.numeric(phi)  # whether to estimate phi
  if(!est.phi){newPhi=phi}else{newPhi=phi=0.5} # start from a low phi values
  
  # early diagnosis of cluster merging 
  early.termination=0
  
  while (i < max.iter && !convg)
  {
    beta0=newbeta0;tau=newtau;Pi=newPi;phi=newPhi
    newPi=colMeans(tau)
    
    # calculate tau_{ik}tau_{jl}
    ZZ=calculateZZ(tau)
    if (min(unlist(lapply(ZZ, max)))<1e-15) {
      early.termination=early.termination+1
      if(early.termination>3){
        break
      }
    }
    
    Xb=xbeta(X,newbeta,p1,length(X)-p1)
    
    # update beta0 in one line
    beta0.top=Reduce("+",lapply(1:length(Y), function(x){exp((1-rho)*(Xb$Xbeta+Xb$Xbetat[[x]]))*Y[[x]]}));#diag(beta0.top)=0
    beta0.bottom=Reduce("+",lapply(1:length(Y), function(x){exp((2-rho)*(Xb$Xbeta+Xb$Xbetat[[x]]))}));#diag(beta0.bottom)=0
    ind=0
    for (ii in 1:K) {
      for (jj in ii:K) {
        ind=ind+1
        newbeta0[ii,jj]=log(sum(ZZ[[ind]]*beta0.top)/sum(ZZ[[ind]]*beta0.bottom)) # 26 Jan
        if(ii!=jj){newbeta0[jj,ii]=newbeta0[ii,jj]}
      }
    }
    
    # update tau with fixed point algorithm
    fp.i=0;fp.convg<-F;newtau.fp=tau
    while (fp.i < fp.iter && !fp.convg) {
      fp.i=fp.i+1;tau.inFP=newtau.fp
      newtau.fp=do.call(rbind, lapply(1:n,function(x){
        tauUpdate.tw(Xb,newPi,newbeta0,tau.inFP,Y,x,rho,newPhi)
      }))
      fp.convg<- (norm(newtau.fp-tau.inFP,"F")<fp.eps)
    }
    newtau=newtau.fp
    
    # estimate phi
    if(est.phi){
      newPhi=phi.solver(Y,newbeta0,Xb,newtau,rho,phi,max.ite=optim.iter,saddle=saddle)
    }
    
    convg <- (norm(beta0-newbeta0,"2")<(K*K*eps))&&(abs(phi-newPhi)<eps)
    i <- i + 1
    
  }
  end.time=Sys.time()
  return(list(beta=newbeta,beta0=newbeta0,newtau=tau,
              community=sapply(1:n, function(x){which.max(newtau[x,])}),
              ite=i,time=end.time-start.time,newrho=rho,newphi=newPhi))
}

# update beta0
beta0.solver <- function(Y,X,beta0,beta,ZZ,rho,phi,max.ite=200,beta0range){
  dat = list(Y=Y,beta=beta,X=X,ZZ=ZZ,rho=rho,phi=phi)
  return(optim(par=0, min.NLKLH.beta0, gr = NULL,
               method ="L-BFGS-B",
               x = dat,
               lower = -beta0range, upper = beta0range,
               control = list(maxit=max.ite), hessian = FALSE)$par)
}
min.NLKLH.beta0 <- function(x, par) {
  tvSum=0
  for(t in 1:length(x$Y)){
    betaXbeta=par+(x$X)*(x$beta[t])
    tvSum=tvSum+x$Y[[t]]*exp((1-x$rho)*betaXbeta)/(1-x$rho)-exp((2-x$rho)*betaXbeta)/(2-x$rho)
  }
  val=x$ZZ*tvSum
  val[lower.tri(val,diag = T)]=0
  return(-sum(val)/(x$phi))
}

beta0.gr <- function(x, par) {
  tvSum=0
  for(t in 1:length(x$Y)){
    betaXbeta=par+(x$X)*(x$beta[t])
    tvSum=tvSum+x$Y[[t]]*exp((1-x$rho)*betaXbeta)-exp((2-x$rho)*betaXbeta)
  }
  val=x$ZZ*tvSum
  val[lower.tri(val,diag = T)]=0
  return(-sum(val)/(x$phi))
}

##########function to calculate Xbeta and Xbeta(t)###########
xbeta <- function(X,beta,p1,p2){
  
  if(p1==0){
    Xbeta=0
  }else{
    Xbeta=Reduce("+",lapply(1:p1, function(i){X[[i]]*beta[[i]]}))
  }
  
    Xbetat=lapply(1:length(beta[[p1+p2]]), function(t){
      Reduce("+",lapply((p1+1):(p1+p2),function(i) {X[[i]]*beta[[i]][t]}))
    })

  return(list(Xbeta=Xbeta,Xbetat=Xbetat))
}

######################## update beta ########################
# fixed-effect covariate
fixed.beta.solver.2steps <- function(Y,X,beta,ZZ,K,rho,p1,max.ite=200){
  dat = list(Y=Y,X=X,ZZ=ZZ,rho=rho,K=K,beta=beta,p1=p1)
  return(optim(par=unlist(beta[1:p1]), min.NLKLH.fixed.beta, gr = NULL,
               method ="L-BFGS-B",
               x = dat,
               lower = -Inf, upper = Inf,
               control = list(maxit=max.ite), hessian = FALSE)$par)
}
min.NLKLH.fixed.beta <- function(x, par) {
  K=x$K
  ind = 0; val=0
  
  b=x$beta
  b[1:x$p1]=as.list(par)
  Xb=xbeta(x$X,b,x$p1,length(x$X)-x$p1)
  
  mat1=Reduce("+",lapply(1:length(x$Y), function(t){
    (x$Y[[t]])*exp((1-x$rho)*( Xb$Xbeta+Xb$Xbetat[[t]]) )
  }))
  diag(mat1)=0
  
  mat2=Reduce("+",lapply(1:length(x$Y), function(t){
    exp((2-x$rho)*( Xb$Xbeta+Xb$Xbetat[[t]]))
  }))
  diag(mat2)=0
  
  for (ii in 1:K) {
    for (jj in ii:K) {
      ind=ind+1
      
      val=val+ifelse(abs(sum(x$ZZ[[ind]]*mat2))<1e-10,0,
                     (sum(x$ZZ[[ind]]*mat1)^(2-x$rho))/(sum(x$ZZ[[ind]]*mat2)^(1-x$rho)))
    }
  }
  return(val)
}

# time-varying covariate effect
tv.beta.solver.2step <- function(Y,X,ZZ,rho,lambda,B,Omega,K,beta,j,p1,max.ite=200){
  dat = list(Y=Y,X=X,ZZ=ZZ,rho=rho,lambda=lambda,B=B,Omega=Omega,K=K,beta=beta,j=j,p1=p1)
  return(optim(par=numeric(length(Y)+4), min.NLKLH.tv.beta, gr = NULL, 
               method ="L-BFGS-B",
               x = dat,
               lower = -Inf, upper = Inf,
               control = list(maxit=max.ite), hessian = FALSE)$par)
}
min.NLKLH.tv.beta <- function(x, par) {
  K=x$K
  ind = 0; mat=0
  
  #########
  beta=x$beta
  beta[[x$j]]=(x$B)%*%matrix(par,,1)
  Xb=xbeta(x$X,beta,x$p1,length(x$X)-x$p1)
  
  mat1=Reduce("+",lapply(1:length(x$Y), function(t){
    (x$Y[[t]])*exp((1-x$rho)*( Xb$Xbeta+Xb$Xbetat[[t]]) )
  }))
  diag(mat1)=0
  
  mat2=Reduce("+",lapply(1:length(x$Y), function(t){
    exp((2-x$rho)*( Xb$Xbeta+Xb$Xbetat[[t]]))
  }))
  diag(mat2)=0
  #########
  
  for (ii in 1:K) {
    for (jj in ii:K) {
      ind=ind+1
      val=ifelse((sum(x$ZZ[[ind]]*mat2))<1e-10,0,(sum(x$ZZ[[ind]]*mat1)^(2-x$rho))/(sum(x$ZZ[[ind]]*mat2)^(1-x$rho)))
      mat=mat+val
    }
  }
  return(-mat/((1-x$rho)*(2-x$rho))+(x$lambda)*matrix(par,1,)%*%(x$Omega)%*%matrix(par,,1))
}

# update phi
phi.solver <- function(Y,beta0,Xb,tau,rho,phi,max.ite=200,saddle=F){
  community=sapply(1:ncol(Y[[1]]), function(x){which.max(tau[x,])})
  dat = list(Y=Y,Xb=Xb,beta0=beta0,community=community,rho=rho,saddle=saddle,p1=p1)
  return(optim(par=phi, est.phi.NLKLH, gr = NULL,
               method ="L-BFGS-B",
               x = dat,
               lower = 0.1, upper = 100, # changed (lower,upper) from (0,inf) to (0.1,1000); 10 Jan
               control = list(maxit=max.ite), hessian = FALSE)$par)
}
est.phi.NLKLH <- function(x, par) {
  Xb=x$Xb
  K=ncol(x$beta0)
  val=0
  for (tt in 1:length(Y)) {
    
    mu=exp((Xb$Xbeta+Xb$Xbetat[[tt]])+x$beta0[x$community,x$community])
    mu=mu[upper.tri(mu,diag = F)]
    y=x$Y[[tt]][upper.tri(x$Y[[tt]],diag = F)]
    if(!(x$saddle)){
      val=val-sum(log(tweedie::dtweedie( y=y, power=x$rho, mu=mu, phi=par)))
    }else{
      val=val-sum(log(tweedie::dtweedie.saddle( y=y, power=x$rho, mu=mu, phi=par)))
    }
    
  }
  
  return(val)
}

# update tau for the node i
tauUpdate.tw <- function(Xb,Pi,beta0,tau,Y,i,rho,phi){
  K=ncol(tau);n=ncol(Y[[1]]);numT=length(Y)
  TAUi=numeric(K)
  for (q in 1:K) {
    a=0
    for (tt in 1:numT) {
      
      for (j in 1:n) {
        if(j!=i){
          for (l in 1:K) {
            if(tau[j,l]!=0){
              b=(ifelse(Y[[tt]][i,j]!=0,
                        Y[[tt]][i,j]*(exp((1-rho)*(beta0[q,l]+Xb$Xbeta[i,j]+Xb$Xbetat[[tt]][i,j] ))/(1-rho) ),0)-
                   exp((2-rho)*(beta0[q,l]+Xb$Xbeta[i,j]+Xb$Xbetat[[tt]][i,j]))/(2-rho))*(tau[j,l])/phi
              a=a+b
            }
          }
        }
      }
      
    }
    TAUi[q]=a
  }
  TAUi=TAUi+log(Pi)
  
  if(Inf%in%abs(TAUi)){ 
    return(tau[i,])
  } else {
    # take the largest value as the "benchmark"
    val=numeric(length(TAUi))
    val[which.max(TAUi)]=1
    for (i in 1:length(TAUi)) {
      if(i!=which.max(TAUi)){
        val[i]=exp(TAUi[i]-TAUi[which.max(TAUi)])
      }
    }
    return(val/ifelse(sum(val)!=0,sum(val),1))
  }
}

# calculate ZZ
calculateZZ <- function(TAU){
  n=nrow(TAU)
  K=ncol(TAU)
  ZZ=vector("list",K*(K+1)/2);ind=0
  for (ii in 1:K) {
    for (jj in ii:K) {
      ind=ind+1
      z=matrix(0,n,n)
      for (i0 in 1:n) {
        for (j0 in 1:n) {
          z[i0,j0]=TAU[i0,ii]*TAU[j0,jj]
        }
      }
      ZZ[[ind]]=z
    }
  }
  return(ZZ)
}