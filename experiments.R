################################################################################
# Title: Algorithms for Determination of Sample Sizes for Bayesian Estimations
#        in Single-Server Markovian Queues
#
# Reference:
# Gomes, E. S., Cruz, F. R. B., & Singh, S. K. (2023).
# Algorithms for Determination of Sample Sizes for Bayesian Estimations
# in Single-Server Markovian Queues.
# American Journal of Mathematical and Management Sciences, 42(4), 307–322.
# https://doi.org/10.1080/01966324.2023.2255316
#
# Authors:
# Eriky S. Gomes
# Frederico R. B. Cruz
#
# Affiliation:
# Department of Statistics
# Universidade Federal de Minas Gerais (UFMG)
#
# Contact:
# eriky-tn@ufmg.br
# fcruz@est.ufmg.br
#
# Copyright (c) 2023 Gomes & Cruz
# Version: v2023
#
# Description:
# R script developed for computational experiments of the article
################################################################################
# auxiliar functions
################################################################################

# maximum likelihood estimator
MLE<-function(smp){
  n<-length(smp)
  y<-sum(smp)
  rhoE<-y/n
  return(rhoE)
}

# gaussian hypergeometric function
GaussHyp<-function(a,b,c,z) {
  Aux<-function(u,a,b,c,z) {
    auxOut<-u^(b-1)*(1-u)^(c-b-1)*(1-u*z)^(-a)
    return(auxOut)}
  output<-integrate(Aux,0,1,a,b,c,z)[[1]]/beta(b,c-b)
  return(output)
}
# average of distribution (dtb) estimator
AverE<-function(smp,dtb,...){
  Aux<-function(p,smp,...){
    auxOut<-p*dtb(p,smp,...)
    return(auxOut)}
  averageE<-integrate(Aux,0.01,0.99,smp,...,subdivisions = 10000)[[1]]
  return(averageE)
}

# maximum of posterior estimator
MaxPost<-function(posterior,...){
  delta<-1e-3                             # accuracy value
  golden<-(sqrt(5)-1)/2                   # inverse of golden number
  iter<-ceiling(log(delta)/log(golden))   # number of iterations
  a<-delta                                # interval inferior limit
  b<-1-delta                              # interval superior limit
  x1<-a+(1-golden)*(b-a)
  x2<-a+golden*(b-a)
  fx1<-posterior(x1,...)
  fx2<-posterior(x2,...)
  for(i in 1:iter){
    if(fx1<fx2){
      a<-x1
      x1<-x2
      x2<-a+golden*(b-a)
      
      fx1<-fx2
      fx2<-posterior(x2,...)
    }
    else{
      b<-x2
      x2<-x1
      x1<-a+(1-golden)*(b-a)
      
      fx2<-fx1
      fx1<-posterior(x1,...)
    }
  }
  xMax<-(a+b)/2
  fxMax<-posterior(xMax,...)
  return(c(xMax,fxMax))
}

# random number function from discretized distribution (dtb)
RandND<-function(rep,dtb,...){
  randN<-numeric(rep)
  discDtb<-function(dtb,...){
    delta<-1e-3
    x<-seq(0.001,0.999,by=delta)
    fx<-numeric(length(x))
    for(i in 1:length(x)){
      fx<-dtb(x,...)
    }
    fx<-fx/sum(fx)
    discDtbOut<-list(x=x,fx=fx)
    return(discDtbOut)
  } # dicretizing distribution
  discP<-discDtb(dtb,...)
  
  # generating randow numbers
  for(i in 1:rep){
    randN[i]<-sample(discP[["x"]],1,prob=discP[["fx"]])
  }
  return(randN)
}

# random number generator by acceptance-rejection method
RandNAR<-function(n,pdf,...){
  tol<-1e-2
  randN<-numeric(n)
  fmax<-MaxPost(pdf,...)[2]
  i=0
  while(i<=n){
    y<-runif(1,min=tol,max=1-tol)
    u<-runif(1,min=tol,max=1-tol)
    if(y<pdf(u,...)/fmax){
      randN[i]=u
      i=i+1}
  }
  return(randN)
}

# posterior inverted function
PostInverted<-function(density,posterior,...){
  tol<-1e-3
  xMax<-MaxPost(posterior,...)[1]
  AuxPosterior<-function(p){
    return(posterior(p,...)-density)
  }
  if(AuxPosterior(tol)*AuxPosterior(xMax)>=0) xLeft<-tol
  else xLeft<-Regula(AuxPosterior,tol,xMax)
  if(AuxPosterior(xMax)*AuxPosterior(1-tol)>=0) xRight<-1-tol
  else xRight<-Regula(AuxPosterior,xMax,1-tol)
  return(c(xLeft,xRight))
}

################################################################################
# beta functions
################################################################################

# beta prior
BPrior<-function(p,a,b){
  num<-function(u,a,b){
    numOut<-u^(a-1)*(1-u)^(b-1)
    return(numOut)
    }
  prob<-num(p,a,b)/beta(a,b)
  return(prob)
}

# beta posterior
BPost<-function(p,smp,a,b){
  n=length(smp)
  y=sum(smp)
  num<-function(u,a,b,n,y){
    numOut<-u^(y+a-1)*(1-u)^(b-1)*(1+u)^(-n-y)
    return(numOut)
    }
  prob<-num(p,a,b,n,y)/beta(y+a,b)/GaussHyp(n+y,y+a,y+a+b,-1)
  return(prob)
}

# beta marginal
BMarg<-function(smp,a,b){
  y<-sum(smp)
  n<-length(smp)
  prob<-beta(y+a,b)*GaussHyp(n+y,y+a,y+a+b,-1)
  return(prob)
}

# beta average estimator
BAver<-function(smp,a,b){
  y=sum(smp)
  n=length(smp)
  average<-exp(log(y+a)+log(GaussHyp(n+y,y+a+1,y+a+b+1,-1))
          -log(GaussHyp(n+y,y+a,y+a+b,-1))-log(y+a+b))
  return(average)
}

################################################################################
# incomplete inverted beta functions
################################################################################

# regularized incomplete beta function
RIB_0.5<-function(a,b){
  Aux<-function(p,a,b){
    auxOut<-p^(a-1)*(1+p)^(-b-a)
  }
  output<-integrate(Aux,0,1,a,b)[[1]]/beta(a,b)
  return(output)
}

# inverted beta prior
IBPrior<-function(p,a,b){
  prob<-p^(a-1)*(1+p)^(-b-a)/beta(a,b)/RIB_0.5(a,b)
  return(prob)
}

# inverted beta posterior
IBPost<-function(p,smp,a,b){
  n<-length(smp)
  y<-sum(smp)
  prob<-p^(y+a-1)*(1+p)^(-n-b-y-a)/beta(y+a,n+b)/RIB_0.5(y+a,n+b)
  return(prob)
}

# inverted beta average
IBAver<-function(smp,a,b){
  n<-length(smp)
  y<-sum(smp)
  average<-beta(y+a+1,n+b-1)*RIB_0.5(y+a+1,n+b-1)/
          beta(y+a,n+b)/RIB_0.5(y+a,n+b)
  return(average)
}

################################################################################
# jeffreys functions
################################################################################

# jeffreys prior
JPrior<-function(p){
  prob<-p^(-1/2)*(1+p)^(-1/2)
  return(prob)
}

# jeffreys posterior
JPost<-function(p,smp){
  n<-length(smp)
  y<-sum(smp)
  prob<-p^(y+1/2-1)*(1+p)^(-(n+y+1/2))/beta(y+1/2,n)/RIB_0.5(y+1/2,n)
  return(prob)
}

# jeffreys average 
JAver<-function(smp){
  n<-length(smp)
  y<-sum(smp)
  average<-(y+1/2)*RIB_0.5(y+3/2,n-1)/(n-1)/RIB_0.5(y+1/2,n)
  return(average)
}

################################################################################
# root-finding algorithms
################################################################################

# # bissection
# Bissection<-function(f,a,b,...){
#   tol<-1e-3
#   k<-0
#   while(abs(b-a)>tol){
#     fa<-f(a,...)
#     fb<-f(b,...)
#     x<-(a+b)/2
#     fx<-f(x,...)
#     ifelse(fa*fx<0,b<-x,a<-x)
#     k<-k+1
#   }
#   root<-(a+b)/2
#   return(root)
# }

# regula-falsi
Regula<-function(f,a,b,...){
  if(f(a)*f(b)>0) return(NaN)
  tol<-1e-3
  # k<-0
  error<-b-a
  x0<-a
  x1<-b
  x<-numeric(1)
  while(abs(error)>tol){
    fa<-f(a,...)
    fb<-f(b,...)
    x<-(a*fb-b*fa)/(fb-fa)
    fx<-f(x,...)
    ifelse(fa*fx<0,b<-x,a<-x)
    lambda<-(x-x1)/(x1-x0)
    error<-lambda/(lambda-1)*(x-x1)
    x0<-x1
    x1<-x
  }
  root<-x
  return(root)
}

################################################################################
# credible interval functions
################################################################################

# equal-tailed credible interval
ETI<-function(posterior,...){
  delta<-1e-3
  alpha<-0.05
  p<-seq(delta,1-delta,by=delta)
  f<-posterior(p,...)
  dArea<-0
  accumulated<-0
  i<-1
  len<-length(p)
  while(accumulated<alpha/2){
    dArea<-delta*(f[i]+f[i+1])/2
    accumulated<-accumulated+dArea
    i<-i+1
  }
  limInf<-p[i-1]
  while(accumulated<1-alpha/2 & i<len){
    dArea<-delta*(f[i]+f[i+1])/2
    accumulated<-accumulated+dArea
    i<-i+1
  }
  limSup<-p[i]
  return(c(limInf,limSup))
}

# highest posterior credible interval
HDI<-function(posterior,...){
  alpha<-0.05
  tol<-1e-3
  fMax<-MaxPost(posterior,...)[2]
  fMin<-tol
  interval<-numeric(2)
  while(fMax-fMin>tol){
    fx<-(fMax+fMin)/2
    interval<-PostInverted(fx,posterior,...)
    #cat("interval=",interval,"\n")
    coverage<-integrate(posterior,interval[1],interval[2],...)[[1]]
    ifelse(coverage>1-alpha,fMin<-fx,fMax<-fx)
  }
  # cat("coverage=",coverage,"\n")
  # cat("length=",interval[2]-interval[1],"\n")
  return(interval)
}

# HDI_length<-function(width,post,...){
#   tol<-1e-2
#   fMax<-MaxPost(post,...)[2]
#   fMin<-tol
#   interval<-numeric(2)
#   length<-numeric(1)
#   while(fMax-fMin>tol){
#     fx<-(fMax+fMin)/2
#     interval<-PostInverted(fx,post,...)
#     length<-interval[2]-interval[1]
#     ifelse(length>width,fMin<-fx,fMax<-fx)
#   }
#   return(interval)
# }

# symmetrical from mean credible interval
CIAver<-function(smp,delta,post,...){
  averageE<-AverE(smp,post,...)
  linf<-max(0,averageE-delta)
  lsup<-min(1,averageE+delta)
  cover<-integrate(post,linf,lsup,smp,...)[[1]]
  return(cover)
}

# # symmetrical from maximum credible interval
# CIMax<-function(smp,delta,post,...){
#   maxE<-MaxPost(post,smp,...)[1]
#   linf<-max(0,maxE-delta)
#   lsup<-min(1,maxE+delta)
#   cover<-integrate(post,linf,lsup,smp,...)[[1]]
#   return(cover)
# }

################################################################################
# sample size criterion functions
################################################################################

# average coverage criterion - CIAver
ACC<-function(delta,priori,post,Averagef,...){
  kappa=0.05                    # significance level
  k=500
  #mx<-numeric(k)
  rhos<-RandNAR(k,priori,...)    # random numbers function
  nMin=2
  nMax=1000
  while(nMax-nMin>1){
    An<-0
    n<-round((nMin+nMax)/2)
    if (n>501){
      return(n=Inf)
    }
    i=1
    while(i<=k){
      x<-rgeom(n,1/(1+rhos[i]))
      y<-sum(x)
      average<-Averagef(x,...)
      #mx[i]<-dMarg(x,...)
      lInf<-max(0,average-delta)
      lSup<-min(1,average+delta)
      if(!is.nan(average)){
        An<-An+integrate(post,lInf,lSup,x,...,subdivisions = 10000)[[1]]#*mx[i]
        i=i+1
      }
    }
    An<-An/k#sum(mx)              # normalization of An
    ifelse(An>=1-kappa,nMax<-n,nMin<-n)
    cat("n=",n,"\n")
  }
  n=nMax
  return(n)
}

# # average coverage criterion function - HDI
# ACC_HDI<-function(width,prior,post,...){
#   kappa=0.05
#   k=1000
#   nMin=2
#   nMax<-900
#   rho<-RandNAR(k,prior,...)
#   while(nMax-nMin>1){
#     An<-0
#     n<-round((nMax+nMin)/2,0)
#     if(n>501) return(n=Inf)
#     i=1
#     while(i<=k){
#       x<-rgeom(n,1/(1+rho[i]))
#       y<-sum(x)
#       interval<-HDI_length(width,post,x,...)
#       #cat("interval=",interval,'\n')
#       linf<-interval[1]
#       lsup<-interval[2]
#       An<-An+integrate(post,linf,lsup,x,...,subdivisions = 10000)[[1]]
#       i=i+1    
#     }
#     An<-An/k
#     ifelse(An>1-kappa,nMax<-n,nMin<-n)
#     cat("n=",n,"\n")
#   }
#   n=nMax
#   return(n)
# }

# average length criterion function
ALC<-function(width,priori,post,CIf,...){
  # works only for significance level of 0.05
  k=500
  rhos<-RandNAR(k,priori,...)
  nMin=2
  nMax=900
  while(nMax-nMin>1){
    Wn<-0
    n<-round((nMax+nMin)/2,0)
    if(n>451){
      return(Inf)
    }
    i<-1
    while(i<=k){
      x<-rgeom(n,1/(1+rhos[i]))
      y<-sum(x)
      ci<-CIf(post,x,...)     # return [lInf,lSup]
      if(!is.finite(ci[1]) | !is.finite(ci[2])){
        next
      }
      ciWidth<-ci[2]-ci[1]
      Wn<-Wn+ciWidth
      i<-i+1
    }
    Wn<-Wn/k                  # mean of credible widths
    if(Wn>width){
      nMin<-n
    } 
    else{ 
      nMax<-n
    }
    cat("n=",n,"\n")
  }
  n<-nMax
  return(n)
}

################################################################################
# monte carlo function
################################################################################

MC<-function(n,p,fEst,...){
  set.seed(24680)
  rep<-1000
  est<-numeric(rep)
  smp<-numeric(n)
  for(i in 1:rep){
    smp<-rgeom(n,1/(1+p))
    est[i]<-fEst(smp,...)
  }
  output<-list(avr=mean(est), var=var(est))
  return(output)
}

################################################################################
# table functions
################################################################################

# abacus tables ----------------------------------------------------------------

# equal-tailed and HDI table function
TabCI<-function(priori,post,size,CIf,...){
  set.seed(13579)
  rep<-1000
  tab<-matrix(nrow=length(size),ncol=2,
              dimnames=list(size,c("length","cover")))
  rho<-RandNAR(rep,priori,...)
  len<-numeric(rep)
  i=1
  while(i<=length(size)){
    cover=0
    j=1
    while(j<=rep){
      smp<-rgeom(size[i],1/(1+rho[j]))
      etEst<-CIf(post,smp,...)
      if(rho[j]>=etEst[1] & rho[j]<=etEst[2]){
        cover=cover+1
      }
      len[j]<-etEst[2]-etEst[1]
      j<-j+1
    }
    cover=cover/rep
    tab[i,]<-c(mean(len),cover)
    cat(round(i/length(size)*100,0),"%\n")
    i<-i+1
  }
  return(tab)
}

# CI from average - beta and inverted beta table function
TabCIAB<-function(delta,prior,post,a,b,sizes){
  set.seed(24680)
  tab<-matrix(nrow=length(sizes),ncol=length(a))
  dimnames(tab)<-list(sizes)
  rep<-1000
  cover<-numeric(rep)
  for(i in 1:length(a)){
    rhos<-RandNAR(rep,prior,a[i],b[i])
    for(j in 1:length(sizes)){
      for(k in 1:rep){
        smp<-rgeom(sizes[j],1/(1+rhos[k]))
        cover[k]<-CIAver(smp,delta,post,a[i],b[i])
      }
      tab[j,i]<-mean(cover)
    }
    cat("\n",round(i/length(a)*100,0),"%")
  }
  return(tab)
}

# CI from average - jeffreys table function
TabCIAJ<-function(delta,prior,post,sizes){
  set.seed(24680)
  tab<-numeric(length(sizes))
  dimnames(sizes)
  rep<-1000
  cover<-numeric(rep)
  for(i in 1:length(sizes)){
    rhos<-RandNAR(rep,prior)
    for(j in 1:rep){
      smp<-rgeom(sizes[i],1/(1+rhos[j]))
      cover[j]<-CIAver(smp,delta,post)
    }
    tab[i]<-mean(cover)
    cat("\n",round(i/length(sizes)*100),"%")
  }
  return(t(tab))
}

# sample size tables -----------------------------------------------------------

# table function for sample size ACC algorithm - beta and inverted 

TabACC<-function(width,prior,post,a,b){
  tab<-matrix(nrow=length(a),ncol=length(width))
  for(i in 1:length(a)){
    for(j in 1:length(width)){
      set.seed(24680)
      tab[i,j]<-ACC_HDI(width[j],prior,post,a[i],b[i])
    }
    cat(round(i/length(a)*100,0),"%\n")
  }
  return(tab)
}

# table function for sample size ALC algorithm - beta and inverted
TabALC<-function(width,prior,post,CIf,a,b){
  tab<-matrix(nrow=length(a),ncol=length(width))
  for(i in 1:length(a)){
    for(j in 1:length(width)){
      set.seed(24680)
      tab[i,j]<-ALC(width[j],prior,post,CIf,a[i],b[i])
    }
    cat(round(i/length(a)*100,0),"%\n")
  }
  return(tab)
}
################################################################################
# generating credible interval abacuses
################################################################################

sizes<-seq(10,400,by=10)

# beta -------------------------------------------------------------------------

a<-c(2.5,1,2,1.5)
b<-c(1.5,1,2,2.5)

# symmetrical CI
BetaSI08<-TabCIAB(0.08,BPrior,BPost,a,b,sizes)
BetaSI10<-TabCIAB(0.10,BPrior,BPost,a,b,sizes)
BetaSI12<-TabCIAB(0.12,BPrior,BPost,a,b,sizes)
BetaSI15<-TabCIAB(0.15,BPrior,BPost,a,b,sizes)

BetaSI08
BetaSI10
BetaSI12
BetaSI15
save(BetaSI08,BetaSI10,BetaSI12,BetaSI15,a,b,file="TabBetaAver.rdata")
load("TabBetaAver.rdata")

# equal-tailed CI
BetaETI1<-TabCIET(BPrior,BPost,sizes,ETI,a[1],b[1])
BetaETI2<-TabCIET(BPrior,BPost,sizes,ETI,a[2],b[2])
BetaETI3<-TabCIET(BPrior,BPost,sizes,ETI,a[3],b[3])
BetaETI4<-TabCIET(BPrior,BPost,sizes,ETI,a[4],b[4])

BetaETI1
BetaETI2
BetaETI3
BetaETI4
save(BetaET1,BetaET2,BetaET3,BetaET4,a,b,file="TabBetaETI.rdata")
#load("TabBetaET.rdata")

# highest density CI
BetaHD1<-TabCI(BPrior,BPost,sizes,HDI,a[1],b[1])
BetaHD2<-TabCI(BPrior,BPost,sizes,HDI,a[2],b[2])
BetaHD3<-TabCI(BPrior,BPost,sizes,HDI,a[3],b[3])
BetaHD4<-TabCI(BPrior,BPost,sizes,HDI,a[4],b[4])

BetaHD1
BetaHD2
BetaHD3
BetaHD4
save(BetaHDI1,BetaHDI2,BetaHDI3,BetaHDI4,a,b,file="TabBetaHD.rdata")
#load("TabBetaHD.rdata")

# inverted beta ----------------------------------------------------------------

a<-c(1.0,1.0,1.0,2.0,2.0,2.0,3.0,3.0,3.0)
b<-c(1.1,2.1,3.0,1.1,2.1,3.0,1.1,2.1,3.0)

# symmetrical CI
IBetaCIA08<-TabCIAB(0.08,IBPrior,IBPost,a,b,sizes)
IBetaCIA10<-TabCIAB(0.10,IBPrior,IBPost,a,b,sizes)
IBetaCIA12<-TabCIAB(0.12,IBPrior,IBPost,a,b,sizes)
IBetaCIA15<-TabCIAB(0.15,IBPrior,IBPost,a,b,sizes)

IBetaCIA08
IBetaCIA10
IBetaCIA12
IBetaCIA15
save(IBetaCIA08,IBetaCIA10,IBetaCIA12,IBetaCIA15,a,b,file="TabIBetaAver.rdata")
load("TabIBetaAver.rdata")

# equal-tailed CI
IBetaET1<-TabCI(IBPrior,IBPost,sizes,ETI,a[1],b[1])
IBetaET2<-TabCI(IBPrior,IBPost,sizes,ETI,a[2],b[2])
IBetaET3<-TabCI(IBPrior,IBPost,sizes,ETI,a[3],b[3])
IBetaET4<-TabCI(IBPrior,IBPost,sizes,ETI,a[4],b[4])
IBetaET5<-TabCI(IBPrior,IBPost,sizes,ETI,a[5],b[5])
IBetaET6<-TabCI(IBPrior,IBPost,sizes,ETI,a[6],b[6])
IBetaET7<-TabCI(IBPrior,IBPost,sizes,ETI,a[7],b[7])
IBetaET8<-TabCI(IBPrior,IBPost,sizes,ETI,a[8],b[8])
IBetaET9<-TabCI(IBPrior,IBPost,sizes,ETI,a[9],b[9])

IBetaET1
IBetaET2
IBetaET3
IBetaET4
IBetaET5
IBetaET6
IBetaET7
IBetaET8
IBetaET9

save(IBetaET1,IBetaET2,IBetaET3,IBetaET4,IBetaET5,IBetaET6,
IBetaET7,IBetaET8,IBetaET9,a,b,file="TabIBetaET.rdata")
load("TabIBetaET.rdata")

# highest density CI
IBetaHD1 <- TabCI(IBPrior,IBPost,sizes,HDI,a[1],b[1])
IBetaHD2 <- TabCI(IBPrior,IBPost,sizes,HDI,a[2],b[2])
IBetaHD3 <- TabCI(IBPrior,IBPost,sizes,HDI,a[3],b[3])
IBetaHD4 <- TabCI(IBPrior,IBPost,sizes,HDI,a[4],b[4])
IBetaHD5 <- TabCI(IBPrior,IBPost,sizes,HDI,a[5],b[5])
IBetaHD6 <- TabCI(IBPrior,IBPost,sizes,HDI,a[6],b[6])
IBetaHD7 <- TabCI(IBPrior,IBPost,sizes,HDI,a[7],b[7])
IBetaHD8 <- TabCI(IBPrior,IBPost,sizes,HDI,a[8],b[8])
IBetaHD9 <- TabCI(IBPrior,IBPost,sizes,HDI,a[9],b[9])

IBetaHD1
IBetaHD2
IBetaHD3
IBetaHD4
IBetaHD5
IBetaHD6
IBetaHD7
IBetaHD8
IBetaHD9

save(IBetaHD1,IBetaHD2,IBetaHD3,IBetaHD4,IBetaHD5,IBetaHD6,
     IBetaHD7,IBetaHD8,IBetaHD9,a,b,file="TabIBetaHDI.rdata")
load("TabIBetaHDI.rdata")

# jeffreys ---------------------------------------------------------------------

# symmetrical CI
JefCIA08<-TabCIAJ(0.08,JPrior,JPost,sizes)
JefCIA10<-TabCIAJ(0.10,JPrior,JPost,sizes)
JefCIA12<-TabCIAJ(0.12,JPrior,JPost,sizes)
JefCIA15<-TabCIAJ(0.15,JPrior,JPost,sizes)

JefCIA08
JefCIA10
JefCIA12
JefCIA15
save(JefCIA08,JefCIA10,JefCIA12,JefCIA15,file="TabJefAver.rdata")
load("TabJefAver.rdata")

# equal-tailed CI
JefET<-TabCI(JPrior,JPost,sizes,ETI)
JefETI

save(JefET,file="TabJefET.rdata")
load("TabJefET.rdata")

# highest density CI
JefHD<-TabCI(JPrior,JPost,sizes,HDI)
JefHDI

save(JefHDI,file="TabJefHD.rdata")
load("TabJefHD.rdata")

################################################################################
# generating sample size tables
################################################################################

# beta -------------------------------------------------------------------------

a<-c(2.5,1,2,1.5)
b<-c(1.5,1,2,2.5)

# symmetrical CI
delta<-c(0.08,0.10,0.12,0.15)
BSampSizeAver<-TabACC(delta,BPrior,BPost,BAver,a,b)
BSampSizeAver

# equal-tailed CI
width<-c(0.16,0.20,0.24,0.30)
BSampSizeET<-TabALC(width,BPrior,BPost,ETI,a,b)
BSampSizeET

# highest density CI - ALC
width<-c(0.16,0.20,0.24,0.30)
BSampSizeHDI<-TabALC(width,BPrior,BPost,HDI,a,b)
BSampSizeHDI
save(BSampSizeHDI,file="BetaSSHDI.rdata")

# highest density CI - ACC
width<-c(0.16,0.20,0.24,0.30)
BSampSizeHDIL<-TabACC(width,BPrior,BPost,a,b)
BSampSizeHDI
save(BSampSizeHDIL,file="BetaSSHDIL.rdata")


save(BSampSizeAver,BSampSizeET,BSampSizeHDI,a,b,delta,width,file="TabBSS.rdata")

# inverted beta ----------------------------------------------------------------

a<-c(1,1,1,2,2,2,3,3,3)
b<-c(1.1,2.1,3,1.1,2.1,3,1.1,2.1,3)

# symmetrical CI
delta<-c(0.08,0.10,0.12,0.15)
ISampSizeAver<-TabACC(delta,IBPrior,IBPost,IBAver,a,b)
ISampSizeAver

# equal-tailed CI
width<-c(0.16,0.20,0.24,0.30)
ISampSizeET<-TabALC(width,IBPrior,IBPost,ETI,a,b)
ISampSizeET

# highest density - ALC
width<-c(0.16,0.20,0.24,0.30)
ISampSizeHDI<-TabALC(width,IBPrior,IBPost,HDI,a,b)
ISampSizeHDI

# # highest density - ACC
# width<-c(0.16,0.20,0.24,0.30)
# ISampSizeHDIL<-TabACC(width,IBPrior,IBPost,a,b)
# ISampSizeHDIL


save(ISampSizeAver,ISampSizeET,ISampSizeHDI,a,b,delta,width,file="TabISS.rdata")
#system('shutdown -s')

# jeffreys ---------------------------------------------------------------------

# symmetrical CI
delta<-c(0.08,0.10,0.12,0.15)
JSampSizeAver<-numeric(length(delta))
for(i in 1:length(delta)){
  set.seed(24680)
  JSampSizeAver[i]<-ACC(delta[i],JPrior,JPost,JAver)
  cat(JSampSizeAver[i],"\n")
}
JSampSizeAver

# equal-tailed CI
width<-c(0.16,0.20,0.24,0.30)
JSampSizeET<-numeric(length(width))
for(i in 1:length(width)){
  set.seed(24680)
  JSampSizeET[i]<-ALC(width[i],JPrior,JPost,ETI)
  cat(JSampSizeET[i],"\n")
}
JSampSizeET

# highest density - ALC
width<-c(0.16,0.20,0.24,0.30)
JSampSizeHDI<-numeric(length(width))
for(i in 1:length(width)){
  set.seed(24680)
  JSampSizeHDI[i]<-ALC(width[i],JPrior,JPost,HDI)
  cat(JSampSizeHDI[i],"\n")
}
JSampSizeHDI

# highest density - ACC
width<-c(0.16,0.20,0.24,0.30)
JSampSizeHDIL<-numeric(length(width))
for(i in 1:length(width)){
  set.seed(24680)
  JSampSizeHDIL[i]<-ACC(width[i],JPrior,JPost)
  cat(JSampSizeHDIL[i],"\n")
}
JSampSizeHDIL

save(JSampSizeAver,JSampSizeET,JSampSizeHDI,a,b,delta,width,file="TabJSS.rdata")
#load('TabJSS.rdata')

################################################################################
# plotting prior distributions
################################################################################

# beta prior
setEPS()
#postscript(paste("FiB.eps",sep=""),width=10.5*0.75,height=8*0.75)
#dev.new(width=10.5*0.75,height=8*0.75)
par(mfrow=c(1,1))
rho<-seq(0.0,1.0,0.05)
plot(rho,BPrior(rho,a=2.5,b=1.5),type="n",xlab=expression(rho),
     ylab=expression(p[b](rho)))
lines(rho,BPrior(rho,a=2.5,b=1.5),type="b",lty=1,pch=1,col="black",lwd=2)
lines(rho,BPrior(rho,a=1.0,b=1.1),type="b",lty=2,pch=2,col="deepskyblue3",lwd=2)
lines(rho,BPrior(rho,a=2.0,b=2.1),type="b",lty=3,pch=3,col="brown3",lwd=2)
lines(rho,BPrior(rho,a=1.5,b=2.5),type="b",lty=4,pch=4,col="green4",lwd=2)
legend("bottom",lty=c(1,2,3,4),pch=c(1,2,3,4),col=c("black","deepskyblue3",
                                                    "brown3","green4"),lwd=2,
       legend=c(
         expression(paste("Beta(",alpha,"=2.5; ",beta,"=1.5)")),
         expression(paste("Beta(",alpha,"=1.0; ",beta,"=1.1)")),
         expression(paste("Beta(",alpha,"=2.0; ",beta,"=2.1)")),
         expression(paste("Beta(",alpha,"=1.5; ",beta,"=2.5)"))))
graphics.off()

# inverted beta prior
setEPS()
#postscript(paste("FiIB.eps",sep=""),width=10.5*0.75,height=8*0.75)
#dev.new(width=10.5*0.75,height=8*0.75)
par(mfrow=c(1,1))
rho<-seq(0.0,1.0,0.05)
plot(rho,IBPrior(rho,a=1.1,b=3.0),type="n",xlab=expression(rho),
     ylab=expression(p[I](rho)),ylim=c(0,3))
lines(rho,IBPrior(rho,a=1.0,b=1.1),type="b",lty=1,pch=1,col="black",lwd=2)
lines(rho,IBPrior(rho,a=1.0,b=3.0),type="b",lty=1,pch=3,col="deepskyblue3",
      lwd=2)
lines(rho,IBPrior(rho,a=2.0,b=3.0),type="b",lty=2,pch=6,col="brown3",lwd=2)
lines(rho,IBPrior(rho,a=3.0,b=1.1),type="b",lty=3,pch=7,col="green4",lwd=2)
lines(rho,IBPrior(rho,a=1.0,b=2.1),type="b",lty=1,pch=2,col="gray50",lwd=2)
lines(rho,IBPrior(rho,a=2.0,b=1.1),type="b",lty=2,pch=4,col="gray50",lwd=2)
lines(rho,IBPrior(rho,a=2.0,b=2.1),type="b",lty=2,pch=5,col="gray50",lwd=2)
lines(rho,IBPrior(rho,a=3.0,b=2.1),type="b",lty=3,pch=8,col="gray50",lwd=2)
lines(rho,IBPrior(rho,a=3.0,b=3.0),type="b",lty=3,pch=9,col="gray50",lwd=2)
legend("top",lty=c(1,1,1,2,2,2,3,3,3),pch=1:9,
       col=c("black","deepskyblue3","brown3","green4","gray50",
             "gray50","gray50","gray50","gray50"),lwd=2,
       legend=c(
         expression(paste("IIB(",alpha,"=1.0; ",beta,"=1.1)")),
         expression(paste("IIB(",alpha,"=1.0; ",beta,"=3.0)")),
         expression(paste("IIB(",alpha,"=2.0; ",beta,"=3.0)")),
         expression(paste("IIB(",alpha,"=3.0; ",beta,"=1.1)")),
         expression(paste("IIB(",alpha,"=1.0; ",beta,"=2.1)")),
         expression(paste("IIB(",alpha,"=2.0; ",beta,"=1.1)")),
         expression(paste("IIB(",alpha,"=2.0; ",beta,"=2.1)")),
         expression(paste("IIB(",alpha,"=3.0; ",beta,"=2.1)")),
         expression(paste("IIB(",alpha,"=3.0; ",beta,"=3.0)"))
       ))
graphics.off()

# jeffreys prior 
setEPS()
#postscript(paste("FiJ.eps",sep=""),width=10.5*0.75,height=8*0.75)
par(mfrow=c(1,1))
rho<-seq(0.0,1.0,0.05)
plot(rho,JPrior(rho),xlab=expression(rho),
     ylab=expression(p[J](rho)),type="b",lty=1,pch=1,col="green4",lwd=2)
legend("top",lty=1,pch=1,col="green4",lwd=2,legend=c("Jeffreys"))
graphics.off()

################################################################################
# tests
################################################################################

# testing equal-tailed
rho<-0.7
smp<-rgeom(500,1/(1+rho))
a<-1.0
b<-1.0
c<-CIET(BPost,smp,a,b)
beta(575+2.5,1.5)
GaussHyp(501+575,575+2.5,575+2.5+1.5,-1)

# testing credible interval from maximum and mean
a=1; b=1; rep<-1000 
rho<-RandNAR(rep,p1,a,b)
coverA<-numeric(rep)
coverM<-numeric(rep)
for(i in 1:rep){
  smp<-rgeom(84,1/(1+rho[i]))
  coverA[i]<-CIAver(smp,0.15,BPost,a,b)
  coverM[i]<-CIMax(smp,0.15,BPost,a,b)
}
mean(coverA)
mean(coverM)

# testing beta average function
a<-1.0
b<-1.0
size<-600
est<-numeric(1000)
for(i in 1:1000){
  smp<-rgeom(size,1/(1+0.95))
  est[i]<-BAver(smp,a,b)
}
length(est[is.finite(est)])

# testing average estimator
smp<-rgeom(20,1/1.5)
AverE(smp,IBPost,1,2)-IBAver(smp,1,2)

# testing ACC.HDI
ACC_HDI(0.20,BPost,BPrior,1,1)

# testing random number function by acceptance-rejection method
RandNAR(10,BPrior,2,2)

# testing maximum of distribution estimator by golden-section
MaxPost(BPrior,2,2)
MaxPost(BPrior,1.5,2.5)

# testing incomplete inverted beta function
RITest<-function(csi,a,b){
  Aux<-function(u,a,b){
    auxOut<-u^(a-1)*(1-u)^(b-1)
  }
  output<-1/beta(a,b)*integrate(Aux,0,csi,a,b)[[1]]
  return(output)
  }
RITest(0.5,1,2)[[1]]-RIB_0.5(1,2)[[1]]
rm(RITest)

# testing average coverage function
ACC(0.15,BPrior,BPost,BAver,2.5,1.5)

# testing ALC algorithm
ALC(0.30,BPrior,BPost,ETI,2.5,1.5)

# comparing algorithms
MeasureTime(ACC(0.15,BPrior,BPost,BAver,2.5,1.5))
MeasureTime(ALC(0.30,BPrior,BPost,ETI,2.5,1.5))

MeasureTime<-function(algorithm){
  start.time<-Sys.time()
  algorithm
  end.time<-Sys.time()
  return(end.time-start.time)
}

# testing Jeffreys prior/posterior and comparison with I.B.
plot(seq(0,1,0.001),JPrior(seq(0,1,0.001)),type='l')
RandNAR(500,JPrior)
RandNAR(500,IBPrior,1,3)
plot(seq(0,1,0.001),BPost(seq(0,1,0.001),rgeom(400,1/(1+1e-3)),2.5,1.5),type='l')
plot(seq(0,1,0.001),JPost(seq(0,1,0.001),rgeom(400,1/(1+1e-3))),type='l')

integrate(JPost,1e-4,1-1e-4,rgeom(400,1/(1+1e-2)))[[1]]

# testing HDI intervals
a<-1; b<-1
rhoE<-numeric(1000)
auxF<-function(rho,x,a,b){
  return(BPost(rho,x,a,b)*rho)
}
for(i in 1:1000){
  x<-rgeom(500,1/(1+0.1))
  rhoE[i]<-integrate(auxF,0,1,x,a,b)[[1]]
}
hist(rhoE)
mean(rhoE)


rhos<-seq(0,1,0.001)
fx<-BPost(rhos,x,a,b)
plot(rhos,fx,type='l')

lengthHDI<-0.3
lengthHDIMC<-numeric(100)
for(i in 1:100){
  x<-rgeom(400,1/(1+0.4))
  intervalHDIL<-HDI_length(lengthHDI,BPost,x,a,b)
  lengthHDIMC[i]<-intervalHDIL[2]-intervalHDIL[1]
  #(intervalHDIL[2]-intervalHDIL[1])-lengthHDI
}
hist(lengthHDIMC)
mean(lengthHDIMC)
(intervaloETI<-ETI(BPost,x,a,b))
(lengthETI<-intervaloETI[2]-intervaloETI[1])
integrate(BPost,intervaloETI[1],intervaloETI[2],x,a,b)[[1]]

(intervalHDI<-HDI(BPost,x,a,b))
(lengthHDI <- intervalHDI[2] - intervalHDI[1])
integrate(BPost,intervalHDI[1],intervalHDI[2],x,a,b)[[1]]



