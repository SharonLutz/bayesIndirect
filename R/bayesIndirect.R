bayesIndirect <-
function(x,k,y,z=NULL,nCov=0,plot0=FALSE,nIts=50000,nBurn=10000,propGamma1Pi=0.7,cx=10,pix=0.5,SEED = 1){

  set.seed(SEED)
  
  library(coda) # Load the coda library

  # Format z to be a matrix, Z
  if(!is.null(z)){
    # format Z to be a matrix (different when only one column)
    if(is.vector(z)){
      Z <- matrix(0,nrow=length(z),ncol=1)
      for(pp in 1:length(z)){Z[pp,]<-z[pp]}
    }
    else{
      Z <- matrix(0,nrow=nrow(z),ncol=ncol(z))
      for(pp in 1:ncol(z)){Z[,pp]<-z[,pp]}
    }
  }else{
    Z <- NULL
  }
  
bIErrorCheck(x,k,y,Z,nCov,plot0,nIts,nBurn,propGamma1Pi,cx,pix)

# All instances of z are changed to Z, (Z is the matrix format from the error check)

# NO covariates
if(nCov==0){
# remove NAs
dataN<-cbind(x,k,y,Z)
dataN<-na.omit(dataN)
x<-dataN[,1]
k<-dataN[,2]
y<-dataN[,3]
y<-y-summary(lm(y~x+k))$coef[1,1] #mean center
}

# 1 covariates
if(nCov==1){
# remove NAs
dataN<-cbind(x,k,y,Z)
dataN<-na.omit(dataN)
x<-dataN[,1]
k<-dataN[,2]
y<-dataN[,3]
z1<-dataN[,4]
y<-y-summary(lm(y~x+k+z1))$coef[1,1] #mean center
}

# 2 covariates
if(nCov==2){
# remove NAs
dataN<-cbind(x,k,y,Z)
dataN<-na.omit(dataN)
x<-dataN[,1]
k<-dataN[,2]
y<-dataN[,3]
z1<-dataN[,4]
z2<-dataN[,5]
y<-y-summary(lm(y~x+k+z1+z2))$coef[1,1] #mean center
}

n<-length(y)
yResid<-summary(lm(y~k))$resid

#proposal parameter for gamma1  for Bernoulli dist
# if gamma1 is NOT mixing well=>change this between 0 and 1
# if gamma1 is mostly 1 than make propGamma1Pi closer to 0 but NOT 0
# if gamma1 is mostly 0 than make propGamma1Pi closer to 1 but NOT 1
#set c1,pi1
c1<-cx 
cz<-10
pi1<-pix

# upper bound of sigmas uniform distribution
A <- 100

###############################
# No covariates
############################### 
if(nCov==0){

# likelihood and priors
logLik <- function(gamma1,beta,sig2){
 return(sum(dnorm(yResid,gamma1*x*beta,sqrt(sig2),log=T)))}

gamma1LogPrior <- function(gamma1){
  return(dbinom(gamma1,1,pi1,log=T))}

#proposal distribution for MH independence sampler
gammaProposal <- function(gamma,propPI){
  return(dbinom(gamma,1,propPI,log=T))}


# closed form conditionals for beta, alpha and sig2
betaSample <- function(gamma1,  sig2){
  v1 <- 1/(sum(x^2)/sig2 + 1/c1)
  m1 <- (sum(x*yResid))/sig2
  if(gamma1==1){bSam<-rnorm(1,v1*m1,sqrt(v1))}
  if(gamma1==0){bSam<-0}
  return(bSam)}


sig2Sample <- function(beta){
  sig2 <- Inf
  while(sig2>A^2){
    sig2 <- 1/rgamma(1,(n-1)/2,sum((yResid-x*beta)^2)/2)
  }
  return(sig2)}

# matrix that will hold posterior samples
outMAT <- matrix(0,nr=nIts,nc=3)
colnames(outMAT) <- c("gamma1", "beta","sig2")

# number of accepted draws of gamma1 for testing MCMC
#should be around 40-45 percent acceptance rate
accGamma1<<- 0

# the sampler
sampler <- function(){
 for(ii in 2:nIts){
 
  # gamma1 - Metropolis Hastings independence proposal
  logPrior <- gamma1LogPrior(gamma1)
  LL <- logLik(gamma1,beta,sig2)
  Jprop<-gammaProposal(gamma1,propGamma1Pi)
  
  #proposal distribution
  gamma1Star <- 1-gamma1
  logPriorStar <- gamma1LogPrior(gamma1Star) 
  LLStar <- logLik(gamma1Star,beta,sig2)
  JpropStar<-gammaProposal(gamma1Star,propGamma1Pi)
  r <- exp(LLStar+logPriorStar-Jprop-LL-logPrior+JpropStar)
  if(runif(1)<r){
    accGamma1 <<- accGamma1+1
    gamma1 <- gamma1Star}

  # Gibbs proposal from closed form conditionals
  beta <- betaSample(gamma1, sig2)
  sig2<-sig2Sample(beta)

outMAT[ii,1] <- gamma1
outMAT[ii,2] <- beta
outMAT[ii,3] <- sig2

  #to keep track of iteration number when testing MCMC
   if(ii%%10000==0){print(paste(ii,"in",nIts,"twice"))}
 }
 return(list(outMAT=outMAT))}  


# randomly generate initial values
    inits <- function(){
    gamma1 <- rbinom(1, 1, 0.5)
    beta <- rnorm(1, 0, sqrt(10*.6))
    sig2 <- runif(1, 0, 20)
    return(c(gamma1, beta, sig2)) }

# 1st run of the mcmc 
outMAT[1,] <- inits()
gamma1 <- outMAT[1,1]
beta <- outMAT[1,2]
sig2 <- outMAT[1,3]
MCMCoutput <- sampler() 
# create an object containing the posterior draws
PosteriorDraws1 <- MCMCoutput$outMAT
# create an object containing the burnt-in draws
BurntInDraws1 <- PosteriorDraws1[nBurn:nIts,]

## 2nd run of the mcmc 
outMAT[1,] <- inits()
gamma1 <- outMAT[1,1]
beta <- outMAT[1,2]
sig2 <- outMAT[1,3]
MCMCoutput <- sampler() 
### create an object containing the posterior draws
PosteriorDraws2 <- MCMCoutput$outMAT
### create an object containing the burnt-in draws
BurntInDraws2 <- PosteriorDraws2[nBurn:nIts,]

}#end of 0 covariates


###############################
# 1 covariates
############################### 

if(nCov==1){

# likelihood and priors
logLik <- function(gamma1,beta,alpha1,sig2){
 return(sum(dnorm(yResid,gamma1*x*beta+alpha1*z1,sqrt(sig2),log=T)))}

gamma1LogPrior <- function(gamma1){
  return(dbinom(gamma1,1,pi1,log=T))}

#proposal distribution for MH independence sampler
gammaProposal <- function(gamma,propPI){
  return(dbinom(gamma,1,propPI,log=T))}


# closed form conditionals for beta, alpha and sig2
betaSample <- function(gamma1, alpha1, sig2){
  v1 <- 1/(sum(x^2)/sig2 + 1/c1)
  m1 <- (sum(x*yResid-alpha1*z1*x))/sig2
  if(gamma1==1){bSam<-rnorm(1,v1*m1,sqrt(v1))}
  if(gamma1==0){bSam<-0}
  return(bSam)}

alpha1Sample <- function(sig2, beta){
  v1 <- 1/(sum(z1^2)/sig2 + 1/cz)
  m1 <- (sum(z1*yResid-beta*x*z1))/sig2
  return(rnorm(1,v1*m1,sqrt(v1)))}

sig2Sample <- function(beta,alpha1){
  sig2 <- Inf
  while(sig2>A^2){
    sig2 <- 1/rgamma(1,(n-1)/2,sum((yResid-x*beta-alpha1*z1)^2)/2)
  }
  return(sig2)}

# matrix that will hold posterior samples
outMAT <- matrix(0,nr=nIts,nc=4)
colnames(outMAT) <- c("gamma1", "beta","alpha1","sig2")

# number of accepted draws of gamma1 for testing MCMC
#should be around 40-45 percent acceptance rate
accGamma1<<- 0

# the sampler
sampler <- function(){
 for(ii in 2:nIts){
 
  # gamma1 - Metropolis Hastings independence proposal
  logPrior <- gamma1LogPrior(gamma1)
  LL <- logLik(gamma1,beta,alpha1,sig2)
  Jprop<-gammaProposal(gamma1,propGamma1Pi)
  
  #proposal distribution
  gamma1Star <- 1-gamma1
  logPriorStar <- gamma1LogPrior(gamma1Star) 
  LLStar <- logLik(gamma1Star,beta,alpha1,sig2)
  JpropStar<-gammaProposal(gamma1Star,propGamma1Pi)
  r <- exp(LLStar+logPriorStar-Jprop-LL-logPrior+JpropStar)
  if(runif(1)<r){
    accGamma1 <<- accGamma1+1
    gamma1 <- gamma1Star}

  # Gibbs proposal from closed form conditionals
  beta <- betaSample(gamma1, alpha1, sig2)
  alpha1 <- alpha1Sample(sig2, beta)
  sig2<-sig2Sample(beta,alpha1)

outMAT[ii,1] <- gamma1
outMAT[ii,2] <- beta
outMAT[ii,3] <- alpha1
outMAT[ii,4] <- sig2

  #to keep track of iteration number when testing MCMC
   if(ii%%10000==0){print(ii)}
 }
 return(list(outMAT=outMAT))}  


# randomly generate initial values
    inits <- function(){
    gamma1 <- rbinom(1, 1, 0.5)
    beta <- rnorm(1, 0, sqrt(10*.6))
    alpha1 <- rnorm(1, 0, sqrt(10*.6))
    sig2 <- runif(1, 0, 20)
    return(c(gamma1, beta, alpha1, sig2)) }

# 1st run of the mcmc 
outMAT[1,] <- inits()
gamma1 <- outMAT[1,1]
beta <- outMAT[1,2]
alpha1 <- outMAT[1,3]
sig2 <- outMAT[1,4]
MCMCoutput <- sampler() 
# create an object containing the posterior draws
PosteriorDraws1 <- MCMCoutput$outMAT
# create an object containing the burnt-in draws
BurntInDraws1 <- PosteriorDraws1[nBurn:nIts,]

## 2nd run of the mcmc 
outMAT[1,] <- inits()
gamma1 <- outMAT[1,1]
beta <- outMAT[1,2]
alpha1 <- outMAT[1,3]
sig2 <- outMAT[1,4]
MCMCoutput <- sampler() 
### create an object containing the posterior draws
PosteriorDraws2 <- MCMCoutput$outMAT
### create an object containing the burnt-in draws
BurntInDraws2 <- PosteriorDraws2[nBurn:nIts,]

}#end of 1 covariate

###############################
# 2 covariates
############################### 

if(nCov==2){

# likelihood and priors
logLik <- function(gamma1,beta,alpha1,alpha2,sig2){
 return(sum(dnorm(yResid,gamma1*x*beta+alpha1*z1+alpha2*z2,sqrt(sig2),log=T)))}

gamma1LogPrior <- function(gamma1){
  return(dbinom(gamma1,1,pi1,log=T))}

#proposal distribution for MH independence sampler
gammaProposal <- function(gamma,propPI){
  return(dbinom(gamma,1,propPI,log=T))}


# closed form conditionals for beta, alpha and sig2
betaSample <- function(gamma1, alpha1, alpha2, sig2){
  v1 <- 1/(sum(x^2)/sig2 + 1/c1)
  m1 <- (sum(x*yResid-alpha1*z1*x-alpha2*z2*x))/sig2
  if(gamma1==1){bSam<-rnorm(1,v1*m1,sqrt(v1))}
  if(gamma1==0){bSam<-0}
  return(bSam)}

alpha1Sample <- function(alpha2, sig2, beta){
  v1 <- 1/(sum(z1^2)/sig2 + 1/cz)
  m1 <- (sum(z1*yResid-beta*x*z1-alpha2*z2*z1))/sig2
  return(rnorm(1,v1*m1,sqrt(v1)))}

alpha2Sample <- function(alpha1, sig2, beta){
  v1 <- 1/(sum(z2^2)/sig2 + 1/cz)
  m1 <- (sum(z2*yResid-beta*x*z2-alpha1*z1*z2))/sig2
  return(rnorm(1,v1*m1,sqrt(v1)))}

sig2Sample <- function(beta,alpha1,alpha2){
  sig2 <- Inf
  while(sig2>A^2){
    sig2 <- 1/rgamma(1,(n-1)/2,sum((yResid-x*beta-alpha1*z1-alpha2*z2)^2)/2)
  }
  return(sig2)}

# matrix that will hold posterior samples
outMAT <- matrix(0,nr=nIts,nc=5)
colnames(outMAT) <- c("gamma1", "beta","alpha1","alpha2","sig2")

# number of accepted draws of gamma1 for testing MCMC
#should be around 40-45 percent acceptance rate
accGamma1<<- 0


# the sampler
sampler <- function(){
 for(ii in 2:nIts){
 
  # gamma1 - Metropolis Hastings independence proposal
  logPrior <- gamma1LogPrior(gamma1)
  LL <- logLik(gamma1,beta,alpha1,alpha2,sig2)
  Jprop<-gammaProposal(gamma1,propGamma1Pi)
  
  #proposal distribution
  gamma1Star <- 1-gamma1
  logPriorStar <- gamma1LogPrior(gamma1Star) 
  LLStar <- logLik(gamma1Star,beta,alpha1,alpha2,sig2)
  JpropStar<-gammaProposal(gamma1Star,propGamma1Pi)
  r <- exp(LLStar+logPriorStar-Jprop-LL-logPrior+JpropStar)
  if(runif(1)<r){
    accGamma1 <<- accGamma1+1
    gamma1 <- gamma1Star}

  # Gibbs proposal from closed form conditionals
  beta <- betaSample(gamma1, alpha1, alpha2, sig2)
  alpha1 <- alpha1Sample(alpha2, sig2, beta)
  alpha2 <- alpha2Sample(alpha1, sig2, beta)
  sig2<-sig2Sample(beta,alpha1,alpha2)

outMAT[ii,1] <- gamma1
outMAT[ii,2] <- beta
outMAT[ii,3] <- alpha1
outMAT[ii,4] <- alpha2
outMAT[ii,5] <- sig2

  #to keep track of iteration number when testing MCMC
   if(ii%%10000==0){print(ii)}
 }
 return(list(outMAT=outMAT))}  


# randomly generate initial values
    inits <- function(){
    gamma1 <- rbinom(1, 1, 0.5)
    beta <- rnorm(1, 0, sqrt(10*.6))
    alpha1 <- rnorm(1, 0, sqrt(10*.6))
    alpha2 <- rnorm(1, 0, sqrt(10*.6))
    sig2 <- runif(1, 0, 20)
    return(c(gamma1, beta, alpha1, alpha2, sig2)) }

# 1st run of the mcmc 
outMAT[1,] <- inits()
gamma1 <- outMAT[1,1]
beta <- outMAT[1,2]
alpha1 <- outMAT[1,3]
alpha2 <- outMAT[1,4]
sig2 <- outMAT[1,5]
MCMCoutput <- sampler() 
# create an object containing the posterior draws
PosteriorDraws1 <- MCMCoutput$outMAT
# create an object containing the burnt-in draws
BurntInDraws1 <- PosteriorDraws1[nBurn:nIts,]

## 2nd run of the mcmc 
outMAT[1,] <- inits()
gamma1 <- outMAT[1,1]
beta <- outMAT[1,2]
alpha1 <- outMAT[1,3]
alpha2 <- outMAT[1,4]
sig2 <- outMAT[1,5]
MCMCoutput <- sampler() 
### create an object containing the posterior draws
PosteriorDraws2 <- MCMCoutput$outMAT
### create an object containing the burnt-in draws
BurntInDraws2 <- PosteriorDraws2[nBurn:nIts,]

}#end of 2 covariates


BurntInSamps<-mcmc.list(mcmc(BurntInDraws1),mcmc(BurntInDraws2))
	
if(plot0==TRUE){
pdf(paste("densityPlot.pdf",sep=""))
plot(density(BurntInDraws2[,2]),main="Posterior Density for Beta",xlab="")
dev.off()
}

summary(BurntInSamps)
}
