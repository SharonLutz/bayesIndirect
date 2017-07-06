\name{bayesIndirect}
\alias{bayesIndirect}
\title{bayesIndirect}

\description{
  bayesIndirect examines the indirect effect of a SNP on the outcome through the mediator in a Bayesian framework with a spike and slab prior.
}
\usage{
bayesIndirect(x, k, y, z = NULL, nCov = 0, plot0 = FALSE, nIts = 50000, nBurn = 10000, propGamma1Pi = 0.7, cx = 10, pix = 0.5, SEED = 1)
}

\arguments{
  \item{x}{
  is the exposure, such as the SNP.
}
  \item{k}{
  is the mediator.
}
  \item{y}{
  is the normally distributed outcome.
}
  \item{z}{
  is the matrix of covariates.
}
  \item{nCov}{
  is the number of covariates. nCov = 0, 1 or 2.
}
  \item{plot0}{
  true if density plot should be outputted or false otherwise.
}
  \item{nIts}{
  is the number of iterations of the MCMC.
}
  \item{nBurn}{
  is the Burn-in period.
}
  \item{propGamma1Pi}{
  is the probability of gamma from the Bernoulli distribution.
}
  \item{cx}{
  is the variance of beta X from the normal distribution.
}
  \item{pix}{
  is the probability that beta x is sampled from the normal distribution.
}
  \item{SEED}{
  is the seed for the SKAT function, default = 1.
}
}

\value{
The posterior mean, standard deviation and quantiles for each variable.
}
\references{
Lutz SM, Hokanson JE, Sharma S, Weiss S, Raby B, Lange C. (2013) On the Integration of Expression Profiles in Genetic Association Studies: A Bayesian Approach to Determine the Path from Gene to Disease. Open Journal of Genetics. 3(3). 
}
\author{
Sharon Lutz, Annie Thwing
}

\section{Warning}{
library(coda) is needed to run this function. 
}

\keyword{ function }
\keyword{ bayesian }
\keyword{ indirect }
