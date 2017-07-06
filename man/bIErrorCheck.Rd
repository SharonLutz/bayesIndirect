\name{bIErrorCheck}
\alias{bIErrorCheck}
\title{bIErrorCheck}
\description{
Performs error checks on the input for the bayesIndirect function.
}
\usage{
bIErrorCheck(x, k, y, Z = NULL, nCov = 0, plot0 = FALSE, nIts = 50000, nBurn = 10000, propGamma1Pi = 0.7, cx = 10, pix = 0.5)
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
  \item{Z}{
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
}
\value{
The function exits with an error message if the error checks are not met.
}
\author{
Sharon Lutz, Annie Thwing
}