\name{gibbsmat}
\alias{gibbsmat}
\alias{gibbsrpt}
\alias{plotvar}
\title{Plot Gibbs sampler output}
\description{
   plot or report on output from Gibbs sampler  
}
\usage{
gibbsmat(x, burnin = 500, R = 2, H = 2, trans = TRUE)
gibbsrpt(x,burnin=500,R=2,H=2,trans=TRUE,plot=TRUE)
plotvar(vec,best,name,...)
}
\arguments{
  \item{x}{results of an MCMC fit}
  \item{burnin}{length of burn-in period for chain}
  \item{R}{number of rookeries}
  \item{H}{number of haplotypes}
  \item{trans}{transform chain values?}
  \item{plot}{plot posterior distributions for each variable?}
  \item{vec}{numeric vector (posterior distribution chain)}
  \item{best}{best-fit/point estimate value}
  \item{name}{name of variable, for x label}
  \item{...}{other arguments to plot}
}
\value{
  sets up an array of plots on the current output device
}
\keyword{misc}
