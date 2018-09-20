\name{tmcmc}
\alias{tmcmc}
\alias{gibbs}
\title{Mixed stock analysis by Markov Chain Monte Carlo}
\description{
   Runs a Gibbs sampler MCMC starting with 95\% contribution from
each source, then combines the chains
}
\usage{
tmcmc(data, n.iter=20000, rseed=1001, n.burnin=floor(n.iter/2),
n.thin=max(1, floor(n.chains * (n.iter - n.burnin)/1000)),
verbose=FALSE, fprior=NULL,
contrib.only=TRUE, rptiter=-1,
outfile=NULL, lang="C",a=NULL,gr=FALSE)
gibbs(sourcesamp, mixsamp, a = 1, startiter, maxiter, startfval = NULL, 
    n.thin = 1, fprior = NULL, rptiter = -1) 
}
\arguments{
  \item{data}{Data: a \code{mixstock.data} object}
  \item{n.iter}{Total length of each chain}
  \item{n.burnin}{Number of burn-in iterations}
  \item{n.thin}{thinning rate.  Must be a positive integer. Set 'n.thin' > 1
    to save memory and computation time if 'n.iter' is large.
    Default is 'max(1, floor(n.chains * (n.iter-n.burnin) /
    1000))' which will only thin if there are at least 2000
    simulations.}
  \item{rseed}{Random-number seed}
  \item{verbose}{Produce lots of output}
  \item{fprior}{Bayesian prior for source contributions}
  \item{contrib.only}{To save memory, store only information about
    contributions from each source and not about
    the estimated marker frequencies in each source}
  \item{rptiter}{How often to issue a progress report.  Negative numbers
    mean no reports}
  \item{outfile}{file to use for output}
  \item{lang}{Run the chain in C or R (for
    debugging/testing purposes only)?}
  \item{a}{prior strength parameter}
  \item{gr}{calculate Gelman-Rubin convergence statistic?}
  \item{sourcesamp}{matrix of marker samples from sources}
  \item{mixsamp}{vector of marker samples from mixed stock}
  \item{startiter}{starting iteration}
  \item{maxiter}{max. number of iterations}
  \item{startfval}{starting source contributions}
}
\value{
   Returns an object of type \code{mixstock.est}
}
\references{Masuda and Pella}
\author{Ben Bolker}
\examples{
  data(bolten98)
  b98c <- markfreq.condense(as.mixstock.data(bolten98))
  t1 <- tmcmc(b98c); t1
}
\keyword{misc}
