\name{mcmc.chainlength.est}
\alias{mcmc.chainlength.est}
\title{Estimate appropriate chain length for mixed stock analysis by MCMC} 
\description{
  Determines an appropriate chain length for MCMC estimation of 
  source contributions to a mixed stock
  by running Raftery and Lewis and Gelman & Rubin
  diagnostics repeatedly until convergence criteria are met
}
\usage{
mcmc.chainlength.est(x, mult=1, inflate=sqrt(2), GR.crit=1.2,
nchains=x$R, verbose=FALSE)
}
\arguments{
  \item{x}{Mixed stock analysis data (a \code{mixstock.data} object or a list containing \code{sourcesamp} and \code{mixsamp}}
  \item{mult}{How many different times to run tests}
  \item{inflate}{How much to increase chain length at every failed
    iteration of Gelman and Rubin}
  \item{GR.crit}{Maximum value for Gelman and Rubin 97.5\% quantile
    in order to declare convergence}
  \item{nchains}{number of separate MCMC chains to run}
  \item{verbose}{print lots of detail while running?}
}
\details{
  If \code{mult} is 1, runs Raftery and Lewis diagnostics on a chain
  starting from equal contributions; if \code{mult} is greater than 1,
  runs them on as many chains as there are sourcees, each starting from
  a 95\% contribution from that source.  Iteratively increases each chain
  length to that suggested by the R&L diagnostic, until all chains
  pass.  Then runs Gelman and Rubin on a set of chains starting from
  each source.  If \code{mult} is greater than 1, it does each step on
  \code{mult} different chains and takes the maximum.
}
\value{
  The maximum chainlength needed to get convergence in all tests
}
\author{Ben Bolker}
\seealso{gibbsC}

\examples{
data(simex)
mcmc.chainlength.est(simex,verbose=TRUE)
}
\keyword{misc}
