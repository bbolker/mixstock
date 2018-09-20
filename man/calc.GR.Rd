\name{calc.GR}
\alias{calc.GR}
\title{Run chains and calculate Gelman and Rubin diagnostics for mixed stock analyses}
\description{
  For a mixed stock analysis data set with \code{R} sources, runs \code{R}
  Gibbs-sampler chains (one starting from a majority of the contribution
  from each source) and calculates Gelman and Rubin diagnostics on the
  combined set of chains.
}
\usage{
calc.GR(data, tot=20000, burn=1, verbose=FALSE, rseed=1001,
           chainfrac=NULL)
}
\arguments{
  \item{data}{A mixed stock analysis data set (either an object of type
    \code{mixstock.data}, or any list containing appropriate
    \code{sourcesamp} and \code{mixsamp} objects)}
  \item{tot}{Total length of each MCMC chain}
  \item{burn}{Burn-in time for MCMC chains}
  \item{verbose}{Produce verbose output?}
  \item{chainfrac}{Fraction of each chain to discard.  By default,
    chainfrac is set to (1-1/R), so that the total length of the
    combined chain is equal to (R*tot/R=tot).}
  \item{rseed}{Seed for random-number generator}
}
\details{
  calls CODA function \code{gelman.diag}
}
\value{
  Returns the diagnostic from the CODA function
}
\author{Ben Bolker}
\keyword{misc}
