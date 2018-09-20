\name{calc.mult.GR}
\alias{calc.mult.GR}
\title{Calculate Gelman and Rubin diagnostic multiple times for mixed stock analyses}
\description{
  Runs Gelman and Rubin diagnostics from CODA multiple times, to get an
  idea of the variation in convergence statistics
}
\usage{
calc.mult.GR(data, n=10, tot=20000, burn=1, verbose=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{Mixed stock analysis data (a list with elements \code{sourcesamp} and \code{mixsamp}}
  \item{n}{Number of replicates to run}
  \item{tot}{Total number of iterates for each chain}
  \item{burn}{Burn-in time for each chain}
  \item{verbose}{Produce verbose output?}
}
\details{
  Runs \code{calc.GR} multiple times, produces a summary table of the
  maximum point estimate and maximum 97.5\% estimate (across all
  variables) for each run with a different random-number seed.
}
\value{
  A matrix with each row giving the random-number seed, max. point
  estimate, max. 97.5\% quantile for each run.
}
\note{The generally accepted criteria for declaring convergence
  according to Gelman and Rubin is that all of the 97.5\% quantiles of
  the estimates of shrink factors are less than 1.2.}							
\author{Ben Bolker}
\seealso{calc.GR}
\examples{
data(simex)
calc.GR(simex,tot=2000)
}
\keyword{misc}
