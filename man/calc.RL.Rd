\name{calc.RL.0}
\alias{calc.RL.0}
\alias{calc.mult.RL}
\alias{RL.max}
\alias{RL.burn}
\title{Use Raftery and Lewis diagnostics to calculate MCMC chain lengths}
\description{
  Uses \code{coda::raftery.diag} (the Raftery and Lewis
  diagnostic) to estimate minimum chain lengths for an MCMC estimate
  for mixed stock analysis.
  Runs R&L iteratively until the criteria are satisfied.
}
\usage{
calc.RL.0(data, startfval, pilot=500, maxit=15, verbose=FALSE,
rseed=1001, debug=FALSE)
calc.mult.RL(data,n=50,debug=FALSE,verbose=FALSE)
RL.max(r)

}
\arguments{
  \item{data}{a \code{mixstock.data} object, or any list with
    \code{sourcesamp} and \code{mixsamp} entries containing source 
    and mixed stock data}
  \item{startfval}{starting value of contribution frequencies for the
    chain, as in \code{gibbs} or \code{gibbsC}: NULL=random start,
    0=equal contributions from all sources, (1..R-1)=95\% contribution
    from one source, with the rest splitting the remainder equally}
  \item{pilot}{Chain length to start with (length of "pilot" run)}
  \item{maxit}{Max. number of iterations of the Raftery and Lewis procedure}
  \item{verbose}{Produce lots of output?}
  \item{rseed}{Random-number seed}
  \item{debug}{produce debugging output?}
  \item{n}{number of different random-number seed chains to try}
  \item{r}{the results of a Raftery and Lewis diagnostic test}
}
\details{
  \code{calc.RL.00} starts by running a Gibbs-sampler chain with the length given by
  \code{pilot}, then repeatedly lengthens the chain until the length is
  greater than that suggested as the total by the Raftery and Lewis
  diagnostic.  (The next suggested step in the procedure is to run
  multiple chains of this length and see whether they pass the Gelman
  and Rubin diagnostic.) \code{calc.mult.RL}
  runs the Raftery and Lewis calculation multiple times, starting
  each chain from a large contribution from each source in turn,
  to see if some starting configurations are slower to converge
 or if there is a lot of variation among chains with different random
 number seeds.
\code{RL.max} picks the expected maximum chain length given
a set of diagnostics; \code{RL.burn} returns the predicted
burn-in required.
}
\value{
for \code{calc.RL.00}:
  \item{current}{Results of the Raftery and Lewis test on the current iteration}
  \item{history}{History of the iterations:}
for \code{calc.mult.RL}, a matrix giving the maximum expected chain
length for each random-number seed/starting point combination
}
\author{Ben Bolker}

\examples{
data(bolten98)
b98c <- markfreq.condense(as.mixstock.data(bolten98))
}
\keyword{misc}
