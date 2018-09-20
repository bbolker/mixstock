\name{gibbsC}
\alias{gibbsC}
\title{Run mixed stock analysis Gibbs sampler, in C}
\description{
  Runs a Gibbs sampler MCMC for mixed stock analysis, calling a
  routine written in C (for code).  Low-level function, called
  by other functions.
}
\usage{
gibbsC(a=1, startiter, maxiter, data, mixsamp=NULL, sourcesamp=NULL,
startfval=NULL, thin=1, fprior=NULL, outfile=FALSE,
outfn="mixstock-gibbs", randseed=1001, rptiter=-1, debug=FALSE,
contrun=FALSE,contrib.start=NULL,sourcefreq.start=NULL)
}
\arguments{
  \item{a}{Prior strength parameter}
  \item{startiter}{Number of iterations to discard (burn-in)}
  \item{maxiter}{Total number of chain steps}
  \item{data}{A \code{mixstock.data} object}
  \item{mixsamp}{Marker sample from mixed populations}
  \item{sourcesamp}{Marker samples from sources}
  \item{startfval}{Where to start the chain: 0=}
  \item{thin}{thinning factor}
  \item{fprior}{Bayesian prior}
  \item{outfile}{send data to an output file?}
  \item{outfn}{name of output file}
  \item{randseed}{random-number seed}
  \item{rptiter}{frequency for sending reports to screen}
  \item{debug}{debug?}
  \item{contrun}{continuation run? set fval and sourcefreq directly}
  \item{contrib.start}{vector of starting contributions}
  \item{sourcefreq.start}{matrix of starting source freqs}
}
\value{
  a numeric matrix containing samples from the chain:
  each row is a vector of estimated contribution frequencies
  from each rookery
}
\note{gibbsC calls C code to generate multinomial deviates
  derived from the randlib.c library
  (version 1.3: currently located at \url{http://lcrocker.github.io/onejoker/randlib/})
  written by Barry W. Brown, James Lovato, Kathy Russell,
  and John Venier, derived in turn from page 559 of: Devroye, Luc,
  Non-Uniform Random Variate Generation.  Springer-Verlag, New York,
  1986.
}
\author{Ben Bolker}
\keyword{misc}
