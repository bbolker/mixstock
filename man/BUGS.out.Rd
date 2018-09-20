\name{BUGS.out}
\alias{BUGS.out}
\title{generate MCMC output in BUGS format}
\description{
  takes an MCMC object produced by \code{\link{tmcmc}} or
  \code{\link{gibbsC}} and sends it to an output file in BUGS output
  format, so that it can be read by CODA and other BUGS post-processors
}
\usage{
BUGS.out(g, file=NULL)
}
\arguments{
  \item{g}{MCMC output as produced by tmcmc}
  \item{file}{output file: if NULL, uses the name of the object}
}
\value{
  Sends output to a file.
}
\author{Ben Bolker}
\keyword{misc}
