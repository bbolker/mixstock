\name{simmixstock2}
\alias{simmixstock2}
\title{Simulate multiple mixed stocks}
\description{
  Simulated multiple mixed stocks
}
\usage{
simmixstock2(sourcefreq, destmat, sourcesize,
 sourcesampsize, mixsampsize,
             nmark, nsource, nmix, rseed = NULL, condense=TRUE)
}
\arguments{
  \item{sourcefreq}{matrix of marker frequencies in the source populations}
  \item{destmat}{matrix of contributions of sources to mixed stocks}
  \item{sourcesize}{vector of relative sizes of sources: default is
    equal sizes for all sources}
  \item{sourcesampsize}{sizes of samples from source populations}
  \item{mixsampsize}{sizes of samples from mixed populations}
  \item{nmark}{number of distinct markers}
  \item{nsource}{number of source populations}
  \item{nmix}{number of mixed populations}
  \item{rseed}{random number seed}
  \item{condense}{automatically use \code{markfreq.condense} to condense
    result?}
}
\details{If \code{sourcefreq} and \code{destmat} are
  specified, computes expected marker frequencies in mixed
  stocks and simulates accordingly.  If \code{sourcefreq}
  and/or \code{destmat} are missing, they are generated
  randomly with uniform probabilities.
  }
\value{
  Returns an object of type \code{mixstock.data}
}
\author{Ben Bolker}
\keyword{misc}
