\name{as.mcmc.bugs}
\alias{as.mcmc.bugs}
\title{convert from R2WinBUGS to CODA format}
\description{
  converts output of R2WinBUGS to standard CODA format
}
\usage{
\method{as.mcmc}{bugs}(x)
}
\arguments{
  \item{x}{an object of class \code{bugs},
    from the \code{bugs} function of R2WinBUGS}
}
\value{
  an item of class \code{mcmc} or \code{mcmc.list}
}
\author{Ben Bolker}
\keyword{misc}
