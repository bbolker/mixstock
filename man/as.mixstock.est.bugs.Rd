\name{as.mixstock.est.bugs}
\alias{as.mixstock.est.bugs}
\alias{as.mixstock.est.mcmc}
\title{Convert BUGS or CODA objects to mixstock estimate object}
\description{
  Converts output of a R2WinBUGS run,
  in bugs format or CODA (mcmc or mcmc.list) format, to a mixed-stock estimate object
}
\usage{
as.mixstock.est.bugs(object,data=NULL,time)
}
\arguments{
  \item{object}{results of a BUGS run (from \code{pm.wbugs} etc.)}
  \item{data}{a \code{mixstock.data} object for filling in source and
    mixture names, etc.}
  \item{time}{time taken to run WinBUGS, in seconds}
}
\value{
  an object of type \code{mixstock.est}
}
\note{
  \code{as.mixstock.est.mcmc} is a stub; for now, you can
  convert R2WinBUGS output to CODA or mixstock, but you
  can't do any other conversions
  }
  \author{Ben Bolker}
\keyword{misc}
