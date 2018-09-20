\name{coef.toptim}
\alias{coef.toptim}
\alias{logLik.toptim}
\title{minimal structure for optim fits}
\description{
  provides coef() and logLik() methods
  for optim fits
}
\usage{
\method{coef}{toptim}(object,\dots)
\method{logLik}{toptim}(object,\dots)
}
\arguments{
  \item{object}{the result of an optim fit}
  \item{\dots}{optional arguments (for generic compatibility)}
}
\details{
  extracts par and -val respectively
}
\value{
  coefficients or value
}
\note{mle is not flexible enough to take vectors as inputs ...}
\keyword{misc}
