\name{numderiv}
\alias{numderiv}
\title{take the numeric derivative of a function}
\description{
  takes the numeric derivative (gradient) of a function by
  finite differencing
}
\usage{
numderiv(p, fn, eps = 1e-05, ...)
}
\arguments{
  \item{p}{parameter vector}
  \item{fn}{function}
  \item{eps}{finite-difference increment}
  \item{\dots}{other arguments to \code{fn}}
}
\value{
  numeric vector representing the approximate gradient of the
  function with respect to the parameters
}
\seealso{\code{\link{numericDeriv}}}
\keyword{misc}
