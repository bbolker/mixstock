\name{blockdiag}
\alias{blockdiag}
\title{Block diagonal matrices}
\description{
   Construct a block-diagonal matrix from a list of submatrices.
}
\usage{
blockdiag(...)
}
\arguments{
  \item{\dots}{a comma-separated list of matrix objects that
will form the submatrices of the block-diagonal matrix.
}}
\value{
   A block-diagonal matrix with dimensions equal to the sum of
the dimensions of the submatrices.
}
\author{Ben Bolker}
\note{If you want to use a previously constructed list of matrices,
use \code{\link{do.call}} (see the examples).
Doesn't preserve names or attributes (yet).}
\examples{
  m1 <- matrix(1,nrow=2,ncol=2)
  m2 <- matrix(2,nrow=1,ncol=3)
  m3 <- matrix(3,nrow=3,ncol=2)
  blockdiag(m1,m2,m3)
  l <- list(m1,m2,m3)
  do.call("blockdiag",l)
}
\keyword{array}
