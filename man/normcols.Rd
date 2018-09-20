\name{normcols}
\alias{normcols}
\title{scale columns to frequencies}
\description{
  divides each column of a matrix by its sum, reducing
  sample values to observed frequencies
}
\usage{
normcols(z)
}
\arguments{
  \item{z}{matrix of sampled numbers (integers)}
}
\value{
  matrix of frequencies
}
\note{This is just a wrapper for \code{scale(z, center = FALSE, scale = colSums(z))}}
\keyword{misc}
