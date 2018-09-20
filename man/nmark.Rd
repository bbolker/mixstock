\name{nmark}
\alias{nmark}
\alias{nmark<-}
\alias{nmix}
\alias{nmix<-}
\alias{nsource}
\alias{nsource<-}
\title{Set and query mixed stock parameters}
\description{
  Sets or queries the number of markers,
  mixed stocks, or sources, for a given mixed-stock
  data set or estimate
}
\usage{
nmark(object)
nmark(object) <- value
nmix(object)
nmix(object) <- value
nsource(object)
nsource(object) <- value
}
\arguments{
  \item{object}{a mixed-stock data object or fit}
  \item{value}{an integer}
}
\value{
  Returns the numeric value or sets the value
}
\author{Ben Bolker}
\keyword{misc}
