\name{expand.bugs.data}
\alias{expand.bugs.data}
\title{convert data to BUGS input format}
\description{
  expands mixed stock data to BUGS input format (indicator variables)
}
\usage{
expand.bugs.data(mixsamp)
}
\arguments{
  \item{mixsamp}{mixed-sample or mixed stock data matrix (numbers from
    each source population)}
}
\value{
  a matrix of indicator (0/1) variables: number of rows = total number of
  individuals, number of columns = number of source populations
}
\author{Ben Bolker}
\examples{
expand.bugs.data(c(A=4,B=2,C=5))
}
\keyword{misc}
