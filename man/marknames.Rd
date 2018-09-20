\name{marknames}
\alias{marknames}
\alias{sourcenames}
\alias{mixnames}
\alias{locnames}
\alias{mixstock.dimnames}
\alias{label.mixstock.data}
\title{Extract/assign names from mixed stock data}
\description{
  Extract names of markers and sources from mixed stock data, or invent
  names for simulated data
}
\usage{
marknames(x)
sourcenames(x)
mixnames(x)
locnames(x)
mixstock.dimnames(H,R)
label.mixstock.data(x,sourcenames=NULL,marknames=NULL,mixnames=NULL)
}
\arguments{
  \item{x}{mixstock.data object}
  \item{H}{number of markers}
  \item{R}{number of sources}
  \item{sourcenames}{vector of source names}
  \item{marknames}{vector of marker names}
  \item{mixnames}{vector of mixed stock names}
}
\details{
  mixstock.dimnames uses Roman numerals for markers and capital letters
  for sources.  label.mixstock.data assigns names to the data if
  the source or mixed components of the data set already have
  names, or use the same rules as mixstock.dimnames.
  \code{locnames} gives a combined vector of source and mixed stock names.
}
\value{
  character vectors of names, or a list with marker and source names
}
\keyword{misc}
