\name{packval}
\alias{packval}
\alias{packval2}
\alias{unpackval}
\alias{unpackval2}
\title{Pack and unpack contribution and marker frequencies}
  \description{
    Packs/unpacks source contributions and source marker frequencies into
  a single vector, after transforming them to an uncorrelated (and
  possibly unbounded) set of parameters.
}
\usage{
packval(f, r, transf="part")
packval2(x,R=2,H=3)
unpackval(p,R=2,H=2,x.orig=NULL,transf="full",input.only=FALSE)
unpackval2(p,R=2,H=2,transf="full",x.orig=NULL)
}
\arguments{
  \item{f}{Source contributions (numeric vector, all between 0 and 1,
    sum equals 1)}
  \item{r}{Source marker frequencies: matrix of marker frequencies
    (row=marker, column=source)}
  \item{x}{a numeric vector containing (elements 1 to R) source
    contributions and (elements R+1 to (R+R*H)) marker frequencies in
    sources}
  \item{p}{a packed/transformed parameter vector}
  \item{transf}{Transform to unbounded variables? ("full","part","none") (See
    \code{\link{p.to.q}} for description/warnings.)}
  \item{R}{number of sources}
  \item{H}{number of markers}
  \item{x.orig}{original data (for extracting source/marker names)}
  \item{input.only}{return only contribution parameters?}
}
\value{
  \code{packval} packs source contributions and marker frequencies
  specified as a separate vector and a matrix: \code{packval2} packs
  them (i.e. transforms them) when they have already been run together
  as a vector. Either produces a numeric vector of the transformed
  values, with length (R-1+R*(H-1)).  \code{unpackval} and
  \code{unpackval2} invert the operation, producing a list
  \item{input.freq}{source contributions}
  \item{source.freq}{marker frequencies in sources}
  or a vector respectively.
}
\author{Ben Bolker}
\seealso{\code{\link{p.to.q}}}
\examples{
data(simex)
sourcefreq <- sweep(simex$sourcesamp,2,apply(simex$sourcesamp,2,sum),"/")
packval(c(0.2,0.8),sourcefreq)
packval(c(0.2,0.8),sourcefreq,transf="full")
}
\keyword{misc}
