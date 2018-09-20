\name{cml}
\alias{MAXMARK}
\alias{MAXSRC}
\alias{cml}
\title{Mixed stock analysis by conditional maximum likelihood}
\description{
  Find the conditional maximum likelihood estimate (assuming marker
  frequencies in the sources are exactly equal to the sample
  frequencies) of the contributions of different
  sources to a mixed stock
}
\usage{
cml(x, start.type="lsolve", fuzz=0, bounds=1e-4,
ndepfac=1000,method="L-BFGS-B",lower=NULL,upper=NULL,
ndeps=NULL,
control=NULL,debug=FALSE,transf="part",grad=NULL,...)
MAXMARK
MAXSRC
}
\arguments{
  \item{x}{a list with elements \code{mixsamp} (a vector of the sampled
  markers in the mixed stock) and \code{sourcesamp} (a matrix,
  with markers in rows and sources in columns, of markers in the
  source samples)}
  \item{grad}{function giving the gradient of the likelhood}
  \item{start.type}{starting values to use: \code{equal} (equal
    contributions from each source); \code{random} (multinomial sample
    with equal probabilities); \code{rand2} (sample out of a transformed
    normal distribution); a number between 1 and the number of sources;
    that source starts with 0.95 contribution and the rest start with
    0.05/(R-1); default \code{lsolve}, the linear solution to the problem}
  \item{fuzz}{min. value (1-min is the max.) for starting contributions}
  \item{bounds}{(bounds,1-bounds) are the lower and upper bounds for mle
    calculations}
  \item{ndepfac}{factor for computing numerical derivatives; numerical
    derivative stepsize is computed as bounds/ndepfac [OBSOLETE with
    gradient function?]}
  \item{method}{optimization method, to be passed to \code{\link{optim}}}
  \item{debug}{produce debugging output?}
  \item{lower}{lower bound}
  \item{upper}{upper bound}
  \item{ndeps}{scaling factor for optimization}
  \item{control}{other control arguments to \code{\link{optim}}}
  \item{transf}{(character) "full": use arctan transform to transform (-Inf,Inf)
    to (0,1) or vice versa; "part": don't; "none"; no transform}
  \item{...}{other arguments to \code{mle} or \code{optim}
    (e.g. \code{hessian=FALSE} to suppress (slow) hessian calculation, etc.)}
}
\details{By default, uses \code{mle} which in turn uses \code{optim}
  with \code{method="L-BFGS-B"} to do bounded optimization}
\note{There is a hard-coded limit of 500 markers and 100 sources;
  contact the package maintainer if this is a problem}
\value{
 an object of class \code{\link{mixstock.est}}, containing the results of
 the fit
}
\author{Ben Bolker}
\examples{
true.freq <- matrix(c(0.65,0.33,0.01,0.01,
                      0.33,0.65,0.01,0.01),ncol=2)
true.contrib <- c(0.9,0.1)
x <- simmixstock0(true.freq,true.contrib,50,100,1004)
cml.est <- cml(x)
cml.est
}
\keyword{misc}


