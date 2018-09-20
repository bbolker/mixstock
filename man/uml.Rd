\name{uml}
\alias{uml}
\alias{uml.ds}
\alias{uml.em}
\title{Mixed stock analysis by unconditional maximum likelihood}
\description{
  Find the unconditional maximum likelihood estimate (jointly estimating
  marker frequencies in sources) of the contributions of different
  sources to a mixed stock, by either a direct-search or an
  expectation-maximization method
}
\usage{
uml(x, method="direct",optmethod="L-BFGS-B",...)
uml.ds(x,grad=NULL,start.type="lsolve",fuzz=0,bounds=1e-4,
ndepfac=1000,method="L-BFGS-B",debug=FALSE,control=NULL,
transf=c("part","full","none"),...)
uml.em(x,prec=1e-8,prior=1)
}
\arguments{
  \item{x}{a list with elements \code{mixsamp} (a vector of the sampled
  markers in the mixed stock) and \code{sourcesamp} (a matrix,
  with markers in rows and sources in columns, of markers in the
  source samples)}
\item{optmethod}{to be passed to optim}
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
  \item{transf}{transformation}
  \item{debug}{produce debugging output?}
  \item{control}{other control arguments to \code{\link{optim}}}
  \item{...}{other arguments to \code{mle} or \code{optim}
    (e.g. \code{hessian=FALSE} to suppress (slow) hessian calculation,
    etc.)}
  \item{prec}{precision for determining convergence of EM algorithm}
  \item{prior}{prior for EM algorithm}
}
\details{\code{uml} uses either a direct-search algorithm or an EM
  algorithm to find the ML estimate}
\value{
  an object of class \code{mixstock.est}, with elements
  \item{fit}{information on the ML fit}
  \item{resample}{bootstrap information, if any}
  \item{data}{original data used for estimate}
  \item{R}{number of sources}
  \item{H}{number of markers}
  \item{contin}{estimation done on transformed proportions?}
  \item{method}{optimization method}
  \item{boot.method}{resampling method}
  \item{boot.data}{raw resampling information}
  \item{gandr.diag}{Gelman-Rubin diagnostic information for MCMC estimates}
  \item{prior}{Prior for MCMC estimates}
  \item{em}{estimation done by EM algorithm?}
}
\author{Ben Bolker}
\examples{
true.freq <- matrix(c(0.65,0.33,0.01,0.01,
                      0.33,0.65,0.01,0.01),ncol=2)
true.contrib <- c(0.9,0.1)
x <- simmixstock0(true.freq,true.contrib,50,100,1004)
uml.est <- uml(x)
uml.est
uml.emest <- uml.em(x)
uml.emest
}
\keyword{misc}
