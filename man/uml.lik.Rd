\name{uml.lik}
\alias{uml.lik}
\alias{uml.grad}
\alias{cml.lik}
\alias{cml.grad}
\alias{dcmat.a}
\title{Likelihoods and gradients of stock analysis models}
\description{
  Calculate negative log likelihoods and gradients
  for unconditional and conditional models
}
\usage{
uml.lik(p, data, transf=c("full","part","none"), verbose=FALSE,
debug=FALSE)
cml.lik(p, sourcefreq, data, transf=c("full","part","none"),
           verbose=FALSE, fulllik=TRUE, debug=FALSE)
cml.grad(p, sourcefreq, data, transf="full",
           verbose=FALSE,fulllik=NULL,debug=FALSE)
uml.grad(p,  data, transf="full", debug=FALSE, verbose=FALSE)
dcmat.a(x,debug=FALSE)
}
\arguments{
  \item{p}{a vector of parameters.}
  \item{data}{a data set in \code{\link{mixstock.data}} format}
  \item{sourcefreq}{source frequencies}
  \item{transf}{how are parameters transformed?}
  \item{verbose}{print messages?}
  \item{debug}{debug?}
  \item{x}{vector of parameters}
  \item{fulllik}{for CML, give likelihood corresponding to source
    samples (test only)?}
}
\details{
  The log likelihood is the log multinomial likelihood of the mixed population
  samples (\code{data$mixsamp}) given the expected frequencies in the
  mixed population, which are computed from the contributions and the
  source marker frequencies, plus the log multinomial likelihoods of
  the samples in each source given the marker frequencies specified
  for each source. \code{dcmat.a} is a utility function for the
  gradient calculations.
}
\value{
  Negative log likelihood, possibly plus a constant corresponding to the
  normalization factor
}
\author{Ben Bolker}
\examples{
data(simex)
rfreq <- normcols(simex$sourcesamp)
tmpf <- function(p) {
  uml.lik(c(p,1-p,rfreq),simex,transf="none")
}
pvec <- seq(0.01,0.99,by=0.01)
plot(pvec,sapply(pvec,tmpf))
}
\keyword{misc}
