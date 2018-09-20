\name{genboot}
\alias{genboot}
\title{Generate bootstrap estimates of mixed stock analyses}
\description{
  Given a data set or simulated data set (in the usual format of a list
  of source samples and samples from the mixed stock), generate
  bootstrap-resampled data sets (either parametrically or
  nonparametrically) and a vector of estimates from the resampled data sets
}
\usage{
genboot(x, method="cml", nboot=1000, rseed=1001, verbose=FALSE,
fuzz=0.001, maxfail=100, pb=TRUE, start.type="lsolve",
param=FALSE,param.match="mean",ndepfac=10000,save.boot=FALSE,
print.boot=FALSE,...)
}
\arguments{
  \item{x}{A mixed stock data set: a list with elements \code{sourcesamp} (a
    matrix, with sources as columns and markers as rows, of samples
  in the sources) and \code{mixsamp} (a vector of markers sampled
  in the mixed stock)}
  \item{method}{"cml" for conditional max. likelihood or "uml" for
    unconditional max. likelihood}
  \item{nboot}{number of bootstrap samples}
  \item{rseed}{random-number seed}
  \item{verbose}{produce verbose output?}
  \item{fuzz}{small value for keeping estimated frequencies away from 0/1}
  \item{maxfail}{number of consecutive fitting tries before entering NAs
  for a given simulation}
  \item{pb}{(logical) draw progress bar?}
  \item{start.type}{starting conditions for fitting procedure (see
    \code{\link{startvec}} and/or \code{\link{startvec0}})}
  \item{param}{Do parametric bootstrapping?}
  \item{param.match}{match mean or mode of resampled frequencies?}
  \item{ndepfac}{scaling factor for CML fit}
  \item{save.boot}{save bootstrap replicates?}
  \item{print.boot}{verbose output from bootstrap?}
  \item{...}{additional arguments to \code{\link{cml}}}
}
\value{
  An object of type \code{mixstock.est} with element \code{resample}
  as a matrix of bootstrap results (where the columns
  are the parameters, the negative log-likelihood of the fit, and a code
  for convergence).
}
\examples{
  data(simex)
  x <- genboot(simex,method="cml",nboot=100)
  r <- x$resample
  r.ok <- r[!is.na(r[,"Convergence"]) & r[,"Convergence"]==0,]
  old.par <- par(no.readonly=TRUE)
  par(mfrow=c(2,2))
  hist(r.ok[,1],main="Contrib. A")
  hist(r.ok[,2],main="Contrib. B")
  hist(r.ok[,3],main="Neg. log likelihoods")
  plot(r.ok[,1],r.ok[,3],
        xlab="Contrib. A",ylab="NLL")
  par(old.par)
}
\keyword{misc}
