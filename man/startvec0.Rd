\name{startvec0}
\alias{startvec0}
\alias{startvec}
\alias{lsolve}
\title{Produce (raw) starting vector of parameters for mixed stock analysis}
\description{
  Provides raw (untransformed) starting vector of source contribution
  parameters for mixed stock analysis.
}
\usage{
startvec0(sourcesamp, mixsamp=NULL, type="equal", sd=1, lmin=0.001)
startvec(sourcesamp, mixsamp=NULL, type="equal", 
         marktype="sample",a=1,cond=FALSE,transf="full",
fuzz=0,sd=1)
lsolve(n, s, tol = 1e-05, warn = FALSE)
}
\arguments{
  \item{sourcesamp}{Marker frequencies in sources, or a list with
    source and mixed samples, or a \code{turtle.data} object.}
  \item{mixsamp}{Marker frequencies in mixed stock}
  \item{type}{Various options for setting starting contributions.
    \code{equal}: equal contributions from all sources. \code{random}:
    random multinomial sample with equal multinomial
    probabilities. \code{rand2}: random sample from transformed normal
    variates with mean 0 and standard deviation \code{sd}.
    A number \code{n} between 1 and the number of sources
    inclusive gives a starting condition with 95\% of the contribution
    from source \code{n} and the other 5\% evenly split between the
    other sources.  The default is to attempt a solution of the linear
    equation (\code{sourcesamp*f=mixsamp}) and use these values as the
    starting contributions.}
  \item{sd}{standard deviations for starting type \code{rand2}}
  \item{lmin}{When doing linear solutions, the boundary values are
    \{\code{lmin}, \code{1-lmin}.}
  \item{marktype}{method for starting marker frequencies:
    \code{sample} uses the observed sample probabilities; \code{random}
    used observed sample probabilities; \code{weighted} does a
    Bayes-weighted start (a la Masuda and Pella)}
  \item{cond}{Conditional likelihood?  (i.e. \code{cond=TRUE} gives just
    the parameters for source contributions, not parameters for source
    marker frequencies)}
  \item{a}{prior strength parameter}
  \item{transf}{transform ("full","part", or "none")}
  \item{fuzz}{fuzz parameter for moving parameters away from the
    boundary}
  \item{n}{source samples (HxR matrix)}
  \item{s}{mixed stock samples (vector)}
  \item{tol}{tolerance for linear fit}
  \item{warn}{warn if numeric problems with solution}
}
\value{
  A parameter vector of the contributions from each of the sources:
  just the raw source contributions in the case of \code{startvec0}, or
  source contributions and possibly marker frequencies (transformed
  or untransformed) in the case of \code{startvec}.
}
\note{\code{lsolve} attempts to get a starting value by solving
  the linear equation (\code{solve(n,s,tol=tol)})}
\author{Ben Bolker}
\examples{
 data(simex)
 startvec0(simex)
 startvec(simex,transf="part")
 startvec(simex)
}
\keyword{misc}
