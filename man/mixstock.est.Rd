\name{mixstock.est}
\alias{mixstock.est}
\alias{print.mixstock.est}
\alias{xyplot.mixstock.est}
\alias{coeftab.mixstock.est}
\title{construct (or print) a mixed stock analysis estimate object}
\description{
  combine a variety of data into a list with class "mixstock.est"
}
\usage{
mixstock.est(fit, resample = NULL,
resample.index=NULL, data = NULL, em = FALSE, sourcesamp =
NULL, mixsamp = NULL, R = NULL, H = NULL, M=1, transf = "full", method =
"unknown", boot.method = "none", boot.data = NULL, GR = NULL, prior =
NULL)
\method{print}{mixstock.est}(x,debug=FALSE,\dots)
\method{xyplot}{mixstock.est}(x,data,\dots)
\method{coeftab}{mixstock.est}(object,clevel=c(0.5,0.95),\dots)
}
\arguments{
  \item{fit}{a mixed stock fit (from cml, uml, etc.)}
  \item{object}{a mixed stock fit}
  \item{resample}{resampling information}
  \item{resample.index}{chain index for resamples}
  \item{data}{original data (ignored in call to \code{xyplot})}
  \item{em}{(logical) was estimation done by EM algorithm?}
  \item{sourcesamp}{source marker samples}
  \item{mixsamp}{mixed stock samples}
  \item{R}{number of sources}
  \item{H}{number of markers}
  \item{M}{number of mixed stocks}
  \item{transf}{transformation of parameters}
  \item{method}{estimation method (cml,uml,mcmc)}
  \item{boot.method}{bootstrap method -- nonpar, parametric, or mcmc}
  \item{boot.data}{bootstrap data}
  \item{GR}{Gelman-Rubin diagnostic results}
  \item{prior}{prior strength}
  \item{x}{a \code{mixstock.est} object}
  \item{debug}{debug?}
  \item{clevel}{confidence regions to plot}
  \item{\dots}{other arguments}
}
\value{
  an object of type \code{mixstock.est}
}
\keyword{misc}

