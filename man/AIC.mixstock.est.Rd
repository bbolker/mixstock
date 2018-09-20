\name{AIC.mixstock.est}
\alias{AIC.mixstock.est}
\alias{boxplot.mixstock.est}
\alias{as.mixstock.est}
\alias{coef.mixstock.est}
\alias{logLik.mixstock.est}
\alias{hist.mixstock.est}
\alias{confint.mixstock.est}
\alias{summary.mixstock.est}
\title{commands for mixed stock analysis estimates}
\description{
  calculate AIC, confidence intervals, etc., for mixed
  stock analysis estimates
}
\usage{
\method{AIC}{mixstock.est}(object, \dots)
\method{boxplot}{mixstock.est}(x,\dots)
\method{coef}{mixstock.est}(object,\dots)
\method{logLik}{mixstock.est}(object,\dots)
\method{confint}{mixstock.est}(object,parm,level=0.95,profile=FALSE,type=c("quantile","credible"),
show.sourcectr=TRUE,show.haps=FALSE,\dots)
\method{summary}{mixstock.est}(object,show.data=TRUE,\dots)
as.mixstock.est(object)
}
\arguments{
  \item{object}{mixed stock analysis estimate}
  \item{x}{mixed stock analysis estimate}
  \item{parm}{parameter}
  \item{level}{confidence level}
  \item{profile}{confidence intervals by likelihood profiling? (stub)}
  \item{type}{for MCMC runs, produce credible intervals or
    quantiles of chains?}
  \item{show.sourcectr}{show confidence intervals for source-centric
    estimates?}
  \item{show.haps}{show confidence limits for marker frequencies?}
  \item{show.data}{print data in summary?}
  \item{\dots}{additional arguments}
}
\value{
  Various summary statistics retrieved from mixed stock analysis estimates
}
\author{Ben Bolker}
\note{AIC and logLik return statistical values etc. ...
  many of these are stubs}
\keyword{misc}
