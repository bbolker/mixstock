\name{loglik2.C}
\alias{loglik2.C}
\title{likelihood calculation, in C}
\description{
  takes a parameter vector corresponding to the (possibly transformed)
  source contributions and source marker frequencies and returns
  the negative log likelihood
}
\usage{
loglik2.C(p, sourcesamp, mixsamp = NULL, verbose = FALSE,
transf = c("full", "part", "none"), full = FALSE, cond = FALSE, debug = FALSE)
}
\arguments{
  \item{p}{parameter vector}
  \item{sourcesamp}{source marker samples}
  \item{mixsamp}{mixed stock samples}
  \item{verbose}{verbose output?}
  \item{transf}{transformation: "full", "part", or "none"}
  \item{full}{full likelihood (including normalization constants?)}
  \item{cond}{conditional likelihood (not including source sample likelihood?}
  \item{debug}{debugging output?}
}
\value{
  numeric value of the negative log-likelihood
}
\keyword{misc}
