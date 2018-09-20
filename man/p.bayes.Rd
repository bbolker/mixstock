\name{p.bayes}
\alias{p.bayes}
\title{Empirical Bayes estimate of prior strength for source marker frequencies}
\description{
  Uses an iterative calculation to estimate an appropriate strength
  for the flat Dirichlet prior probability for
  source marker frequencies
}
\usage{
p.bayes(sourcesamp, bold = 1, cut = 1e-04)
}
\arguments{
  \item{sourcesamp}{matrix of source samples}
  \item{bold}{starting value for iterations (??)}
  \item{cut}{cutoff value for iterations}
}
\details{
  See Masuda and Pella's paper
}
\value{
  prior strength parameter (a)
}
\keyword{misc}
