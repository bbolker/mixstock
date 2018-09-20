\name{q.to.p}
\alias{q.to.p}
\alias{p.to.q}
\title{Real-to-multifrequency transformation}
\description{
  Transform a vector of n real-valued variables in (-Inf,Inf) [or (0,1)]
  to a vector of n+1 variables in (0,1) that sum to 1, or vice versa.
}
\usage{
q.to.p(q, transf="full")
p.to.q(p, transf="full")
}
\arguments{
  \item{q}{Unconstrained/transformed values:
    vector of n numeric values in (-Inf,Inf) [if \code{transf="full"}]
  or (0,1) [if \code{transf="part"}]}
  \item{p}{Vector of n+1 numeric values in (0,1) that sum to 1}
  \item{transf}{(character) "full": use arctan transform to transform (-Inf,Inf)
    to (0,1) or vice versa; "part": don't; "none"; no transform}
}
\details{
  Essentially, this is a
  transformation from an unconstrained set of variables to a bounded,
  constrained set of variables.  If contin is TRUE, an arctan
  transformation (\code{v <-> atan(v)/pi+0.5}) is used to transform
  (-Inf,Inf) to (0,1) or vice versa.  In either case, the correlated
  set of variables (which sum to 1) is transformed to an unconstrained
  set by taking each variable to be a remainder: \code{x[1]=x[1]},
  \code{x[2]=x[2]/(1-x[1])}, and so forth.
}
\value{
  Vector of transformed values.
}
\author{Ben Bolker}
\note{This transformation is designed to deal with the problems of
  bounded optimization and constraints.  It actually behaves quite badly
  because small values are transformed to large negative values, messing
  up the uniform scaling of the parameters.  Now that the bounded
  optimization of \code{optim} has improved, \code{contin="full"} may not
  be a good idea.  It's not clear whether the other transformation
  (remainders) is better or worse than just optimizing on the first
  (n-1) components and assuming that the last frequency equals one minus
  the sum of the rest.
}
\examples{
  p.to.q(c(0.3,0.3,0.4))
  p.to.q(c(0.3,0.3,0.4),transf="part")
  q.to.p(c(-4,3))
  q.to.p(c(0,0))
  q.to.p(c(0.5,0.5),transf="part")
}
\keyword{misc}
