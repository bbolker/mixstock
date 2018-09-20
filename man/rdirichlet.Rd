\name{rdirichlet}
\alias{rdirichlet}
\title{Dirichlet deviates}
\description{
  Produces random deviates for
  the Dirichlet distribution, the multivariate analogue of the beta
  distribution.
}
\usage{
rdirichlet(n, alpha)
}
\arguments{
  \item{n}{Number of random deviates to generate.}
  \item{alpha}{Vector of shape parameters.}
}
\details{
  The Dirichlet distribution is a multivariate distribution that describes
  distributions of frequencies.  If a multinomial sample \eqn{n_1, \ldots n_N}
  is taken from a population, the estimated distribution of 
  frequencies is Dirichlet(\eqn{n_1+1, \ldots n_2+N}).  If \eqn{\{g_i\}} is
  a set of Gamma deviates with shape parameters \eqn{\{\alpha_i\}}, then
  \eqn{\{g_i\}/sum(\{g_i\})} is a Dirichlet deviate.

}
\value{
  A matrix of random deviates
  (each in a distinct row) from \code{rdirichlet}
}
\author{Ben Bolker}
\examples{
  rdirichlet(5,c(7,4,4))
 }

\keyword{distribution}
