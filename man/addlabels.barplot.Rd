\name{addlabels.barplot}
\alias{addlabels.barplot}
\title{add labels to a barplot}
\description{
  Adds labels (at specified heights) to an existing barplot
}
\usage{
addlabels.barplot(x, vals, names, min = 0.1, cols = par("fg"), horiz = FALSE, ...)
}
\arguments{
  \item{x}{x positions (may be derived from a call to barplot)}
  \item{vals}{heights of labels}
  \item{names}{label text}
  \item{min}{minimum size for adding labels}
  \item{cols}{colors of label text}
  \item{horiz}{horizontal barplot}
  \item{\dots}{additional arguments to text()}
}
\value{
  none.
}
\author{Ben Bolker}
\examples{
set.seed(1001)
bvals <- matrix(runif(12),nrow=3)
b <- barplot(bvals)
addlabels.barplot(b,bvals,LETTERS[1:12])
}
\keyword{aplot}
