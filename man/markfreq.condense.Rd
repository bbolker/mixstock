\name{markfreq.condense}
\alias{markfreq.condense}
\title{condense marker frequency list}
\description{
  "condenses" standard marker-frequency-list data (marker samples
  from sources and mixed stocks) by removing unrepresented
  markers, lumping markers where that removes no information, etc.}
\usage{
markfreq.condense(sourcesamp=NULL, mixsamp=NULL, debug=FALSE, exclude.nomix=FALSE)
}
\arguments{
  \item{sourcesamp}{Turtle data in standard form: list of sourcesamp (see below)
    and mixsamp OR matrix of source genotype samples (sources=columns,
    markers=rows)}
  \item{mixsamp}{Vector of mixed stock marker samples}
  \item{debug}{Print out debugging information?}
  \item{exclude.nomix}{Exclude markers not found in the mixed stock?}
}
\details{
  Criteria:
  \itemize{
    \item any marker found in the mixed stock but not in any of the sources is
    deleted from both: it provides no further information about the source contributions
    (although it does contribute to goodness-of-fit tests)
  \item any set of markers found in only one source (and possibly
    also in the mixed stock, although not necessarily) is lumped
    together
  \item "exclude.nomix" determines whether to exclude markers not found
    in the mixed stock
  }
}
\value{
  "Standard" mixed stock analysis data set, suitably condensed.
  \item{sourcesamp}{source samples}
  \item{mixsamp}{mixed stock samples}
  \item{err}{was there an error in attempting to lump?}
}
\author{Ben Bolker}
\examples{
data(lahanas98raw,package="mixstock")
lahan98c <- markfreq.condense(as.mixstock.data(lahanas98raw))
}
\keyword{misc}
