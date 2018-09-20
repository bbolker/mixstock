\name{simmixstock0}
\alias{simmixstock0}
\title{Simulate a sample from sources and mixed stock}
\description{
  Given information on marker frequencies in sources, the
  contributions of each source to a mixed stock, and sample
  sizes, generates multinomial samples of markers from sources and
  the mixed stock.
}
\usage{
simmixstock0(sourcefreq, input.freq, sourcesampsize, mixsampsize, rseed=NULL)
}
\arguments{
  \item{sourcefreq}{Matrix (markers=rows, source=cols) of true marker
    frequencies in sources}
  \item{input.freq}{True contributions from each source to the mixed stock}
  \item{sourcesampsize}{Sample size from (each) source}
  \item{mixsampsize}{Sample size from mixed stock}
  \item{rseed}{Random-number seed}
}
\value{
Multinomial sample.
  \item{mixsamp}{Sample from mixed stock (vector)}
  \item{sourcesamp}{Sample from sources (matrix, markers=rows, source=cols)}
  ...
}
\author{Ben Bolker}

\examples{
true.freq <- matrix(c(0.65,0.33,0.01,0.01,
                      0.33,0.65,0.01,0.01),ncol=2)
true.contrib <- c(0.9,0.1)
x <- simmixstock0(true.freq,true.contrib,sourcesampsize=50,mixsampsize=100,rseed=1004)
}
\keyword{misc}
