\name{plot.mixstock.est}
\alias{plot.mixstock.est}
\title{plot mixed stock analysis estimates}
\description{
  plots a mixed stock analysis estimate as a point plot
  with 95\% confidence limites
}
\usage{
\method{plot}{mixstock.est}(x, plot.freqs = FALSE,
   sourcectr=FALSE, contrib.lab = "Estimated source contributions",
   sourcefreq.lab = "Estimated source marker freqs",
   markcolors = rainbow(x$H),
   alength=0.25,
   aunits="inches",
   abbrev, level=0.95,
   axes=TRUE, \dots)
}
\arguments{
  \item{x}{mixstock estimate}
  \item{plot.freqs}{plot marker frequency estimates?}
  \item{sourcectr}{plot source-centric estimates?}
  \item{contrib.lab}{label for source contribution plot}
  \item{sourcefreq.lab}{label for marker frequency plot}
  \item{markcolors}{colors corresponding to markers}
  \item{alength}{length of error bar ends in many-to-many plots}
  \item{aunits}{units of error bar end lengths in many-to-many plots}
  \item{abbrev}{abbreviate names on horizontal axis in many-to-many
    plots? Default if TRUE is 3 characters, but may also be an integer
    specifying the number of characters}
  \item{level}{confidence level}
  \item{axes}{(logical) plot axes?}
  \item{\dots}{other arguments to \code{\link{barplot}}}
}
\note{
  The \code{\dots} argument can contain a great number of
  optional arguments to \code{barplot}: see \code{\link{barplot}}
  and \code{xyplot} (in the \code{lattice} package).  Among many others, one can specify
  (for example) \code{scale=list(x=list(cex=0.6),y=list(log=TRUE))}
  for smaller labels on the horizontal axis and a log scale on
  the vertical axis; \code{layout=c(3,5)} to change the number
  of rows and columns of panels; or \code{as.table=TRUE} to change
  the ordering of the panel to top-to-bottom, left-to-right.
  (One can also use \code{abbreviate=TRUE} in the \code{scale}
  list, but the \code{abbrev} argument above seems to be more powerful.)
  }
\value{
  Produces a plot on the current graphics device.
}
\examples{
  data(simex)
  u1 = uml(simex)
  plot(u1)
  plot(u1,plot.freqs=TRUE)
}
\keyword{misc}
