\name{simex}
\title{A small sample mixed stock analysis data set}
\alias{simex}
\usage{data(simex)}
\description{This data set is just a small example, generated with four
  markers (haplotypes) and two sources (rookeries: four markers later condensed with
  \code{markfreq.condense} to three), using true marker frequencies
  \{0.65,0.33,0.01,0.1\} in rookery A and \{0.33,0.65,0.01,0.1\} in
  rookery B; and true contributions of 90\% from rookery A and 10\% from
  rookery B.  In other words, it was generated with the
  commands
  \code{true.freq <- matrix(c(0.65, 0.33, 0.01, 0.01, 0.33, 0.65, 0.01,
    0.01), ncol = 2)}
  \code{true.contrib <- c(0.9, 0.1)}
  \code{simmixstock0(true.freq, true.contrib, sourcesampsize = 50, 
    mixsampsize = 100, rseed = 1004)}
  }
\format{an object, called \code{simex}, of class \code{mixstock.data}}
\keyword{datasets}

