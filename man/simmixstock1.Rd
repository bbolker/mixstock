\name{simmixstock1}
\alias{simmixstock1}
\alias{sim.mark.freq}
\title{Simulate marker frequencies and distributions in a mixed stock}
\description{
  Functions for simulating marker frequency distributions and samples
  in source and mixed populations
}
\usage{
simmixstock1(sampsize = NULL, true.freq = matrix(c(0.65, 0.33, 0.01, 
    0.01, 0.33, 0.65, 0.01, 0.01), ncol = 2), true.contrib = c(0.9, 
    0.1), boot = FALSE, param = FALSE, data = NULL, rseed = 1004, 
    nboot = 1000, chainlen = NULL, ests = c("cmlboot.nonpar", 
        "cmlboot.par", "umlboot.nonpar", "umlboot.par", "mcmc"), 
    verbose = FALSE, contrib.only = FALSE)
sim.mark.freq(H,R,g.mark,g.source)
}
\arguments{
  \item{sampsize}{total sampsize: half from mixed population,
    (1/(2R)) from each source}
  \item{true.freq}{matrix of marker frequencies in sources}
  \item{true.contrib}{contributions from each source to source
    population}
  \item{boot}{bootstrap existing data?}
  \item{param}{parametric bootstrap?}
  \item{data}{original data set to bootstrap}
  \item{nboot}{number of bootstrap samples}
  \item{chainlen}{chain length for MCMC}
  \item{ests}{list of estimates to produce (parametric or nonparametric
    bootstrap for CML or UML estimation, MCMC)}
  \item{H}{number of markers}
  \item{R}{number of sources}
  \item{g.mark}{geometric distribution parameter for marker frequency}
  \item{g.source}{geometric distribution parameter for source
    contribution}
  \item{rseed}{random number seed}
  \item{contrib.only}{save only source contributions in MCMC chain
    results?}
  \item{verbose}{verbose output?}
}
\value{
  \code{sim.mark.freq} just returns an HxR matrix of  marker
  \code{simmixstock1} returns a list with a \code{\link{genboot}} result
  for each type of estimate requested;
}
\author{Ben Bolker}
\seealso{\code{\link{genboot}}, \code{\link{mysumvec}}}
\keyword{misc}
