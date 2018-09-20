\name{runsims}
\alias{runsims}
\title{Run mixed stock simulations}
\description{
   Run multiple simulations of a mixed stock systems
   with specified marker frequencies and source contributions, running
   multiple estimations (bootstrap samples or MCMC chains) for each
   simulation
}
\usage{
runsims(sim.n=10, mc.n=10, totsamp=200, which="all",
true.freq=matrix(c(0.65, 0.31, 0.01, 0.01, 0.01, 0.01,
               0.31,0.65,0.01,0.01,0.01,0.01),ncol=2),
       true.contrib=c(0.9, 0.1), est="MCMC", verbose=FALSE,
fuzz=0.001, nboot=1000, bootrpt=20, minmarks=3)
}
\arguments{
  \item{sim.n}{Number of simulations to run}
  \item{mc.n}{Number of bootstrap/MCMC chains to run for each simulation}
  \item{totsamp}{total sample size (to be distributed half in 
    mixed stock and half, evenly, among sources)}
  \item{which}{which markers to use: "common", "rare", or "all"}
  \item{true.freq}{matrix of true marker frequencies in the sources
  (column=source, row=marker)}
  \item{true.contrib}{vector of true contributions of sources to the
    mixed stock}
  \item{est}{estimation method: "MCMC" (Markov Chain Monte Carlo), "cml"
    (conditional max. likelihood), or "uml" (unconditional max. likelihood)}
  \item{verbose}{produce lots of debugging output?}
  \item{fuzz}{"fuzz" parameters for (e.g.) keeping estimated values away
  from 0/1}
  \item{nboot}{Number of bootstrap samples/length of MCMC chain}
  \item{bootrpt}{Frequency for reporting on the progress of bootstrap code}
  \item{minmarks}{Minimum number of markers to allow for a simulation}
}
\value{
  Array of the results from all simulations and bootstrap samples
}
\examples{
## mild kluge to drop unneeded dimensions
  runsims(sim.n=1,mc.n=1,nboot=100,est="cml")[,,]
  }
\keyword{misc}
