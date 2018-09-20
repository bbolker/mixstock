\name{get.bot}
\alias{get.bot}
\alias{get.frq}
\alias{get.ctl}
\alias{get.bse}
\alias{get.mix}
\alias{put.ctl}
\alias{put.bse}
\alias{put.mix}
\alias{put.mp}
\alias{get.mp.input}
\title{Masuda-Pella interchange functions}
  \description{
    Read and write control files for Masuda and Pella's program
} 

\usage{
get.bot(fn)
get.frq(fn)
get.ctl(fn)
get.bse(fn)
put.ctl(fn,which=1,tot,ranseed=c(899271,4480026,90092812),
                 thin=1,title="auto-gen R input",sourcenames,fprior,startval,H)
put.bse(fn,sourcesamp)
put.mix(fn,mixsamp)
put.mp(data,tot=25000,title="Input file auto-gen from R",
                fn="test",ranseed=c(899271,4480026,90092812))
get.mp.input(fn,which=1)
}
\arguments{
  \item{fn}{file name (without extension)}
  \item{which}{numeric value of control file}
  \item{data}{data in mixstock.data format (list of sourcesamp and mixsamp)}
  \item{tot}{total number of Gibbs steps}
  \item{ranseed}{numeric vector for random-number seed}
  \item{thin}{thinning factor}
  \item{title}{title for input files}
  \item{sourcenames}{vector of source names}
  \item{fprior}{prior on source contributions (default=c(1/R,R)}
  \item{startval}{starting values for source contributions}
  \item{H}{number of markers}
  \item{sourcesamp}{matrix (HxR) of marker samples in sources}
  \item{mixsamp}{vector of marker sample in mixed population}
}
\value{
  produces output files or returns a list
  of
  \item{sourcesamp}{source samples (matrix)}
  \item{mixsamp}{pooled samples (vector)}
  \item{a}{prior strength}
  \item{startiter}{starting iteration}
  \item{maxiter}{maximum iteration}
  \item{startfval}{starting fvalue}
  \item{thin}{thinning factor}
  \item{fprior}{source prior}
  \item{rptiter}{how often to print output to the screen}
}
\note{put.mp and get.mp.input are portmanteau functions}
\author{Ben Bolker}
\keyword{misc}
