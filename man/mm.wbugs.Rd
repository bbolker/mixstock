\name{mm.wbugs}
\alias{mm.wbugs}
\alias{write.TO.bugscode}
\title{Run many-to-many model via WinBUGS (or JAGS)}
\description{
  Sets up the many-to-many model and passes it to WinBUGS via R2WinBUGS
  (or R2jags)}
\usage{
mm.wbugs(x, sourcesize,
n.iter=20000,  n.burnin=floor(n.iter/2),
n.chains=x$R,
n.thin=max(1, floor(n.chains * (n.iter - n.burnin)/1000)),
files.only = FALSE,
inittype=c("dispersed","random"),bugs.code=c("TO","BB"),
                     returntype=c("mixstock","coda","bugs"),
pkg=c("WinBUGS","JAGS"),
mixprior=1,
which.init,
debug=FALSE,
working.directory,...)
write.TO.bugscode(fn,MIX)
}
\arguments{
  \item{x}{a mixstock data object}
  \item{sourcesize}{Relative sizes of sources}
  \item{n.iter}{Total length of each chain}
  \item{n.burnin}{Number of burn-in iterations}
  \item{n.chains}{Number of chains (default, number of sources)}
  \item{n.thin}{thinning rate.  Must be a positive integer. Set 'n.thin' > 1
    to save memory and computation time if 'n.iter' is large.
    Default is 'max(1, floor(n.chains * (n.iter-n.burnin) /
    1000))' which will only thin if there are at least 2000
    simulations.}
  \item{files.only}{(unimplemented) don't run WinBUGS, just produce
    input files}
  \item{inittype}{"dispersed" or "random" depending on how you
    want multiple chains to be initialized}
  \item{bugs.code}{"TO" or "BB" depending on whether you
    want old-style (Toshi Okuyama=TO) or more compact but
    possibly slower (Ben Bolker=BB) code}
  \item{mixprior}{Dirichlet prior for contributions to mixed
    stocks.  Should be either (1) a a single numeric value
    which will be replicated the appropriate number of times;
    (2) a vector of length (nmix+1) [one extra for the unknown stock],
    identical for all sources -- or (3) a matrix with (nmix+1)
    columns and a row for each source}
  \item{which.init}{for "dispersed" start with fewer chains than
    sources, which sources should be used as the dominant sources
    in the chains? (default is a random sample without replacement from
    the list of sources)}
  \item{debug}{run BUGS in debug mode?  (i.e. don't
    exit and go back to R automatically)}
  \item{\dots}{other arguments to \code{bugs}}
  \item{fn}{file name to write BUGS code to}
  \item{MIX}{number of mixed stocks}
  \item{returntype}{return value as a \code{mixstock.est} object,
    a CODA object, or a BUGS object?}
  \item{pkg}{which package to use for back-end calculations}
  \item{working.directory}{working directory for BUGS calculations}
}
\value{
  results of WinBUGS run, as a \code{mixstock.est} object
  by default: type varies according to \code{returntype}.
  \code{write.TO.code} produces a BUGS model file.
}
\note{
  For diagnostic purposes, it may be worth running the
  code initially with \code{returntype="bugs"}
  and using \code{as.mcmc.bugs} and \code{as.mixstock.est.bugs}
  to convert the result to either CODA format or mixstock
  format.

  pkg="JAGS" is still experimental.
}
\author{Ben Bolker}
\keyword{misc}
