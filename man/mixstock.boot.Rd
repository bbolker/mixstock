\name{mixstock.boot}
\alias{mixstock.boot}
\title{Bootstrap samples of mixed stock analysis data}
\description{
Create a bootstrap (multinomial) sample from a given set of 
marker data (source and mixed population genotypes),
either parametric or non-parametric
}
\usage{
mixstock.boot(x, param=FALSE, condense=TRUE,save.freq=FALSE,
param.match="mean")
}
\arguments{
  \item{x}{"Standard" mixstock data, a list with \code{mixsamp} (vector
    of samples in the mixed population) and \code{sourcesamp} (matrix of
    samples in the sources, rows=markers, columns=sources)}
  \item{param}{parametric bootstrapping or not?}
  \item{condense}{use \code{markfreq.condense} on the results?}
  \item{save.freq}{save frequencies?}
  \item{param.match}{match mean or mode of distribution when parametric bootstrapping?}
}
\details{
Nonparametric bootstrapping just resamples the observed data from the
mixed population and from each source with replacement; this is
equivalent to taking a multinomial sample with the probabilities equal
to the observed sample frequencies.  Parametric bootstrapping takes the
observed samples and resamples the probabilities themselves from a
Dirichlet distribution, then takes a multinomial sample.
}
\value{
  A bootstrapped data set, in the same format as the input data: i.e.,
  \item{mixsamp}{samples in mixed population}
  \item{sourcesamp}{samples in sources}
  ...
}
\author{Ben Bolker}

\examples{
true.freq <- matrix(c(0.65,0.33,0.01,0.01,
                       0.33,0.65,0.01,0.01),ncol=2)
true.contrib <- c(0.9,0.1)
x <- simmixstock0(true.freq,true.contrib,50,100,1004)
nboot <- 1000
boot.results.par <- matrix(NA,ncol=12,nrow=nboot)
boot.results.npar <- matrix(NA,ncol=12,nrow=nboot)
for (i in 1:nboot) {
  x.par  <- mixstock.boot(x,param=TRUE,condense=FALSE)
  x.npar <- mixstock.boot(x,condense=FALSE)
  boot.results.par[i,] <- c(x.par$sourcesamp,x.par$mixsamp)
  boot.results.npar[i,] <- c(x.npar$sourcesamp,x.npar$mixsamp)
}
summary(boot.results.par[,7:8])
summary(boot.results.npar[,7:8])
par(mfrow=c(1,2))
hist(boot.results.par[,7])
hist(boot.results.npar[,7])
}
\keyword{misc}
