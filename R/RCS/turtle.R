.First.lib <- function(lib,pkg) {
  library.dynam("turtle",pkg,lib)
  require(bbfuns)
  require(bbdists)
  require(bbgraphics)
  require(mleprof)
  require(coda)
#  require(boa)
}

as.turtle.par <- function(object) {
  class(object) <- "turtle.par"
  object
}

as.turtle.data <- function(object) {
  result <- list()
  if (is.data.frame(object) || is.matrix(object)) {
    ## data frame/matrix, convert to list
    result$R <- ncol(object)-1
    result$H <- nrow(object)
    result$rooksamp <- object[,1:result$R]
    result$poolsamp <- object[,ncol(object)]
  } else {
    result <- object
    if (is.null(object$R))
      result$R <- ncol(object$rooksamp)
  if (is.null(object$H))
    result$H <- nrow(object$rooksamp)
  }
  result <- label.turtle.data(result)
  class(result) <- "turtle.data"
  result
}

rooknames <- function(x) {
  dimnames(x$rooksamp)[[2]]
}

hapnames <- function(x) {
  dimnames(x$rooksamp)[[1]]
}

samp2freq <- function(x)
  sweep(x,2,apply(x,2,sum),"/")

as.turtle.est <- function(object) {
  class(object) <- c("turtle.est",class(object))
  object
}

turtle.est <- function(fit,resample=NULL,data=NULL,em=FALSE,
                       rooksamp=NULL,poolsamp=NULL,R=NULL,H=NULL,
                       contin=FALSE,method="unknown",
                       boot.method="none",
                       boot.data=NULL,gandr=NULL,prior=NULL) {
  if (is.null(rooksamp))
    if (!is.null(data))
      rooksamp <- data$rooksamp
  if (is.null(poolsamp))
    if (!is.null(data))
      poolsamp <- data$poolsamp
  if (is.null(R) | is.null(H))
    if (is.null(rooksamp) & is.null(data))
      stop("must provide either # of rookeries and haplotypes or data")
    else {
      R <- ncol(rooksamp)
      H <- nrow(rooksamp)
    }
  object <- list(fit=fit,resample=resample,
                 data=list(rooksamp=rooksamp,poolsamp=poolsamp,R=R,H=H),
                 R=R,H=H,contin=contin,method=method,boot.method=boot.method,
                 boot.data=boot.data,gandr.diag=gandr,prior=prior,em=em)
  class(object) <- "turtle.est"
  object
}

print.turtle.data <- function(x,...) {
  cat(x$R,"rookeries,",x$H,"distinct haplotypes\n")
  cat("Sample data:\n")
  x.dat <- cbind(x$rooksamp,x$poolsamp)
  dimnames(x.dat) <- list(hap=dimnames(x.dat)[[1]],
                       rook=c(dimnames(x$rooksamp)[[2]],"mixed"))
  print(x.dat)
}

## print parameters for a turtle estimate
print.turtle.par <- function(x,...) {
}

## print info for a turtle estimate
## want to print -logLik as well?
print.turtle.est <- function(x,debug=FALSE,...) {
  if (x$method=="mcmc") {
    if (debug) cat("MCMC\n")
    vals <- x$fit
  }
  else if (x$em==TRUE && x$method == "uml"){
    vals <- list(input.freq=x$fit$input.freq,rook.freq=x$fit$rook.freq)
  }     
  else if (!is.null(class(x$fit)) && class(x$fit)=="turtle.est") {
    if (debug) cat("calculating coef(x$fit)\n")
    vals <- coef(x$fit)
  }
  else { ## mle
    if (debug) cat("calculating xval\n")
    vals <- xval(coef(x$fit),R=x$R,H=x$H,x.orig=x$data,contin=x$contin,
                 input.only=(x$method=="cml"))
  }
  cat("Estimated input contributions:\n")
  print(vals$input.freq)
  cat("\nEstimated rookery frequencies:\n")
  if (x$method != "cml")
    print(vals$rook.freq)
  else cat("(cml: no estimate)\n")
  cat("\nmethod:",x$method,"\n")
  if (x$method=="mcmc")
    cat("prior strength:",x$prior,"\n")
}

coef.turtle.est <- function(object,...) {
  if (object$method=="mcmc")
    object$fit
  else
    xval(coef(object$fit),R=object$data$R,H=object$data$H,x.orig=object$data,contin=object$contin,
         input.only=(object$method=="cml"))
}

summary.turtle.est <- function(object,...) {
  print.turtle.data(object$data)
  print.turtle.est(object)
  if (!is.null(object$resample)) {
    cat("\nResampling method:",object$boot.method,"\n")
    print(summary(object$resample))
  }
  invisible(list(data=object$data,fit=object$fit,resample.sum=mysum(object$resample)))
}

plot.turtle.data <- function(x,prop=TRUE,legend=TRUE,
                             colors=rainbow(x$H),leg.space=0.3,
                             leg.ncol=3,pool.off=0.5,
                             ...) {
  H <- nrow(x$rooksamp)
  R <- ncol(x$rooksamp)
  vals <- cbind(as.matrix(x$rooksamp),x$poolsamp)
  dimnames(vals)[[2]][R+1] <- "Mixed"
  hapnames <- dimnames(vals)[[1]]
  if (prop) vals <- sweep(vals,2,apply(vals,2,sum),"/")
  if (!legend) leg.space <- 0
  y.ht <- 1/(1-leg.space)
  b <- barplot(vals,ylim=c(0,y.ht),axes=FALSE,col=colors,
               space=c(rep(0.2,R),pool.off),...)
  axis(side=2,at=seq(0,1,by=0.2))
  if (legend) {
    legend(c(0,H),c(1+(y.ht-1)/10,y.ht),
           hapnames,fill=colors,ncol=leg.ncol)
  }
}


plot.turtle.est <- function(x,legend=FALSE,
                            label=TRUE,min.contrib=0.01,
                            colors=rainbow(x$data$R),
                            hapcolors=rainbow(x$data$H),
                            leg.space=0.5,
                            leg.ncol=3,pool.off=0.5,
                            plot.freqs=FALSE,
                            contrib.lab="Estimated rookery contributions",
                            rookfreq.lab="Estimated rookery haplotype freqs",
                            cex=1,
                            ...) {
  if (is.null(plot.freqs))
    plot.freqs <- (x$method=="uml")
  if (x$method=="mcmc")
    vals <- x$fit
  else if (!is.null(class(x$fit)) && class(x$fit)=="turtle.est") {
    vals <- coef(x$fit)
  }
  else { ## mle
    vals <- xval(coef(x$fit),R=x$data$R,H=x$data$H,
                  x.orig=x$data,contin=x$contin,
                 input.only=(x$method=="cml"))
  }
  ifreq <- vals$input.freq
  w.min <- rep(TRUE,length(ifreq))
  if (min.contrib>0)
    w.min <- ifreq>min.contrib
  ifreq <- ifreq[w.min]
  R <- length(ifreq)
  if (!legend) leg.space <- 0
  x.ht <- 1/(1-leg.space)
  b <- barplot(matrix(ifreq),col=colors[w.min],xlim=c(0,x.ht),
               main=contrib.lab,cex=cex,ylim=c(0,1),...)
  if (label)
    barplot.txt(b,ifreq,names(ifreq),min=0.01,cex=cex)
  if (legend) {
    legend(c(0,H),c(1+(y.ht-1)/10,y.ht),
           hapnames,fill=colors,ncol=leg.ncol)
  }
  if (plot.freqs) {
    if (x$method=="cml")
      warning("CML estimate, no rookery frequencies estimated")
    else
      b <- barplot(coef(x)$rook.freq,col=hapcolors,
                   main=rookfreq.lab,cex=cex,ylim=c(0,1))
  }
}

intervals.turtle.est <- function(x,profile=FALSE,all.par=FALSE,
                                 alpha=0.05) {
  R <- x$R; H <- x$H
  npar <- ifelse(all.par,(R+H*R),R)
  if (!profile && !is.null(x$resample))
    apply(x$resample[,1:npar],2,quantile,c(alpha/2,1-alpha/2),na.rm=TRUE)
  else
    warning("can't do profile yet")
}

deviance.turtle.est <- function(object,...) {
  2*-logLik(object$fit)
}

aic.turtle.est <- function(object,...) {
  R <- object$R
  H <- object$H
  if (object$method=="cml") npar <- R-1
  else npar <- R-1 + (H-1)*R
  deviance(x)+2*npar
}

mcmc.chainlength.est <- function(x,mult=1,inflate=sqrt(2),
                                 gandr.crit=1.2,
                                 nchains=x$R,
                                 verbose=FALSE) {
  if (verbose) cat("Calculating Raftery and Lewis diagnostic\n")
  if (nchains==1) {
    randl.est <- calc.randl.0(x,startfval=0,verbose=verbose)
    chainlen <- randl.max(randl.est)
  }
  else {
    randl.est <- calc.mult.randl(x,n=mult,verbose=verbose)
    chainlen <- max(randl.est[,"total"],na.rm=TRUE)
  }
  ok <- FALSE
  if (verbose) cat("R&L estimated chain length:",chainlen,"\n")
  if (verbose) cat("Calculating Gelman and Rubin diagnostic\n")
  while (!ok) {
    if (mult==1) {
      gandr.est <- calc.gandr(x,tot=chainlen,verbose=verbose)
      ## allow shrink factor to work for newer and older versions of CODA
      psrf <- ifelse(is.null(gandr.est$confshrink),
                        gandr.est$psrf[,2],
                        gandr.est$confshrink[,2])
      ## consider using "multivariate" version??
      ## apply transformation to deal with non-normality?
      gandr.val <- max(psrf,na.rm=TRUE)
    }
    else {
      gandr.est <- calc.mult.gandr(x,n=mult,tot=chainlen,verbose=verbose)
      gandr.val <- max(gandr.est[,3],na.rm=TRUE)
    }
    ok <- gandr.val<=gandr.crit
    if (!ok) {
        chainlen <- ceiling(chainlen*inflate)
        if (verbose)
	   cat("Gelman and Rubin failed: increasing chain length to",
              chainlen,"\n")
      }
  }  ## loop until G&R passes
  chainlen
}
  
q.to.p.R <- function(q,contin=TRUE) {
  # take a vector of (n-1) numbers in (-Inf,Inf) and produce a
  # vector of n numbers in (0,1) that sum to 1
  n <- length(q)+1
  p <- numeric(n)
  if (contin) q <- atan(q)/pi+0.5
  rem <- 1
  for (i in 1:(n-1)) {
    p[i] <- q[i]*rem
    rem <- rem-p[i]
  }
  p[n] <- rem
  p
}

q.to.p <- function(q,contin=TRUE) {
  if (any(is.na(q)))
    rep(NA,length(q)+1)
  else
    .C("q_to_p2",as.double(c(q,0)),
       as.integer(length(q)+1),as.integer(contin))[[1]]
}

## could do more efficiently with (1-cumsum(p))?

## running into trouble with the "branch cut"
p.to.q <- function(p,contin=TRUE) {
  n <- length(p)-1
  rem <- c(1,1-cumsum(p))
  v <- p[1:n]/rem[1:n]
  v[is.na(v)] <- 0
  v[v>1] <- 1  ## clean up branch cut problems
  v[v<0] <- 0
  if (contin) v <- tan((v-0.5)*pi)
  v
}

## log-likelihood of a multinomial mixture
## p: parameter vector.
##    First (R-1) elements are frequency contributions from different sources, expressed in
##    transformed (-infty to infty) coordinates, vector of length (R-1)
##    Last (R:(R*(H+1)-1)) elements are frequencies in the source pools, expressed in
##    transformed coordinates, really an H*(R-1) matrix (by column)
## R: number of sources (rookeries)
## H: number of types (haplotypes)
## rooksamp: sampled proportions in sources (rookery haplotypes), vector (R*H)
## poolsamp: sampled proportions in pool (pelagic pop haplotypes), vector (H)
## second try: need to constrain all freqs to 0-1 and sum(freq)=1
## new parameterization q[] for probability vector from 1..n:
##  p[j] = (atan(q[j])*pi+0.5)*(1-sum(p[1:(j-1)])), p[0] defined as 1
trans.par <- function(p,rooksamp,poolsamp=NULL,contin=TRUE) {
  if (is.list(rooksamp) & is.null(poolsamp)) {
    poolsamp <- rooksamp$poolsamp
    rooksamp <- rooksamp$rooksamp
  }
  H <- length(poolsamp)
  R <- length(as.matrix(rooksamp))/H  ## just in case not in a matrix?
  f <- q.to.p(p[1:(R-1)],contin=contin) # transform to probabilities
  h <- matrix(p[-(1:(R-1))],nrow=H-1)
  h <- apply(h,2,q.to.p,contin=contin) # transform
  pool.freq <- as.vector(h %*% f)  # expected frequency in pool
  list(pool.freq=pool.freq,h=h,R=R,H=H)
}

trans.Rpar <- function(p,rooksamp,poolsamp,h,contin=TRUE) {
  H <- length(poolsamp)
  R <- length(as.matrix(rooksamp))/H
  f <- q.to.p(p,contin=contin) # transform input parameters to probabilities
  h <- apply(rooksamp,2,function(z)z/sum(z))
  pool.freq <- as.vector(h %*% f)  # expected frequency in pool
  list(pool.freq=pool.freq,h=h,R=R,H=H)
}

multilik <- function(prob,samp) {
#  cat(prob," ",samp,"\n")
  samp <- as.matrix(samp)
  lik <- samp*log(prob)
### handle zero probs: say we have a term (p^n) in the likelihood.
###  p==0 and n==0; everything's OK (set log-lik to 0)
###  p==0 and n>0; likelihood=0 (set log-lik to -Inf)
  lik[prob==0 & samp==0] <- 0
  lik[prob==0 & samp>0]  <- -Inf
  lik
}

cumcount <- 0

## log-likelihood calling C function
loglik2C <- function(p,rooksamp,poolsamp=NULL,
                     do.grad=FALSE,semiverb=FALSE,verbose=FALSE,
                     contin=TRUE,debug=FALSE,raw=FALSE) {
  if (is.list(rooksamp) && is.null(poolsamp)) {
    poolsamp <- rooksamp$poolsamp
    rooksamp <- rooksamp$rooksamp
  }
  rooksamp <- as.vector(as.matrix(rooksamp))
  lik <- 0
  cumcount <- 0
  H <- length(poolsamp)
  R <- length(rooksamp)/H
  if (!raw)
  result <- .C("loglik2wrap",
               as.double(lik),
               as.double(p),
               as.integer(R),
               as.integer(H),
               as.double(poolsamp),
               as.double(rooksamp),
               as.integer(cumcount),
               as.integer(contin),
               as.integer(debug))[[1]]
  else
  result <- .C("loglik3wrap",
               as.double(lik),
               as.double(p),
               as.integer(R),
               as.integer(H),
               as.double(poolsamp),
               as.double(rooksamp),
               as.integer(cumcount),
               as.integer(debug))[[1]]
  if (verbose) cat("Parms=",p,"\nRook=",rooksamp,"\nMix=",poolsamp,
                   "\ncontin=",contin,"\nLik=",
                   result,"\n")
  return(result)
}

truncvals <- function(x,lower,upper,penfac) {
  penalty <- 0 
  if (any(x<lower | x>upper)) {
    penalty <- penfac*sum(ifelse(x<lower,(lower-x)^2,ifelse(x>upper,(x-upper)^2,0)))
    x <- pmin(upper,pmax(lower,x))
  }
  list(x=x,penalty=penalty)
}

## log likelihood with transformations
loglik2 <- function(p,rooksamp,poolsamp=NULL,
                    do.grad=FALSE,semiverb=FALSE,verbose=FALSE,
                    contin=TRUE,penfac=1000,lower=0,upper=1,fuzz=1e-4) {
  if (is.list(rooksamp) & is.null(poolsamp)) {
    poolsamp <- rooksamp$poolsamp
    rooksamp <- as.vector(as.matrix(rooksamp$rooksamp))
  }
  cumcount <<- cumcount+1
  tr <- trans.par(p,rooksamp,poolsamp,contin=contin)
  trunc <- truncvals(tr$pool.freq,fuzz,1-fuzz,penfac)
  tr$pool.freq <- trunc$x
  penalty <- trunc$penalty
  trunc <- truncvals(tr$h,fuzz,1-fuzz,penfac)
  tr$h <- trunc$x
  penalty <- penalty+trunc$penalty
  lik1 <- sum(multilik(tr$pool.freq,poolsamp)) # log-likelihood of pool samples
  lik2 <- sum(multilik(as.vector(tr$h),as.matrix(rooksamp))) # log-likelihood of source samples
  totlik <- -(lik1+lik2)+penalty
  if (verbose) {
    cat("Raw parameters: ",p,"\n")
    cat("Transformed parameters:\n")
    print(xval(p,R=tr$R,H=tr$H,contin=contin))
    cat("Mixed frequencies: ",tr$pool.freq,"\n")
    cat("Mixed samples: ",poolsamp,"\n")
    cat("Mixed nll=",multilik(tr$pool.freq,poolsamp),"\n")
    cat("sum=",lik1,"\n")
    cat("Mixed nll=",lik1,"\n")
    cat("Source nll=",apply(matrix(multilik(tr$h,rooksamp),nrow=tr$H),2,sum),"\n")
    cat("Total source nll=",sum(multilik(tr$h,rooksamp)),"\n")
    cat("Penalty=",penalty,"\n")
  }
  ### these are bogus!  have to deal with transformation of variables
  ### ugh.  argument for dealing with constraints in a different way??
  if (do.grad) {
    stop("do.grad isn't working")
    all.prob <- c(tr$pool.freq,tr$h)
    all.samp <- c(poolsamp,rooksamp)
    grad <- -all.prob/all.samp
    grad[all.prob==0 & all.samp==0] <- 0
    attr(totlik,"gradient") <- grad
  }
  if (semiverb) cat("loglik2: ",format(totlik,digits=15),"\n")
  totlik
}

loglik2P <- function(p,rooksamp,poolsamp,verbose=FALSE) {
  H <- length(poolsamp)
  R <- length(as.matrix(rooksamp))/H
  loglik2(packval2(p,R,H),rooksamp,poolsamp)
}
  
loglik2R <- function(p,rooksamp,poolsamp,n.freq,verbose=FALSE) {
  ## restricted/conditional maximization -- fix source frequencies (usually to MLE estimates from
  ##  pool samples), get log-likelihood as a function of input freqs only
  tr <- trans.Rpar(p,rooksamp,poolsamp,n.freq)
  lik1 <- sum(multilik(tr$pool.freq,poolsamp)) # log-likelihood of pool samples
  lik2 <- sum(multilik(tr$h,rooksamp)) # log-likelihood of source samples
  totlik <- -(lik1+lik2)
  if (verbose) {
    cat("Raw parameters: ",p,"\n")
    cat("Transformed parameters:\n")
    print(xval(p,R=tr$R,H=tr$H,contin=contin))
    cat("Mixed frequencies: ",tr$pool.freq,"\n")
    cat("Mixed samples: ",poolsamp,"\n")
    cat("Mixed nll=",multilik(tr$pool.freq,poolsamp),"\n")
    cat("sum=",lik1,"\n")
    cat("Source nll=",apply(matrix(multilik(tr$h,rooksamp),nrow=tr$H),2,sum),"\n")
  }
  totlik  
}

## pack input frequencies (1*R) and rookery frequencies (H*R)
# into a transformed
## parameter vector
## N.B. changed index to 2: rookeries as columns
packval <- function(f,r,contin=FALSE) {
  c(p.to.q(f,contin),apply(r,2,p.to.q,contin=contin))
}

packval2 <- function(x,R=2,H=3) { # inverse of xval2
  if (length(x) != R+R*H) stop("Unequal frequencies in packval: did you forget R and H?")
  c(p.to.q(x[1:R]),apply(matrix(x[(R+1):length(x)],nrow=H),2,p.to.q))
}

xval <- function(p,R=2,H=2,x.orig=NULL,contin=TRUE,input.only=FALSE) {
  # unpack a parameter vector (transformed frequencies (R-1) of pool and (H*(R-1)) of sources)
  if (!input.only & length(p) != (R-1)+R*(H-1))
    stop("Unequal frequencies in xval: did you forget R and H?")
  input.freq <- q.to.p(p[1:(R-1)],contin=contin)
  if (!input.only)
    rook.freq <- zapsmall(apply(matrix(p[-(1:(R-1))],
                                       nrow=(H-1)),
                                2,q.to.p,contin=contin))
  else rook.freq <- NULL
  if (!is.null(x.orig)) {
    names(input.freq) <- dimnames(x.orig$rooksamp)[[2]]
    if (!input.only) dimnames(rook.freq) <- dimnames(x.orig$rooksamp)
  }
  list(input.freq=input.freq,
       rook.freq=rook.freq)
}

xvalR <- function(p,contin=TRUE) {
  # unpack a parameter vector for conditional (restricted) likelihood
  R <- length(p)
  list(input.freq=q.to.p(p))
}

xval2 <- function(p,R=2,H=2,contin=TRUE,x.orig=NULL) {
  # unpacked parameter vector as a vector
  u <- unlist(xval(p,R,H,contin=contin,x.orig=x.orig))
  names(u)[(R+1):length(u)] <- as.vector(outer(dimnames(x.orig$rooksamp)[[1]],
                                     dimnames(x.orig$rooksamp)[[2]],
                                     function(x,y)paste("rookfreq",y,x,sep=".")))
  u
}

gibbsmat <-  function(x,burnin=500,R=2,H=2,trans=TRUE) {
  v <- x$retvals[burnin:nrow(x$retvals),]
  v <- v[v[,"accept"]==1,]  # take just acceptances
  np <- (H-1)*(R+1)
  if (trans)
    t(apply(v[,1:np],1,xval2,R=R,H=H))
  else
    v[,1:np]
}

gibbsrpt <- function(x,burnin=500,R=2,H=2,trans=TRUE,plot=TRUE) {
  ## assume x is metropolis output
  outtab <- matrix()
  v2 <- gibbsmat(x,burnin,R,H,trans)
  np <- (H-1)*(R+1)
  np2 <- H*(R+1)
  if (trans)
    v3 <- xval2(x$est)
  else
    v3 <- x$est
  if (plot)
  par(mfrow=rep(ceiling(sqrt(np2)),2))
  for (i in 1:np2)
    plotvar(v2[,i],v3[i],dimnames(v2)[[2]][i])
}

plotvar <- function(vec,best,name,...) {
  hist(vec,freq=FALSE,xlab=name,...)
  abline(v=quantile(vec,c(0.05,0.95)),col="red",lty=2)
  abline(v=best,col="blue")
}

lsolve <- function(n,s,tol=1e-5,warn=FALSE) {
  n <- apply(n,2,function(x)x/sum(x))
  R <- ncol(n)
  s <- as.vector(s/sum(s))
  sing <- FALSE
  if (nrow(n)==ncol(n))
    sing <- any(abs(eigen(n)$values)<tol)
  if (sing) {
    warning("singular matrix, returning equal contribs")
    rep(1/R,R)
  } else {
    options(show.error.messages=FALSE)
    m <- try(solve(n,s,tol=tol))
    options(show.error.messages=TRUE)
    if (!is.null(class(m))) {
      if (class(m)=="try-error") {
        if (warn) warning("solve failed, returning equal contribs")
        rep(1/R,R)
      } else {
        if (warn) warning("solve failed (??), returning equal contribs")
        rep(1/R,R)
      }
    }
    else m
  }
}

startvec0 <- function(rooksamp,poolsamp=NULL,type="equal",sd=1,lmin=1e-3) {
  if (is.null(poolsamp) && is.list(rooksamp)) {
      poolsamp <- rooksamp$poolsamp
      rooksamp <- rooksamp$rooksamp
    }
  H <- length(poolsamp)  # number of haplotypes
  R <- length(as.matrix(rooksamp))/H # number of rookeries
  if (type=="equal") {
    # set input contribs all equal
    f <- rep(1/R,R)
  }
  else if (type=="random") { ## random sample with equal multinom probs
    f <- rmulti(1,sum(rooksamp),rep(1/R,R))
  }
  else if (type=="rand2") { ## random sample, equally distributed
    f <- q.to.p(rnorm(R-1,mean=0,sd=sd))
  }
  else if (is.numeric(type) & type>=1 & type<=length(poolsamp)) {
    f <- rep(0.05/(R-1),R)
    f[type] <- 0.95
  } 
  else {
    ## set input frequencies to linear solution
    #  rookery frequencies to MLEs
    f <- lsolve(rooksamp,poolsamp)
  }
  f[f<=0] <- lmin
  f[f>=1] <- 1-lmin
  f/sum(f)
}

startvec <- function(rooksamp,poolsamp=NULL,type="equal",
                     haptype="sample",a=1,cond=FALSE,
                     transf="full",fuzz=0,sd=1) {
  if (is.list(rooksamp) & is.null(poolsamp)) {
    poolsamp <- rooksamp$poolsamp
    rooksamp <- rooksamp$rooksamp
  }
  # rooksamp: H*R matrix of samples from rookeries
  # s: H vector of samples from pool
  H <- length(poolsamp)  # number of haplotypes
  R <- length(as.matrix(rooksamp))/H # number of rookeries
  f <- startvec0(rooksamp,poolsamp,type=type,lmin=fuzz,sd=sd)
  f <- switch(transf,full=p.to.q(f),part=p.to.q(f,contin=FALSE),
              none=f)
  if (!cond) {
    if (haptype=="sample")  ## take observed sample probs
      xfreq <- apply(rooksamp,2,function(x)x/sum(x))
    else if (haptype=="random") ## random probs
      xfreq <- apply(rnorm((H-1)*R,mean=0,sd=sd),2,q.to.p)
    else if (haptype=="weighted") {  ## Bayes weighted start
      harm.n <- 1/mean(1/apply(rooksamp,2,sum))
      ybar <- apply(apply(rooksamp,2,function(z)z/sum(z)),1,mean)
      rprior <- a*sqrt(harm.n)*matrix(rep(ybar,R),ncol=R)
      xfreq <- apply(rprior,2,function(x)x/sum(x))
    }
    xfreq <- apply(xfreq,2,
                   function(x)switch(transf,
                                     full=p.to.q(x),
                                     part=p.to.q(x,contin=FALSE),
                                     none=x))
    r <- c(f,xfreq)
    if (fuzz>0 & transf!="full") {
      r[r<=0] <- fuzz
      r[r>=1] <- 1-fuzz
    }
    r
  }
  else f
}

## starting vector: restricted (conditional) likelihood
startvecR <- function(rooksamp,poolsamp,type="equal")
  startvec(rooksamp,poolsamp,type=type,cond=TRUE)

## log likelihood, no transformation
loglik3 <- function(p,rooksamp,poolsamp=NULL,do.grad=FALSE,verbose=FALSE) {
  if (is.list(rooksamp) & is.null(poolsamp)) {
    poolsamp <- rooksamp$poolsamp
    rooksamp <- as.vector(as.matrix(rooksamp$rooksamp))
  }
  if (is.list(p))
    if (any(p>1)) {
      warning("prob>1! returning 1000")
      return(1000*max(p)^2)
    }
  cumcount <<- cumcount+1
  H <- length(poolsamp)
  R <- length(as.matrix(rooksamp))/H
  f <- p[1:(R-1)]  # extract mixed-pool frequencies
  f <- c(f,1-sum(f))
  h <- matrix(p[-(1:(R-1))],nrow=H-1) # extract rook freq
  h <- rbind(h,1-apply(h,2,sum))       # complete
  if (verbose) {
    print(h)
    print(f)
  }
  pool.freq <- as.vector(h %*% f)  # expected frequency in pool
  lik1 <- sum(multilik(pool.freq,poolsamp)) # log-likelihood of pool samples
  lik2 <- sum(multilik(h,rooksamp)) # log-likelihood of source samples
  totlik <- -(lik1+lik2)
  if (verbose) {
    cat("Raw parameters: ",p,"\n")
    cat("Mixed frequencies: ",pool.freq,"\n")
    cat("Mixed samples: ",poolsamp,"\n")
    cat("Mixed nll=",multilik(pool.freq,poolsamp),"\n")
    cat("sum=",lik1,"\n")
    cat("Mixed nll=",lik1,"\n")
    cat("Source nll=",apply(matrix(multilik(h,rooksamp),nrow=H),2,sum),"\n")
  }
  totlik
}

## utility function for staggered labels (on multiple lines)
stagger.labs <- function(side,at=NULL,labels=NULL, levels=2, n=NULL, ...) {
  par.old <- par(no.readonly=TRUE)
  mgp.old <- par("mgp")
  if (side==1 | side==3)
    axp <- par("xaxp")
  else
    axp <- par("yaxp")
  if (!is.null(labels))
    n <- length(labels)
  else if (!is.null(at))
    n <- length(at)
  else if (is.null(n)) n <- axp[3]
  if (is.null(at)) at <- pretty(axp[1:2],3)
#  axis(side,labels=FALSE,...)
  par(tcl=0) # shut off tick marks
  for (i in 1:levels) {
    par(mgp=c(mgp.old[1],mgp.old[2]+i-1,mgp.old[3]))
    w <- seq(i,length(at),by=levels)
    axis(side,at=at[w],labels=labels[w],...)
  }
  par(par.old)
}

## utility function for normalizing a vector
normfun <-  function(z)z/sum(z)

#source("misc.R")

gibbsC <- function(a=1,startiter,maxiter,data,poolsamp=NULL,
                   rooksamp=NULL,startfval=NULL,thin=1,
                   fprior=NULL,outfile=FALSE,outfn="turtle-gibbs",
                   randseed=1001,
                   rptiter=-1,
                   debug=FALSE,
                   contrun=FALSE,contrib.start=NULL,
                   rookfreq.start=NULL) {
  ## must match defs in turtle-gibbs.c (or at least be no larger)
  MAXHAP <- 100
  MAXROOK <- 100
  if (is.null(poolsamp) & is.null(rooksamp))
    if (is.null(data))
      stop("must provide data!")
    else {
      poolsamp <- data$poolsamp
      rooksamp <- data$rooksamp
    }
  R <- ncol(rooksamp)
  H <- nrow(rooksamp)
  if (H>MAXHAP) stop("# haplotypes > MAXHAP")
  if (R>MAXROOK) stop("# rookeries > MAXROOK")
  rooksamp <- as.matrix(rooksamp)
  rooksum <- apply(rooksamp,1,sum)
  if (is.null(contrib.start)) contrib.start <- numeric(R)
  if (is.null(rookfreq.start)) rookfreq.start <- numeric(R*H)
  if (any(rooksum==0))
    stop("gibbsC will crash with haplotypes absent from all rookeries")
  ## cat("allocating results vector ...\n")
  if (is.null(fprior)) fprior <- -1
  ## cat("starting gibbswrap ...\n")
  if (debug)
    cat("H:",H,
        "\nR:",R,
        "\na:",a,
        "\nstartiter:",startiter,
        "\nmaxiter:",maxiter,
        "\npoolsamp:",poolsamp,
        "\nrooksamp:",as.vector(rooksamp),
        "\nstartfval:",startfval,
        "\nthin:",thin,
        "\nfprior:",fprior,
        "\nreslen:",(maxiter-startiter)/thin*(R+H*R),
        "\nrandseed:",randseed,
        "\nrptiter:",rptiter,
        "\noutfn:",outfn,"\n")
  r <- .C("gibbswrap",as.integer(H),as.integer(R),
          as.double(a),as.integer(startiter),
          as.integer(maxiter),as.integer(poolsamp),
          as.integer(as.vector(rooksamp)),
          as.integer(startfval),as.integer(thin),
          as.double(fprior),
          as.double(numeric((maxiter-startiter)/thin*(R+H*R))),
          as.integer(outfile),
          as.character(outfn),
          ## as.character(randphrase),
          as.integer(randseed),
          as.integer(rptiter),
          as.integer(contrun),
          as.double(contrib.start),
          as.double(rookfreq.start))
  ##  cat("gibbsC finished\n")
  ##  results <- matrix(r[[11]],nrow=(maxiter-startiter)/thin,ncol=R+H*R)
  ## results  
  tot <- (maxiter-startiter)/thin
  x <- matrix(r[[11]],tot,R+H*R,FALSE)
  ## add default dimnames, if necessary: haplotypes by roman numeral,
  ## rookeries by letter
  if (!is.null(dimnames(rooksamp)[[2]]))
    rooknames <- dimnames(rooksamp)[[2]]
  else
    rooknames <- LETTERS[1:R]
  if (!is.null(dimnames(rooksamp)[[1]]))
    hapnames <- dimnames(rooksamp)[[1]]
  else
    hapnames <- dec.to.roman(1:H)
  dimnames(x) <- list(1:tot,c(paste("contrib",rooknames,sep="."),
                outer(hapnames,rooknames,
                      function(x,y)paste("frq",y,x,sep="."))))
  x
}

label.turtle.data <- function(x,rooknames=NULL,hapnames=NULL) {
  if (is.null(rooknames))
    if (!is.null(dimnames(x$rooksamp)[[2]]))
      rooknames <- dimnames(x$rooksamp)[[2]]
    else
      rooknames <- LETTERS[1:x$R]
  if (is.null(hapnames))
    if (!is.null(dimnames(x$rooksamp)[[1]]))
      hapnames <- dimnames(x$rooksamp)[[1]]
    else if (!is.null(names(x$poolsamp)))
      hapnames <- names(x$poolsamp)
    else
      hapnames <- dec.to.roman(1:x$H)
  dimnames(x$rooksamp) <- list(haplotype=hapnames,rookery=rooknames)
  names(x$poolsamp) <- hapnames
  x
}

"p.bayes" <-
function(rooksamp,bold=1,cut=0.0001){
  aa1_apply(rooksamp,2,function(z) z/sum(z))
  n_apply(rooksamp,2,sum)
  lambda_apply(aa1,1,mean)    #same as ybar
  bnew_bold; bold_10
  while(abs(bnew-bold) > cut){
  bold_bnew
  num_sum(n^2*(1-apply(aa1^2,2,sum))/(n+bold)^3)
  denom_sum(n^2*apply((aa1-lambda)^2,2,sum)/(n+bold)^3)
  bnew_num/denom
}
bnew
}

gibbs <- function(rooksamp,poolsamp,a=1,startiter,maxiter,startfval=NULL,thin=1,
                  fprior=NULL,rptiter=-1) {
  if (any(apply(rooksamp,1,sum)==0))
    stop("Can't do Gibbs with all-missing loci ...")
  R <- ncol(rooksamp)
  H <- nrow(rooksamp)
  totpool <- sum(poolsamp)
  tmpmat <- matrix(nrow=H,ncol=R)
  ## calculate prior according to Pella and Masuda from harmonic mean: a scales the strength
  harm.n <- 1/mean(1/apply(rooksamp,2,sum))
  ybar <- apply(apply(rooksamp,2,function(z)z/sum(z)),1,mean)
  rookprior <- a*sqrt(harm.n)*matrix(rep(ybar,R),ncol=R)
  if (is.null(fprior)) fprior <- rep(1/R,R)  ## default prior for contributions is EQUAL contrib from all rooks
  results <- matrix(nrow=(maxiter-startiter)/thin,ncol=R+H*R)
  ## FIXME: not robust to missing dimnames?
  if (!is.null(dimnames(rooksamp)))
    dimnames(results) <- list(NULL,
                              c(paste("contrib",colnames(rooksamp),sep="."),
                                outer(1:H,1:R,function(x,y)paste("rookhap",colnames(rooksamp)[y],
                                                                 rownames(rooksamp)[x],sep="."))))
  ## rdirichlet fails with freq 0! -- but HAGs not represented in the
  ## source pops should be excluded anyway, unless we're doing fancy stuff
  ## with "unobserved" sources ...
  # set initial rookery freqs
  rookfreqval <- apply(rookprior,2,function(z)rdirichlet(1,z))
  if (is.null(startfval))  ## use random start
    fval <- rdirichlet(1,fprior)
  else if (length(startfval)==1) {
    if (startfval==0)  ## equal-contribution start
      fval <- matrix(rep(1/R,R),ncol=R)
    else if (startfval <= R) {  ## start with 95% in one rookery, the rest evenly divided
      fval <- matrix(rep(0.05/(R-1),R),ncol=R)
      fval[startfval] <- 0.95
    }
    else stop("startfval must be between 0 and R")
  }
  ## pool contribs (f): vector, length R
  ## rook haplotype freqs (h): matrix, H rows (haplotypes) x R cols (rookeries)
  ## pool freqs (pool.freq): h %*% f, vector, length R
  ## "val" indicates realized (Gibbs-sampler) value as opposed to Dirichlet params
  for (j in 1:maxiter) {
    pool.freq <- as.vector(rookfreqval %*% t(fval))  ## expected frequency in pool
    ## probability that an individual with hap H (row) comes from rookery R (column);
    ## use R's "columns first" rule to calculate
    w <- t(apply(rookfreqval,1,"*",fval))/pool.freq
    w[as.logical(match(w[,1],NA,nomatch=0)),] <- fprior
    ## take multinomial samples of each type ...
    for (i in 1:H) {
      tmpmat[i,] <- rmulti(1,poolsamp[i],w[i,])
    }   
    ## get posteriors for p (pool contribs, f) and Q (rook freqs, rookfreq)
    ## supposing we're doing things the easy way:
    ## posterior of p =  (pool sample plus any priors if desired)
    rookcontrib <- apply(tmpmat,2,sum)
    fval <- rdirichlet(1,rookcontrib+fprior)
    ## posterior of Q = (rookery sample + pool sample (known) + priors)
    rookfreqval <- apply(rooksamp+tmpmat+rookprior,2,function(z)rdirichlet(1,z))
    if (j>startiter) {
      w <- j-startiter
      if (w %% thin == 0)
        results[w/thin,] <- c(fval,as.vector(rookfreqval))
    }
  }
  results
}


BUGS.out <- function(g,file=NULL) {
  if (is.null(file)) file <- deparse(substitute(g))
  ## take a gibbs run (nxp dataframe) and output to files in BUGS/CODA-compatible format
  indfn <- paste(file,".ind",sep="")
  outfn <- paste(file,".out",sep="")
  cat("",file=indfn)  ## zero indfn
  cat("",file=outfn)  ## zero outfn
  p <- ncol(g)
  n <- nrow(g)
  vnames <- dimnames(g)[[2]]
  if (is.null(vnames))
    vnames <- paste("var",1:p,sep="")
  for (i in 1:p) {
    cat(vnames[i]," ",(i-1)*n+1,i*n,"\n",file=indfn,append=TRUE)
    write.table(cbind(1:n,g[,i]),file=outfn,append=TRUE,col.names=FALSE,
                row.names=FALSE,quote=FALSE)
  }
}

## rotate a vector
rotvec _ function(x,i) {
  n <- length(x)
  v <- ((1:n)+i+n) %% n
  v[v==0] <- n
  x[v]
}       

## simulate pop. frequencies
## PROBLEM: after normalization, this doesn't preserve
## total haplotype frequencies
sim.hap.freq <- function(H,R,g.hap,g.rook) {
  x <- matrix(NA,nrow=H,ncol=R)
  tothaps <- g.hap^(1:H)
  for (i in 1:H) {
    x[i,] <- rotvec(tothaps[i]^(1:R),i-1)
  }
  x <- apply(x,2,normfun)
  dimnames(x) <- turtle.dimnames(H,R)
  x
}

turtle.dimnames <- function(H,R)
  list(haplotype=dec.to.roman(1:H),rookery=LETTERS[1:R])

## plot haplotype frequencies
hapfreq.plot <- function(x,maxcols=8,colors=NULL) {
  H <- nrow(x)
  R <- ncol(x)
  if (all(x==0 || x>1)) {
    warning("matrix looks like a sample matrix, not a frequency matrix: normalizing")
    x <- apply(x,2,normfun)
  }
  if (is.null(colors))
    colors <- topo.colors(H)
  barplot(as.matrix(x),names=dimnames(x)[[2]],col=colors,ylim=c(0,1.3),axes=FALSE,
          xlab="Rookery/population",ylab="Haplotype frequency")
  axis(side=2,at=seq(0,1,by=0.2))
  legend(0,1.3,dimnames(x)[[1]],ncol=min(H,maxcols),fill=colors)
}

## Lump all haps in a given category (e.g. found only in one rookery)?
## Criteria:
##  1. any haplotype found in the pooled sample but not in any of the rookeries is
##     deleted from both: it provides no further information about the source contributions
##     (although it does contribute to goodness-of-fit tests)
##  2. any set of haplotypes found in only one rookery (and possibly also in the pooled
##     population, although not necessarily) is lumped together
## "exclude.nopool" determines whether to exclude haplotypes not found in pooled pop.
hapfreq.condense <- function(rooksamp=NULL,poolsamp=NULL,debug=FALSE,
                             exclude.nopool=FALSE) {
  err <- FALSE
  if (!is.null(rooksamp) & is.null(poolsamp)) {
    poolsamp <- rooksamp$poolsamp
    rooksamp <- rooksamp$rooksamp
  }
  exclude <- numeric(0)  ## list of haps to exclude
  lump <- list()         ## list of haps to lump together
  for (i in 1:nrow(rooksamp)) {
    if (debug) cat(i,"\n")
    if (is.na(match(i,c(lump,recursive=TRUE)))) {
      ## this hap not already lumped in with another
      ##      if (sum(rooksamp[i,]>0)==0 && poolsamp[i]>0) exclude <- c(exclude,i)
      if (any(is.na(rooksamp)))
        stop(paste("NA in rooksamp:",rooksamp,collapse=""))
      if (!any(rooksamp[i,]>0)) {
        exclude <- c(exclude,i)   ## not found in any rookeries
        if (debug) cat("hap ",i," found in no rookeries\n")
      }
      else if (exclude.nopool & poolsamp[i]==0) { ## not found in pooled pop
        exclude <- c(exclude,i)
        if (debug) cat("hap ",i," not found in pooled sample: exclude\n")
      }
      else
        if (sum(rooksamp[i,]>0)==1) { ## hap present in only one rookery
          if (debug) cat("hap ",i," only in 1 rookery: try lumping\n")
          tmplump <- i
          if (i<nrow(rooksamp))
            for (j in ((i+1):nrow(rooksamp)))
              if (sum(rooksamp[j,]>0)==1 &&
                  which(rooksamp[i,]>0)==which(rooksamp[j,]>0)) {
                ## present in same rookery
                if (debug) cat("hap ",j," also only in same rookery: lump\n")
                tmplump <- c(tmplump,j)
              }
          ## add to "lump" list
          if (length(tmplump)>1) lump <- c(lump,list(tmplump))
          ## if more than one hap like this exists, lump them
        }
    }
  }
  if (length(lump)>0)  {
    if (debug) cat("lumping\n")
    for (i in 1:length(lump)) {
      rooksamp[lump[[i]][1],] <- apply(rooksamp[lump[[i]],],2,sum)   ## add up sample numbers for haps to be lumped
    poolsamp[lump[[i]][1]] <- sum(poolsamp[lump[[i]]])   ## add up sample numbers for haps to be lumped
    dimnames(rooksamp)[[1]][lump[[i]][1]] <- paste(dimnames(rooksamp)[[1]][lump[[i]]],collapse="/") ## combine hap IDs
      exclude <- c(exclude,lump[[i]][-1],recursive=TRUE)        ## throw away rows
    }
  }
  if (length(exclude)>0) {
    if (debug) cat("excluding ",exclude,"\n")
    rooksamp <- rooksamp[-exclude,]
    poolsamp <- poolsamp[-exclude]
  }
  ## eliminate rookeries with no informative haps left
  rooksamp <- rooksamp[,apply(rooksamp,2,sum)>0]
  if (length(poolsamp)<2) {
    warning("Not enough haplotypes left")
    err <- TRUE
  }
  res <- as.turtle.data(list(rooksamp=rooksamp,poolsamp=poolsamp,err=err))
  res
}

## run Heidelberger & Welch tests, report if all passed
diag.check <- function(x,verbose=FALSE) {
   mc <- mcmc(x)
   d <- heidel.diag(mc)
   if (verbose) print(d)
   (all(d[,"stest"]==1 & d[,"htest"]==1))
}

# amoebaC <- function(H,R,poolsamp,rooksamp,startparm=NULL,disp=0.05,tol=1e-8,maxit=100000,
#                     chkfreq=-1,rptfreq=0,cond=FALSE,debuglevel=0) {
#   amcode <- 0
#   if (is.null(startparm)) 
#     startparm <- startvec(rooksamp,poolsamp)
# value <- loglik2(startparm,rooksamp,poolsamp)
#   if (debuglevel>20) {
#     print(startparm)
#     print(value)
#   }
#   retval <- .C("amoebawrap",as.integer(cond),as.integer(H),as.integer(R),
#                as.integer(poolsamp),as.integer(as.vector(as.matrix(rooksamp))),
#                as.double(startparm),as.double(value),as.integer(chkfreq),
#                as.integer(amcode),as.double(disp),
#                as.double(tol),as.integer(maxit),as.integer(rptfreq),as.integer(debuglevel))
#   #  print(retval[[6]])
#   list(par=retval[[6]],val=retval[[7]],code=retval[[9]])
#   ## code 1: < tol (probably OK)
#   ## code 2: too many iterations
#   ## code 3: stuck?
#}

"calc.randl.00" <-
function (data, startfval = 0, pilot = 500, maxit = 15, verbose = FALSE, 
    rseed = 1001, debug = FALSE) 
{
    if (debug == TRUE) 
        verbose <- TRUE
    rooksamp <- data$rooksamp
    poolsamp <- data$poolsamp
    H <- nrow(rooksamp)
    R <- ncol(rooksamp)
    res <- matrix(ncol = 2, nrow = maxit + 1)
    if (debug) 
        cat("Beginning gibbsC\n")
    g0 <- gibbsC(start = 0, thin = 1, maxiter = pilot, poolsamp = 
poolsamp, 
        rooksamp = rooksamp, startfval = startfval, randseed = rseed)
    contrib.end <- g0[pilot,1:R]
    rookfreq.end <- matrix(g0[pilot,(R+1):(R*H+R)],H,R)
    if (debug) 
        print(summary(g0))
    randl.mat <- raftery.diag(g0, 0.975, 0.02, 0.95, 0.001)$resmatrix
    dimnames(randl.mat)[[2]] <- c("Burn-in", "Total", "Lower bound", 
        "Dependence factor")
    if (debug) 
        print(randl.mat)
    if (any(is.na(randl.mat))) 
        warning("NaNs detected in R&L (ignored)")
    randl.parms <- apply(randl.mat, 2, max, na.rm = TRUE)
    if (debug) 
        print(randl.parms)
    if (randl.parms["Lower bound"] > pilot) 
        stop("error: pilot run too short")
    curlen <- pilot
    it <- 1
    if (debug) 
        cat(curlen, randl.parms["Total"], it, maxit, "\n")
    while (curlen < randl.parms["Total"] && it < maxit) {
        if (verbose) 
            cat("it", it, "Current: ", curlen, " Total: ", 
randl.parms["Total"], 
                "\n")
        res[it, ] <- c(curlen, randl.parms["Total"])
        curlen <- randl.parms["Total"]
        add.num <- curlen - pilot
        g1 <- gibbsC(start = 0, thin = 1, maxiter = add.num, 
            poolsamp = poolsamp, rooksamp = rooksamp, startfval = 
startfval, 
            randseed = rseed + it,contrib.start = contrib.end, 
            rookfreq.start = rookfreq.end)
        g1 <- rbind(g0,g1)
        contrib.end <- g1[curlen,1:R]
        rookfreq.end <- matrix(g1[curlen,(R+1):(R*H+R)],H,R)
        pilot <- curlen
        if (debug) 
            print(summary(g1))
        randl.mat <- raftery.diag(g1, 0.975, 0.02, 0.95, 0.001)$resmatrix
        dimnames(randl.mat)[[2]] <- c("Burn-in", "Total", "Lower bound", 
            "Dependence factor")
        if (debug) 
            print(randl.mat)
        randl.parms <- apply(randl.mat, 2, max)
        if (debug) 
            print(randl.parms)
        it <- it + 1
        g0 <- g1
    }
    if (verbose) 
        cat("Current: ", curlen, " Total: ", randl.parms["Total"], 
            "\n")
    res[it, ] <- c(curlen, randl.parms["Total"])
    res <- res[!is.na(res[, 1]), , drop = FALSE]
    dimnames(res) <- list(iteration = 1:nrow(res), c("Current", 
        "Suggested"))
    list(current = randl.mat, history = res)
}

"calc.randl.0" <- 
function (data, startfval = 0, pilot = 500, maxit = 15, verbose = FALSE, 
    rseed = 1001, debug = FALSE, CODA = TRUE){
  if(CODA==TRUE){
    calc.randl.00(data, startfval = startfval, pilot = pilot, maxit = 
maxit, verbose = verbose, 
    rseed = rseed, debug = debug)}
  else{
    calc.randl.11(data, startfval = startfval, pilot = pilot, maxit = 
maxit, verbose = verbose, 
    rseed = rseed, debug = debug)}
}

calc.randl.11 <- function(data,startfval=0,
                         pilot=500,maxit=15,
                         verbose=FALSE,rseed=1001,debug=FALSE) {
  ## calculating Raftery and Lewis diagnostics for data set poolsamp/rooksamp:
  ## run R&L until the chain passes
  ##  require(coda)
  if (debug==TRUE) verbose <- TRUE
  rooksamp <- data$rooksamp
  poolsamp <- data$poolsamp
  H <- nrow(rooksamp)
  R <- ncol(rooksamp)
  res <- matrix(ncol=2,nrow=maxit+1)  ## results matrix
  if (!exists(".boa.par")) {
##  following shouldn't be necessary, since package now requires
##    boa library ...
##   if (file.exists("boa.r")) source("boa.r")  ## have to make more general!
    boa.start()
  }
  ## calculate R&L diagnostics
  if (debug) cat("Beginning gibbsC\n")
  g0 <- gibbsC(start=0,thin=1,
             maxiter=pilot,poolsamp=poolsamp,
             rooksamp=rooksamp,startfval=startfval,randseed=rseed)
  if (debug) print(summary(g0))
  randl.mat <- boa.randl(g0,0.975,0.02,0.95,0.001)
  if (debug) print(randl.mat)
  ## FIXME: there shouldn't *be* an NaN at this point!!
  if (any(is.na(randl.mat)))
    warning("NaNs detected in R&L (ignored)")
  randl.parms <- apply(randl.mat,2,max,na.rm=TRUE)
  if (debug) print(randl.parms)
  if (randl.parms["Lower Bound"]>pilot)
    stop("error: pilot run too short")
  curlen <- pilot
  it <- 1
  if (debug) cat(curlen,randl.parms["Total"],it,maxit,"\n")
  while (curlen<randl.parms["Total"] && it<maxit) {
    if (verbose) cat("it",it,"Current: ",curlen," Total: ",randl.parms["Total"],"\n")
    res[it,] <- c(curlen,randl.parms["Total"])
    curlen <- randl.parms["Total"]
    g1 <- gibbsC(start=0,thin=1,
             maxiter=randl.parms["Total"],
                 poolsamp=poolsamp,
                 rooksamp=rooksamp,startfval=startfval,randseed=rseed+it)
    if (debug) print(summary(g1))
    randl.mat <- boa.randl(g1,0.975,0.02,0.95,0.001)
    if (debug) print(randl.mat)
    randl.parms <- apply(randl.mat,2,max)
    if (debug) print(randl.parms)
    it <- it+1
  }
  if (verbose) cat("Current: ",curlen," Total: ",randl.parms["Total"],"\n")
  res[it,] <- c(curlen,randl.parms["Total"])
  res <- res[!is.na(res[,1]),,drop=FALSE]
#  print(str(res))
#  print(dim(res))
#  print(res)
  dimnames(res) <- list(iteration=1:nrow(res),c("Current","Suggested"))
  list(current=randl.mat,history=res)
}

randl.max <-
function (r) 
{
    if (is.list(r) && (is.null(names(r)) || names(r)[1] != "current")) 
        max(sapply(r, randl.max))
    else {
        Q <- ncol(r$current)
        tot <- r$current[, "Total"]
        if(Q==4){thin <- 1}
        else{thin <- r$current[, "Thin"]}
        return(max(c((tot/thin) * max(thin), max(r$history))))
    }
}

randl.burn <- function(r) {
  ## return burn-in from a R&L diagnostic
  if (is.list(r) && (is.null(names(r)) || names(r)[1] != "current"))
    max(sapply(r,randl.burn))
  else
    max(r$current[,"Burn-in"])
}

calc.mult.randl <- function(data,n=50,debug=FALSE,verbose=FALSE) {
  ## do Raftery and Lewis calculation multiple times, starting
  ## each chain from a large contribution from each rookery in turn
  if (debug) print(data)
  rooksamp <- data$rooksamp
  poolsamp <- data$poolsamp
  R <- ncol(rooksamp)
  resmat <- matrix(nrow=n,ncol=R+2)
  for (i in 1:n) {
    if (verbose) cat("Mult R&L: iteration:",i,"\n")
    rseed <- 1000+i
    set.seed(rseed)
    r <- list()
    if (debug) cat(i,"\n")
    r[[1]] <- calc.randl.0(data,startfval=1,debug=debug,verbose=verbose)
    for (j in 2:R) {
      if (verbose) cat("Mult R&L: chain",j,"\n")
      r[[j]] <- calc.randl.0(data,startfval=j,
                             ## DON'T use pilot from previous runs
                             ## pilot=randl.max(r[[1]]),
                             verbose=verbose,debug=debug)
    }
    tot <- randl.max(r)
    resmat[i,] <- c(rseed,sapply(r,randl.max),tot)
  }
  dimnames(resmat) <- list(iteration=1:nrow(resmat),
                           c("RNG seed",paste("chain",1:R),"total"))
  resmat
}
  
calc.gandr <- function(data,tot=20000,burn=1,verbose=FALSE,rseed=1001,
                       chainfrac=NULL) {
  require(coda)
  poolsamp <- data$poolsamp
  rooksamp <- data$rooksamp
  R <- ncol(rooksamp)
  if (!is.null(rseed)) set.seed(rseed)
  if (is.null(chainfrac)) chainfrac <- 1-1/R
  chain.start <- round(chainfrac*tot)
  g <- lapply(1:R,
        	function(z) {
	      if (verbose) cat("G&R: running chain",z,"of",R,"\n");
              gibbsC(start=burn+chain.start,
                     thin=1,
                     maxiter=tot+burn,
                     poolsamp=poolsamp,
                     rooksamp=rooksamp,
                     randseed=rseed,
                     startfval=z)})
  tmc <- mcmc.list(lapply(g,mcmc))
  if (verbose) cat("Starting G&R diagnostic ...\n")
  g <- gelman.diag(tmc,transform=TRUE)
  if (verbose) cat("                        ... done\n")
  g
}

calc.mult.gandr <- function(data,n=10,tot=20000,burn=1,verbose=FALSE) {
  require(coda)   ## CODA instead of boa (why??)
  poolsamp <- data$poolsamp
  rooksamp <- data$rooksamp
  R <- ncol(rooksamp)
  gresmat <- matrix(nrow=n,ncol=R+1)
  for (i in 1:n) {
    cat(i,"\n")
    rseed <- 1000+i
    g <- calc.gandr(data,tot,burn,verbose,rseed)
    gresmat[i,] <- c(rseed,apply(g$confshrink,2,max))
  }
  dimnames(gresmat) <- list(iteration=1:n,
                            c("RNG seed","Max point est.","Max 97.5% quantile"))
  gresmat
}

plot.turtle.mcmc <- function(x,...) {
  boxplot(as.data.frame(do.call("rbind",x$mcmc),row.names=NULL)[,1:x$data$R],...)
}

plot.turtle.mcmc <- function(x,...) {
  boxplot(as.data.frame(do.call("rbind",x$mcmc),row.names=NULL)[,1:x$data$R],...)
}

boxplot.turtle.est <- function(x,...) {
  R <- arglist[[1]]$data$R
  if (!all(lapply(arglist,function(z)z$data$R)==R))
    stop("Different numbers of rookeries in estimates")
  ## get statistics from each type: summary??
}

hist.turtle.est <- function(x,log="",...) {
  if (is.null(x$resample))
     stop("object doesn't contain resampling information: can't do histogram")
  if (log=="x") {
    v <- log10(x$resample+min(x$resample[x$resample>0]))
  }
  else
    v <- x$resample
  invisible(sapply(1:x$R,function(i)hist(v[,i],freq=FALSE,
                                         xlab="Contribution",
                                         main=paste("Contribution from",
                                           dimnames(x$data$rooksamp)[[2]][i]),
                                         ...)))
}

tmcmc <- function(data,tot=20000,rseed=1001,burn=100,
                  chainfrac=NULL,verbose=FALSE,fprior=NULL,
                  contrib.only=TRUE,rptiter=-1,
                  outfile=NULL,thin=1,
                  lang="C",a=NULL,gr=FALSE) {
  ## run MCMC chain
  ## FIXME: implement thinning???
  ## FIXME: return calculated value of a (prior strength)
  rooksamp <- data$rooksamp
  poolsamp <- data$poolsamp
  H <- length(poolsamp)
  R <- ncol(rooksamp)
  gandr <- NA
  if (is.null(a)) {
    ## set a values according to p.bayes calculation
    harm.n <- sqrt(1/mean(1/apply(rooksamp,2,sum)))
    beta <- p.bayes(rooksamp)
    a <- beta/harm.n
  }
  rooknames <- dimnames(rooksamp)[[2]]
  set.seed(rseed)
  if (is.null(chainfrac)) chainfrac <- 1-1/R
  chain.start <- round(chainfrac*tot)
  add <- (tot-chain.start)*thin - (tot-chain.start)
  if (contrib.only) maxpar <- R
  else maxpar <- R+R*H
  g <- lapply(1:R,function(z) {
    if (verbose) cat("chain",z,"of",R,"\n");
    if (lang=="C") {
      if (is.null(outfile)) {
        outfile <- FALSE
        outfn <- ""
      }
      else {
        outfn <- paste(outfile,".",z,sep="")
        outfile <- TRUE
      }
      gibbsC(start=burn+chain.start,
             thin=thin,
             maxiter=tot+burn+add,
             poolsamp=poolsamp,
             outfile=outfile,outfn=outfn,fprior=fprior,
             rooksamp=rooksamp,rptiter=rptiter,
             startfval=z,a=a)[,1:maxpar]
    }
    else gibbs(start=burn+chain.start,
           thin=thin,
           maxiter=tot+burn+add,
           poolsamp=poolsamp,
           rooksamp=rooksamp,rptiter=rptiter,
           startfval=z,a=a)[,1:maxpar]})
    if(gr==T){ tmc <- mcmc.list(lapply(g, mcmc))
               gandr.est <- gelman.diag(tmc, transform = TRUE)
               psrf <- ifelse(is.null(gandr.est$confshrink),gandr.est$psrf[, 2],gandr.est$confshrink[, 2])
               gandr.val <- max(psrf, na.rm = TRUE)
               gandr <- gandr.val 
    }
  if (verbose) cat("Chains done: trying to combine them\n")
  allchains <- as.data.frame(do.call("rbind",g),row.names=NULL)
  allchains.sum <- apply(allchains,2,mean)
  if (!contrib.only) {
    rook.freq <- matrix(allchains.sum[(R+1):length(allchains.sum)],ncol=R)
    dimnames(rook.freq) <- dimnames(data$rooksamp)
  }
  else
    rook.freq <- NULL
  turtle.est(fit=list(input.freq=allchains.sum[1:R],
               rook.freq=rook.freq),
             resample=allchains,
             data=data,
             method="mcmc",
             boot.method="mcmc",
             contin=FALSE,gandr=gandr,
             prior=a)
}

mysumvec <- function(x,names=NULL) {
  if (is.null(names))
    if (!is.null(names(x)))
      names <- names(x)
    else if (!is.null(dimnames(x)[[2]]))
      names <- dimnames(x)[[2]]
    else
      names <- rep("",ncol(x))
  y <- c(apply(x,2,mean,na.rm=TRUE),apply(x,2,median,na.rm=TRUE),
         apply(x,2,sd,na.rm=TRUE),apply(x,2,quantile,0.025,na.rm=TRUE),
         apply(x,2,quantile,0.05,na.rm=TRUE),apply(x,2,quantile,0.95,na.rm=TRUE),
         apply(x,2,quantile,0.975,na.rm=TRUE))
  names(y) <- t(outer(c("mean","median","sd","Q02.5","Q05","Q95","Q97.5"),
                      names,paste,sep="."))
  y
}

mysum <- function(x,names=NULL) {
    if (is.null(names))
    if (!is.null(names(x)))
      names <- names(x)
    else if (!is.null(dimnames(x)[[2]]))
      names <- dimnames(x)[[2]]
    else
      names <- NULL
    y <- data.frame(apply(x,2,mean,na.rm=TRUE),apply(x,2,median,na.rm=TRUE),
                  apply(x,2,sd,na.rm=TRUE),
                  apply(x,2,quantile,0.025,na.rm=TRUE),
                  apply(x,2,quantile,0.05,na.rm=TRUE),
                  apply(x,2,quantile,0.95,na.rm=TRUE),
                  apply(x,2,quantile,0.975,na.rm=TRUE))
    dimnames(y) <- list(names,
                      c("mean","median","sd","Q02.5","Q05","Q95","Q97.5"))
  y
}

runsims <- function(sim.n=10,mc.n=10,totsamp=200,which="all",
                    true.freq = matrix(c(0.65,0.31,0.01,0.01,0.01,0.01,
                      0.31,0.65,0.01,0.01,0.01,0.01),ncol=2),
                    true.contrib = c(0.9,0.1),
                    est="MCMC",verbose=FALSE,fuzz=1e-3,
                    nboot=1000,bootrpt=20,
                    minhaps=3) {
  statnames <- c("mean","median","sd","Q02.5","Q05","Q95","Q97.5")
  nstats <- length(statnames)
  nrook <- length(true.contrib)
  rooknames <- LETTERS[1:nrook]
  labnames <- rooknames
  resarray <- array(dim=c(sim.n,mc.n,nrook*nstats))
  dimnames(resarray) <- list(sim=as.character(1:sim.n),
                             mcchain=as.character(1:mc.n),
                             stats=t(outer(statnames,
                               labnames,paste,sep=".")))
  rare.only <- FALSE
  common.only <- FALSE
  if (which=="rare") {
    rare.only <- TRUE
  }
  else if (which=="common") {
    common.only <- TRUE
  }
  R <- ncol(true.freq)
  H <- nrow(true.freq)
  if (verbose) cat(R," rookeries ",totsamp," total samples\n")
## assume half of sample in pooled population, rest evenly divided betw. rookeries
  for (i in 1:sim.n) {
    cat("runsim: sim ",i,"\n")
    ok <- FALSE
    failct <- 0
    set.seed(1000+i)
    while (!ok) {
      sim0 <- simturtle0(true.freq,true.contrib,totsamp/(2*R),totsamp/2,NULL)
      if (verbose) print(sim0)
      sim <- hapfreq.condense(sim0)
      if (verbose) print(sim)
      if (sim$err==FALSE)
        ok <- (nrow(sim$rooksamp)>minhaps)  ## more than one rare haplotype in sample?
      if (!ok) {
        failct <- failct+1
      }
      else {
        cat(failct," tries for good sample \n")
        if (common.only) {
          sim$poolsamp <- sim$poolsamp[1:2]
          sim$rooksamp <- sim$rooksamp[1:2,]
        } else
        if (rare.only) {
          sim$poolsamp <- sim$poolsamp[-(1:2)]
          sim$rooksamp <- sim$rooksamp[-(1:2),]
        }
      }
    } ## while !ok
#    print(sim)   
    for (j in 1:mc.n) {
      cat("   runsim: boot/MCMC run ",j,"\n")
      if (est=="MCMC")
        x <- tmcmc(sim,tot=nboot,rseed=1000+sim.n*i+j)$resample[,1:R]
      else if (est=="cml")
        x <- genboot(sim,type="cml",rseed=1000+sim.n*i+j,fuzz=fuzz,nboot=nboot,
                     rpt=bootrpt)$resample[,1:R]
      else if (est=="uml")
        x <- genboot(sim,type="uml",
                     rseed=1000+2*sim.n*i+j,fuzz=fuzz,nboot=nboot,
                     rpt=bootrpt)$resample[,1:R]
      
      else stop(paste("Unknown est type ",est))
      resarray[i,j,] <- mysumvec(x,names=labnames)
    }
  }
  resarray
}

cml <- function(x,start.type="lsolve",fuzz=1e-3,bounds=1e-4,ndepfac=1000,
                method="L-BFGS-B",lower=NULL,upper=NULL,ndeps=NULL,
                control=NULL,
                debug=FALSE,
                grad=cml.grad,
                oldstyle=FALSE,...) {
  ## find contribs (by CML) that give
  ## assume x is a list with elements rooksamp, poolsamp
  R <- ncol(x$rooksamp)
  H <- nrow(x$rooksamp)
  poolfreq  <- x$poolsamp/sum(x$poolsamp)
  rookfreq <- sweep(x$rooksamp,2,apply(x$rooksamp,2,sum),"/")
  trookfreq <- apply(rookfreq,2,p.to.q,contin=FALSE)
  start <- startvec(x,type=start.type,transf="part",fuzz=fuzz,
                    cond=TRUE)
#  start <- startvec0(x$rooksamp,x$poolsamp,type=start.type,lmin=fuzz)
#  start <- pmax(fuzz,pmin(1-fuzz,start))
  if(debug) {
    cat("Starting values:\n")
    print(start)
    cat("Starting log likelihood:\n")
    print(cml.lik(start,rookfreq=rookfreq,data=x))
  }
  npar <- length(start)
  if (is.null(lower)) lower <- rep(bounds,npar)
  if (is.null(upper)) upper <- rep(1-bounds,npar)
  ndeps <- rep(bounds/ndepfac,npar)
  if (!is.null(control)) {
    if (is.null(control$"ndeps"))
      control <- c(control,ndeps=ndeps)
  }
  else control <- list(ndeps=ndeps)
  if (oldstyle) ## old likelihood function, no gradient
  m <- try(mle(start,cml.lik,
               rookfreq=rookfreq,data=x,
               method=method,
#               contin=FALSE,
               lower=lower,
               upper=upper,
               control=control,debug=debug,...))
  else
    m <- try(mle(start,cml.lik2,gr=cml.grad,
                 rookfreq=rookfreq,
                 trookfreq=trookfreq,data=x,
                 method=method,
                 ##               contin=FALSE,
                 lower=lower,
                 upper=upper,
                 control=control,debug=debug,...))
  if ((!is.null(class(m))) && (class(m)=="try-error"))
    m <- list(coefficients=rep(NA,R-1),value=NA,convergence=NA)
  turtle.est(m,data=x,R=R,H=H,method="cml",contin=FALSE)
}

## Jacobian of the parameter transformation
## deriv of c[row] wrt p[column]
rem <- function(x)
  1 - c(0,cumsum(x)[1:(length(x)-1)])

rem <- function(x) {
  .C("rem",x,as.integer(length(x)))[[1]]
}

## numeric derivatives
numderiv <- function(p,fn,eps=1e-5,...) {
  n <- length(p)
  ret <- numeric(n)
  f0 <- fn(p,...)
  for (i in 1:n) {
    v <- p
    v[i] <- p[i]+eps
    tmp <- (fn(v,...)-f0)/eps
    if (is.na(tmp)) {
      v[i] <- p[i]-eps
      tmp <- (fn(v,...)-f0)/(-eps)
    }
    ret[i] <- tmp
  }
  ret
}

numjacob <- function (p, fn, eps = 1e-05, ...) {
    n <- length(p)
    f0 <- fn(p, ...)
    m <- length(f0)
    ret <- matrix(nrow=n,ncol=m)
    for (i in 1:n) {
        v <- p
        v[i] <- p[i] + eps
        tmp <- (fn(v, ...) - f0)/eps
        if (is.na(tmp)) {
            v[i] <- p[i] - eps
            tmp <- (fn(v, ...) - f0)/-eps
        }
        ret[i,] <- tmp
    }
    ret
}


dcmat.a.R <- function(x) {
  len <- length(x)
  remmat <- sapply(1:len,function(z)rem(q.to.p(x[-z],contin=FALSE)))
  ret <- matrix(0,nrow=len,ncol=len)
  m <- row(ret)[lower.tri(ret)]
  n <- col(ret)[lower.tri(ret)]
  ret[lower.tri(ret)] <- x[m]*-diag(remmat[m-1,n])
  diag(ret) <- rem(q.to.p(x,contin=FALSE))[1:len]
  ret <- rbind(ret,-apply(ret,2,sum))
  ret
}  

dcmat.a <- function(x,debug=FALSE) {
  n <- length(x)
  matrix(.C("dcmat",as.double(x),as.double(numeric(n*(n+1))),as.integer(n),
            as.integer(debug))[[2]],
         nrow=n+1)
}
     
cml.grad <- function(p,rookfreq,data,trookfreq=NULL,
                     debug=FALSE,verbose=FALSE) {
  c0 <- q.to.p(p,contin=FALSE)
  M <- data$poolsamp  ## true samples
  m <- as.matrix(rookfreq) %*% c0  ## calc freqs in pool
  j1 <- t(as.matrix(rookfreq)) %*% as.vector(M/m)
  r <-   -as.vector(t(dcmat.a(p)) %*% j1)
  if (debug) {
    cat("cml.grad",r,"\n")
    tstnd <- numderiv(p,cml.lik,rookfreq=rookfreq,data=data)
    cat("cml.grad (numeric):",tstnd,"\n")
    devs <- r-tstnd
    reldev <- (r-tstnd)/tstnd
    cat("cml.grad, devs:",max(abs(devs),na.rm=TRUE),
        max(abs(reldev),na.rm=TRUE),"\n")
  }
  r
}

uml.grad <- function(p,data,
                     rooksamp=NULL,contin=FALSE, ## kluges
                     debug=FALSE,verbose=FALSE) {
  R <- data$R
  H <- data$H
  x <- xval(p,R=R,H=H,contin=FALSE)  ## unpack
  c0 <- x$input.freq
  ##  rookfreq <- as.matrix(x$rook.freq)
  ## M <- data$poolsamp  ## true samples
  m <- x$rook.freq %*% c0  ## calc freqs in pool
  mvec <- as.vector(data$poolsamp/m)
  mvec[!is.finite(mvec)] <- 0
  ## gradient on input contributions
  j1 <- t(x$rook.freq) %*% mvec
  ## new ## j1 <- t(mvec %*% x$rook.freq)
  g1 <- as.vector(-t(dcmat.a(p[1:(R-1)])) %*% j1)
  ## likelihood on mixed pool from rook freqs
  ## need to keep matrix orientation consistent:
  ##  e.g., rookeries=cols,haplotypes=rows
  fmat <- as.matrix(data$rooksamp/x$rook.freq)
  fmat[!is.finite(fmat)] <- 0 ## 0 est, 0 sample
  j2 <- outer(mvec,c0)+fmat
  pmat <- matrix(p[-(1:(R-1))],nrow=H-1)
  #### old code: alternative is, maybe, marginally faster
  ##  tmat <- list()
  ##  for (i in 1:R)
  ##    tmat[[i]] <- -t(dcmat.a(pmat[,i]))
  ##  tmat <- do.call("blockdiag",tmat)
  tmat <- do.call("blockdiag",
                  lapply(as.list(1:R),function(z)-t(dcmat.a(pmat[,z]))))
  g2 <- as.vector(tmat %*% as.vector(j2))
  r <- c(g1,g2)
  if (debug) cat("uml.grad:",r,"\n")
  if (debug) {
    tstnd <- numderiv(p,uml.lik,data=data)
    cat("uml.grad (numeric):",tstnd,"\n")
    devs <- r-tstnd
    reldev <- (r-tstnd)/tstnd
    cat("uml.grad, devs:",max(abs(devs),na.rm=TRUE),
        max(abs(reldev),na.rm=TRUE),"\n")
  }
  r
}

cml.lik2 <- function(p,trookfreq,
                     rookfreq=NULL,  ##kluge
                     data,verbose=FALSE,debug=FALSE) {
  l <- loglik2C(c(p,trookfreq),
                rooksamp=data$rooksamp,
                poolsamp=data$poolsamp,
                verbose=verbose,contin=FALSE)
  if (debug) cat("cml.lik:",p,l,"\n")
  l
}

cml.lik <- function(p,rookfreq,data,verbose=FALSE,debug=FALSE) {
  l <- loglik2C(packval(q.to.p(p,contin=FALSE),rookfreq),
               rooksamp=data$rooksamp,
                poolsamp=data$poolsamp,
               verbose=verbose,contin=FALSE)
  if (debug) cat("cml.lik:",p,l,"\n")
  l
}

uml <- function(x, method="direct",...) {
    switch(method,
               direct = uml.ds(x,...),
               EM =, em = uml.em(x,...))
}

uml.em <- function(x,prec=1e-8,prior=1){
      R <- x$R
      H <- x$H
      rookname <- colnames(x$rooksamp)
      rooksamp <- x$rooksamp
      poolsamp <- x$poolsamp
      poolsum  <- sum(poolsamp)
      fval <- matrix(rep(1/R, R),ncol=R)
      rookfreqval <- apply(rooksamp+prior, 2, function(z) z/sum(z))
      fold <- rep(10,R)
      while(sum(abs(fval-fold)) > prec){
           fold <- fval
           pool.freq <- as.vector(rookfreqval %*% t(fval))
           w <- t(apply(rookfreqval, 1, "*", fval))/pool.freq
           w[as.logical(match(w[, 1], NaN, nomatch = 0)), ] <- rep(0,R)
           rval <- w*matrix(rep(poolsamp,R),H,R)
           fval <- matrix(apply(rval, 2, sum)/poolsum,ncol=R)
           rookfreqval <- apply(rval+rooksamp, 2, function(z) z/sum(z))
      }
      fval <- as.vector(fval)
      names(fval) <- rookname
      m <- list(input.freq=fval,rook.freq=rookfreqval)
      turtle.est(m, data = x, R = R, H = H, method = "uml", em = TRUE)
}

uml.ds <- function(x,
                grad=uml.grad,
                start.type="lsolve",fuzz=1e-3,
                bounds=1e-4,ndepfac=1000,method="L-BFGS-B",
                oldstyle=FALSE,
                debug=FALSE,
                control=NULL,
                ...) {
##  R <- ncol(x$rooksamp)
##  H <- nrow(x$rooksamp)
  R <- x$R
  H <- x$H
  start <- startvec(x,type=start.type,transf="part",fuzz=fuzz)
  defndeps <- rep(bounds/ndepfac,length(start))
  if (!is.null(control)) {
    if (is.null(control$"ndeps"))
      control <- c(control,ndeps=defndeps)
  }
  else control <-   list(ndeps=defndeps)
  if (oldstyle) grad <- NULL
  m <- try(mle(start,
           fn=uml.lik,
           gr=grad,
##           fn=loglik2C,
##           rooksamp=x,
##           contin=FALSE,
           data=x,
           debug=debug,
           method=method,
           lower=rep(bounds,length(start)),
           upper=rep(1-bounds,length(start)),
           control=control,...))
  if ((!is.null(class(m))) && (class(m) == "try-error")){
  m <- list(coefficients = rep(NA, ((R - 1) + R * (H - 1))),
            value = NA, convergence = NA)}
  turtle.est(m,data=x,R=R,H=H,method="uml")
}

uml.lik <- function(p,data,verbose=FALSE,debug=FALSE,raw=FALSE) {
     r <- loglik2C(p,rooksamp=data$rooksamp,
           poolsamp=data$poolsamp,verbose=verbose,contin=FALSE,
           raw=raw)
  if (debug) cat("uml.lik:",r,"\n")
  r
}

turtle.boot <- function(x,param=FALSE,condense=TRUE,save.freq=FALSE,
                        param.match="mean") {
  ## bootstrap rookery samples and haplotype samples
  ## sim-turtle based on observed freqs
  if (!param) {
    poolfreq  <- x$poolsamp/sum(x$poolsamp)
    rookfreq <- sweep(x$rooksamp,2,apply(x$rooksamp,2,sum),"/")
  }
  else {
    ## should observed sample frequencies match the mean, or the mode (==MLE), of
    ## the resampled frequencies????
    if (param.match=="mean")
      alpha.plus <- 0
    else if (param.match=="mode")
      alpha.plus <- 1
    else
      stop("param.match not equal to mean or mode")
    poolfreq <- rdirichlet(1,x$poolsamp+alpha.plus)
    rookfreq <- apply(x$rooksamp,2,function(x)rdirichlet(1,x+alpha.plus))
  }
  rooksampsize <- apply(x$rooksamp,2,sum)
  poolsamp <- as.vector(rmulti(1,sum(x$poolsamp),poolfreq))
  rooksamp <- matrix(nrow=nrow(x$rooksamp),ncol=ncol(x$rooksamp))
  for (i in 1:ncol(rookfreq))
    rooksamp[,i] <- rmulti(1,rooksampsize[i],rookfreq[,i])
#  rooksamp <- matrix(rmulti(1,sum(x$rooksamp),as.vector(as.matrix(rookfreq))),nrow=nrow(x$rooksamp))
#  rooksamp <- matrix(rmulti(2,sum(x$rooksamp),as.vector(as.matrix(rookfreq))),nrow=nrow(x$rooksamp))
  names(poolsamp) <- names(x$poolsamp)
  dimnames(rooksamp) <- dimnames(x$rooksamp)
  data <- as.turtle.data(list(rooksamp=rooksamp,poolsamp=poolsamp))
  if (condense)
    data <- hapfreq.condense(data)
  if (param && save.freq) {
    attr(data,"rookfreq") <- rookfreq
    attr(data,"poolfreq") <- poolfreq
  }

  data
}


genboot <- function(x,type="cml",nboot=1000,rseed=1001,
                    verbose=FALSE,fuzz=1e-3,
                    maxfail=100,rpt=20,
                    start.type="lsolve",
                    param=FALSE,
                    param.match="mean",
                    ndepfac=10000,
                    save.boot=FALSE,
                    print.boot=FALSE,...) {
  set.seed(rseed)
  H <- nrow(x$rooksamp)
  R <- ncol(x$rooksamp)
  hapnames <- dimnames(x$rooksamp)[[1]]
  rooknames <- dimnames(x$rooksamp)[[2]]
  fitnames <- c(paste("contrib",rooknames,sep="."))
  if (type=="cml") {
    basefit <- cml(x,fuzz=fuzz,ndepfac=ndepfac,...)
    npar <- R-1
    noutpar <- R
  } else {
    basefit <- uml(x,start.type=start.type,fuzz=fuzz)
    npar <- (R-1)+R*(H-1)
    noutpar <- R+R*H
    fitnames <- c(fitnames,
                  paste("rookfreq",
                        outer(rooknames,
                              hapnames,paste,sep="."),sep="."))
  }
  nbootres <- noutpar+2
  boot.results <- matrix(nrow=nboot,ncol=nbootres)
  if (save.boot) boot.data <- list()
  else boot.data <- NULL
  for (i in 1:nboot) {
    if (i %% rpt == 0) cat("genboot ",i,"\n")
    ok <- FALSE
    failct <- 0
    while (!ok && failct<maxfail) {
      ## sample until at least 2 OK haplotypes in pooled sample
      x0 <- turtle.boot(x,param=param,param.match=param.match)
      ok <- (length(x0$poolsamp)>1)
      failct <- failct+1
    }
    if (print.boot) cat(as.vector(x0$poolsamp),
                        as.vector(as.matrix(x0$rooksamp)),"\n")
    if (failct<maxfail) {
      ok <- FALSE
      failct <- 0
 #    while (!ok & failct<maxfail) {
        ## why keep trying here??
      if (type=="cml")
        fit <- cml(x0)
      else
        fit <- uml(x0)
      tmpfit <- unlist(coef(fit))
      ok <- all(!is.na(tmpfit)) & all(tmpfit>=0 & tmpfit<=1)
      if (is.na(ok)) ok <- FALSE  ## double-check
      failct <- failct+1
      if (!ok) {
        cat("Failure at",i,":\n")
        print(x0)
      }
#      }  ## while (!ok & failct<maxfail)
    } ## if failct<maxfail
    if (failct>=maxfail | !ok) {
      cat("Reporting NAs\n")
      boot.results[i,] <- as.numeric(rep(NA,nbootres))
    }
    else {
      boot.results[i,] <- c(pad.NA(tmpfit,noutpar),
                            -logLik(fit$fit),fit$fit$convergence)
      if (save.boot) boot.data[[i]] <- x0
    }
  } 
  dimnames(boot.results) <- list(NULL,c(fitnames,"NLL","Convergence"))
  turtle.est(fit=basefit,data=x,
             resample=boot.results,
             method=type,
             boot.method=ifelse(param,
               paste("parametric",param.match,sep="."),
               "nonparametric"),
             boot.data=boot.data)
}

pad.NA <- function(x,len) {
  y <- as.numeric(rep(NA,len))
  y[1:length(x)] <- x
  y
}

resarr.plot <- function(r,np=2,start=1,gap=0.05,fix.names=FALSE,...) {
##  cat("resarr.plot: np=",np,"\n")
  if (fix.names) {
##    cat("fixing names\n")
    r <- fix.names(r)
  }
  nsim <- length(dimnames(r)$sim)
  nmc <- length(dimnames(r)$mc)
  plotCI(x=rep(start:(nsim+start-1),rep(nmc,nsim)),
         y=t(r[,,"median.A"]),ali=t(r[,,"Q05.A"]),aui=t(r[,,"Q95.A"]),
                add=TRUE,gap=gap,...)
}


## generate likelihood profile for rookery contributions
## by a mild kluge: rearrange data to make the rookery
## of interest come first in the data set, then run
## standard profile calculation
turtle.prof <- function(w,data) {
  if (!is.numeric(w))
    w <- which(w==rooknames(data))
  ## change order
  hap0 <- hapnames(data)
  rook0 <- rooknames(data)
  data$rooksamp <- cbind(data$rooksamp[,w],data$rooksamp[,-w])
  dimnames(data$rooksamp) <- list(hap=hap0,rook=
                                  c(rook0[w],rook0[-w]))
  u0 <- uml(data)
  pr <- profile(u0$fit,prange=c(0,1),fun=uml.lik,data=data)
  pr
}

rmnom <- function(size,prob) {
  res <- numeric(length(prob))
  result <- .C("rmnom",as.integer(size),as.double(prob),
               as.integer(length(prob)),as.integer(res))
  result[[4]]
}

rdirich <- function(shape) {
  res <- numeric(length(shape))
  result <- .C("rdirichwrap",as.numeric(shape),
               as.integer(length(shape)),
               as.numeric(res))
  result[[3]]
}

