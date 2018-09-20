## .First.lib <- function(lib,pkg) {
##  library.dynam("mixstock",pkg,lib)
## }

## hard-coded limits in C functions (!!)
MAXMARK <- 500
MAXSRC  <- 100

## blockdiag, rdirichlet, dec.to.roman, addlabels.barplot extracted from
##  bbmisc package -- avoid extra dependencies

##  function to plot labels on a (single, stacked) barplot
addlabels.barplot <- function(x,vals,names,min=0.1,cols=par("fg"),
                        horiz=FALSE,...) {
  if (is.matrix(vals)) {
    x <- rep(x,rep(nrow(vals),length(x)))
    names <- rep(names,ncol(vals))
  }
  else {
    x <- rep(x,length(vals))
    names <- rep(names,length(vals))
  }
  cols <- rep(cols,length.out=length(vals))
  if (is.vector(vals))
      yvals <- (c(0,cumsum(vals)[-length(vals)])+cumsum(vals))/2
  else
    yvals <- apply(vals,2,function(z)(c(0,cumsum(z)[-length(z)])+cumsum(z))/2)
  for (i in 1:length(vals))
    if (vals[i]>min)
      if (horiz)
        text(yvals[i],x[i],names[i],col=cols[i],...)
      else
        text(x[i],yvals[i],names[i],col=cols[i],...)
}

blockdiag <- function(...) {
  args <- list(...)
  nc <- sapply(args,ncol)
  cumnc <- cumsum(nc)
  ##  nr <- sapply(args,nrow)
  ## NR <- sum(nr)
  NC <- sum(nc)
  rowfun <- function(m,zbefore,zafter) {
    cbind(matrix(0,ncol=zbefore,nrow=nrow(m)),m,
          matrix(0,ncol=zafter,nrow=nrow(m)))
  }
  ret <- rowfun(args[[1]],0,NC-ncol(args[[1]]))
  for (i in 2:length(args)) {
    ret <- rbind(ret,rowfun(args[[i]],cumnc[i-1],NC-cumnc[i]))
  }
  ret
}  

rdirichlet<-function(n,alpha)
  ## pick n random deviates from the Dirichlet function with shape parameters a
{
  l  <- length(alpha)
  x  <- matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE)
  sm <- x%*%rep(1,l)
  x/as.vector(sm)
}

## temporary class for optim objects
coef.toptim <- function(object,...) {
  object$par
}

logLik.toptim <- function(object,...) {
  -object$val
}

mixstock.data <- function(sourcesamp,mixsamp,nsource,nmix,sourcesize) {
  if (missing(mixsamp) && (!missing(nsource) | !missing(nmix))) {
    if (!missing(nsource)) nmix <- ncol(sourcesamp)-nsource else
    if (!missing(nmix)) nsource <- ncol(sourcesamp)-nmix else
    stop("must specify mixsamp, nsource, or nmix")
    mixsamp <- sourcesamp[,-(1:nsource)]
    sourcesamp <- sourcesamp[,1:nsource]
  }
  result <- list()
  result$R <- ncol(sourcesamp)
  result$M <- ncol(mixsamp)
  result$H <- nrow(sourcesamp)
  result$sourcesamp <- sourcesamp
  result$mixsamp <- mixsamp
  if (!missing(sourcesize)) result$sourcesize <- sourcesize
  class(result) <- "mixstock.data"
  result
}

as.mixstock.data <- function(object,nmix=1,sourcesize) {
  result <- list()
  if (is.data.frame(object) || is.matrix(object)) {
    ## data frame/matrix, convert to list
    result$R <- ncol(object)-nmix
    result$M <- nmix
    if (!missing(sourcesize)) {
      if (is.numeric(sourcesize)) {
        result$sourcesize <- result$sourcesize
      } else if (is.character(sourcesize)) {
        sourcesizerow <- switch(sourcesize,first=1,last=nrow(sourcesize))
        sourcesize <- object[sourcesizerow,1:result$R]
        object <- object[-sourcesizerow,]
      }
    }
    result$H <- nrow(object)
    result$sourcesamp <- object[,1:result$R]
    result$mixsamp <- object[,(result$R+1):ncol(object)]
  } else {
    result <- object
    if (is.null(object$R))
      result$R <- ncol(object$sourcesamp)
    if (is.null(object$H))
      result$H <- nrow(object$sourcesamp)
  }
  result <- label.mixstock.data(result,
                                sourcenames=sourcenames(result),
                                marknames=marknames(result),
                                mixnames=mixnames(result))
  class(result) <- "mixstock.data"
  result
}

sourcenames <- function(x) {
  dimnames(x$sourcesamp)[[2]]
}

marknames <- function(x) {
  dimnames(x$sourcesamp)[[1]]
}

mixnames <- function(x) {
  m = if (is.null(x$mixsamp)) NULL else dimnames(x$mixsamp)[[2]]
  if (length(m)==0) m="Mixed"
  m
}

locnames <- function(x) {
  c(sourcenames(x),mixnames(x))
}

as.mixstock.est <- function(object) {
  class(object) <- c("mixstock.est",class(object))
  object
}

mixstock.est <- function(fit,resample=NULL,
                         resample.index=NULL,
                         data=NULL,em=FALSE,
                         sourcesamp=NULL,mixsamp=NULL,R=NULL,H=NULL,M=1,
                         transf="full",method="unknown",
                         boot.method="none",
                         boot.data=NULL,GR=NULL,prior=NULL) {
  if (is.null(sourcesamp))
    if (!is.null(data))
      sourcesamp <- data$sourcesamp
  if (is.null(mixsamp))
    if (!is.null(data))
      mixsamp <- data$mixsamp
  if (is.null(R) | is.null(H))
    if (is.null(sourcesamp) & is.null(data))
      stop("must provide either # of sources and markers or data")
    else {
      R <- ncol(sourcesamp)
      H <- nrow(sourcesamp)
    }
  object <- list(fit=fit,resample=resample,resample.index=resample.index,
                 ## used only for one-to-many analyses
                 data=list(sourcesamp=sourcesamp,mixsamp=mixsamp,R=R,H=H,M=1),
                 R=R,H=H,M=1,
                 transf=transf,method=method,boot.method=boot.method,
                 boot.data=boot.data,GR.diag=GR,prior=prior,
                 em=em)
  class(object) <- "mixstock.est"
  object
}

print.mixstock.data <- function(x,...) {
  cat(x$R,"sources,",
      x$M,"mixed stock(s),",
      x$H,"distinct markers\n")
  cat("Sample data:\n")
  x.dat <- cbind(x$sourcesamp,x$mixsamp)
  ##mixnames <- if (!is.matrix(x$mixsamp)) "mixed" else colnames(x$mixsamp)
  ##dimnames(x.dat) <- list(mark=rownames(x$sourcesamp),
  ## source=c(colnames(x$sourcesamp),mixnames))
  print(x.dat)
}

print.mixstock.est <- function(x,debug=FALSE,...) {
  if (x$method=="mcmc") {
    if (debug) cat("MCMC\n")
    vals <- x$fit
  }
  else if (x$em==TRUE && x$method == "uml"){
    vals <- list(input.freq=x$fit$input.freq,source.freq=x$fit$source.freq)
  }     
  else if (!is.null(class(x$fit)) && class(x$fit)=="mixstock.est") {
    if (debug) cat("calculating coef(x$fit)\n")
    vals <- coef(x$fit)
  }
  else { ## mle
    if (debug) cat("calculating unpackval\n")
    vals <- unpackval(coef(x$fit),R=x$R,H=x$H,x.orig=x$data,transf=x$transf,
                 input.only=(x$method=="cml"))
  }
  cat("Estimated input contributions:\n")
  print(vals$input.freq)
  if (!is.null(vals$sourcectr.freq)) {
    cat("\nEstimated input contributions (source-centric):\n")
    print(vals$sourcectr.freq)
  }
  cat("\nEstimated marker frequencies in sources:\n")
  if (x$method != "cml")
    print(vals$source.freq)
  else cat("(cml: no estimate)\n")
  cat("\nmethod:",x$method,"\n")
  if (x$method=="mcmc")
    cat("prior strength:",x$prior,"\n")
}

coef.mixstock.est <- function(object,...) {
  if (object$method=="mcmc")
    object$fit
  else
    unpackval(coef(object$fit),R=object$data$R,H=object$data$H,x.orig=object$data,transf=object$transf,
         input.only=(object$method=="cml"))
}

summary.mixstock.est <- function(object,show.data=TRUE,...) {
  if (show.data) {
    if (is.null(object$data)) {
      cat("No data stored\n")
    } else {
      print.mixstock.data(object$data)
    }
  }
  ## print.mixstock.est(object)
  ## if (!is.null(object$resample)) {
  ## cat("\nResampling method:",object$boot.method,"\n")
  ## print(summary(object$resample))
  ##  }
  cobj <- confint(object,...)
  fit <- object$fit[!sapply(object$fit,is.null)] ## drop NULLs
  all <- cbind(unlist(lapply(fit,t)),cobj)
  rownames(all) <- rownames(cobj)
  cat("\nEstimates:\n")
  nm <- nmix(object)
  ns <- nsource(object)
  ## nh <- nmark(object)
  mm <- (nm>1)
  if (mm)
    cat("\nMixed-stock-centric:\n")
  print(all[1:(nm*ns),])
  if (mm) {
    cat("\nSource-centric:\n")
    print(all[(nm*ns+1):(nm*ns+(nm+1)*ns),])
  }
  if (!mm & nrow(all)>(nm*ns)) {
    cat("\nSource marker frequency estimates:\n")
    print(all[-(1:(nm*ns+1)),])
  }
  invisible(list(data=object$data,fit=object$fit,resample.sum=mysum(object$resample)))
}

plot.mixstock.data <- function(x,prop=TRUE,legend=TRUE,
                               colors=rainbow(x$H),
                               leg.space=0.3,
                               leg.ncol,leg.cex=1,
                               mix.off=0.5,
                               stacklabels=FALSE,
                               sampsize=FALSE,
                               horiz=TRUE,
                               vlab="Haplotype frequency",
                               ...) {
  ## H <- nrow(x$sourcesamp)
  R <- ncol(x$sourcesamp)
  ## combine all values into a single matrix
  vals <- cbind(as.matrix(x$sourcesamp),x$mixsamp)
  ssize = colSums(vals)
  if (length(dim(x$mixsamp))<2) {  ## not matrix or data frame
    dimnames(vals)[[2]][R+1] <- "Mixed"
    M <- 1
  } else M <- ncol(x$mixsamp)
  if (prop) vals <- sweep(vals,2,ssize,"/")
  if (legend) {
    if (horiz) {
      layout(matrix(1:2,nrow=1,ncol=2),widths=c(1,leg.space))
    } else {
      layout(matrix(2:1,nrow=2,ncol=1),heights=c(leg.space,1))
    }
    ## if (!legend) leg.space <- 0
    ## y.ht <- 1/(1-leg.space)
    ## leg.bot <- 1.1
    ## if (!prop) {
    ## y.ht <- y.ht*max(vals)
    ## leg.bot <- leg.bot*max(vals)
  }
  if (horiz) {
    op <- par(las=1)
    on.exit(par(op))
    b <- barplot(as.matrix(vals), ## xlim=c(0,y.ht),
                 axes=FALSE,col=colors,
                 space=c(rep(0.2,R),mix.off,rep(0.2,M-1)),
                 horiz=TRUE,
                 names.arg=rep("",R+M),
                 xlab=vlab,...)
  } else {
    b <- barplot(as.matrix(vals), ## ylim=c(0,y.ht),
                 axes=FALSE,col=colors,
                 space=c(rep(0.2,R),mix.off,rep(0.2,M-1)),
                 horiz=FALSE,
                 names.arg=rep("",R+M),
                 ylab=vlab,...)
  }
  b.axis = if (horiz) 2 else 1
  f.axis = if (horiz) 1 else 2
  if (stacklabels) {
    staxlab(side=b.axis,at=b,labels=locnames(x))
  } else mtext(side=b.axis,line=par("mgp")[2],at=b,locnames(x))
  axis(side=f.axis,at=seq(0,1,by=0.2))
  if (missing(leg.ncol)) {
    leg.ncol <- if (horiz) 1 else 3
  }
  if (sampsize) {
    mtext(side = if (!horiz) 3 else 4,         
          at=b,line=1,ssize)
  }
  if (legend) {
    op <- par(xpd=NA,mar=rep(0,4),xaxs="i",yaxs="i")
    on.exit(par(op))
    plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
    L=legend(0:1,0:1, ## if (!horiz) "topleft" else "right",
           marknames(x),fill=colors,ncol=leg.ncol,cex=leg.cex,
      bty="n")
    if ((horiz && any(L$text$x+strwidth(marknames(x))>1)) ||
        (!horiz && any(L$text$y-strheight(marknames(x))<0)))
      warning("some legend text may be truncated: increase leg.space?")
  }
}

plot.mixstock.est <- function(x,
                              plot.freqs=FALSE,
                              sourcectr=FALSE,
                              contrib.lab="Estimated source contributions",
                              sourcefreq.lab="Estimated source marker freqs",
                              markcolors=rainbow(x$H),
                              alength=0.25,
                              aunits="inches",
                              abbrev,
                              level=0.95,
                              axes=TRUE,
                              ## scales=list(),
                              ## lattice.args=list(),
                              ## layout=NULL,
                              ...) {
  mm <- nmix(x)>1
  if (mm) {
    rnames <- colnames(x$fit$input.freq)
    mnames <- rownames(x$fit$input.freq)
    nx <- dim(x$resamplist$theta)[1]  ## number of samples
    nm <- nmix(x)
    nr <- nsource(x)
    ## little bit of a kluge to work with lattice --
    ##  has to repeat functionality from confint()
    if (!missing(abbrev)) {
        if (is.logical(abbrev)) {
          if (abbrev) abbn <- 3
        } else {
          abbn <- abbrev
          abbrev <- TRUE
        }
      } else abbrev <- FALSE

    if (!sourcectr) {
      if (abbrev) rnames <- abbreviate(rnames,abbn)
      ld2 <- data.frame(x=c(x$resamplist$theta),
                        mix=factor(rep(rep(mnames,nr),each=nx)),
                        source=factor(rep(rnames,each=nm*nx)))
      f <- formula(x~source|mix)
    } else {
      mnames <- c(mnames,"Unk")
      if (abbrev) mnames <- abbreviate(mnames,abbn)
      ld2 <- data.frame(x=c(x$resamplist$div),
                        source=factor(rep(rep(rnames,nm+1),each=nx)),
                        mix=factor(rep(mnames,each=nr*nx)))
      f <- formula(x~mix|source)
    }
    cipanel <- function(x,y,...) {
      L <- split(y,x)
      xvec <- 1:length(L)
      m <- sapply(L,mean)
      ci <- sapply(L,quantile,c((1-level)/2,(1+level)/2))
      lpoints(xvec,m,pch=16)
      larrows(xvec,m,xvec,ci[1,],angle=90,length=alength,units=aunits,)
      larrows(xvec,m,xvec,ci[2,],angle=90,length=alength,units=aunits)
    }
    bwplot(f,data=ld2,ylab="Contribution",
           panel=cipanel,...)
  } else {
    if (x$method=="mcmc")
      vals <- x$fit
    else
      vals <- coef(x)
    ifreq <- vals$input.freq
    ci <- confint(x,level=level)
    if (is.null(ci)) {
      plot(1:x$R,ifreq,
           axes=FALSE,
           xlab="Source",
           ylab=contrib.lab,...)
    } else {
      suppressWarnings(plotCI(1:x$R,y=ifreq,li=ci[,1],ui=ci[,2],
                              axes=FALSE,
                              xlab="Source",
                              ylab=contrib.lab,...))
    }
    if (axes) {
        axis(side=2)
        axis(side=1,at=1:x$R,labels=names(ifreq))
    }
    if (plot.freqs) {
      if (x$method=="cml")
        warning("CML estimate, no source frequencies estimated")
      else
        barplot(coef(x)$source.freq,col=markcolors,
                main=sourcefreq.lab,ylim=c(0,1))
    }
  }
}

confint.mixstock.est <- function(object,parm,level=0.95,
                                 profile=FALSE,
                                 type=c("quantile","credible"),
                                 show.sourcectr=TRUE,
                                 show.haps=FALSE,...) {
  type <- match.arg(type)
  R <- nsource(object)
  ## H <- nmark(object)
  M <- nmix(object)
  mm <- (M>1)
  alpha <- 1-level
  if (!profile)
    if (!is.null(object$resample)) {
      ##      if (length(dim(object$resample))==2) {
      ci <- switch(type,
                   quantile=t(apply(object$resample,2,
                     quantile,c(alpha/2,1-alpha/2),na.rm=TRUE)),
                   credible=HPDinterval(as.mcmc(object$resample)))
      if (!mm && !show.haps) ci <- ci[1:R,]
      if (mm && !show.sourcectr) ci <- ci[1:(R*M),]
      ##       } else {
      ##      return(aperm(apply(object$resample,c(2,3),
      ##                         quantile,c(alpha/2,1-alpha/2),na.rm=TRUE),c(2,3,1)))
      ##  }
      return(ci)
    } else return(NULL)
  if (profile) {
    warning("can't do profile yet")
    return(NULL)
  }
}

logLik.mixstock.est <- function(object,...) {
  logLik(object$fit)
}

AIC.mixstock.est <- function(object,...) {
  R <- object$R
  H <- object$H
  if (object$method=="cml") npar <- R-1
  else npar <- R-1 + (H-1)*R
  -2*logLik(object)+2*npar
}

mcmc.chainlength.est <- function(x,mult=1,inflate=sqrt(2),
                                 GR.crit=1.2,
                                 nchains=x$R,
                                 verbose=FALSE) {
  if (verbose) cat("Calculating Raftery and Lewis diagnostic\n")
  if (nchains==1) {
    RL.est <- calc.RL.0(x,startfval=0,verbose=verbose)
    chainlen <- RL.max(RL.est)
  }
  else {
    RL.est <- calc.mult.RL(x,n=mult,verbose=verbose)
    chainlen <- max(RL.est[,"total"],na.rm=TRUE)
  }
  ok <- FALSE
  if (verbose) cat("R&L estimated chain length:",chainlen,"\n")
  if (verbose) cat("Calculating Gelman and Rubin diagnostic\n")
  while (!ok) {
    if (mult==1) {
      GR.est <- calc.GR(x,tot=chainlen,verbose=verbose)
      ## allow shrink factor to work for newer and older versions of CODA
      psrf <- ifelse(is.null(GR.est$confshrink),
                        GR.est$psrf[,2],
                        GR.est$confshrink[,2])
      ## consider using "multivariate" version??
      ## apply transformation to deal with non-normality?
      GR.val <- max(psrf,na.rm=TRUE)
    }
    else {
      GR.est <- calc.mult.GR(x,n=mult,tot=chainlen,verbose=verbose)
      GR.val <- max(GR.est[,3],na.rm=TRUE)
    }
    ok <- GR.val<=GR.crit
    if (!ok) {
        chainlen <- ceiling(chainlen*inflate)
        if (verbose)
	   cat("Gelman and Rubin failed: increasing chain length to",
              chainlen,"\n")
      }
  }  ## loop until G&R passes
  chainlen
}
  

q.to.p <- function(q,transf="full") {
  if (length(q)>MAXSRC*MAXMARK) stop(paste("exceeded hard-coded size limits:",
                                            "(total vector length=",MAXSRC*MAXMARK,")",sep=""))
  contin = (transf=="full")
  if (any(is.na(q)))
    rep(NA,length(q)+1)
  else
    .C("q_to_p2",as.double(c(q,0)),
       as.integer(length(q)+1),as.integer(contin),PACKAGE="mixstock")[[1]]
}

## running into trouble with the "branch cut"
p.to.q <- function(p,transf="full") {
  if (transf=="none") return(p)
  n <- length(p)-1
  rem <- c(1,1-cumsum(p))
  v <- p[1:n]/rem[1:n]
  v[is.na(v)] <- 0
  v[v>1] <- 1  ## clean up branch cut problems
  v[v<0] <- 0
  if (transf=="full") v <- tan((v-0.5)*pi)
  v
}

## log-likelihood calling C function
loglik2.C <- function(p,
                      sourcesamp,
                      mixsamp=NULL,
                      verbose=FALSE,
                      transf=c("full","part","none"),
                      full=FALSE,
                      cond=FALSE,
                      debug=FALSE) {
  transf <- match.arg(transf)
  if (is.list(sourcesamp) && is.null(mixsamp)) {
    mixsamp <- sourcesamp$mixsamp
    sourcesamp <- sourcesamp$sourcesamp
  }
  transf <- match.arg(transf)
  contin <- (transf=="full")
  transf <- (transf!="none")
  sourcesamp <- as.vector(as.matrix(sourcesamp))
  lik <- 0
  cumcount <- 0
  H <- length(mixsamp)
  R <- length(sourcesamp)/H
  if (R>MAXSRC || H > MAXMARK) stop(paste("exceeded hard-coded size limits:",
                                    " (R=",MAXSRC,", H=",MAXMARK,")",sep=""))
  result <- .C("loglik2wrap",
               as.double(lik),
               as.double(p),
               as.integer(R),
               as.integer(H),
               as.integer(mixsamp),
               as.integer(sourcesamp),
               as.integer(cumcount),
               as.integer(contin),
               as.integer(transf),
               as.integer(full),
               as.integer(cond),
               as.integer(debug),PACKAGE="mixstock")[[1]]
  if (verbose) cat("Parms=",p,"\nSource=",sourcesamp,"\nMix=",mixsamp,
                   "\ncontin=",contin,
                   "\ntransf=",transf,
                   "\nfull=",full,
                   "\nLik=",result,"\n")
  return(result)
}

## pack input frequencies (1*R) and source frequencies (H*R)
## into a transformed parameter vector
## N.B. changed index to 2: sources as columns
packval <- function(f,r,transf="part") {
  c(p.to.q(f,transf=transf),apply(r,2,p.to.q,transf=transf))
}


## inverse of unpackval2
packval2 <- function(x,R=2,H=3) { 
  if (length(x) != R+R*H) stop("Unequal frequencies in packval: did you forget R and H?")
  c(p.to.q(x[1:R]),apply(matrix(x[(R+1):length(x)],nrow=H),2,p.to.q))
}

## unpack a parameter vector (transformed frequencies (R-1) of mix and (H*(R-1)) of sources)
unpackval <-   function(p,R=2,H=2,x.orig=NULL,transf="full",input.only=FALSE) {
  if (transf=="none") {
    input.freq = p[1:R]
    source.freq = matrix(p[(R+1):length(R)],nrow=H)
  } else {
    if (!input.only & length(p) != (R-1)+R*(H-1))
      stop("Unequal frequencies in unpackval: did you forget R and H?")
    input.freq <- q.to.p(p[1:(R-1)],transf=transf)
    if (!input.only)
      source.freq <- zapsmall(apply(matrix(p[-(1:(R-1))],
                                         nrow=(H-1)),
                                  2,q.to.p,transf=transf))
    else source.freq <- NULL
  }
  if (!is.null(x.orig)) {
    names(input.freq) <- dimnames(x.orig$sourcesamp)[[2]]
    if (!input.only) dimnames(source.freq) <- dimnames(x.orig$sourcesamp)
  }
  list(input.freq=input.freq,
       source.freq=source.freq)
}


## unpacked parameter vector as a vector
unpackval2 <- function(p,R=2,H=2,transf="full",x.orig=NULL) {
u <- unlist(unpackval(p,R,H,transf=transf,x.orig=x.orig))
  names(u)[(R+1):length(u)] <- as.vector(outer(dimnames(x.orig$sourcesamp)[[1]],
                                     dimnames(x.orig$sourcesamp)[[2]],
                                     function(x,y)paste("sourcefreq",y,x,sep=".")))
  u
}

gibbsmat <-  function(x,burnin=500,R=2,H=2,trans=TRUE) {
  v <- x$retvals[burnin:nrow(x$retvals),]
  v <- v[v[,"accept"]==1,]  # take just acceptances
  np <- (H-1)*(R+1)
  if (trans)
    t(apply(v[,1:np],1,unpackval2,R=R,H=H))
  else
    v[,1:np]
}

gibbsrpt <- function(x,burnin=500,R=2,H=2,trans=TRUE,plot=TRUE) {
  ## assume x is metropolis output
  ## outtab <- matrix()
  v2 <- gibbsmat(x,burnin,R,H,trans)
  ## np <- (H-1)*(R+1)
  np2 <- H*(R+1)
  if (trans)
    v3 <- unpackval2(x$est)
  else
    v3 <- x$est
  if (plot) {
    par(mfrow=rep(ceiling(sqrt(np2)),2))
    for (i in 1:np2)
      plotvar(v2[,i],v3[i],dimnames(v2)[[2]][i])
  }
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
    return(rep(1/R,R))
  }
  m <- try(qr.solve(n,s,tol=tol),silent=TRUE)
  if (class(m)=="try-error") {
    if (warn) warning("solve failed, returning equal contribs")
    return(rep(1/R,R))
  }
  m
}

startvec0 <- function(sourcesamp,mixsamp=NULL,type="equal",sd=1,lmin=1e-3) {
if (is.null(mixsamp) && is.list(sourcesamp)) {
      mixsamp <- sourcesamp$mixsamp
      sourcesamp <- sourcesamp$sourcesamp
    }
  H <- length(mixsamp)  # number of markers
  R <- length(as.matrix(sourcesamp))/H # number of sources
  if (type=="equal") {
    # set input contribs all equal
    f <- rep(1/R,R)
  }
  else if (type=="random") { ## random sample with equal multinom probs
    f <- rmultinom(1,size=sum(sourcesamp),prob=rep(1/R,R))
  }
  else if (type=="rand2") { ## random sample, equally distributed
    f <- q.to.p(rnorm(R-1,mean=0,sd=sd))
  }
  else if (is.numeric(type) & type>=1 & type<=length(mixsamp)) {
    f <- rep(0.05/(R-1),R)
    f[type] <- 0.95
  } 
  else {
    ## set input frequencies to linear solution
    #  source frequencies to MLEs
    f <- lsolve(sourcesamp,mixsamp)
  }
  f[f<=0] <- lmin
  f[f>=1] <- 1-lmin
  f/sum(f)
}

startvec <- function(sourcesamp,mixsamp=NULL,type="equal",
                     marktype="sample",a=1,cond=FALSE,
                     transf="full",fuzz=0,sd=1) {
  if (is.list(sourcesamp) & is.null(mixsamp)) {
    mixsamp <- sourcesamp$mixsamp
    sourcesamp <- sourcesamp$sourcesamp
  }
  # sourcesamp: H*R matrix of samples from sources
  # s: H vector of samples from mix
  H <- length(mixsamp)  # number of markers
  R <- length(as.matrix(sourcesamp))/H # number of sources
  f <- startvec0(sourcesamp,mixsamp,type=type,lmin=fuzz,sd=sd)
  f <- p.to.q(f,transf)
  if (!cond) {
    if (marktype=="sample")  ## take observed sample probs
      xfreq <- apply(sourcesamp,2,function(x)x/sum(x))
    else if (marktype=="random") ## random probs
      xfreq <- apply(rnorm((H-1)*R,mean=0,sd=sd),2,q.to.p)
    else if (marktype=="weighted") {  ## Bayes weighted start
      harm.n <- 1/mean(1/apply(sourcesamp,2,sum))
      ybar <- apply(apply(sourcesamp,2,function(z)z/sum(z)),1,mean)
      rprior <- a*sqrt(harm.n)*matrix(rep(ybar,R),ncol=R)
      xfreq <- apply(rprior,2,function(x)x/sum(x))
    }
    xfreq <- apply(xfreq,2,p.to.q,transf)
    r <- c(f,xfreq)
    if (fuzz>0 & transf!="full") {
      r[r<=0] <- fuzz
      r[r>=1] <- 1-fuzz
    }
    r
  }
  else f
}



## utility function for scaling column sums to 1
normcols <- function(z) {
  scale(z,center=FALSE,scale=colSums(z))
}

gibbsC <- function(a=1,startiter=0, maxiter,data,mixsamp=NULL,
                   sourcesamp=NULL,startfval=NULL,thin=1,
                   fprior=NULL,outfile=FALSE,outfn="mixstock-gibbs",
                   randseed=1001,
                   rptiter=-1,
                   debug=FALSE,
                   contrun=FALSE,contrib.start=NULL,
                   sourcefreq.start=NULL) {
  ## must match defs in mixstock-gibbs.c (or at least be no larger)
  if (is.null(mixsamp) & is.null(sourcesamp))
    if (is.null(data))
      stop("must provide data!")
    else {
      mixsamp <- data$mixsamp
      sourcesamp <- data$sourcesamp
    }
  R <- ncol(sourcesamp)
  H <- nrow(sourcesamp)
  if (H>MAXMARK) stop(paste("# markers > MAXMARK(=",MAXMARK,")",sep=""))
  if (R>MAXSRC) stop(paste("# sources > MAXSOURCE (=",MAXSRC,")",sep=""))
  sourcesamp <- as.matrix(sourcesamp)
  sourcesum <- apply(sourcesamp,1,sum)
  if (thin>1 && (maxiter-startiter %% thin >0)) {
    maxiter=(((maxiter-startiter) %/% thin+1) *thin)+startiter
    ## warning("maxiter adjusted upwards for thin>1")
  }
  if (is.null(contrib.start)) contrib.start <- numeric(R)
  if (is.null(sourcefreq.start)) sourcefreq.start <- numeric(R*H)
  if (any(sourcesum==0))
    stop("gibbsC will crash with markers absent from all sources")
  ## cat("allocating results vector ...\n")
  if (is.null(fprior)) fprior <- -1
  ## cat("starting gibbswrap ...\n")
  if (debug)
    cat("H:",H,
        "\nR:",R,
        "\na:",a,
        "\nstartiter:",startiter,
        "\nmaxiter:",maxiter,
        "\nmixsamp:",mixsamp,
        "\nsourcesamp:",as.vector(sourcesamp),
        "\nstartfval:",startfval,
        "\nthin:",thin,
        "\nfprior:",fprior,
        "\nreslen:",(maxiter-startiter)/thin*(R+H*R),
        "\nrandseed:",randseed,
        "\nrptiter:",rptiter,
        "\noutfn:",outfn,"\n")
  if (R>MAXSRC || H > MAXMARK) stop(paste("exceeded hard-coded size limits:",
                                          " (R=",MAXSRC,", H=",MAXMARK,")",sep=""))
  r <- .C("gibbswrap",
          as.integer(H),
          as.integer(R),
          as.double(a),
          as.integer(startiter),
          as.integer(maxiter),
          as.integer(mixsamp),
          as.integer(as.vector(sourcesamp)),
          as.integer(startfval),
          as.integer(thin),
          as.double(fprior),
          as.double(numeric((maxiter-startiter)/thin*(R+H*R))),
          as.integer(outfile),
          as.character(outfn),
          ## as.character(randphrase),
          as.integer(randseed),
          as.integer(rptiter),
          as.integer(contrun),
          as.double(contrib.start),
          as.double(sourcefreq.start),PACKAGE="mixstock")
##  cat("gibbsC finished\n")
  ##  results <- matrix(r[[11]],nrow=(maxiter-startiter)/thin,ncol=R+H*R)
  ## results  
  tot <- floor((maxiter-startiter)/thin)
  x <- matrix(r[[11]],nrow=tot,ncol=R+H*R,byrow=FALSE)
  ## add default dimnames, if necessary: markers by roman numeral,
  ## sources by letter
  if (!is.null(dimnames(sourcesamp)[[2]]))
    sourcenames <- dimnames(sourcesamp)[[2]]
  else
    sourcenames <- paste("R",1:R,sep="")
if (!is.null(dimnames(sourcesamp)[[1]]))
    marknames <- dimnames(sourcesamp)[[1]]
  else
    marknames <- paste("H",1:H,sep="")
  dimnames(x) <- list(1:tot,c(paste("contrib",sourcenames,sep="."),
                outer(marknames,sourcenames,
                      function(x,y)paste("frq",y,x,sep="."))))
  x
}

label.mixstock.data <- function(x,sourcenames=NULL,marknames=NULL,mixnames=NULL) {
  if (is.null(sourcenames)) 
    if (!is.null(colnames(x$sourcesamp))) {
      sourcenames <- colnames(x$sourcesamp)
    } else sourcenames <- LETTERS[1:x$R]
  if (is.null(mixnames)) {
    M <- if (is.null(x$M)) 1 else x$M
    if (!is.null(colnames(x$mixsamp))) {
      mixnames <- colnames(x$mixsamp)
    } else mixnames <- paste("M",1:M)
  }
  if (is.null(marknames))
    if (!is.null(rownames(x$sourcesamp))) {
      marknames <- rownames(x$sourcesamp)
    } else if (!is.null(names(x$mixsamp))) {
      marknames <- names(x$mixsamp)
    } else marknames <-  paste("H",1:x$H,sep="")
  dimnames(x$sourcesamp) <- list(marker=marknames,source=sourcenames)
  if (is.null(dim(x$mixsamp))) {
    names(x$mixsamp) <- marknames
  } else dimnames(x$mixsamp) <- list(marker=marknames,mixstock=mixnames)
  x
}

"p.bayes" <- function(sourcesamp,bold=1,cut=0.0001){
  aa1<-apply(sourcesamp,2,function(z) z/sum(z))
  n<-apply(sourcesamp,2,sum)
  lambda<-apply(aa1,1,mean)    #same as ybar
  bnew<-bold; bold<-10
  while(abs(bnew-bold) > cut){
    bold<-bnew
    num<-sum(n^2*(1-apply(aa1^2,2,sum))/(n+bold)^3)
    denom<-sum(n^2*apply((aa1-lambda)^2,2,sum)/(n+bold)^3)
    bnew<-num/denom
  }
  bnew
}

gibbs <- function(sourcesamp,mixsamp,a=1,startiter,maxiter,startfval=NULL,n.thin=1,
                  fprior=NULL,rptiter=-1) {
  if (any(apply(sourcesamp,1,sum)==0))
    stop("Can't do Gibbs with all-missing loci ...")
  R <- ncol(sourcesamp)
  H <- nrow(sourcesamp)
  ## totmix <- sum(mixsamp)
  tmpmat <- matrix(nrow=H,ncol=R)
  ## calculate prior according to Pella and Masuda from harmonic mean: a scales the strength
  harm.n <- 1/mean(1/apply(sourcesamp,2,sum))
  ybar <- apply(apply(sourcesamp,2,function(z)z/sum(z)),1,mean)
  sourceprior <- a*sqrt(harm.n)*matrix(rep(ybar,R),ncol=R)
  if (is.null(fprior)) fprior <- rep(1/R,R)  ## default prior for contributions is EQUAL contrib from all sources
  results <- matrix(nrow=(maxiter-startiter)/n.thin,ncol=R+H*R)
  ## FIXME: not robust to missing dimnames?
  if (!is.null(dimnames(sourcesamp)))
    dimnames(results) <- list(NULL,
                              c(paste("contrib",colnames(sourcesamp),sep="."),
                                outer(1:H,1:R,function(x,y)paste("sourcemark",
                                                                 colnames(sourcesamp)[y],
                                                                 rownames(sourcesamp)[x],sep="."))))
  ## rdirichlet fails with freq 0! -- but HAGs not represented in the
  ## source pops should be excluded anyway, unless we're doing fancy stuff
  ## with "unobserved" sources ...
  # set initial source freqs
  sourcefreqval <- apply(sourceprior,2,function(z)rdirichlet(1,z))
  if (is.null(startfval))  ## use random start
    fval <- rdirichlet(1,fprior)
  else if (length(startfval)==1) {
    if (startfval==0)  ## equal-contribution start
      fval <- matrix(rep(1/R,R),ncol=R)
    else if (startfval <= R) {  ## start with 95% in one source, the rest evenly divided
      fval <- matrix(rep(0.05/(R-1),R),ncol=R)
      fval[startfval] <- 0.95
    }
    else stop("startfval must be between 0 and R")
  }
  ## mix contribs (f): vector, length R
  ## source marker freqs (h): matrix, H rows (markers) x R cols (sources)
  ## mix freqs (mix.freq): h %*% f, vector, length R
  ## "val" indicates realized (Gibbs-sampler) value as opposed to Dirichlet params
  for (j in 1:maxiter) {
    mix.freq <- as.vector(sourcefreqval %*% t(fval))  ## expected frequency in mix
    ## probability that an individual with mark H (row) comes from source R (column);
    ## use R's "columns first" rule to calculate
    w <- t(apply(sourcefreqval,1,"*",fval))/mix.freq
    w[as.logical(match(w[,1],NA,nomatch=0)),] <- fprior
    ## take multinomial samples of each type ...
    for (i in 1:H) {
      tmpmat[i,] <- rmultinom(1,mixsamp[i],w[i,])
    }   
    ## get posteriors for p (mix contribs, f) and Q (source freqs, sourcefreq)
    ## supposing we're doing things the easy way:
    ## posterior of p =  (mix sample plus any priors if desired)
    sourcecontrib <- apply(tmpmat,2,sum)
    fval <- rdirichlet(1,sourcecontrib+fprior)
    ## posterior of Q = (source sample + mix sample (known) + priors)
    sourcefreqval <- apply(sourcesamp+tmpmat+sourceprior,2,function(z)rdirichlet(1,z))
    if (j>startiter) {
      w <- j-startiter
      if (w %% n.thin == 0)
        results[w/n.thin,] <- c(fval,as.vector(sourcefreqval))
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

## bbfuns no longer loaded
## rotate a vector
.rotvec <- function(x,i) {
  n <- length(x)
  v <- ((1:n)+i+n) %% n
  v[v==0] <- n
  x[v]
}       

## simulate pop. frequencies
## PROBLEM: after normalization, this doesn't preserve
## total marker frequencies
sim.mark.freq <- function(H,R,g.mark,g.source) {
  x <- matrix(NA,nrow=H,ncol=R)
  totmarks <- g.mark^(1:H)
  for (i in 1:H) {
    x[i,] <- .rotvec(totmarks[i]^(1:R),i-1)
  }
  x <- normcols(x) ## check!!
  dimnames(x) <- mixstock.dimnames(H,R)
  x
}

mixstock.dimnames <- function(H,R) {
  list(marker=paste("H",1:H,sep=""),source=paste("R",1:R,sep=""))
}

## Lump all marks in a given category (e.g. found only in one source)?
## Criteria:
##  1. any marker found in the mixed stock but not in any of the sources is
##     deleted from both: it provides no further information about the source contributions
##     (although it does contribute to goodness-of-fit tests)
##  2. any set of markers found in only one source (and possibly also in the mixed
##     population, although not necessarily) is lumped together
## "exclude.nomix" determines whether to exclude markers not found in mixed pop.
markfreq.condense <- function(sourcesamp=NULL,mixsamp=NULL,debug=FALSE,
                             exclude.nomix=FALSE) {
  err <- FALSE
  if (!is.null(sourcesamp) & is.null(mixsamp)) {
    mixsamp <- sourcesamp$mixsamp
    sourcesamp <- sourcesamp$sourcesamp
  }
  mixsamp <- as.matrix(mixsamp)
  exclude <- numeric(0)  ## list of marks to exclude
  lump <- list()         ## list of marks to lump together
  for (i in 1:nrow(sourcesamp)) {
    if (debug) cat(i,"\n")
    if (is.na(match(i,c(lump,recursive=TRUE)))) {
      ## this mark not already lumped in with another
      ##      if (sum(sourcesamp[i,]>0)==0 && mixsamp[i]>0) exclude <- c(exclude,i)
      if (any(is.na(sourcesamp)))
        stop(paste("NA in sourcesamp:",sourcesamp,collapse=""))
      if (!any(sourcesamp[i,]>0)) {
        exclude <- c(exclude,i)   ## not found in any sources
        if (debug) cat("mark ",i," found in no sources\n")
      }
      else if (exclude.nomix & all(mixsamp[i,]==0)) { ## not found in mixed pop
        exclude <- c(exclude,i)
        if (debug) cat("mark ",i," not found in mixed stock(s): exclude\n")
      }
      else
        if (sum(sourcesamp[i,]>0)==1) { ## mark present in only one source
          if (debug) cat("mark ",i," only in 1 source: try lumping\n")
          tmplump <- i
          if (i<nrow(sourcesamp))
            for (j in ((i+1):nrow(sourcesamp)))
              if (sum(sourcesamp[j,]>0)==1 &&
                  which(sourcesamp[i,]>0)==which(sourcesamp[j,]>0)) {
                ## present in same source
                if (debug) cat("mark ",j," also only in same source: lump\n")
                tmplump <- c(tmplump,j)
              }
          ## add to "lump" list
          if (length(tmplump)>1) lump <- c(lump,list(tmplump))
          ## if more than one mark like this exists, lump them
        }
    }
  }
  if (length(lump)>0)  {
    if (debug) cat("lumping\n")
    for (i in 1:length(lump)) {
      sourcesamp[lump[[i]][1],] <- colSums(sourcesamp[lump[[i]],])   ## add up sample numbers for haps to be lumped
      mixsamp[lump[[i]][1],] <- colSums(mixsamp[lump[[i]],,drop=FALSE])   ## add up sample numbers for marks to be lumped
      rownames(sourcesamp)[lump[[i]][1]] <- paste(rownames(sourcesamp)[lump[[i]]],collapse="/") ## combine mark IDs
      rownames(mixsamp)[lump[[i]][1]] <- paste(rownames(mixsamp)[lump[[i]]],collapse="/") ## combine mark IDs
      exclude <- c(exclude,lump[[i]][-1],recursive=TRUE)        ## throw away rows
    }
  }
  if (length(exclude)>0) {
    if (debug) cat("excluding ",exclude,"\n")
    sourcesamp <- sourcesamp[-exclude,]
    mixsamp <- mixsamp[-exclude,]
  }
  ## eliminate sources with no informative marks left
  sourcesamp <- sourcesamp[,apply(sourcesamp,2,sum)>0]
  if (length(mixsamp)<2) {
    warning("Not enough markers left")
    err <- TRUE
  }
  res <- as.mixstock.data(list(sourcesamp=sourcesamp,mixsamp=mixsamp,err=err))
  res
}

## run Heidelberger & Welch tests, report if all passed
hwdiag.check <- function(x,verbose=FALSE) {
  mc <- mcmc(x)
  d <- heidel.diag(mc)
  if (verbose) print(d)
  (all(d[,"stest"]==1 & d[,"htest"]==1))
}

calc.RL.0 <- function (data, startfval = 0, pilot = 500, maxit = 15, verbose = FALSE, 
    rseed = 1001, debug = FALSE)
{
    if (debug == TRUE) 
        verbose <- TRUE
    sourcesamp <- data$sourcesamp
    mixsamp <- data$mixsamp
    H <- nrow(sourcesamp)
    R <- ncol(sourcesamp)
    res <- matrix(ncol = 2, nrow = maxit + 1)
    if (debug) 
        cat("Beginning gibbsC\n")
    g0 <- gibbsC(startiter = 0, thin = 1, maxiter = pilot, mixsamp = 
                 mixsamp, 
                 sourcesamp = sourcesamp, startfval = startfval, randseed = rseed)
    contrib.end <- g0[pilot,1:R]
    sourcefreq.end <- matrix(g0[pilot,(R+1):(R*H+R)],H,R)
    if (debug) 
        print(summary(g0))
    RL.mat <- raftery.diag(g0, 0.975, 0.02, 0.95, 0.001)$resmatrix
    dimnames(RL.mat)[[2]] <- c("Burn-in", "Total", "Lower bound", 
        "Dependence factor")
    if (debug) 
        print(RL.mat)
    if (any(is.na(RL.mat))) 
        warning("NaNs detected in R&L (ignored)")
    RL.parms <- apply(RL.mat, 2, max, na.rm = TRUE)
    if (debug) 
        print(RL.parms)
    if (RL.parms["Lower bound"] > pilot) 
        stop("error: pilot run too short")
    curlen <- pilot
    it <- 1
    if (debug) 
        cat(curlen, RL.parms["Total"], it, maxit, "\n")
    while (curlen < RL.parms["Total"] && it < maxit) {
        if (verbose) 
            cat("it", it, "Current: ", curlen, " Total: ", 
RL.parms["Total"], 
                "\n")
        res[it, ] <- c(curlen, RL.parms["Total"])
        curlen <- RL.parms["Total"]
        add.num <- curlen - pilot
        g1 <- gibbsC(startiter = 0, thin = 1, maxiter = add.num, 
            mixsamp = mixsamp, sourcesamp = sourcesamp, startfval = 
startfval, 
            randseed = rseed + it,contrib.start = contrib.end, 
            sourcefreq.start = sourcefreq.end)
        g1 <- rbind(g0,g1)
        contrib.end <- g1[curlen,1:R]
        sourcefreq.end <- matrix(g1[curlen,(R+1):(R*H+R)],H,R)
        pilot <- curlen
        if (debug) 
            print(summary(g1))
        RL.mat <- raftery.diag(g1, 0.975, 0.02, 0.95, 0.001)$resmatrix
        dimnames(RL.mat)[[2]] <- c("Burn-in", "Total", "Lower bound", 
            "Dependence factor")
        if (debug) 
            print(RL.mat)
        RL.parms <- apply(RL.mat, 2, max)
        if (debug) 
            print(RL.parms)
        it <- it + 1
        g0 <- g1
    }
    if (verbose) 
        cat("Current: ", curlen, " Total: ", RL.parms["Total"], 
            "\n")
    res[it, ] <- c(curlen, RL.parms["Total"])
    res <- res[!is.na(res[, 1]), , drop = FALSE]
    dimnames(res) <- list(iteration = 1:nrow(res), c("Current", 
        "Suggested"))
    list(current = RL.mat, history = res)
}

RL.max <-
function (r) 
{
    if (is.list(r) && (is.null(names(r)) || names(r)[1] != "current")) 
        max(sapply(r, RL.max))
    else {
        Q <- ncol(r$current)
        tot <- r$current[, "Total"]
        if(Q==4){thin <- 1}
        else{thin <- r$current[, "Thin"]}
        return(max(c((tot/thin) * max(thin), max(r$history))))
    }
}

RL.burn <- function(r) {
  ## return burn-in from a R&L diagnostic
  if (is.list(r) && (is.null(names(r)) || names(r)[1] != "current"))
    max(sapply(r,RL.burn))
  else
    max(r$current[,"Burn-in"])
}

calc.mult.RL <- function(data,n=50,debug=FALSE,verbose=FALSE) {
  ## do Raftery and Lewis calculation multiple times, starting
  ## each chain from a large contribution from each source in turn
  if (debug) print(data)
  sourcesamp <- data$sourcesamp
  ## mixsamp <- data$mixsamp
  R <- ncol(sourcesamp)
  resmat <- matrix(nrow=n,ncol=R+2)
  for (i in 1:n) {
    if (verbose) cat("Mult R&L: iteration:",i,"\n")
    rseed <- 1000+i
    set.seed(rseed)
    r <- list()
    if (debug) cat(i,"\n")
    r[[1]] <- calc.RL.0(data,startfval=1,debug=debug,verbose=verbose)
    for (j in 2:R) {
      if (verbose) cat("Mult R&L: chain",j,"\n")
      r[[j]] <- calc.RL.0(data,startfval=j,
                             ## DON'T use pilot from previous runs
                             ## pilot=RL.max(r[[1]]),
                             verbose=verbose,debug=debug)
    }
    tot <- RL.max(r)
    resmat[i,] <- c(rseed,sapply(r,RL.max),tot)
  }
  dimnames(resmat) <- list(iteration=1:nrow(resmat),
                           c("RNG seed",paste("chain",1:R),"total"))
  resmat
}
  
calc.GR <- function(data,tot=20000,burn=1,verbose=FALSE,rseed=1001,
                       chainfrac=NULL) {
  mixsamp <- data$mixsamp
  sourcesamp <- data$sourcesamp
  R <- ncol(sourcesamp)
  if (!is.null(rseed)) set.seed(rseed)
  if (is.null(chainfrac)) chainfrac <- 1-1/R
  chain.start <- round(chainfrac*tot)
  g <- lapply(1:R,
        	function(z) {
	      if (verbose) cat("G&R: running chain",z,"of",R,"\n");
              gibbsC(startiter=burn+chain.start,
                     thin=1,
                     maxiter=tot+burn,
                     mixsamp=mixsamp,
                     sourcesamp=sourcesamp,
                     randseed=rseed,
                     startfval=z)})
  tmc <- mcmc.list(lapply(g,mcmc))
  if (verbose) cat("Starting G&R diagnostic ...\n")
  g <- gelman.diag(tmc,transform=TRUE)
  if (verbose) cat("                        ... done\n")
  g
}

calc.mult.GR <- function(data,n=10,tot=20000,burn=1,verbose=FALSE) {
  ## mixsamp <- data$mixsamp
  sourcesamp <- data$sourcesamp
  R <- ncol(sourcesamp)
  gresmat <- matrix(nrow=n,ncol=R+1)
  for (i in 1:n) {
    cat(i,"\n")
    rseed <- 1000+i
    g <- calc.GR(data,tot,burn,verbose,rseed)
    gresmat[i,] <- c(rseed,apply(g$confshrink,2,max))
  }
  dimnames(gresmat) <- list(iteration=1:n,
                            c("RNG seed","Max point est.","Max 97.5% quantile"))
  gresmat
}

boxplot.mixstock.est <- function(x,...) {
  arglist <- list(...)
  R <- arglist[[1]]$data$R
  if (!all(lapply(arglist,function(z)z$data$R)==R))
    stop("Different numbers of sources in estimates")
  ## get statistics from each type: summary??
}

hist.mixstock.est <- function(x,log="",...) {
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
                                           dimnames(x$data$sourcesamp)[[2]][i]),
                                         ...)))
}

tmcmc <- function(data,n.iter=20000,rseed=1001,n.burnin=floor(n.iter/2),
                  n.thin=max(1,floor(n.chains*(n.iter-n.burnin)/1000)),
                  verbose=FALSE,fprior=NULL,
                  contrib.only=TRUE,rptiter=-1,
                  outfile=NULL,
                  lang="C",a=NULL,gr=FALSE) {
  ## run MCMC chain
  ## FIXME: implement thinning???
  ## FIXME: return calculated value of a (prior strength)
  sourcesamp <- data$sourcesamp
  mixsamp <- data$mixsamp
  H <- length(mixsamp)
  R <- ncol(sourcesamp)
  n.chains <- R
  GR <- NA
  if (is.null(a)) {
    ## set a values according to p.bayes calculation
    harm.n <- sqrt(1/mean(1/apply(sourcesamp,2,sum)))
    beta <- p.bayes(sourcesamp)
    a <- beta/harm.n
  }
  ## sourcenames <- dimnames(sourcesamp)[[2]]
  set.seed(rseed)
  if (contrib.only) maxpar <- R
  else maxpar <- R+R*H
  chainfun <- function(z) {
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
      gibbsC(startiter=n.burnin,
             thin=n.thin,
             maxiter=n.iter,
             mixsamp=mixsamp,
             outfile=outfile,outfn=outfn,fprior=fprior,
             sourcesamp=sourcesamp,rptiter=rptiter,
             startfval=z,a=a)[,1:maxpar]
    }
    else gibbs(startiter=n.burnin,
               ## thin=1, ## thin not implemented?
               maxiter=n.iter,
               mixsamp=mixsamp,
               sourcesamp=sourcesamp,rptiter=rptiter,
               startfval=z,a=a)[,1:maxpar]
  }
  g <- lapply(1:R,chainfun)
  if(gr){ tmc <- mcmc.list(lapply(g, mcmc))
          GR.est <- gelman.diag(tmc, transform = TRUE)
          psrf <- ifelse(is.null(GR.est$confshrink),GR.est$psrf[, 2],GR.est$confshrink[, 2])
          GR.val <- max(psrf, na.rm = TRUE)
          GR <- GR.val 
        }
  if (verbose) cat("Chains done: trying to combine them\n")
  allchains <- as.data.frame(do.call("rbind",g),row.names=NULL)
  allchains.sum <- colMeans(allchains)
  if (!contrib.only) {
    source.freq <- matrix(allchains.sum[(R+1):length(allchains.sum)],ncol=R)
    dimnames(source.freq) <- dimnames(data$sourcesamp)
  }
  else
    source.freq <- NULL
  mixstock.est(fit=list(input.freq=allchains.sum[1:R],
                 source.freq=source.freq),
                 resample=allchains,
                 resample.index=rep(1:length(g),sapply(g,nrow)),
                 data=data,
                 method="mcmc",
                 boot.method="mcmc",
                 transf="part",GR=GR,
               prior=a)
}

mysumvec <- function(x,names=NULL,levels=c(0.95,0.9)) {
  if (is.null(names))
    if (!is.null(names(x)))
      names <- names(x)
    else if (!is.null(dimnames(x)[[2]]))
      names <- dimnames(x)[[2]]
    else
      names <- rep("",ncol(x))
  Qvec <- c((1-levels[1])/2,(1-levels[2])/2,(1+levels[2])/2,(1+levels[1])/2)
  y <- c(apply(x,2,mean,na.rm=TRUE),apply(x,2,median,na.rm=TRUE),
         apply(x,2,sd,na.rm=TRUE),
         apply(x,2,quantile,Qvec[1],na.rm=TRUE),
         apply(x,2,quantile,Qvec[2],na.rm=TRUE),
         apply(x,2,quantile,Qvec[3],na.rm=TRUE),
         apply(x,2,quantile,Qvec[4],na.rm=TRUE))
  names(y) <- t(outer(c("mean","median","sd",paste("Q",Qvec,sep="")),
                        names,paste,sep="."))
  y
}

mysum <- function(x,names=NULL) {
    if (is.null(names))
    if (!is.null(names(x)))
      names <- names(x)
    else if (!is.null(dimnames(x)[[2]]))
      names <- dimnames(x)[[2]]
    else names <- 1:ncol(x)
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
                    minmarks=3) {
  statnames <- c("mean","median","sd","Q02.5","Q05","Q95","Q97.5")
  nstats <- length(statnames)
  nsource <- length(true.contrib)
  sourcenames <- LETTERS[1:nsource]
  labnames <- sourcenames
  resarray <- array(dim=c(sim.n,mc.n,nsource*nstats))
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
  ## H <- nrow(true.freq)
  if (verbose) cat(R," sources ",totsamp," total samples\n")
## assume half of sample in mixed population, rest evenly divided betw. sources
  for (i in 1:sim.n) {
    cat("runsim: sim ",i,"\n")
    ok <- FALSE
    failct <- 0
    set.seed(1000+i)
    while (!ok) {
      sim0 <- simmixstock0(true.freq,true.contrib,totsamp/(2*R),totsamp/2,NULL)
      if (verbose) print(sim0)
      sim <- markfreq.condense(sim0)
      if (verbose) print(sim)
      if (sim$err==FALSE)
        ok <- (nrow(sim$sourcesamp)>minmarks)  ## more than one rare marker in sample?
      if (!ok) {
        failct <- failct+1
      }
      else {
        cat(failct," tries for good sample \n")
        if (common.only) {
          sim$mixsamp <- sim$mixsamp[1:2]
          sim$sourcesamp <- sim$sourcesamp[1:2,]
        } else
        if (rare.only) {
          sim$mixsamp <- sim$mixsamp[-(1:2)]
          sim$sourcesamp <- sim$sourcesamp[-(1:2),]
        }
      }
    } ## while !ok
#    print(sim)   
    for (j in 1:mc.n) {
      cat("   runsim: boot/MCMC run ",j,"\n")
      if (est=="MCMC")
        x <- tmcmc(sim,n.iter=nboot,rseed=1000+sim.n*i+j)$resample[,1:R]
      else if (est=="cml")
        x <- genboot(sim,method="cml",rseed=1000+sim.n*i+j,fuzz=fuzz,nboot=nboot,
                     rpt=bootrpt)$resample[,1:R]
      else if (est=="uml")
        x <- genboot(sim,method="uml",
                     rseed=1000+2*sim.n*i+j,fuzz=fuzz,nboot=nboot,
                     rpt=bootrpt)$resample[,1:R]
      
      else stop(paste("Unknown est type ",est))
      resarray[i,j,] <- mysumvec(x,names=labnames)
    }
  }
  resarray
}

cml <- function(x,start.type="lsolve",
                fuzz=0,bounds=1e-4,
                ndepfac=1000,
                method="L-BFGS-B",
                lower=NULL,
                upper=NULL,
                ndeps=NULL,
                control=NULL,
                debug=FALSE,
                transf="part",
                grad=NULL,  ## was cml.grad
                ...) {
  ## find contribs (by CML) that give MLE
  ## assume x is a list with elements sourcesamp, mixsamp
  R <- ncol(x$sourcesamp)
  H <- nrow(x$sourcesamp)
  ## mixfreq  <- x$mixsamp/sum(x$mixsamp)
  sourcefreq <- normcols(x$sourcesamp)
  tsourcefreq <- apply(sourcefreq,2,p.to.q,transf=transf)
  if (H==2) tsourcefreq <- matrix(tsourcefreq,nrow=1)
  if (is.character(start.type)) {
    start <- startvec(x,type=start.type,transf=transf,fuzz=fuzz,
                      cond=TRUE)
  } else start <- start.type
  ##  start <- startvec0(x$sourcesamp,x$mixsamp,type=start.type,lmin=fuzz)
  ##  start <- pmax(fuzz,pmin(1-fuzz,start))
  if(debug) {
    cat("Starting values:\n")
    print(start)
    cat("Starting log likelihood:\n")
    print(cml.lik(start,sourcefreq=tsourcefreq,data=x,transf=transf))
  }
  npar <- length(start)
  if (is.null(lower)) lower <- rep(bounds,npar)
  if (is.null(upper)) upper <- rep(1-bounds,npar)
  ndeps <- rep(bounds/ndepfac,npar)
  start <- pmin(pmax(start,lower),upper)
  if (!is.null(control)) {
    if (is.null(control$"ndeps"))
      control <- c(control,ndeps=ndeps)
  }
  else control <- list(ndeps=ndeps)
  m <- try(optim(par=start,
                 fn=cml.lik,
                 gr=grad,
                 sourcefreq=tsourcefreq,
                 data=x,
                 method=method,
                 transf=transf,
                 lower=lower,
                 upper=upper,
                 control=control,debug=debug,...),silent=TRUE)
  if ((!is.null(class(m))) && (class(m)=="try-error"))
    m <- list(coefficients=rep(NA,R-1),value=NA,convergence=NA)
  else class(m) <- "toptim"
  mixstock.est(m,data=x,R=R,H=H,M=1,method="cml",transf=transf)
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

## FIXME: remind myself what this does and what the limits should be
dcmat.a <- function(x,debug=FALSE) {
  n <- length(x)
  if (n>MAXMARK) stop(paste("exceeded hard-coded size limits:",
                                           "(R=",MAXMARK,")",sep=""))
  matrix(.C("dcmat",as.double(x),as.double(numeric(n*(n+1))),as.integer(n),
            as.integer(debug),PACKAGE="mixstock")[[2]],
         nrow=n+1)
}

## BUGGY??
cml.grad <- function(p,
                     sourcefreq,
                     data,
                     transf="full",
                     verbose=FALSE,
                     fulllik=NULL,
                     debug=FALSE) {
  ## sourcefreq should be _transformed_ source frequencies
  R <- ncol(sourcefreq)
  H <- nrow(sourcefreq)+1
  c0 <- unpackval(c(p,sourcefreq),transf=transf,R=R,H=H)
  M <- data$mixsamp  ## true samples
  m <- as.matrix(c0$source.freq) %*% c0$input.freq  ## calc freqs in mix
  j1 <- t(as.matrix(c0$source.freq)) %*% as.vector(M/m)
  r <-   -as.vector(t(dcmat.a(p)) %*% j1)
  if (debug) {
    cat("cml.grad",r,"\n")
    tstnd <- numderiv(p,cml.lik,
                      sourcefreq=sourcefreq,
                      transf=transf,
                      data=data)
    cat("cml.grad (numeric):",tstnd,"\n")
    devs <- r-tstnd
    reldev <- (r-tstnd)/tstnd
    cat("cml.grad, devs:",max(abs(devs),na.rm=TRUE),
        max(abs(reldev),na.rm=TRUE),"\n")
  }
  r
}

uml.grad <- function(p,data,
                     transf="full", ## kluges
                     debug=FALSE,verbose=FALSE) {
  R <- data$R
  H <- data$H
  x <- unpackval(p,R=R,H=H,transf=transf)  ## unpack
  c0 <- x$input.freq
  m <- x$source.freq %*% c0  ## calc freqs in mix
  mvec <- as.vector(data$mixsamp/m)
  mvec[!is.finite(mvec)] <- 0
  ## gradient on input contributions
  j1 <- t(x$source.freq) %*% mvec
  ## new ## j1 <- t(mvec %*% x$source.freq)
  g1 <- as.vector(-t(dcmat.a(p[1:(R-1)])) %*% j1)
  ## likelihood on mixed pool from source freqs
  ## need to keep matrix orientation consistent:
  ##  e.g., sources=cols,markers=rows
  fmat <- as.matrix(data$sourcesamp/x$source.freq)
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
    tstnd <- numderiv(p,
                      uml.lik,
                      data=data,
                      transf=transf)
    cat("uml.grad (numeric):",tstnd,"\n")
    devs <- r-tstnd
    reldev <- (r-tstnd)/tstnd
    cat("uml.grad, devs:",max(abs(devs),na.rm=TRUE),
        max(abs(reldev),na.rm=TRUE),"\n")
  }
  r
}

cml.lik <- function(p,
                    sourcefreq,
                    data,
                    transf=c("full","part","none"),
                    verbose=FALSE,
                    fulllik=TRUE,
                    debug=FALSE) {
  transf <- match.arg(transf)
  p <- c(p,sourcefreq)
  l <- loglik2.C(p=p,
                 sourcesamp=data$sourcesamp,
                 mixsamp=data$mixsamp,
                 cond=!fulllik,
                 verbose=verbose,
                 transf=transf)
  if (debug) cat("cml.lik:",p,l,"\n")
  l
}

uml <- function(x, method="direct",optmethod="L-BFGS-B",...) {
  switch(method,
         direct = uml.ds(x,method=optmethod,...),
         em = uml.em(x,...))
}

uml.em <- function(x,prec=1e-8,prior=1) {
  R <- x$R
  H <- x$H
  sourcename <- colnames(x$sourcesamp)
  sourcesamp <- x$sourcesamp
  mixsamp <- x$mixsamp
  mixsum  <- sum(mixsamp)
  fval <- matrix(rep(1/R, R),ncol=R)
  sourcefreqval <- apply(sourcesamp+prior, 2, function(z) z/sum(z))
  fold <- rep(10,R)
  while(sum(abs(fval-fold)) > prec){
    fold <- fval
    mix.freq <- as.vector(sourcefreqval %*% t(fval))
    w <- t(apply(sourcefreqval, 1, "*", fval))/mix.freq
    w[as.logical(match(w[, 1], NaN, nomatch = 0)), ] <- rep(0,R)
    rval <- w*matrix(rep(mixsamp,R),H,R)
    fval <- matrix(apply(rval, 2, sum)/mixsum,ncol=R)
    sourcefreqval <- apply(rval+sourcesamp, 2, function(z) z/sum(z))
  }
  fval <- as.vector(fval)
  names(fval) <- sourcename
  m <- list(input.freq=fval,source.freq=sourcefreqval)
  mixstock.est(m, data = x, R=R, H=H, M=1,method = "uml",
             transf="none", em = TRUE)
}

## UML by direct search
uml.ds <- function(x,
                   grad=NULL, ## uml.grad,
                   start.type="lsolve",
                   fuzz=0,
                   bounds=1e-4,
                   ndepfac=1000,
                   method="L-BFGS-B",
                   debug=FALSE,
                   control=NULL,
                   transf=c("part","full","none"),
                   ...) {
  R <- x$R
  H <- x$H
  transf <- match.arg(transf)
  start <- startvec(x,type=start.type,transf=transf,fuzz=fuzz)
  defndeps <- rep(bounds/ndepfac,length(start))
  if (!is.null(control)) {
    if (is.null(control$"ndeps"))
      control <- c(control,ndeps=defndeps)
  }
  else control <-   list(ndeps=defndeps)
  if (method=="L-BFGS-B") {
    lower=rep(bounds,length(start))
    upper=rep(1-bounds,length(start))
    if (transf=="full")
      warning("using bounds/L-BFGS-B and full transformation is probably wrong")
    m <- try(optim(par=start,
                   fn=uml.lik,
                   gr=grad,
                   data=x,
                   transf=transf,
                   debug=debug,
                   method=method,
                   lower=lower,upper=upper,
                   control=control,...),silent=TRUE)
  } else {  ## ugly! without bounds for non-L-BFGS-B code
    m <- try(optim(par=start,
                   fn=uml.lik,
                   gr=grad,
                   data=x,
                   transf=transf,
                   debug=debug,
                   method=method,
                   control=control,...),silent=TRUE)

  }
  if ((!is.null(class(m))) && (class(m) == "try-error")){
    m <- list(coefficients = rep(NA, ((R - 1) + R * (H - 1))),
              value = NA, convergence = NA)}
  class(m) <- "toptim"
  mixstock.est(m,data=x,R=R,H=H,M=1,transf=transf,method="uml")
}

uml.lik <- function(p,data,transf=c("full","part","none"),
                    verbose=FALSE,debug=FALSE) {
  transf <- match.arg(transf)
  if (debug) cat("args:",p,"\n")
  r <- loglik2.C(p,
                 sourcesamp=data$sourcesamp,
                 mixsamp=data$mixsamp,
                 transf=transf,
                 verbose=verbose)
  if (debug) cat("uml.lik:",r,"\n")
  r
}

mixstock.boot <- function(x,param=FALSE,condense=TRUE,save.freq=FALSE,
                        param.match="mean") {
  ## bootstrap source samples and marker samples
  ## sim-mixstock based on observed freqs
  if (!param) {
    mixfreq  <- x$mixsamp/sum(x$mixsamp)
    sourcefreq <- normcols(x$sourcesamp)
  } else {
    ## should observed sample frequencies match the mean, or the mode (==MLE), of
    ## the resampled frequencies????
    ## check ???
    if (param.match=="mean")
      alpha.plus <- 1
    else if (param.match=="mode")
      alpha.plus <- 0
    else
      stop("param.match not equal to mean or mode")
    if (any(x$sourcesamp==0) && alpha.plus==0)
      stop("can't sample dirichlet with zero values and zero prior")
    mixfreq <- rdirichlet(1,x$mixsamp+alpha.plus)
    sourcefreq <- apply(x$sourcesamp,2,function(x)rdirichlet(1,x+alpha.plus))
  }
  sourcesampsize <- apply(x$sourcesamp,2,sum)
  mixsamp <- as.vector(rmultinom(1,sum(x$mixsamp),mixfreq))
  sourcesamp <- matrix(nrow=nrow(x$sourcesamp),ncol=ncol(x$sourcesamp))
  for (i in 1:ncol(sourcefreq))
    sourcesamp[,i] <- rmultinom(1,sourcesampsize[i],sourcefreq[,i])
#  sourcesamp <- matrix(rmulti(1,sum(x$sourcesamp),as.vector(as.matrix(sourcefreq))),nrow=nrow(x$sourcesamp))
#  sourcesamp <- matrix(rmulti(2,sum(x$sourcesamp),as.vector(as.matrix(sourcefreq))),nrow=nrow(x$sourcesamp))
  names(mixsamp) <- names(x$mixsamp)
  dimnames(sourcesamp) <- dimnames(x$sourcesamp)
  data <- as.mixstock.data(list(sourcesamp=sourcesamp,mixsamp=mixsamp))
  if (condense)
    data <- markfreq.condense(data)
  if (param && save.freq) {
    attr(data,"sourcefreq") <- sourcefreq
    attr(data,"mixfreq") <- mixfreq
  }
  data
}


genboot <- function(x,method="cml",
                    nboot=1000,rseed=1001,
                    verbose=FALSE,fuzz=1e-3,
                    maxfail=100,pb=TRUE,
                    start.type="lsolve",
                    param=FALSE,
                    param.match="mean",
                    ndepfac=10000,
                    save.boot=FALSE,
                    print.boot=FALSE,...) {
  set.seed(rseed)
  H <- nrow(x$sourcesamp)
  R <- ncol(x$sourcesamp)
  marknames <- dimnames(x$sourcesamp)[[1]]
  sourcenames <- dimnames(x$sourcesamp)[[2]]
  fitnames <- c(paste("contrib",sourcenames,sep="."))
  if (method=="cml") {
    basefit = cml(x,fuzz=fuzz,ndepfac=ndepfac,...)
    ## npar <- R-1
    noutpar <- R
  } else {
    basefit <- uml(x,start.type=start.type,fuzz=fuzz)
    ## npar <- (R-1)+R*(H-1)
    noutpar <- R+R*H
    fitnames <- c(fitnames,
                  paste("sourcefreq",
                        outer(sourcenames,
                              marknames,paste,sep="."),sep="."))
  }
  nbootres <- noutpar+2
  boot.results <- matrix(nrow=nboot,ncol=nbootres)
  if (save.boot) boot.data <- list()
  else boot.data <- NULL
  if (pb) tpb <- txtProgressBar(max=nboot)
  for (i in 1:nboot) {
    ## if (i %% rpt == 0) cat("genboot ",i,"\n")
    if (pb) setTxtProgressBar(tpb,i)
    ok <- FALSE
    failct <- 0
    while (!ok && failct<maxfail) {
      ## sample until at least 2 OK markers in mixed stock
      x0 <- mixstock.boot(x,param=param,param.match=param.match)
      ok <- (length(x0$mixsamp)>1)
      failct <- failct+1
    }
    if (print.boot) cat(as.vector(x0$mixsamp),
                        as.vector(as.matrix(x0$sourcesamp)),"\n")
    if (failct<maxfail) {
      ok <- FALSE
      failct <- 0
 #    while (!ok & failct<maxfail) {
        ## why keep trying here??
      if (method=="cml") {
        fit <- cml(x0)
      } else {
        fit <- uml(x0)
      }
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
      boot.results[i,] <- c(c(tmpfit,rep(NA,noutpar-length(tmpfit))),
                            -logLik(fit$fit),fit$fit$convergence)
      if (save.boot) boot.data[[i]] <- x0
    }
  }
  if (pb) close(tpb)
  dimnames(boot.results) <- list(NULL,c(fitnames,"NLL","Convergence"))
  mixstock.est(fit=basefit$fit,data=x,
             resample=boot.results,
             method=method,
             transf=basefit$transf,
             boot.method=ifelse(param,
               paste("parametric",param.match,sep="."),
               "nonparametric"),
             boot.data=boot.data)
}


#######################

## expand mixsamp data
expand.bugs.data <- function (mixsamp) {
  if (!is.null(dim(mixsamp))) {  ## matrix or data frame
    ## recursive call to expand.bugs.data
    Tm <- colSums(mixsamp)
    maxn <- max(Tm)
    if (min(Tm)==maxn) {
      xlist <- lapply(split(mixsamp,col(mixsamp)),expand.bugs.data) ## kluge
    } else {
      xlist <- apply(mixsamp,2,expand.bugs.data)
      xlist <- lapply(xlist,
                      function(x) {if (nrow(x)==maxn) x else rbind(x,matrix(0,nrow=maxn-nrow(x),ncol=ncol(x)))})
    }
    xarr <- abind(xlist,along=3)
    return(aperm(xarr,c(3,1,2)))  ## rearrange to (m,n,h)
  }
  if (!is.null(names(mixsamp))) {
    names(mixsamp) <- gsub("/", "", names(mixsamp))
  } else names(mixsamp) <- paste("H",1:length(mixsamp),sep="")
  x1 <- factor(rep(names(mixsamp), mixsamp), levels = names(mixsamp)); x1 ## codetools kluge 
  x2 <- model.matrix(~x1 - 1)
  attr(x2, "assign") <- attr(x2, "contrasts") <- NULL
  x2
}


## FIXME: common thin/burn/tot interpretation for tmcmc and pm.wbugs
## FIXME: don't know how to set seed for WinBUGS?
pm.wbugs <- function(x,
                     n.iter=20000,n.burnin=floor(n.iter/2),
                     n.chains=x$R,
                     n.thin=max(1,floor(n.chains*(n.iter-n.burnin)/1000)),
                     ...) {
  pm.bugscode <- system.file(package = "mixstock", "bugs", "pellamasuda.bug")
  ## variables repeated at end of line are codetools kluges
  sourcesamp <- t(x$sourcesamp); sourcesamp
  mixsamp <- expand.bugs.data(x$mixsamp); mixsamp
  R = ncol(x$sourcesamp)
  harm.n <- sqrt(1/mean(1/apply(x$sourcesamp,2,sum)))
  beta <- p.bayes(x$sourcesamp)
  a <- beta/harm.n
  ybar <- rowMeans(apply(x$sourcesamp,2,function(z)z/sum(z)))
  sourceprior <- t(a*sqrt(harm.n)*matrix(rep(ybar,R),ncol=R)); sourceprior
  contprior <- rep(1/R,R); contprior
  sourceprop.obs <- t(normcols(x$sourcesamp))
  T <- colSums(x$sourcesamp); T
  H <- x$H; H
  N <- sum(x$mixsamp)
  data <- list("contprior","sourceprior","mixsamp","sourcesamp","T","R","H","N")
  parameters <- c("theta")
  basetheta <-  rep(0.05/(R-1),R)
  inits <- list()
  for (i in 1:R) {
    fval <- basetheta
    fval[i] <- 0.95
    Zrand <- sample(1:R,size=N,replace=TRUE,prob=fval)
    inits[[i]] <- list(Z=Zrand,pi=sourceprop.obs)
  }
  ## convert from tmcmc to BUGS defaults
  b1 = bugs(data,inits,parameters,model.file=pm.bugscode,
    n.chains=R,n.iter=n.iter,
    n.burnin=n.burnin,n.thin=n.thin,...)
  b1
}

## write many-to-many model in original (T. Okuyama) format
write.TO.bugscode = function(fn,MIX) {
  cat("model {\n",file=fn)
  for (m in 1:MIX) {
    bstr = paste("for(i in 1:(Tm[",m,"])) {\n",
      "   mixsamp",m,"[i,1:H] ~ dmulti(pi[Z",m,"[i],],1);\n",
      "   Z",m,"[i] ~ dcat(theta[",m,",1:R]);\n",
      "}\n",sep="")
    cat(bstr,file=fn,append=TRUE)
  }
  block2 = c("# model for pi (the multinomial parameter)",
    "for(i in 1:R){sourcesamp[i,1:H] ~ dmulti(pi[i,1:H],T[i])}",
    ## "for(i in 1:R){T[i] <- sum(sourcesamp[i,])}",
    "for(j in 1:R){pi[j,1:H] ~ ddirch(beta[j,])",
    "for(k in 1:H){beta[j,k] <- fp[k]}}",
    "",
    "for(j in 1:R){",
    "for(k in 1:(MIX+1)) {",
    "	div[j,k] <- delta[j,k]/sum(delta[j,])",
    "	delta[j,k] ~ dgamma(dp[j,k],1)",
    "}",
    "}",
    "")
  cat(block2,file=fn,append=TRUE,sep="\n")
  cat("for(i in 1:R){",file=fn,append=TRUE)
  for (m in 1:MIX) {
    bstr = paste("   DERIV[",m,",i] <- div[i,",m,"]*sourcesize[i]",sep="")
    cat(bstr,file=fn,append=TRUE,sep="\n")
  }
  block3 = c("}",
    "",
    "for(i in 1:MIX){",
    "   mixsize[i] <- sum(DERIV[i,])",
    "   rmixsize[i] <- mixsize[i]/sum(mixsize[])",
    "}",
    "",
    "for(j in 1:MIX){",
    "   for(i in 1:R){",
    "      theta[j,i] <- DERIV[j,i]/sum(DERIV[j,])",
    "   }",
    "}",
    "")
##    "for(i in 1:(MIX+1)){dp[i] <- 1}")  ## now set as data
  cat(block3,file=fn,append=TRUE,sep="\n")
  cat("}",file=fn,append=TRUE,sep="\n")
}

mm.wbugs <- function(x,
                     sourcesize,
                     n.iter=20000,n.burnin=floor(n.iter/2),
                     n.chains=x$R,
                     n.thin=max(1,floor(n.chains*(n.iter-n.burnin)/1000)),
                     files.only=FALSE,
                     inittype=c("dispersed","random"),
                     bugs.code=c("TO","BB"),
                     returntype=c("mixstock","coda","bugs"),
                     pkg=c("WinBUGS","JAGS"),
                     mixprior=1,
                     which.init,debug=FALSE,
                     working.directory,...) {
  if (missing(working.directory)) working.directory <- tempdir()
  pkg <- match.arg(pkg)
  inittype <- match.arg(inittype)
  returntype <- match.arg(returntype)
  bugs.code <- match.arg(bugs.code)
  sourcesamp <- t(x$sourcesamp)
  R <- ncol(x$sourcesamp)
  H <- nrow(x$sourcesamp); H ## codetools kluge
  MIX <- ncol(x$mixsamp)
  mm.bugscode <- paste(tempfile(paste("mm",bugs.code,"bugs",sep="_"),
                                tmpdir=working.directory),"txt",sep=".")
  if (bugs.code=="TO") {
    write.TO.bugscode(mm.bugscode,MIX=MIX)
    mixsamplist <- as.list(paste("mixsamp",1:MIX,sep=""))
    for (i in 1:MIX) {
      assign(paste("mixsamp",i,sep=""),expand.bugs.data(x$mixsamp[,i]))
    }
  } else {
    file.copy(system.file(package = "mixstock", "bugs", "manymany.txt"),
              mm.bugscode)
    mixsamplist <- list("mixsamp")
    mixsamp <- expand.bugs.data(x$mixsamp)
    mixsamp ## codetools kluge
  }
  ## n.chains <- R
  harm.n <- sqrt(1/mean(1/colSums(sourcesamp)))
  beta <- p.bayes(sourcesamp)
  a <- beta/harm.n
  ybar <- colMeans(normcols(sourcesamp))
  sourceprior <- t(a * sqrt(harm.n) * matrix(rep(ybar, R), 
                                             ncol = R))
  ## contprior <- rep(1/R, R)
  sourceprop.obs <- t(normcols(sourcesamp))
  if (missing(sourcesize)) sourcesize <- rep(1,R)
  attr(sourceprop.obs,"scaled:scale") <- NULL
  T <- rowSums(sourcesamp); T ## codetools kluge
  Tm <- colSums(x$mixsamp)
  maxTm <- max(Tm)
  fp <- sourceprior[1,]  ## haplotype freq prior
  dp <- mixprior
  ## following could be done more compactly but maybe clearer this way
  if (length(dp)==1) {  ## scalar: expand
    dp <- matrix(dp,nrow=R,ncol=MIX+1)
  } else if (length(dp)==MIX+1) {  ## single row: expand
    dp <- matrix(dp,nrow=R,ncol=MIX+1,byrow=TRUE)
  } else if (!all(dim(dp)==c(R,MIX+1))) stop("prior for source-to-mixed contributions is the wrong shape")
  ## dp <- rep(dp,length.out=MIX+1) ## replicate mixed prior if necessary
  c(fp,dp) ## codetools kluge
  data <- c(list("Tm",
                 "sourcesamp","sourcesize",
                 "R","H","MIX","fp","dp","T"),
            mixsamplist)
  if (bugs.code=="BB" && pkg=="JAGS") data <- c(data,"maxTm")
  parameters <- c("theta","div")
  ## initial values:
  if (inittype=="random") {
    inits <- NULL
  } else {
    inits <- list()
    ## we run into problems trying to
    ## do, e.g., ... ?
    basetheta <- rep(0.05/(R - 1), R)
    if (missing(which.init)) {
      if (n.chains==R) {
        which.init <- 1:R
      } else which.init <- sort(sample(1:R,size=n.chains,replace=FALSE))
    }
    for (i in 1:n.chains) {
      w <- which.init[i]
      if (w>0) {
        ## set initial values for hap freqs in rookeries (= observed prop),
        ##   ##  origins (=95% in chain i)
        fval <- basetheta
        fval[w] <- 0.95
      } else {
        stop("stub: should set chain to UML solution ...\n")
      }
      Zrand <- matrix(sample(1:R, size = MIX*max(Tm), replace = TRUE, prob = fval),
                      nrow=MIX)
      if (bugs.code!="TO") {
        inits[[i]] <- list(Z=Zrand)
      } else {
        initlist <- list()
        for (j in 1:MIX) {
          Zname <- paste("Z",j,sep="")
          initlist<- c(initlist,list(Zrand[j,1:Tm[j]]))
          names(initlist)[j] <- Zname
        }
        inits[[i]] <- initlist
      } ## TO code
    } ## loop over number of chains
      ## pi = t(sourceprop.obs),
      ## TO DO: try adjusting hap freq start, or using one or the other
  } ## inittype != random
  if (pkg=="WinBUGS") {
    elapsed <- system.time(b1 <- bugs(data,inits,parameters,mm.bugscode,
                                      n.chains=n.chains,
                                      n.iter=n.iter,n.burnin=n.burnin,
                                      n.thin=n.thin,
                                      working.directory=working.directory,
                                      debug=debug,...))[3]
  } else {
    ## browser()
    elapsed <- system.time(b1 <- jags(data,inits,parameters,mm.bugscode,
                                      n.chains=n.chains,
                                      n.iter=n.iter,n.burnin=n.burnin,
                                      n.thin=n.thin,
                                      working.directory=working.directory,
                                      ...)$BUGSoutput)[3]
  }
  return(switch(returntype,
                mixstock=as.mixstock.est.bugs(b1,data=x,time=elapsed),
                coda=as.mcmc.bugs(b1),
                bugs=b1))
}

## WinBUGS code with hierarchical size covariates (not many-to-many)
sourcesize.wbugs <- function(x,
                             n.iter=20000,n.burnin=floor(n.iter/2),
                             n.chains=x$R,
                             n.thin=max(1,floor(n.chains*(n.iter-n.burnin)/1000)),
                             ...) {
  stop("stub")
  ## sourcesize.bugscode <- system.file(package = "mixstock", "bugs", "sourcesize.bug")
}

as.mcmc.bugs <- function(x) {
  if (x$n.chains>1) {
    z <- list()
    for (i in 1:x$n.chains) {
      z[[i]] <- mcmc(x$sims.array[,i,],start=1,thin=x$n.thin)
    }
    class(z) <- "mcmc.list"
  } else {
    z <- mcmc(x$sims.matrix,start=1,thin=x$n.thin)
  } 
  return(z)
}
## higher-level functions for mixstock analysis

## simulate mixstocks
simmixstock0 <-
  function(sourcefreq,input.freq,sourcesampsize,mixsampsize,rseed=NULL) {
  if (!is.null(rseed))
    set.seed(rseed)
  sourcesampsize <- round(sourcesampsize,0)
  mixsampsize <- round(mixsampsize,0)
  H <- nrow(sourcefreq)
  R <- ncol(sourcefreq)
 ## simulate (multinomial samples from sources and mixed pop)
  sourcesamp.r <- apply(sourcefreq,2,function(p)rmultinom(1,sourcesampsize,p))
  mix.freq <- as.vector(sourcefreq %*% input.freq)
  mixsamp.r <- as.vector(rmultinom(1,mixsampsize,mix.freq))
  dimnames(sourcesamp.r) <- mixstock.dimnames(H,R)
  names(mixsamp.r) <- dimnames(sourcesamp.r)[[1]]
  return(as.mixstock.data(list(mixsamp=mixsamp.r,sourcesamp=sourcesamp.r)))
}

simmixstock1 <-
  ## general-purpose simulation & analysis
function(sampsize=NULL,
         true.freq=matrix(c(0.65, 0.33, 0.01, 0.01, 0.33, 0.65, 0.01, 0.01),
           ncol = 2),
         true.contrib=c(0.9, 0.1),
         boot=FALSE,param=FALSE,
         data=NULL, rseed=1004, nboot=1000, chainlen=NULL,
         ests=c("cmlboot.nonpar",
           "cmlboot.par","umlboot.nonpar","umlboot.par","mcmc"),
         verbose=FALSE,
   contrib.only=FALSE) {
  ## generate "data"
  if (!is.null(data)) {
    x <- data
    if (!is.null(sampsize)) {
      datasampsize <- sum(x$sourcesamp)+sum(x$mixsamp)
      x$sourcesamp <- round(x$sourcesamp*sampsize/datasampsize)
      x$mixsamp <- round(x$mixsamp*sampsize/datasampsize)
      if (boot)
        x <- mixstock.boot(x,param=param)
    }
  } else {
    R <- ncol(true.freq)
    mixsampsize <- sampsize/2
    sourcesampsize <- sampsize/(2*R)
    x <- simmixstock0(true.freq, true.contrib, sourcesampsize = sourcesampsize, 
                 mixsampsize=mixsampsize, rseed = rseed)
  }
  x <- markfreq.condense(x)
  if (verbose) {
    cat("Data for analysis:\n")
    print(x)
  }
  results <- list()
  if ("cml" %in% ests) {
    if (verbose) cat("cml\n")
    results <- c(results,cml=list(cml(x)))
  }
  if ("uml" %in% ests) {
    if (verbose) cat("uml\n")
    results <- c(results,uml=list(uml(x)))
  }
  if ("cmlboot.nonpar" %in% ests) {
    if (verbose) cat("cmlboot.nonpar\n")
    results <- c(results,cmlboot.nonpar=list(genboot(x,type="cml",nboot=nboot,param=FALSE)))
  }
  if ("cmlboot.par" %in% ests) {
    if (verbose) cat("cmlboot.par\n")
    results <- c(results,cmlboot.par=list(genboot(x,type="cml",nboot=nboot,param=TRUE)))
  }
  if ("umlboot.nonpar" %in% ests) {
    if (verbose) cat("umlboot.nonpar\n")
    results <- c(results,umlboot.nonpar=list(genboot(x,type="uml",nboot=nboot,param=FALSE)))
  }
  if ("umlboot.par" %in% ests) {
    if (verbose) cat("umlboot.par\n")
    results <- c(results,umlboot.par=list(genboot(x,type="uml",nboot=nboot,param=TRUE)))
  }
  if ("mcmc" %in% ests) {
    if (verbose) cat("mcmc\n")
    if (is.null(chainlen))
      chainlen <- mcmc.chainlength.est(x,mult=5)
    results <- c(results,mcmc=list(tmcmc(x,n.iter=chainlen,
   contrib.only=contrib.only)))
  }
  results
}

simmixstock2 <- function(sourcefreq,       ## haplotype freqs in sources: dim (H,S)
                         destmat,          ## frequencies of sources in mixed: dim(S,M)
                         sourcesize,       ## relative source sizes: length S
                         sourcesampsize,   ## 
                         mixsampsize,
                         nmark, nsource, nmix,
                         rseed = NULL,
                         condense = TRUE) 
{
  nmat <- function(x) { sweep(x,2,colSums(x),"/") }
  if (!is.null(rseed)) 
    set.seed(rseed)
  sourcesampsize <- round(sourcesampsize)
  if (missing(nsource)) nsource <- length(sourcesampsize)
  mixsampsize <- round(mixsampsize)
  if (missing(nmix)) nmix <- length(mixsampsize)
  if (missing(nmark)) {
    if (missing(sourcefreq))
      stop("must specify source  frequencies or number of marker classes")
    nmark <- nrow(sourcefreq)
  }
  if (missing(sourcefreq)) sourcefreq <- nmat(matrix(runif(nmark*nsource),nrow=nmark,ncol=nsource))
  if (missing(destmat)) destmat <- nmat(matrix(runif(nmix*nsource),ncol=nmix,nrow=nsource))
  if (missing(sourcesize)) sourcesize <- rep(1,nsource)
  sourcesamp.r <- apply(sourcefreq, 2, rmultinom,
                        n=1, size=sourcesampsize)
  mix.freq <- nmat(sweep(sourcefreq,2,sourcesize,"*") %*% destmat)
  ## ?? did this work before ??
  ## mix.freq <- nmat(t(t(sourcefreq)*sourcesize) %*% t(destmat))
  mixsamp.r <-  apply(mix.freq, 2, rmultinom,
                      n=1, size=mixsampsize)
  dimnames(sourcesamp.r) <- mixstock.dimnames(nmark, nsource)
  dimnames(mixsamp.r) <- list(marker=dimnames(sourcesamp.r)[[1]],mix=paste("M",1:nmix,sep=""))
  m <- mixstock.data(sourcesamp=sourcesamp.r,mixsamp=mixsamp.r)
  if (condense) m <- markfreq.condense(m)
  return(m)
}

## code for reading from/writing to Masuda&Pella control files

get.bot <- function(fn) {
   fn <- paste(fn,".BOT",sep="")
   x <- read.table(fn)
   x
}


get.frq <- function(fn) {
   fn <- paste(fn,".FRQ",sep="")
   x <- read.table(fn)
   x
}

get.ctl <- function(fn) {
   fn <- paste(fn,".CTL",sep="")
   fnconn <- file(fn,open="rt")
   x <- readLines(fnconn,n=19)  ## read file as character strings
   title <- x[1]
   length <- as.numeric(x[8])
   R <- as.numeric(x[9])
   ranseed <- as.numeric(x[11:13])
## is there a way to do this with read.table?
   tab <- scan(fnconn,what=list(double(1),double(1),"",double(1)))
   list(title=title,tot=length,R=R,ranseed=ranseed,fprior=tab[[2]],
         sourcename=tab[[3]],startval=tab[[4]])
}

get.bse <- function(fn) {
 fn <- paste(fn,".BSE",sep="")
 tab <- read.table(fn)
 sourcesamp <- t(tab[,-(1:3)])
 list(sourcesamp=sourcesamp)
}

get.mix <- function(fn) {
  fn <- paste(fn,".MIX",sep="")
  tab <- read.table(fn)
  mixsamp <- apply(tab,2,sum)
  list(mixsamp=mixsamp)
}

put.ctl <- function(fn,which=1,tot,ranseed=c(899271,4480026,90092812),
                 thin=1,title="auto-gen R input",sourcenames,fprior,startval,H) {
   ctl.fn <- paste(fn,which,".CTL",sep="")
   sink(ctl.fn)
   cat(title,"\n")
   sapply(c("BSE","MIX"),
        function(ext) cat(paste(fn,".",ext,"\n",sep="")))
   sapply(c("SUM","BOT","FRQ","B01"),
        function(ext) cat(paste(fn,which,".",ext,"\n",sep="")))
   cat(tot,"\n")
   R <- length(sourcenames)
   cat(R,"\n")
   cat(1,"\n")  ## number of characters -- always 1
   cat(ranseed,sep="\n")
   cat(thin,"\n") ## composition thin
   cat(thin,"\n") ## frequency thin
   cat(paste("(",paste(rep(c("I1","X"),H),collapse=","),")\n",sep=""))
   cat(paste("(I3,2X,I1,3X,",  ## number and num. of characters
                 paste(rep(c("I3","X"),H+1),collapse=","),")\n",sep=""))
                 ## total and 
   cat("T T T T\n")
   cat(formatC(1,width=3),
       formatC(H,width=3),
       formatC("F",width=3,flag="+"),  ## Hardy-Weinberg,
       "  ",
       formatC("mtDNA",width=18,flag="+"),"\n",sep="")
   for (i in 1:R) {
     ## FORTRAN format specified: I3,1X,F7.6,2X,A18,1X,F7.6
     ## can't quite match FORTRAN FX.Y format: digits=n always
     ## gives a minimum width of n+2 (decimal point+leading zero)
      cat(formatC(i,width=3),
          " ",
          formatC(fprior[i],width=6,digits=5,format="f"),
          "  ",
          formatC(sourcenames[i],width=18),
          " ",
          formatC(startval[i],width=6,digits=5,format="f"),
          "\n",
          sep="")
    }
  sink()
}

put.bse <- function(fn,sourcesamp) {
  fn <- paste(fn,".BSE",sep="")
  sink(fn)
  sapply(1:ncol(sourcesamp),
       function(i){
          cat(formatC(i,format="d",width=3),
              "  ",
              1,
              "   ",
              formatC(sum(sourcesamp[,i]),width=3,flag="-")," ",
              paste(formatC(sourcesamp[,i],width=3,flag="-"),collapse=" "),
              " ",  ## extra space necessary if "X" present at end of format line in CTL
              "\n",
              sep="")})
   sink()
}

put.mix <- function(fn,mixsamp) {
  fn <- paste(fn,".MIX",sep="")
  tmpvec <- function(which,H) {
     v <- rep(0,H)
     v[which] <- 1
     paste(formatC(v,format="d",width=1),collapse=" ")
  }
  sink(fn)
  H <- length(mixsamp)
  for (i in 1:H)
    if (mixsamp[i]>0)
      for(j in 1:mixsamp[i])
        cat(tmpvec(i,H),"\n",sep="")
  sink()
}

put.mp <- function(data,tot=25000,title="Input file auto-gen from R",
                fn="test",ranseed=c(899271,4480026,90092812)) {
  put.mix(fn,data$mixsamp)
  put.bse(fn,data$sourcesamp)
  R <- data$R
  fprior <- rep(1/R,R)
  sourcenames <- dimnames(data$sourcesamp)[[2]]
  for (r in 1:R) {
    startval <- rep(0.05/(R-1),R)
    startval[r] <- 0.95
    put.ctl(fn,title=paste(title," [source ",r,": ",
      sourcenames[r],"]",sep=""),which=r,tot=tot,ranseed=ranseed,
                 thin=1,sourcenames=sourcenames,fprior=fprior,startval=startval,
                 H=data$H) 
   }
}

get.mp.input <- function(fn,which=1) {
  ctl <- get.ctl(paste(fn,which,sep=""))
  bse <- get.bse(fn)
  mix <- get.mix(fn)
  sourcesamp <- bse$sourcesamp
  dimnames(sourcesamp) <- list(NULL,ctl$sourcename)
  list(sourcesamp=sourcesamp, mixsamp=mix$mixsamp, a = 1, startiter=1, maxiter=ctl$tot, 
       startfval = ctl$startval, thin = 1, fprior = ctl$fprior, rptiter = -1)
}

as.mixstock.est.bugs <- function(object,data=NULL,time) {
  X <- list()
  mm <- (!is.null(object$sims.list$div)) ## ?? is this adequate?
  X$fit <- list(input.freq=colMeans(object$sims.list$theta),source.freq=NULL)
  if (!is.null(object$sims.list$div)) {
    X$fit <- c(X$fit,list(sourcectr.freq=colMeans(object$sims.list$div)))
  }
  X$resample <- object$sims.matrix[,!colnames(object$sims.matrix) %in% "deviance"]
  ##                                   1:(ncol(object$sims.matrix)-1)]
  X$resamplist <- object$sims.list
  X$R <- ncol(X$fit$input.freq)
  X$M <- nrow(X$fit$input.freq)
  X$data <- data
  X$method <- "mcmc"
  X$boot.method <- "mcmc"
  if (!missing(time)) X$time <- time
  class(X) <- "mixstock.est"
  X$H <- NA
  if (!is.null(data)) {
    X$H <- data$H
    if (mm) {
      ## need data to get source and mix names
      mnames <- colnames(data$mixsamp)
      rnames <- colnames(data$sourcesamp)
      ## fixed bug: names were backwards
      dimnames(X$fit$input.freq) <- list(mnames,rnames)
      dimnames(X$fit$sourcectr.freq) <- list(rnames,c(mnames,"Unknown"))
      dimnames(X$resample) <- list(NULL,
                                   ## n.b. sims.matrix is ordered ALPHABETICALLY
                                   c(outer(rnames,c(mnames,"Unk"),paste,sep="."),  ## 'div'
                                     outer(mnames,rnames,paste,sep=".")))          ## 'theta'
      dimnames(X$resamplist$theta) <- list(NULL,mnames,rnames)
      dimnames(X$resamplist$div) <- list(NULL,rnames,c(mnames,"Unknown"))
    }
  }
  X
}

as.mixstock.est.mcmc <- function(object,data=NULL) {
  stop("stub")
  X <- list()
  X$fit <- list(input.freq=object,source.freq=NULL)
  X$resample <- object$sims.list$theta
  X$R <- length(X$fit$input.freq)
  X$data <- data
  X$method <- "mcmc"
  X$boot.method <- "mcmc"
  class(X) <- "mixstock.est"
  X$H <- NA
  if (!is.null(data)) {
    X$H <- data$H
    if (is.matrix(X$fit$input.freq)) {
      mnames <- colnames(data$mixsamp)
      rnames <- colnames(data$sourcesamp)
      ## fixed bug: names were backwards
      dimnames(X$resample) <- list(NULL,mnames,rnames)
      dimnames(X$fit$input.freq) <- list(mnames,rnames)
    }
  }
  X
}
  

nsource <- function(object) {
  object$R
}
nmix <- function(object) {
  if (is.null(object$M)) 1 else object$M
}
nmark <- function(object) {
  object$H
}
"nsource<-" <- function(object,value) {
  object$R <- value
  object
}
"nmix<-" <- function(object,value) {
  object$M <- value
  object
}
"nmark<-" <- function(object,value) {
  object$H <- value
  object
}

##  NOTE: S3 method for coefplot2 package, but we don't define
## default to avoid dependence on coefplot2
coeftab.mixstock.est <- function(object,clevel=c(0.5,0.95),...)  {
  ## loc should be a function with an na.rm argument (mean,median)
  ## FIXME:
  pvec <- c((1-clevel)/2,(1+clevel)/2)
  est<- coef(object)$input.freq
  sd <- rep(NA,length(est))
  Qvec <- matrix(NA,nrow=length(est),ncol=2*length(clevel),
                 dimnames=list(names(est),as.character(sort(paste(pvec,"%",sep="")))))
  if (!is.null(object$resample)) {
    rr <- object$resample
    rr <- rr[,!colnames(rr) %in% c("NLL","Convergence")]
    ## est <- apply(rr,2,loc,na.rm=TRUE)
    sd <- apply(rr,2,sd,na.rm=TRUE)
    Qvec <-  t(apply(rr,2,quantile,probs=pvec))
    Qvec <- Qvec[,order(pvec)]
  }
  cc <- cbind(Estimate=est,"Std. Error"=sd,Qvec)
  rownames(cc) <- gsub("^contrib\\.","",rownames(cc))
  class(cc) <- "coeftab"
  cc
}

## designed for MCMC only
xyplot.mixstock.est <- function(x,data,...) {
  ## convert allchains back into an mcmc list
  m <- as.mcmc.list(lapply(split.data.frame(x$resample,x$resample.index),as.mcmc))
  xyplot(m,...)
}

as.data.frame.mixstock.est <- function(x,row.names, optional, ...) {
    if (!missing(row.names)) warning("'row.names' argument ignored")
    if (!missing(optional)) warning("'optional' argument ignored")
    ## unpack from $resample
    mm <- (!is.null(x$resamplist$div)) ## ?? is this adequate?
    if (!mm) stop("as.data.frame is only implemented for many-to-many fits")
    coefs <- ldply(coef(x),function(x) if (is.null(x)) NULL else melt(x))
    ci <- confint(x)
    ci2 <- transform(as.data.frame(ci),
                    Var1=gsub("\\..+$","",rownames(ci)),
                    Var2=gsub("Unk","Unknown",gsub("^.+\\.","",rownames(ci))))
    setNames(merge(coefs,ci2),c("to","from","type","est","lwr","upr"))
}
