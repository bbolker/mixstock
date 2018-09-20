library(mixstock)
Z <- simmixstock2(nsource=4,nmark=5,nmix=3,
                 sourcesize=c(4,2,1,1),
                 sourcesampsize=rep(25,4),
                 mixsampsize=rep(30,3),rseed=1001)               

if (FALSE) {
  ## argh, don't actually try to run this on r-forge, requires JAGS
Zfit0 <- mm.wbugs(Z,sourcesize=c(4,2,1,1),pkg="JAGS",n.iter=200)
Zfit1 <- mm.wbugs(Z,sourcesize=c(4,2,1,1),pkg="JAGS",n.iter=200,
                                    mixprior=0.5)
Zfit2 <- mm.wbugs(Z,sourcesize=c(4,2,1,1),pkg="JAGS",n.iter=200,
                                    mixprior=c(3,1,1,1))
mp <- matrix(sample(1:5,replace=TRUE,size=16),nrow=4)
Zfit3 <- mm.wbugs(Z,sourcesize=c(4,2,1,1),pkg="JAGS",n.iter=200,
                                    mixprior=mp)

}
