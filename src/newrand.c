#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <R_ext/PrtUtil.h>
#define MAXNCAT 500L

void genmul(int n, double *p, int ncat,int *ix, int round)
/***************************************************************
  GENerate an observation from the MULtinomial distribution
  arguments: 
     n Number of events that will be classified into one of
           the categories 1..ncat
     p Vector of probabilities.  p[i] is the probability that
       an event will be classified into category i. p[i]
        must be [0,1], sum(p[i])<=1.0.
        Only the first (ncat-1) p[i] must be defined
           since p[ncat-1] is 1.0 minus the sum of the first
           NCAT-1 P(i).
     ncat Number of categories.  Length of P and IX.
     ix  Observation from multinomial distribution.  All IX(i)
            will be nonnegative and their sum will be N.
     round round sums>1 down to 1?

     Algorithm from page 559 of: Devroye, Luc, Non-Uniform Random Variate
Generation.  Springer-Verlag, New York, 1986.

modified for use inside R from randlib.c version 1.3
written by Barry W. Brown, James Lovato, Kathy Russell, and John Venier
http://odin.mdacc.tmc.edu/anonftp/
 
*******************************************************************/
{
  double prob,ptot,sum;
  int i,icat,ntot;
  double sum_tol=1.0e-4, sum_max=0.99999999;
  double pextra;
  int noisy=0;

  /* input checking */
  if(n < 0) error("n < 0");
  if(ncat <= 1) error("ncat <= 1");
  ptot = 0.0;
  for(i=0; i<ncat-1; i++) {
    if(*(p+i) < 0.0) error("some P(i) < 0");
    if(*(p+i) > 1.0) {
      for (i=0; i<ncat; i++)
	Rprintf("%f ",p[i]);
      Rprintf("\n");
      error("some P(i) > 1");
    }
    ptot += *(p+i);
  }
  /* if sum of probabilities is > sum_max but < sum_max+sum_tol,
   * divide by (sum_max+sum_tol) to bring the total down 
   * is there a better solution??? */
  if (round==1) {
    pextra = ptot-sum_max;
    if(pextra>0.0) {
      if (pextra<sum_tol) { /* rescale probs */
	if (noisy) {
	  warning("Reducing GENMUL sum (%g)\n",pextra);
	}
	for (i=0; i<ncat-1; i++) {  /* rescale probs */
	  *(p+i) /= (ptot/sum_max);
	}
      } else
	error("sum of P(i)>1 [too big to round]");
    } /* if sum_max < ptot < sum_max+sum_tol */
  } else {
    if (ptot > sum_max) {
      Rprintf("Sum of P(i)=%f (1-sum(P(i))=%g)\n",ptot,1.0-ptot);
      for (i=0; i<ncat; i++)
	Rprintf("%f ",p[i]);
      Rprintf("\n");
      error("sum of P(i) > 1");
    }
  }

/*  Initialize variables  */
    ntot = n;
    sum = 1.0;
    for(i=0; i<ncat; i++) ix[i] = 0;
    GetRNGstate();

/*    Generate the observation  */
    for(icat=0; icat<ncat-1; icat++) {
      prob = *(p+icat)/sum;
      *(ix+icat) = rbinom(ntot,prob);
      ntot -= *(ix+icat);
      if(ntot <= 0) return;
      sum -= *(p+icat);
    }
    *(ix+ncat-1) = ntot;

    PutRNGstate();
    return;
}


void rdirich(double *shape, int ncat, double *results) {
  double *gam;
  static int debug=-1;
  int i;
  double tot;
  static double minval=1.0e-200;

  if (ncat>MAXNCAT) {
    Rprintf("ERROR: number of categories (%d) > max (%ld) in rdirich\n",
	    ncat,MAXNCAT);
    /* exit(1); */  /* FIXME, should throw error */
  }
  for (i=0; i<ncat; i++)
    if (shape[i]<0) {
      error("shape parameter <0: %d, %g\n",i,shape[i]);
    }

  gam = (double *)Calloc(ncat,double);
  GetRNGstate();
  if (debug>0) Rprintf("RDIRICH:\n");
  for (i=0, tot=0.0; i<ncat; i++) {
    if (debug>0) Rprintf("%d %f\n",i,shape[i]);
    gam[i]=rgamma(shape[i],1.0);
    /* Rprintf("%g %g\n",gam[i],shape[i]); */
    tot += gam[i];
  }

  for (i=0; i<ncat; i++) {
    results[i] = gam[i]/tot;
    /* can't stop rgamma from returning zero for small shape parameters;
     * setting minimum on gam works less well */
    if (results[i]<minval)
      results[i]=minval;
  }
  PutRNGstate();
  Free(gam);
  return;
}  

void frdirich(float *shape, int ncat, float *results) {
  double *dshape, *dresults;
  static float fminval = 1e-30; /* ?? don't know exactly what
				 * the minimum floating-point number is,
				 * or how to make it portable;
				 * on RH 7.2/Intel it's approx. 1e-38 */
  int i;
  
  dshape = (double *)Calloc(ncat,double);
  dresults = (double *)Calloc(ncat,double);
  for (i=0; i<ncat; i++)
    dshape[i] = (double)shape[i];
  rdirich(dshape,ncat,dresults);
  for (i=0; i<ncat; i++) {
    results[i] = (float)dresults[i];
    if (results[i]<fminval)
      results[i]=fminval;
  }
  Free(dshape);
  Free(dresults);
  return;
}  

void fgenmul(int n, float *p, int ncat,int *ix, int round) {
  double *dp;
  int i;
  
  dp = (double *)Calloc(ncat,double);
  for (i=0; i<ncat; i++)
    dp[i] = (double)p[i];
  genmul(n,dp,ncat,ix,round);
  Free(dp);
  return;
}

void rmnom(int *size, double *prob, int *len, int*result) {
  genmul(*size, prob, *len, result, FALSE);
  return;
}

void rdirichwrap(double *shape, int *ncat, double *results) {
  rdirich(shape,*ncat,results);
  return;
}
