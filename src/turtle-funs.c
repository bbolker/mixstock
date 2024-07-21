#include <math.h>
#include <stdio.h>
#include <R_ext/PrtUtil.h>
#include <Rmath.h>
#include "utils.h"
#define MAXMARK 500
#define MAXMARK1 501
#define MAXSRC  100
#define MAXSRC1 101
/* MAXDIM was 100 */
#define MAXVECDIM 50000  /* product */

void rem(double * x, int * n) {
    /* compute 'fraction remaining' */
    /* argh! this should be easy ... fixed stupid stupid bugs */
  int i;
  double tmp1,tmp2;
  
  for (i=0, tmp1=1.0; i<*n; i++) {
    tmp2=x[i];
    x[i]=tmp1;
    tmp1 -= tmp2;
  }
}

  
void q_to_p(double * q, int n, int contin) {
  /* take a vector of (n-1) numbers in (-Inf,Inf) and produce a */
  /* vector of n numbers in (0,1) that sum to 1: arctangent of each
   * number is the fraction of remaining stuff that goes to that
   * category.
   * contin determines whether we have to transform back from (-Inf,Inf)
   * or just from (0,1) */
  double v[MAXVECDIM];
  double rem;
  int i;

  if (contin==1) {
    /*  fprintf(stderr,"q_to_p: transforming variables\n"); */
    for (i=0; i<(n-1); i++)
      v[i] = atan(q[i])/M_PI+0.5;
  } else 
    /* more efficient way to do this? */
    for (i=0; i<(n-1); i++)
      v[i] = q[i];
  rem = 1.0;
  for (i=0; i<(n-1); i++) {
    q[i] = v[i]*rem;
    rem -= q[i];
  }
  q[n-1] = rem;
}

void q_to_p2(double * q, int *n, int *contin) {
  /* take a vector of (n-1) numbers in (-Inf,Inf) and produce a */
  /* vector of n numbers in (0,1) that sum to 1: arctangent of each
   * number is the fraction of remaining stuff that goes to that
   * category */
  /* double v[MAXDIM]; */
  double rem;
  int i;

  /* Rprintf("q_to_p2: n=%d\n",*n); */
  if (*contin==1) {
    /*  fprintf(stderr,"q_to_p: transforming variables\n"); */
    for (i=0; i<(*n-1); i++)
      q[i] = atan(q[i])/M_PI+0.5;
  }
  rem = 1.0;
  for (i=0; i<(*n-1); i++) {
    /*    tmp = q[i]; */
    q[i] *= rem;
    rem -= q[i];
  }
  q[*n-1] = rem;
}

void dcmat(double *p, double *xmat, int *len, int *debug) {
    int r,c,k,nrow,ncol;
    int contin=0;
    /* how big must this be? how big can it be? */
    double remmat[MAXMARK][MAXMARK1];
    double tmpsum;
    nrow = *len+1;
    ncol = *len;
    /* Rprintf("dcmat\n"); */
    /* Rprintf("dcmat: len=%d\n",*len); */
    /* fill in remmat with x[-z] */
    for (r=0; r<(nrow-1); r++) {
	for (c=0, k=0; c<(ncol/*-1 */); c++, k++) {
	    if (k==r) k++; /* skip self */
	    remmat[r][c] = p[k];
	}
	/* translate remmat from transformed probs to probs */
	/* Rprintf("calling q_to_p2: %d\n",*len); */
	q_to_p2(remmat[r],len,&contin);
	/* translate to 'fraction remaining' */
	rem(remmat[r],len);
    }
    if (*debug==1) {
	for (r=0; r<(nrow-1); r++) {
	    for (c=0; c<ncol; c++)
		REprintf("%1.3f ",remmat[r][c]);
	    REprintf("\n");
	}
    }
    /* fill in lower triangle */
    for (r=1; r<(nrow-1); r++)
	for (c=0; c<r; c++) {
	    /*      xmat[r+c*(nrow)] = -p[r]*remmat[r-1][c]; */
	    /* have to transpose remmat here */
	    xmat[r+c*(nrow)] = -p[r]*remmat[c][r-1];
	}
    if (*debug==1) {
	for (r=0; r<(nrow-1); r++) {
	    for (c=0; c<ncol; c++)
		REprintf("%1.3f ",xmat[r+c*nrow]);
	    REprintf("\n");
	}
    }
    /* fill in diagonal */
    for (c=0; c<ncol; c++)
	remmat[0][c] = p[c];
    /* Rprintf("calling q_to_p2 2: %d\n",*len); */
    q_to_p2(remmat[0],len,&contin);
    rem(remmat[0],len);
    for (r=0; r<(nrow-1); r++)
	xmat[r+r*(nrow)] = remmat[0][r];
    /* fill in last row */
    for (c=0; c<ncol; c++) {
	for (tmpsum=0, r=0; r<nrow-1; r++)
	    tmpsum-=xmat[r+c*(nrow)];
	xmat[(nrow-1)+c*(nrow)] = tmpsum;
    }
}
  
void p_to_q(double * p, int n) {
  double rem[MAXSRC];
  double v[MAXSRC];
  int i;

  rem[0] = 1.0;
  for (i=1; i<n; i++)
    rem[i]=rem[i-1]-p[i-1];
  for (i=0; i<n; i++) {
    if (rem[i]==0.0)
      v[i]=0;
    else
      v[i]=p[i]/rem[i];
    /* branch cut problems */
    if (v[i]>1.0) v[i]=1.0;
    if (v[i]<0.0) v[i]=0.0;
  }
  for (i=0; i<n; i++)
    p[i]=tan((v[i]-0.5)*M_PI);
}
  
/* # log-likelihood of a multinomial mixture */
/* # p: parameter vector. */
/* #    First (R-1) elements are frequency contributions from different sources, expressed in */
/* #    transformed (-infty to infty) coordinates, vector of length (R-1) */
/* #    Last (R:(R*(H+1)-1)) elements are frequencies in the source pools, expressed in */
/* #    transformed coordinates, really an H*(R-1) matrix (by column) */
/* # R: number of sources (rookeries) */
/* # H: number of types (haplotypes) */
/* # n.samp: sampled proportions in sources (rookery haplotypes), vector (R*H) */
/* # s.samp: sampled proportions in pool (pelagic pop haplotypes), vector (H) */
/* # second try: need to constrain all freqs to 0-1 and sum(freq)=1 */
/* # new parameterization q[] for probability vector from 1..n: */
/* #  p[j] = (atan(q[j])*pi+0.5)*(1-sum(p[1:(j-1)])), p[0] defined as 1 */


void trans_par(double poolfreq[],
	       double h[MAXSRC][MAXMARK],
	       double *prmvec,
	       int R, int H,
	       int contin,
	       int transf) {
  double f[MAXVECDIM];
  int i,j,k;

  if (transf==1) {
    /*  transform contribution parameters to probabilities */
    for (i=0; i<(R-1); i++) f[i]=prmvec[i];
    q_to_p(f,R,contin);                      
    k=R-1;
  } else {
    for (i=0; i<R; i++) f[i]=prmvec[i];
    k=R;
  }
  for (i=0; i<R; i++) {
    if (transf==1) { /* transform rookery freq parameters to proportions */
      for (j=0; j<(H-1); j++)
	h[i][j]=prmvec[k++];
      q_to_p(h[i],H,contin);
    } else {
      for (j=0; j<H; j++)
	h[i][j]=prmvec[k++];
    }
  }
  for (i=0; i<H; i++) {
    poolfreq[i]=0.0;
    for (j=0; j<R; j++)
      poolfreq[i]+=h[j][i]*f[j];
  }
}

/* unused */
void invtrans_par(double poolfreq[], 
		  double h[MAXSRC][MAXMARK], 
		  double *prmvec, int R, int H,
		  int contin) {
  /* translate parameters from */
  double f[MAXSRC];
  int i,j,k;

  for (i=0; i<(R-1); i++)
     f[i]=prmvec[i];
  q_to_p(f,R,contin);                      /*  transform to probabilities */
  k=R-1;
  for (i=0; i<R; i++)
    for (j=0; j<(H-1); j++)
      h[i][j]=prmvec[k++];
  for (i=0; i<R; i++)
    q_to_p(h[i],H,contin);
  for (i=0; i<H; i++) {
    poolfreq[i]=0.0;
    for (j=0; j<R; j++)
      poolfreq[i]+=h[j][i]*f[j];
  }
}

/* trans.Rpar = function(p,n.samp,s.samp,h) { */
/*   H = length(s.samp) */
/*   R = length(n.samp)/H */
/*   f = q.to.p(p)  *//*  transform to probabilities */
/*   pool.freq = as.vector(h %*% f)   *//*  expected frequency in pool */
/*   list(pool.freq=pool.freq,h=h,R=R,H=H) */
/* } */

double multilik(double * prob, int * samp, int n, 
		int full, int debug) {
  int i;
  double lik;
  int tot;
  
  lik=0.0;
  for (i=0; i<n; i++) {
    if (!(prob[i]==0 && samp[i]==0)) {
      lik -= (double)samp[i]*log(prob[i]);
      if (debug==1)
	Rprintf("mlik: %d %f %d %f\n",i,lik,samp[i],log(prob[i]));
    }
  }
  if (full) {
    for (i=0,tot=0; i<n; i++) {
      lik += lgammafn((double)samp[i]+1.0);
      tot += samp[i];
    }
    lik -= lgammafn((double)tot+1.0);
  }
  if (debug==1)
    Rprintf("full mlik: %f\n",lik);
  return(lik);
}

/* obsolete??*/
/* calculate UML likelihood from untransformed contrib/haplotype
  parameters */
void loglik3wrap(double *lik,
		 double *p, 
		 int *R, 
		 int *H, 
		 int * s_samp, 
		 int * n_samp,
		 int * cumcount, 
		 int * debug) {
  double pool_freq[MAXMARK];
  double contrib[MAXSRC];
  /* double h1[MAXDIM][MAXDIM]; */
  double *hvec;

  int r,h,k;
  double m1,m2;

  for (k=0; k<*R; k++)  /* copy first R items in param vector */
    contrib[k] = p[k];
  k=*R;  /* redundant */
  for (h=0; h<*H; h++) pool_freq[h]=0;
  for (r=0; r<*R; r++) {
     for (h=0; h<*H; h++) {
        pool_freq[h] += contrib[r]*p[k];
        k++;
     }
  }
  hvec = p+*R;  /* set the rookery frequencies vector */
  m1= multilik(pool_freq,s_samp,*H,0,*debug);
  /* need to calculate m2 one rookery at a time to
   * get combinatorial coefficient right */
  m2=0;
  for (r=0; r<(*R); r++) {
    /* step through rookeries */
    m2+= multilik(hvec,n_samp,(*H),0,*debug);
    hvec+=(*H);
    n_samp+=(*H);
  }
  if (*debug==1) Rprintf("LL2W: %f %f\n",m1,m2);
  *lik=m1+m2;
}


void loglik2wrap(double *lik,  /* likelihood (return value) */
		 double *p,    /* parameter vector */
		 int *R, int *H,  /* numbers of rookeries and haplotypes */
		 int * s_samp, /* pooled sample */ 
		 int * n_samp, /* rookery sample */
		 int * cumcount,  /* total number of evaluations */
		 int * contin,    /* transform from -Inf/Inf */
		 int * transf,    /* transform? */
		 int * full,      /* return full likelihood? */
		 int * cond,      /* conditional likelihood? */
		 int * debug) {
  double pool_freq[MAXMARK];
  double h1[MAXSRC][MAXMARK];
  double hvec[MAXVECDIM];

  int r,h,k;
  double m1,m2;

   *cumcount++;
   if (*debug==1) Rprintf("loglik2wrap: contin=%d, transf=%d\n",*contin,*transf);
   trans_par(pool_freq,h1,p,*R,*H,*contin,*transf); /* transform/compute pooled frequencies */
   k=0;
   /* flatten rookery haplotype parameters into a vector */
   for (r=0; r<*R; r++) {
     for (h=0; h<*H; h++) {
       hvec[k++]=h1[r][h];
     }
   }
   if (*debug==1) {
     for (k=0; k<(*R)*(*H); k++)
       Rprintf("LL2W: %f %d\n",hvec[k],n_samp[k]);
     for (k=0; k<(*H); k++)
       Rprintf("LL2W: %f %d\n",pool_freq[k],s_samp[k]);
   }
   m1= multilik(pool_freq,s_samp,*H,*full,*debug);
   m2 = 0;
   if (*cond!=1) m2= multilik(hvec,n_samp,(*R)*(*H),*full,*debug);
   if (*debug==1) Rprintf("LL2W: %f %f\n",m1,m2);
   *lik=m1+m2;
}


