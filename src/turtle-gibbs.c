#include <stdio.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include <R.h>
#include <R_ext/PrtUtil.h>

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define MAXNCAT 500L
#define MAXHAP  MAXNCAT
#define MAXROOK 100L

typedef float  frvec[MAXROOK];
typedef float  fhvec[MAXHAP];
typedef int    irvec[MAXROOK];
typedef int    ihvec[MAXHAP];
typedef float  fhrmat[MAXHAP][MAXROOK];
typedef float  frhmat[MAXROOK][MAXHAP];
typedef int    ihrmat[MAXHAP][MAXROOK];
typedef int    irhmat[MAXROOK][MAXHAP];

/* 25 July 2000: revamp memory allocation etc. to avoid double-allocating:
 * single, etc. etc. */
extern void frdirich(float *prob, int ncat, float *results);
extern void fgenmul(int n, float *p, int ncat, int *ix, int round);

void gibbswrap(int *H, int *R, double *a, int *startiter, int *maxiter,
	       int * poolsamp, int * rooksamp, int *startfval, int *thin,
	       double *dprior, double *results,
	       int *outfile,
 	       char **outfn,
	       int *ranseed,
	       int *rptiter,
               int *contrun,
               double * startcontrib,
               double * startrookfreq);

int gibbs(int H, int R, float a, long startiter, long maxiter,
	  int * poolsamp, irhmat rooksamp, int startfval, int thin,
	  frvec fprior, double *resmat,
	  int outfile, char *outfn, int from_R, int rptiter,
          int contrun, frvec scontrib, frhmat srookfreq);

void gibbswrap(int *H, int *R, double *a, int *startiter, int *maxiter,
	       int * poolsamp, int * rooksamp, int *startfval, int *thin,
	       double *dprior, double *results,
	       int  *outfile,
 	       char **outfn,
	       int *ranseed,
	       int *rptiter,
               int *contrun,
               double * startcontrib,
               double * startrookfreq) {

  /* get everything in double format from R */

  frvec fprior;
  frvec scontrib;
  irhmat rookmat;
  frhmat srookfreq;
  int h,r;
  int npts,nvars;
  /* char *outfn="gibbsout"; */

  /* Rprintf("%d %d %d %f %f\n",*ranseed,*rptiter,*contrun,startcontrib[0],startrookfreq[0]); */

  /*  if (*outfile==1) {
   * Rprintf("Printing output to file %s.{BOT,FRQ}\n",outfn[0]);
   *   }
   */

    for (r=0; r<(*R); r++) {
      scontrib[r] = (float)startcontrib[r];
      for (h=0; h<(*H); h++)  {
        rookmat[r][h] = rooksamp[(*H)*r+h];
        srookfreq[r][h] = (float)startrookfreq[(*H)*r+h];
      }
    }

  /*   if (dprior[0]<0) */
  /*     fprior=NULL; */
  /*   else { */
  /*     Rprintf("allocating fprior\n"); */
  /*     fprior=farray(*R); */
    for (r=0; r<(*R); r++)
      fprior[r]=(float)dprior[r];
    /*   } */

  gibbs(*H,*R,(float)*a,(long)*startiter,(long)*maxiter,
	poolsamp,rookmat,*startfval,*thin,fprior,results,
	*outfile,*outfn,TRUE,*rptiter,*contrun,
	scontrib,srookfreq);
  npts =(*maxiter-*startiter)/(*thin);
  nvars = (*R)+(*H)*(*R);
  /*   for (i=0; i<npts; i++) */
  /*     for (j=0; j<nvars; j++) */
  /*       results[i+j*npts]=(double)(resmat[i][j]); */
      
  /* Rprintf("%f %f %f %f\n",
   * resmat[0][0],resmat[0][1],resmat[1][0],resmat[1][1]); */
}

int gibbs(int H, int R, float a, long startiter, long maxiter,
	  int * poolsamp, irhmat rooksamp, int startfval, int thin,
	  frvec fprior, double * results,
	  int outfile, char * outfn, int from_R, int rptiter,
          int contrun, frvec startcontrib, frhmat startrookfreq) {

  int r,h,w0,w1;
  float sum;
  float totpool;
  float harmn;
  long it;
  char tmpbuf[1000];
  FILE * botfile=NULL, * frqfile=NULL;
  int print_out=FALSE;
  int save_out=TRUE;
  int badflag;
  
  long npts;

  frhmat rookprior,  /* prior (Dirichlet) parameters for rookery haplotype freqs */
    rookfreqval,     /* current Gibbs sample of rookery haplotype freqs */
    rooktot2;        /* total Gibbs sample+prior+real sample rook-hap estimates */
  irvec rooktot;     /* total sample size for each rookery */
/*   frvec fval;         *//* current Gibbs sample for rookery contributions */
  float fval[MAXROOK]; /* test of funny display in gdb */
  fhvec ybar;        /* average overall haplotype freqs in rookery samples */
  fhvec poolfreq;    /* current computed Gibbs sample for pool haplotype freqs */
  fhrmat w;          /* temporary (transposed) Gibbs prob. that indiv with hap H comes from rookery R */
  ihrmat tmpmat;     /* Gibbs sample of numbers of each hap from each rookery */
  frvec rookcontrib, /* total (Gibbs) contribution of each rookery */
    rooktotp;        /* Gibbs contr + prior for rookery contribs */


  if (outfile>0) {
     print_out=TRUE;
     save_out=FALSE;
     /* open output files */
     strncpy(tmpbuf,outfn,100);
     strncat(tmpbuf,".bot",4);
     botfile = fopen(tmpbuf,"w"); 
     strncpy(tmpbuf,outfn,100);
     strncat(tmpbuf,".frq",4);
     frqfile = fopen(tmpbuf,"w"); 
  }
     
  if (H>MAXHAP) {
    Rprintf("# haplotypes (%d) exceeds maximum (%d)\n",H,MAXHAP);
    return(1);
  }
  if (R>MAXROOK) {
    Rprintf("# rookeries (%d) exceeds maximum (%d)\n",R,MAXROOK);
    return(1);
  }
  for (r=0; r<R; r++) {
    for (rooktot[r]=0, h=0; h<H; h++)
      rooktot[r] += rooksamp[r][h];
    if (rooktot[r]==0) {
      Rprintf("Can't do Gibbs with all-missing loci ...\n");
      return(2);
    }
  }
  for (h=0, totpool=0.0; h<H; h++)
    totpool += poolsamp[h];
  /*  calculate prior according to Pella and Masuda from harmonic mean:
      a parameter scales the strength of the prior */
  for (r=0, sum=0; r<R; r++) {
    sum += (float)1/rooktot[r];
  }
  harmn = 1/(sum/R);
  for (h=0; h<H; h++) {
    for (r=0,ybar[h]=0.0; r<R; r++)
      ybar[h] += (float)rooksamp[r][h]/rooktot[r];
    ybar[h] /= R;
  }
  for(h=0; h<H; h++)
    for (r=0; r<R; r++)
      rookprior[r][h] =  a*sqrt(harmn)*ybar[h];
/*   if (fprior==NULL) { */
/*     Rprintf("allocating fprior\n"); */
/*     fprior = farray(R); */
    /** default prior for contributions is EQUAL contrib from all rooks **/
  if (fprior[0]<0)
    for (r=0; r<R; r++)
      fprior[r]=(float)1/R;
  /* } */
  /* allocate results matrix if necessary */
  npts = (maxiter-startiter)/thin;
  if (!from_R && outfile<=0) {
    results = (double *)calloc(npts*(R+H*R),sizeof(double));
    Rprintf("Allocating space: from_R=%d\n",from_R);
  }
    /*     results = fmatrix((maxiter-startiter)/thin,R+H*R); */
    /*  dimnames(results) <- list(NULL,
	c(paste("contrib",dimnames(rooksamp)$rookery,sep="."),
	outer(1:H,1:R,function(x,y)paste("rookhap",dimnames(rooksamp)$rookery[y],
	dimnames(rooksamp)$haplotype[x],sep="."))))
    */
  if (contrun==1) { /* continuation run; set fval and rookfreq directly */
     Rprintf("continuation run\n");
     for (r=0; r<R; r++) {
       fval[r]=startcontrib[r];
       for (h=0; h<H; h++)
          rookfreqval[r][h]=startrookfreq[r][h];
    }
  } else {
    /* set initial rookery freqs */
    for (r=0; r<R; r++)
      frdirich(&(rookprior[r][0]),H,&(rookfreqval[r][0]));
    /*      Rprintf("\nrookprior:\n"); */
    /*      for (r=0; r<R; r++) { */
    /*        for (h=0; h<H; h++) */
    /*  	Rprintf("%f ",rookprior[r][h]); */
    /*        Rprintf("\n"); */
    /*      } */
    /*      Rprintf("\nINITIAL rookfreqval:\n"); */
    /*      for (r=0; r<R; r++) { */
    /*        for (h=0; h<H; h++) */
    /*  	Rprintf("%g ",rookfreqval[r][h]); */
    /*      Rprintf("\n"); */
    /* set initial contributions */
    if (startfval<0)  /* ## use random start */
      frdirich(&(fprior[0]),R,&(fval[0]));
    else if (startfval==0)  /* equal-contribution start */
      for (r=0; r<R; r++)
	fval[r]=(float)1/R;
    else if (startfval<=R) {  /* start with 95% in one rookery,
			       * the rest evenly divided */
      for (r=0; r<R; r++) {
	if (r==startfval)
	  fval[r]=0.95;
	else
	  fval[r]=0.05/(R-1);
      }
    }
    else {
      Rprintf("startfval must be between 0 and R\n");
      return(3);
    }
  } /* if not contrun */
    /* pool contribs (f): vector, length R
     * rook haplotype freqs (h): matrix,
     *      H rows (haplotypes) x R cols (rookeries)
     * pool freqs (pool.freq): h %*% f, vector, length R
     * "val" indicates realized (Gibbs-sampler) value as
     * opposed to Dirichlet params
     */
    for (it=0; it<maxiter; it++) {
      if (rptiter>0 && it % rptiter==0)
	Rprintf("it. %ld\n",it);
      for (h=0; h<H; h++) poolfreq[h]=0.0;
      for (h=0; h<H; h++)
	/* total expected frequency in mixed population  */
	for (r=0; r<R; r++)
	  poolfreq[h]+=rookfreqval[r][h]*fval[r];
      /* probability that an individual with hap H (row) comes
       * from rookery R (column);
       * use R's "columns first" rule to calculate */
      badflag=FALSE;
      for (h=0; h<H; h++)
	for (r=0; r<R; r++) {
	  if (rookfreqval[r][h]==0 && poolfreq[h]==0) {
	    Rprintf("flagged w[%d,%d]\n",r,h);
	    w[h][r] = 0.0; badflag=TRUE;
	  }
	  else
	    w[h][r] = rookfreqval[r][h]*fval[r]/poolfreq[h];
	}
      /*     ## take multinomial samples of each type ... */
      for (h=0; h<H; h++)
	/* 	genmul(poolsamp[h],&(w[h][0]),H,&(tmpmat[h][0])); */
	/*** DUMB DUMB DUMB DUMB! fixed bug -- was "H" instead of "R" for
	 *** number of categories  ... ***/
	fgenmul(poolsamp[h],&(w[h][0]),R,&(tmpmat[h][0]),TRUE);
	/*  get posteriors for p (pool contribs, f) and Q
	 * (rook freqs, rookfreq) 
	 * supposing we're doing things the easy way: 
	 * posterior of p =  (pool sample plus any priors if desired)
	 * rookcontrib <- apply(tmpmat,2,sum) */
      for (h=0; h<H; h++)
	for (r=0; r<R; r++) 
	  if (tmpmat[h][r]<0) {
	    Rprintf("flagged tmpmat[%d,%d]\n",h,r);
	    badflag=TRUE;
	  }
      if (badflag) {
	Rprintf("\npoolsamp: ");
	for (h=0; h<H; h++)
	      Rprintf("%d ",poolsamp[h]);
	Rprintf("\nw:\n");
	for (h=0; h<H; h++) {
	  for (r=0; r<R; r++)
	    Rprintf("%f ",w[h][r]);
	  Rprintf("\n");
	}
	Rprintf("\ntmpmat:\n");
	for (h=0; h<H; h++) {
	  for (r=0; r<R; r++)
	    Rprintf("%d ",tmpmat[h][r]);
	  Rprintf("\n");
	}
	Rprintf("\nrookfreqval:\n");
	for (r=0; r<R; r++) {
	  for (h=0; h<H; h++)
	    Rprintf("%f ",rookfreqval[r][h]);
	  Rprintf("\n");
	}
	Rprintf("\nfval:\n");
	for (r=0; r<R; r++)
	  Rprintf("%f ",fval[r]);
	Rprintf("\npoolfreq:\n");
	for (h=0; h<H; h++)
	  Rprintf("%f ",poolfreq[h]);
	Rprintf("\n");
      }
      for (r=0; r<R; r++) rookcontrib[r]=0.0;
      for (h=0; h<H; h++)
	for (r=0; r<R; r++)
	  rookcontrib[r] += tmpmat[h][r];
      for (r=0; r<R; r++) {
	rooktotp[r] = rookcontrib[r]+fprior[r];
	if (rooktotp[r]<0) Rprintf("rooktotp<0\n");
      }
      frdirich(&(rooktotp[0]),R,&(fval[0]));
      /* posterior of Q = (rookery sample + pool sample (known) + priors) */
      for (r=0; r<R; r++) {
	for (h=0; h<H; h++)
	  rooktot2[r][h] = rooksamp[r][h]+tmpmat[h][r]+rookprior[r][h];
	frdirich(&(rooktot2[r][0]),H,&(rookfreqval[r][0]));
      }
      if (it>=startiter) {
	w0 = it-startiter;
	if (w0 % thin == 0) {
	  w1 = w0 / thin;
          if (print_out) {
  	     for (r=0; r<R; r++)
	        fprintf(botfile,"%1.10g ",fval[r]);
             fprintf(botfile,"\n");
	      for (r=0; r<R; r++)
	        for (h=0; h<H; h++)
	           fprintf(frqfile,"%1.10g ",rookfreqval[r][h]);
             fprintf(frqfile,"\n");
          } /* print output */
          if (save_out) {
  	     for (r=0; r<R; r++)
	        results[w1+r*npts] = (double)fval[r];
	      for (r=0; r<R; r++)
	        for (h=0; h<H; h++)
	           results[w1+(R+r*H+h)*npts] = (double)rookfreqval[r][h];
	  } /* save output */
	}  /* if thin */
      }  /* if beyond burn-in */
    } /* main Gibbs iteration loop */
    /* if (outfile) {
       Rprintf("check BUGS output file -- not tested since changes\n");
       BUGSout(results,H,R,(maxiter-startiter)/thin,outfn);
       } */
    if (print_out) {
      fclose(botfile);
      fclose(frqfile);
    }
    return(0);
}


