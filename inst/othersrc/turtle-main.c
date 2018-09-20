int main(int argc, char ** argv) {
  int H, R;
  float a;
  long startiter, maxiter;
  int rptiter=1000;
  int poolsamp[MAXHAP];
  irhmat rooksamp;
  int startfval;
  int thin;
  frvec fprior,scontrib;
  frhmat srookfreq;
  double *resmat=0;
  int tmp; /* to prevent warnings on unused result from scanf */

  int h,r;

  /* no sanity checks: should check maxiter<startiter, (maxiter-startiter)/thin is
   * an integer, startfval <= R, etc., etc. */
  printf("Enter # haps, # rooks:\n");
  tmp=scanf("%d %d",&H,&R);
  /*   rooksamp=i2matrix(R,H); */
  printf("Enter hap samples from pooled pop (%d)\n",H);
  for (h=0; h<H; h++)
      tmp=scanf("%d",poolsamp+h);
  for (r=0; r<R; r++) {
    printf("Hap samples from rookery %d:\n",r+1);
    for (h=0; h<H; h++)
	tmp = scanf("%d",rooksamp[r]+h);
  }
  printf("startfval (0 for equal contribs, 1<=n<=R for biased contrib):\n");
  tmp = scanf("%d",&startfval);
  printf("burn-in, total, thinning factor:\n");
  tmp = scanf("%ld %ld %d",&startiter,&maxiter,&thin);
  a=1.0;
  gibbs(H,R,a,startiter,maxiter,poolsamp,rooksamp,startfval,thin,fprior,resmat,
	1,"turtle-gibbs",FALSE,rptiter,FALSE,scontrib,srookfreq);
  exit(0);
}

void BUGSout(double * g, int H, int R, long tot, char *fn);

void BUGSout(double * g, int H, int R, long tot, char *fn) {
  /*   ## take a gibbs run (nxp dataframe) and output to files in BUGS/CODA-compatible format */
  FILE *ind, *out;
  int p;
  int i;
  long j;

  ind=tryfile(fn,".ind",1,1);
  out=tryfile(fn,".out",1,2);

  p = R+H*R; /* number of columns */

  for (i=0; i<p; i++) {
    fprintf(ind,"var%d %ld %ld\n",i+1,i*tot+1,(i+1)*tot);
    for (j=0; j<tot; j++)
      fprintf(out,"%ld %f\n",j+1,g[j+i*tot]);
  }
  return;
}
