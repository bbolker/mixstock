#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <R_ext/Print.h>
#include "machdefs.h"
#include "utils.h"

#ifdef EXTRAMATH
/* round double to nearest integer; part of std Unix libraries */
int irint(double x) {
     double y; 

    y = (double)((int)x);
  if (fabs(frac(y))>0.5)
	y += bbsign(y);
	return(y);
  /* return((int)rint(x)); */
}

double log2(double x) {
   return(log(x)/M_LN2);
}

double aint(double x) {
   return(x<0 ? floor(x)-1.0 : floor(x));
}

#endif

void FatalError(char * msg, int ecode) {
    REprintf("%s\n",msg);
    /* exit(ecode); */
}

void WarnError(char * msg, int ecode) {
    REprintf("%d: %s\n",ecode,msg);
}

unsigned gettok( char * instr, unsigned maxlen, float * tmpcor )
{
  /* passed an arbitrary-length string, a maximum length, and a
   *  pointer to double, scans as many tokens as it can (up to
   *  maxlen) and fills in the array, returning the number it got
   */
  int lastcol = 0;
  int ok;
  char * tmptok;

  tmptok = strtok(instr," ,\t");
  if (tmptok != NULL)
	if (sscanf(tmptok,"%f",tmpcor))
	{
	  lastcol=1;
	  do
	    {
	      tmpcor++;
	      tmptok = strtok((char *)NULL," ,\t");
		  if (tmptok != NULL)
		ok = sscanf(tmptok,"%f",tmpcor);
	      else ok=0;
	      if (ok>0) lastcol++;
		}
		  while (lastcol<maxlen && ok);
	  return(lastcol);
	} /* if first one was OK */
    else return(0); /* got first token but it wasn't good */
  else return(0);  /* failed to find anything */
  
}
/*********************************************************/
void * salloc(size_t s)
{
  /* "safe alloc" */
  char * t=0;
  static float tot_alloc=0.0;

  if (s>0)
    {
      t = (char *)malloc(s);
      if (t==NULL) {
	  REprintf("Sorry: out of memory\n");
	  /* fflush(stderr); */
	  REprintf("(tried to allocate %ld bytes, total allocation so far=%1.0f)\n",s,tot_alloc);
	  /* fflush(stderr); */
	  /* exit(12); */
      }
      tot_alloc += (float)s;
    }
  return((void *)t);
}
/*********************************************************/
void * scalloc(size_t s, unsigned n)
{
  /* "safe calloc" */
  char * t=0;

  if ((s*n)>0) {
	  t = (char *)calloc(s,n);
	  if (t==NULL) { REprintf("Sorry: out of memory\n");
	      /* exit(12); */
		   };
	  };
	  return((void *)t);
}

/*********************************************************/
void skipcomm(FILE * infil, char * inlin, int len) {
  int commline;
  do {
    if (fgets(inlin,len,infil)==NULL) {
      inlin[0]=(char)0;
      commline=FALSE;
    } else
    commline=(inlin[strspn(inlin," \t")]=='#');
  } while (commline);
}
/*********************************************************/
FILE * tryfile(char * descrip, char * suff, int outp, int excode) {
  char fn[STRLEN];
  char emsg[STRLEN];
  char icode[2];
  char ioname[7];
  

  FILE * tmpfil;
  strcpy(fn,descrip);
  strcat(fn,suff);
  
  if (outp) {
    strcpy(icode,"w"); strcpy(ioname,"output"); } else {
      strcpy(icode,"r"); strcpy(ioname,"input"); };
   
    if ((tmpfil=fopen(fn,icode))==NULL) {
      sprintf(emsg,"Can't open %s file %s\n",ioname,fn);
      FatalError(emsg,excode);
      /* exit(0); */ /* make compiler happy */
      return(NULL); /* make Borland C compiler happy */
    } else
      return(tmpfil);

}
/*********************************************************/
