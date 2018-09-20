/* generically useful utility definitions etc. */
#include <stdlib.h>  /* for size_t def. */
#  ifdef PROTOS_MISSING
   int printf(char * format, ...);
   int fprintf( FILE * stream, char * format, ...);
   int scanf(char * format, ...);
   int fscanf(FILE * stream,char * format, ...);
   int sscanf(char * s,char * format, ...);
   int fclose(FILE * stream);
   int fflush(FILE * stream);
   int unlink(char * name);
#  endif
#  ifndef sun
#     define M_PI 	   3.14159265358979323846
#     define M_SQRT2   1.41421356237309504880
#     ifndef M_LN2
#        define M_LN2     0.693147181 /* more? */
#     endif
#  endif
#  ifdef EXTRAMATH
/* these are defined in the SunOS 4 math libraries */
      double aint(double);
      int irint(double);
       double log2(double);
#     ifndef M_LN2
#         define M_LN2     0.693147181 /* more? */
#     endif
#     define M_PI 	   3.14159265358979323846
#  endif
#  ifdef SOLARIS2
#      define index(x,y) strchr(x,y)
#  endif

#ifndef SQR
#   define SQR(x) ( (x) * (x) )
#endif

#define	TWOPI	   6.28318530717958647692
#define SQRT_TWOPI 2.50662827463100050242
#define PISQ       9.8696044010893586191

#ifndef TINY
#  define TINY 1e-8 
#endif
#ifndef BIGNUM
#  define BIGNUM    1.0e20
#endif
#ifndef TRUE
#   define TRUE 1
#endif
#ifndef FALSE
#   define FALSE (!TRUE)
#endif
#define STRLEN 256    /* a reasonable "default" string length */
#define BIGSTR 65535
#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
/* avoid conflicts, particularly with John Saponara's cview.
 *  Pray that anyone who has defined MIN and MAX before me
 *    has done it sensibly!!!
 */
#ifndef MIN
#  define MIN(x,y)  ((x)<(y) ? (x) : (y))
#endif
#ifndef MAX
#  define MAX(x,y)  ((x)>(y) ? (x) : (y))
#endif
#define bbsign(x) ((x)/fabs(x))
/* consider a macro: this only works for float/double! */
#define frac(x) ((x)-(int)(x))

typedef short bool;
unsigned gettok( char * instr, unsigned maxlen, float * tmpcor );
void FatalError(char * msg, int ecode);
void WarnError(char * msg, int ecode);
void skipcomm(FILE * infil, char * inlin, int len);
FILE * tryfile(char * descrip, char * suff, int outp, int excode);

/**************************************************
*  safe procedures for allocating memory:
*   copied from Kay & Kummerfeld, _C Programming
*   in a UNIX Environment_ */
void * salloc(size_t s);
void * scalloc(size_t s, unsigned n);
#define talloc(type)          (type *)salloc(sizeof(type))
#define tsalloc(type,size)    (type *)salloc((unsigned)(sizeof(type)*size))

typedef struct param {
  char descrip[256];   /* description */
  char name[4];      /* short name */
  int type;        /* 0=scalar, 1=vector, 2=matrix */
  float *value;    /* pointer to value */
  float defval;      /* default value */
} param;

double complex_modarg_to_reim(double *x, double *y);
double complex_mod(double x, double y);
double complex_arg(double x, double y);
