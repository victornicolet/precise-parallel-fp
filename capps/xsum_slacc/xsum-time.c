/* TIME MEASUREMENT PROGRAM FOR FUNCTIONS FOR EXACT SUMMATION. */

/* Written by Radford M. Neal, 2015 */

/* Program to perform timing tests of exact and non-exact summation methods.

   Run with a command of the form:

       xsum-time N M R [ data ... ] [ inf ]

   Here, N is the size of the vectors that are summed, M is the number
   of such vectors, and R is the number of times these M vectors of size
   N are summed, by each method.

   After running these tests for summing vector elements, the sum of the
   squares is also computed, and then the dot product of each vector
   with a corresponding vector in a different set of M vectors.

   The data in each of the vectors is as specified after R, with any
   remaining elements filled in randomly, differently for each of the
   M vectors.  This random data sums to zero, since the second half
   of it is the mirror of the first half, negated (with zero in the 
   middle if the number of elements is odd).  In addition, if "inf"
   is the last argument, every eigth element is set to infinity (so
   the result will also be infinity).

   This program can be compiled to compare with Zhu & Hayes' iFastSum
   and ExactSum methods, if the source for those is present.  The -DZHU
   option must be given to the compiler for this to happend.

   This program can be compiled to compare with 128-bit floating point
   arithmetic using the gcc __float128 type.  The -DFLOAT128 option
   must be given to the compiler for this to happen.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "xsum.h"

#ifdef ZHU
#include "ExactSum.h"
#endif

#define START_CLOCK ( start_clock = clock() )
#define END_CLOCK ( clock_dur = (double)(clock()-start_clock)/CLOCKS_PER_SEC )

int different (double a, double b)
{ 
  return isnan(a) != isnan(b) || !isnan(a) && !isnan(b) && a != b;
}

int main (int argc, char **argv)
{
  clock_t start_clock;
  double clock_dur;
  xsum_flt result_s, result_l, result_d, result_dno, result_k;
# ifdef ZHU
  xsum_flt result_iFast, result_Online;
  ExactSum mysum;
  xsum_flt *ai;
# endif
# ifdef FLOAT128
  xsum_flt result_128;
# endif
  int used_small, used_large;
  xsum_small_accumulator sacc;
  xsum_large_accumulator lacc;
  xsum_flt *a, *a2;
  int N, M, R;
  int i, j, k;
  int ndata, n;
  char **data;
  int inf;

  inf = 0;
  if (argc>0 && strcmp(argv[argc-1],"inf")==0)
  { inf = 1;
    argc -= 1;
  }
  ndata = argc-4;

  if (argc<4 || (N=atoi(argv[1]))<1 || ndata>N
             || (M=atoi(argv[2]))<1 || (R=atoi(argv[3]))<1)
  { fprintf(stderr,"Usage: xsum-time N M R [ data ... ] [ inf ]\n");
    exit(1);
  }

  data = argv+4;

  /* On an Intel machine, set the 387 FPU to do double rounding, in order
     to get correct IEEE 64-bit floating point.  (Only relevant if SSE2
     instructions not used instead of FPU.)  This will disable higher
     precision for long double, however! */

# ifdef DOUBLE
  { unsigned int mode = 0x27f;
    __asm__ ("fldcw %0" : : "m" (*&mode));
  }
# endif

  printf("\nTIMING TESTS\n\n");
  printf("N = %d, M = %d, R = %d",N,M,R);
  if (ndata > 0 || inf)
  { printf("  data:");
    for (j = 0; j < ndata; j++) printf(" %s",data[j]);
    if (inf) printf(" inf"); 
  }
  printf("\n");
 
  a = (xsum_flt *) calloc (2*N*M, sizeof *a);  /* Put a and a2 into one block */
  a2 = a + N*M;        /* to suppress possible variation in cache performance */

# ifdef ZHU
  ai = (xsum_flt *) calloc (N, sizeof *ai);
# endif

  n = N - ndata;
  for (i = 0; i < M; i++)
  { for (j = 0; j < ndata; j++) a[i*N+j] = atof(data[j]);
    if (n > 0)
    { int64_t rnd;
      a[ndata+i*N+n/2] = 0;
      rnd = 2345678;
      for (j = 0; j < n/2; j++)
      { rnd = (rnd*8192) % 67101323;/* from "Multiplicative congruential pseudo-
                                   random number generators", Downham&Roberts */
        a[ndata+i*N+j] = exp (30 * (double)rnd / 67101323);
        rnd = (rnd*8192) % 67101323;
        a[ndata+i*N+j] *= (double)rnd / 67101323;
        a[ndata+i*N+n-1-j] = -a[ndata+i*N+j];
        if (inf && j%8 == 0) 
        { a[ndata+i*N+n-1-j] = a[ndata+i*N+j] = 1.0/0.0;
        }
      }
    }
  }

  for (i = 0; i < N*M; i++) a2[i] = a[i];

  printf("\nVECTOR SUM\n\n");

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { xsum_small_init(&sacc);
      xsum_small_addv(&sacc,a+k*N,N);
      result_s = xsum_small_round(&sacc);
    }
  }
  END_CLOCK;
  printf("Small accumulator:  %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_s, clock_dur, clock_dur*1e9/R/N/M);

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { xsum_large_init(&lacc);
      xsum_large_addv(&lacc,a+k*N,N);
      result_l = xsum_large_round(&lacc);
    }
  }
  END_CLOCK;
  printf("Large accumulator:  %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_l, clock_dur, clock_dur*1e9/R/N/M);

# ifdef ZHU

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { memcpy (ai, a+k*N, N * sizeof *ai);
      result_iFast = mysum.iFastSum(ai-1,N);
    }
  }
  END_CLOCK;
  printf("iFastSum result:    %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_iFast, clock_dur, clock_dur*1e9/R/N/M);

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { result_Online = mysum.OnlineExactSum(a+k*N-1,N);
    }
  }
  END_CLOCK;
  printf("OnlineExact result: %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_Online, clock_dur, clock_dur*1e9/R/N/M);

# endif

# ifdef FLOAT128

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { result_128 = xsum_sum_float128(a+k*N,N);
    }
  }
  END_CLOCK;
  printf("Float 128 result:   %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_128, clock_dur, clock_dur*1e9/R/N/M);

# endif

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { result_k = xsum_sum_kahan(a+k*N,N);
    }
  }
  END_CLOCK;
  printf("Kahan sum result:   %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_k, clock_dur, clock_dur*1e9/R/N/M);

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { result_d = xsum_sum_double(a+k*N,N);
    }
  }
  END_CLOCK;
  printf("Double result:      %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_d, clock_dur, clock_dur*1e9/R/N/M);

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { result_dno = xsum_sum_double_not_ordered(a+k*N,N);
    }
  }
  END_CLOCK;
  printf("Double not ordered: %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_dno, clock_dur, clock_dur*1e9/R/N/M);

  if (ndata <= 2)
  { double correct;
    int i;
    if (inf && n/2 > 0)
    { correct = 1.0/0.0;
    }
    else
    { correct = 0.0;
      for (i = 0; i < ndata; i++) correct += a[i];
    }
    if (different(result_s,correct))
    { printf("RESULT USING SMALL ACCUMULATOR IS WRONG\n");
    }
    if (different(result_l,correct))
    { printf("RESULT USING LARGE ACCUMULATOR IS WRONG\n");
    }
  }

  if (result_s != (double) xsum_small_round(&sacc))
  { printf ("RESULT DIFFERS AFTER ROUNDING SMALL ACCUMULATOR TWICE\n");
  }

  if (result_l != (double) xsum_large_round(&lacc))
  { printf ("RESULT DIFFERS AFTER ROUNDING LARGE ACCUMULATOR TWICE\n");
  }

  if (result_s != result_l)
  { printf("RESULTS DIFFER FOR SMALL AND LARGE ACCUMULATORS\n");
  }

# ifdef ZHU

  if (result_iFast != result_l)
  { printf("RESULTS DIFFER FOR IFASTSUM AND LARGE ACCUMULATOR\n");
  }

  if (result_Online != result_l)
  { printf("RESULTS DIFFER FOR ONLINEEXACT AND LARGE ACCUMULATOR\n");
  }

# endif

  if (N <= 2 && result_s != result_d)
  { printf("RESULTS DIFFER FOR DOUBLE AND SMALL ACCUMULATOR\n");
  }

  if (N <= 2 && result_l != result_d)
  { printf("RESULTS DIFFER FOR DOUBLE AND LARGE ACCUMULATOR\n");
  }

  xsum_small_init(&sacc); /* Get used counts half-way, before cancels to zero */
  xsum_small_addv(&sacc,a,N/2);
  used_small = xsum_small_chunks_used(&sacc);

  used_large = xsum_large_chunks_used(&lacc);

  printf("\n");
  printf("Small accumulator chunks used: %d\n", used_small);
  printf("Large accumulator chunks used: %d\n", used_large);

  printf("\nVECTOR NORM\n\n");

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { xsum_small_init(&sacc);
      xsum_small_add_sqnorm(&sacc,a+k*N,N);
      result_s = xsum_small_round(&sacc);
    }
  }
  END_CLOCK;
  printf("Small accumulator:  %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_s, clock_dur, clock_dur*1e9/R/N/M);

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { xsum_large_init(&lacc);
      xsum_large_add_sqnorm(&lacc,a+k*N,N);
      result_l = xsum_large_round(&lacc);
    }
  }
  END_CLOCK;
  printf("Large accumulator:  %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_l, clock_dur, clock_dur*1e9/R/N/M);

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { result_d = xsum_sqnorm_double(a+k*N,N);
    }
  }
  END_CLOCK;
  printf("Double result:      %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_d, clock_dur, clock_dur*1e9/R/N/M);

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { result_dno = xsum_sqnorm_double_not_ordered(a+k*N,N);
    }
  }
  END_CLOCK;
  printf("Double not ordered: %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_dno, clock_dur, clock_dur*1e9/R/N/M);

  if (result_s != (double) xsum_small_round(&sacc))
  { printf ("RESULT DIFFERS AFTER ROUNDING SMALL ACCUMULATOR TWICE\n");
  }

  if (result_l != (double) xsum_large_round(&lacc))
  { printf ("RESULT DIFFERS AFTER ROUNDING LARGE ACCUMULATOR TWICE\n");
  }

  if (result_s != result_l)
  { printf("RESULTS DIFFER FOR SMALL AND LARGE ACCUMULATORS\n");
  }

  if (N <= 2 && result_s != result_d)
  { printf("RESULTS DIFFER FOR DOUBLE AND SMALL ACCUMULATOR\n");
  }

  if (N <= 2 && result_l != result_d)
  { printf("RESULTS DIFFER FOR DOUBLE AND LARGE ACCUMULATOR\n");
  }

  xsum_small_init(&sacc); /* Get used counts half-way, before cancels to zero */
  xsum_small_addv(&sacc,a,N/2);
  used_small = xsum_small_chunks_used(&sacc);

  used_large = xsum_large_chunks_used(&lacc);

  printf("\n");
  printf("Small accumulator chunks used: %d\n", used_small);
  printf("Large accumulator chunks used: %d\n", used_large);

  printf("\nVECTOR DOT PRODUCT\n\n");

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { xsum_small_init(&sacc);
      xsum_small_add_dot(&sacc,a+k*N,a2+k*N,N);
      result_s = xsum_small_round(&sacc);
    }
  }
  END_CLOCK;
  printf("Small accumulator:  %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_s, clock_dur, clock_dur*1e9/R/N/M);

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { xsum_large_init(&lacc);
      xsum_large_add_dot(&lacc,a+k*N,a2+k*N,N);
      result_l = xsum_large_round(&lacc);
    }
  }
  END_CLOCK;
  printf("Large accumulator:  %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_l, clock_dur, clock_dur*1e9/R/N/M);

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { result_d = xsum_dot_double(a+k*N,a2+k*N,N);
    }
  }
  END_CLOCK;
  printf("Double result:      %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_d, clock_dur, clock_dur*1e9/R/N/M);

  START_CLOCK;
  for (j = 0; j < R; j++)
  { for (k = 0; k < M; k++)
    { result_dno = xsum_dot_double_not_ordered(a+k*N,a2+k*N,N);
    }
  }
  END_CLOCK;
  printf("Double not ordered: %+.16le  time: %7.3lf s, %7.2lf ns/term\n",
         result_dno, clock_dur, clock_dur*1e9/R/N/M);

  if (result_s != (double) xsum_small_round(&sacc))
  { printf ("RESULT DIFFERS AFTER ROUNDING SMALL ACCUMULATOR TWICE\n");
  }

  if (result_l != (double) xsum_large_round(&lacc))
  { printf ("RESULT DIFFERS AFTER ROUNDING LARGE ACCUMULATOR TWICE\n");
  }

  if (result_s != result_l)
  { printf("RESULTS DIFFER FOR SMALL AND LARGE ACCUMULATORS\n");
  }

  if (N <= 2 && result_s != result_d)
  { printf("RESULTS DIFFER FOR DOUBLE AND SMALL ACCUMULATOR\n");
  }

  if (N <= 2 && result_l != result_d)
  { printf("RESULTS DIFFER FOR DOUBLE AND LARGE ACCUMULATOR\n");
  }

  xsum_small_init(&sacc); /* Get used counts half-way, before cancels to zero */
  xsum_small_addv(&sacc,a,N/2);
  used_small = xsum_small_chunks_used(&sacc);

  used_large = xsum_large_chunks_used(&lacc);

  printf("\n");
  printf("Small accumulator chunks used: %d\n", used_small);
  printf("Large accumulator chunks used: %d\n", used_large);

  printf("\n");

  return 0;
}
