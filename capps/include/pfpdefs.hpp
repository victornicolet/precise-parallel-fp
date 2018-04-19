//
// Created by nicolet on 28/10/17.
//
#include <stdlib.h>
#include <climits>
#include <cstdio>

#ifndef PRECISE_PARALLEL_FP_PFPDEFS_H
#define PRECISE_PARALLEL_FP_PFPDEFS_H

// Using haswell tsc for timing
#ifndef PFP_TIME
#define PFP_TIME(call,start,mem)\
start = pfp_rdtsc();\
call;\
mem = (double) (pfp_rdtsc() - start);
#endif

#ifndef PPF_WTIME
#define PFP_WTIME(call,start,mem,start2,mem2)\
start2 = omp_get_wtime();\
start = pfp_rdtsc();\
call;\
mem = (double) (pfp_rdtsc() - start);\
mem2 = (double) (omp_get_wtime() - start2);
#endif


#include "omp.h"
#include <cstdint>

typedef long msize_t;
static double rand_dbl_limit = 10.0;

static inline void frand_init(double* a, msize_t n){
    for(msize_t i = 0; i < n; i++){
        a[i] = (double)rand()/(RAND_MAX/rand_dbl_limit);
    }
}

static inline double dabs(double a){
    return (a < 0) ? -a : a;
}

struct _mts_ {double sum; double mts;};

void propagate_error(double* input, double error, msize_t begin, msize_t asize);
void println_darray(double *a, msize_t from, msize_t to);
double myparallelmax(double* data, int ndata);
_mts_ myparallelmts(double* data, int ndata);
_mts_ custom_reduce_mts(double* data, int ndata);

inline uint64_t pfp_rdtsc()
{
    uint32_t hi, lo;
    asm volatile ("rdtsc" : "=a"(lo), "=d"(hi));
    return lo | ((uint64_t)hi << 32);
}

#endif //PRECISE_PARALLEL_FP_PFPDEFS_H
