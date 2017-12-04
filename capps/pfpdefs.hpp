//
// Created by nicolet on 28/10/17.
//
#include <stdlib.h>
#include <climits>
#include <cstdio>

#ifndef PRECISE_PARALLEL_FP_PFPDEFS_H
#define PRECISE_PARALLEL_FP_PFPDEFS_H
#include "omp.h"

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

#endif //PRECISE_PARALLEL_FP_PFPDEFS_H
