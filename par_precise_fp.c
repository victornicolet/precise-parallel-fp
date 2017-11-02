//
// Created by nicolet on 28/10/17.
//
#include "pfpdefs.h"
#include "par_precise_fp.h"
#include "stdbool.h"
#include "omp.h"

double dsum(double *a, msize_t n){
    double sum = 0;
    for(msize_t i = 0; i < n; i++) {
        sum += a[i];
    }
    return sum;
}

double dsum_par(double *a, msize_t n){

    double sum = 0;
    msize_t i = 0;
#pragma omp parallel for private(i) schedule(static) reduction(+:sum)
    for (i = 0; i < n; i++) {
        sum = sum + a[i];
    }

    return sum;
}

void dpxsum(double *a, double* prefixes, msize_t size){
    double psum = 0;
    for (msize_t i = 0; i < size; ++i) {
        psum += a[i];
        prefixes[i] = psum;
    }
}

// Compute a prefix sum in prefixes
void dpxsum_par(double* a, double* prefixes, msize_t size){
    {
        double *z, *x = prefixes;
        int nthr;

#pragma omp parallel
        {
            msize_t i;
#pragma omp single
            {
                nthr = omp_get_num_threads();
                z = malloc((size_t) (sizeof(double) * (nthr + 1)));
                z[0] = 0;
            }
            int tid = omp_get_thread_num();

            double sum = 0;
#pragma omp for schedule(static)
            for(i=0; i<size; i++) {
                sum += a[i];
                x[i] = sum;
            }
            z[tid+1] = sum;
#pragma omp barrier

            int offset = 0;
            for(i=0; i<(tid+1); i++) {
                offset += z[i];
            }

#pragma omp for schedule(static)
            for(i=0; i<size; i++) {
                x[i] += offset;
            }
        }
        free(z);
    }
}

void dprecise_parallel(double* input, double* prefix, double* err, double* errprefix,
                       double precision, msize_t n, int *cnt) {
    _Bool terminate = true;
    (*cnt)++;
    // Step 2
    double secSum;
    msize_t i = 0;

    err[0] = 0;

#pragma omp parallel for shared(terminate) schedule(static) private(i)
    for(i=1; i < n; i++) {
        secSum = input[i] + prefix[i - 1];
        err[i] = prefix[i] - secSum;
        if (dabs(err[i]) > precision)
            terminate = false;
    }
//    err'[i] = prefix[i] - prefix[i-1] - input[i]

    if (!terminate) {
        dpxsum(err, errprefix, n);

#pragma omp parallel for simd
        for (i=1; i<n; i++) {
            prefix[i] = prefix[i] - errprefix[i];
        }
//        prefix'[i] = prefix[i] - SUM(0,i)[err]
//                   = prefix[i] - SUM(j = 1, .. i)[prefix(j) - input(j) - prefix(j-1)]
//                   = prefix[i] - (prefix[i] - prefix[0]) - SUM(j = 1, .. i)[input(j)]
//                   = prefix[i] - (prefix[i] - SUM(j = 0, .. i)[input(j)])
        dprecise_parallel(input, prefix, err, errprefix, precision, n, cnt);
    } else {
        return;
    }
}