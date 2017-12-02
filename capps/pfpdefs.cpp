//
// Created by nicolet on 31/10/17.
//
#include <stdio.h>
#include "pfpdefs.hpp"

void dclean(double *a, msize_t n){
    for (msize_t i = 0; i < n; ++i) {
        a[i] = 0.0;
    }
}

void propagate_error(double* a, double err, msize_t k, msize_t n) {
    for(msize_t i = k; i < n; i++){
        a[i] += err;
    }
}

void println_darray(double *a, msize_t from, msize_t to){
    if(from <= to - 1) {
        printf("[");
        for (msize_t i = from; i < to - 1; i++) {
            printf("%f, ", a[i]);
        }
        printf("%f]\n", a[to - 1]);
    } else {
        printf("[]\n");
    }
}