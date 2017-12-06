//
// Created by nicolet on 28/10/17.
//

#ifndef PRECISE_PARALLEL_FP_PAR_PRECISE_FP_H
#define PRECISE_PARALLEL_FP_PAR_PRECISE_FP_H
#include "pfpdefs.hpp"

double dsum(double*, msize_t);
double dsum_par(double*, msize_t);
void dpxsum(double*, double*, msize_t);
void dpxsum_par(double*, double*, msize_t);
void dprecise_parallel(double* input,
                       double* prefix,
                       double* error,
                       double* errorprefix,
                       double precision,
                       msize_t size,
                       int *counter);

#define f64factor (1 << 21) + 1

template<typename T> inline T split(T a, T &s) {
    T c, x;
    c = f64factor * a;
    x = c - (c - a);
    s = a -x;
    return x;
}

template<typename T> inline T Dekker2Mul(T a,T b, T &s){
    T x, a1, a2, b1, b2;
    x = a * b;
    a1 = split(a, a2); // a = a1 + a2 and a1, a2 non overlapping
    b1 = split(b, b2); // a = a1 + a2 and a1, a2 non overlapping
    s = a2 * b2 - (((x - a1 * b1) - a2 * b1) - a1 * b2);
    return x;
}

#if INSTRSET > 7
template<typename T> inline T FMA2Mul(T a, T b, T &s){
    T x;
    x = a * b;
    s = fma(a, b, -x);
    return x;
}
#endif

template<typename T> inline T twomul(T a, T b, T &s) {
#if INSTRSET > 7
    FMA2Mul(a, b, s);
#else
    Dekker2Mul(a, b, s);
#endif
}

#endif //PRECISE_PARALLEL_FP_PAR_PRECISE_FP_H
