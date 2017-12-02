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

#endif //PRECISE_PARALLEL_FP_PAR_PRECISE_FP_H
