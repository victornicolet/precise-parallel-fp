//
// Created by nicolet on 05/12/17.
//

#ifndef PRECISE_PARALLEL_FP_TEST_POLY_H
#define PRECISE_PARALLEL_FP_TEST_POLY_H

#include "common.hpp"
#include "pfpdefs.hpp"
#include "blas1.hpp"


void m_test_poly(int argc, char** argv);
__mts sequential_poly(int, double*, double);
__mts inexact_parallel_poly(int, double*, double);

#endif //PRECISE_PARALLEL_FP_TEST_MTS_H
