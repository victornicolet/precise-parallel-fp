//
// Created by nicolet on 05/12/17.
//

#ifndef PRECISE_PARALLEL_FP_TEST_MTS_H
#define PRECISE_PARALLEL_FP_TEST_MTS_H

#include "common.hpp"
#include "pfpdefs.hpp"
#include "blas1.hpp"

/* This function times and compares the different implementations of mts.
 * It takes an optional argument (passed from commad line to main), which defines the size N of the array to be generated, by doing N = 2^arg (default is N = 2^20) */ 
void m_test_mts(int argc, char** argv);

#endif //PRECISE_PARALLEL_FP_TEST_MTS_H
