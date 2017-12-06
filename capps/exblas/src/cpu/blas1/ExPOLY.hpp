/*
 *  Copyright (c) 2016 Inria and University Pierre and Marie Curie
 *  All rights reserved.
 */

/**
 *  \file cpu/blas1/ExPOLY.hpp
 *  \brief Provides a set of summation routines
 *
 *  \authors
 *    Developers : \n
 *        Roman Iakymchuk  -- roman.iakymchuk@lip6.fr \n
 *        Sylvain Collange -- sylvain.collange@inria.fr \n
 */

#ifndef EXPOLY_HPP_
#define EXPOLY_HPP_

#include "superaccumulator.hpp"
#include "ExPOLY_FPE.hpp"
#define TBB_PREVIEW_DETERMINISTIC_REDUCE 1
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>
#include <omp.h>
#include <blas1.hpp>

#ifdef EXBLAS_MPI
    #include <mpi.h>
#endif
#include "common.hpp"


/**
 * \class TBBlongsum
 * \ingroup ExPOLY
 * \brief This class is meant to be used in our multi-level reproducible and 
 *  accurate algorithm with superaccumulators only
 */
class TBBlongPoly {
    double* a; /**< a real vector to sum */
public:
    double factor;
    Superaccumulator acc; /**< supperaccumulator */

    /**
     * The main function that performs summation of the vector's elelements into the 
     * superaccumulator
     */
    void operator()(tbb::blocked_range<size_t> const & r) {
        double f, tmp;
        f = factor;
        for(size_t i = r.begin(); i != r.end(); i ++) {
            f *= factor;
            tmp = a[i] * factor;
            if (tmp > 0.) acc.Accumulate(tmp);
        }
        factor = f;
    }

    /** 
     * Construction that uses another object of TBBlongsum for initialization
     * \param x a TBBlongsum instance
     */
    TBBlongPoly(TBBlongPoly & x, tbb::split) : factor(x.factor), a(x.a), acc(e_bits, f_bits) {}

    /** 
     * Joins two superaccumulators of two different instances
     * \param y a TBBlongsum instance
     */
    void join(TBBlongPoly & y) {
        double rres = y.acc.Round() * factor;
        acc.Accumulate(rres);
        factor *= y.factor;
    }

    /** 
     * Construction that initiates a real vector to sum and a supperacccumulator
     * \param a a real vector
     * \param initfactor argument of the polynomial
     */
    TBBlongPoly(double a[], double initfactor) :
        factor(initfactor), a(a), acc(e_bits, f_bits)
    {}
};


/**
 * \ingroup ExPOLY
 * \brief Parallel summation computes the sum of elements of a real vector with our 
 *     multi-level reproducible and accurate algorithm that solely relies upon superaccumulators
 *
 * \param N vector size
 * \param a vector
 * \param inca specifies the increment for the elements of a
 * \param offset specifies position in the vector to start with 
 * TODO: not done for offset
 * \return Contains the reproducible and accurate sum of elements of a real vector
 */
__mts ExPOLYSuperacc(int N, double *a, double factor, int inca, int offset);

/**
 * \ingroup ExPOLY
 * \brief Parallel summation computes the sum of elements of a real vector with our 
 *     multi-level reproducible and accurate algorithm that relies upon 
 *     floating-point expansions of size CACHE and superaccumulators when needed
 *
 * \param N vector size
 * \param a vector
 * \param inca specifies the increment for the elements of a
 * \param offset specifies position in the vector to start with 
 * TODO: not done for inca and offset
 * \return Contains the reproducible and accurate sum of elements of a real vector
 */
template<typename CACHE> __mts ExPOLYFPE(int N, double *a, double factor, int inca, int offset);

#endif // EXSUM_HPP_
