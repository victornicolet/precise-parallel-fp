/*
 *  Copyright (c) 2016 Inria and University Pierre and Marie Curie
 *  All rights reserved.
 */

/**
 *  \file cpu/blas1/ExSUM.hpp
 *  \brief Provides a set of summation routines
 *
 *  \authors
 *    Developers : \n
 *        Roman Iakymchuk  -- roman.iakymchuk@lip6.fr \n
 *        Sylvain Collange -- sylvain.collange@inria.fr \n
 */

#ifndef EXSUM_HPP_
#define EXSUM_HPP_

#include "superaccumulator.hpp"
#include "ExMTS_FPE.hpp"
#define TBB_PREVIEW_DETERMINISTIC_REDUCE 1
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>
#include <omp.h>
#include <algorithm>

#ifdef EXBLAS_MPI
#include <mpi.h>
#endif
#include "common.hpp"
#include "blas1.hpp"

static inline __mts join(__mts const l, __mts const r){
    return {r.sum + l.sum, std::max(l.mts, l.mts + r.sum)};
}


/**
 * \class TBBlongsum
 * \ingroup ExSUM
 * \brief This class is meant to be used in our multi-level reproducible and 
 *  accurate algorithm with superaccumulators only
 */
class TBBlongmts { // TODO figure out what can be done to do mts instead of sum
    double* a; /**< a real vector to sum */
public:
    Superaccumulator acc, mtsacc; /**< supperaccumulator */

    /**
     * The main function that performs summation of the vector's elelements into the 
     * superaccumulator
     */
    void operator()(tbb::blocked_range<size_t> const & r) {
        for(size_t i = r.begin(); i != r.end(); i += r.grainsize()) {
            acc.Accumulate(a[i]);
            mtsacc.Accumulate(a[i]);
        }
    }

    /** 
     * Construction that uses another object of TBBlongsum for initialization
     * \param x a TBBlongsum instance
     */
    TBBlongmts(TBBlongmts & x, tbb::split) : a(x.a), acc(e_bits, f_bits), mtsacc(e_bits, f_bits) {}

    /** 
     * Joins two superaccumulators of two different instances
     * \param y a TBBlongsum instance
     */
    void join(TBBlongmts & y) {
        double mtsl, mtsr;
        mtsr = y.mtsacc.Round();
        mtsacc.Accumulate(y.acc);
        acc.Accumulate(y.acc);
        mtsl = mtsacc.Round();
        if(mtsr > mtsl){
            mtsacc = y.mtsacc;
        }
    }

    /** 
     * Construction that initiates a real vector to sum and a supperacccumulator
     * \param a a real vector
     */
    TBBlongmts(double a[]) :
            a(a), acc(e_bits, f_bits), mtsacc(e_bits, f_bits)
    {}
};


/**
 * \ingroup ExSUM
 * \brief Parallel summation computes the sum of elements of a real vector with our 
 *     multi-level reproducible and accurate algorithm that solely relies upon superaccumulators
 *
 * \param N vector size
 * \param a vector
 * \return Contains the reproducible and accurate sum of elements of a real vector
 */
__mts ExMTSSuperacc(int N, double *a);

/**
 * \ingroup ExSUM
 * \brief Parallel summation computes the sum of elements of a real vector with our 
 *     multi-level reproducible and accurate algorithm that relies upon 
 *     floating-point expansions of size CACHE and superaccumulators when needed
 *
 * \param N vector size
 * \param a vector
 * \return Contains the reproducible and accurate sum of elements of a real vector
 */
template<typename CACHE> __mts ExMTSFPE(int N, double *a);

#endif // EXSUM_HPP_
