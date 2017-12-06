/*
 *  Copyright (c) 2016 Inria and University Pierre and Marie Curie
 *  All rights reserved.
 */

/**
 *  \file cpu/blas1/ExSUM.FPE.hpp
 *  \brief Provides a set of routines concerning floating-point expansions
 *
 *  \authors
 *    Developers : \n
 *        Roman Iakymchuk  -- roman.iakymchuk@lip6.fr \n
 *        Sylvain Collange -- sylvain.collange@inria.fr \n
 */
#ifndef EXMTS_FPE_HPP_
#define EXMTS_FPE_HPP_

#include "superaccumulator.hpp"


struct dmts {
    double asum;
    double mts;
};

/**
 * \struct FPExpansionTraits
 * \ingroup ExMTS
 * \brief This struct is meant ot specify optimization or other technique used
 */

template<bool EX=false, bool FLUSHHI=false, bool H2SUM=false, bool CRF=false, bool CSWAP=false, bool B2SUM=true, bool SORT=false, bool VICT=false>
struct FPExpansionTraits
{
    static bool constexpr EarlyExit = EX;
    static bool constexpr FlushHi = FLUSHHI;
    static bool constexpr Horz2Sum = H2SUM;
    static bool constexpr CheckRangeFirst = CRF;
    static bool constexpr ConditionalSwap = CSWAP;
    static bool constexpr Biased2Sum = B2SUM;
    static bool constexpr Sort = SORT;
    static bool constexpr Victimcache = VICT;
};

/**
 * \struct FPExpansionVect
 * \ingroup ExMTS
 * \brief This struct is meant to introduce functionality for working with
 *  floating-point expansions in conjuction with superaccumulators
 */
template<int N, typename TRAITS=FPExpansionTraits<false,false> >
struct FPExpansionVectM1
{
    /**
     * Constructor
     * \param sa superaccumulator
     */
    FPExpansionVectM1(Superaccumulator & sa, Superaccumulator & sa2);

    /**
     * This function accumulates value x to the floating-point expansion
     * \param x input value
     */
    void XAccumulate(double x);

    /**
     * This function is used to flush the floating-point expansion to the superaccumulator
     */
    void Flush();
private:
    // Store the sum
    Superaccumulator & sum_superacc;
    double sum[N] __attribute__((aligned(32)));
    // Store the mts
    Superaccumulator & mts_superacc;
    double mts[N] __attribute__((aligned(32)));
};

template<int N, typename TRAITS>
FPExpansionVectM1<N,TRAITS>::FPExpansionVectM1(Superaccumulator & sa, Superaccumulator & sa2) :
    sum_superacc(sa),
    mts_superacc(sa2)
{
    std::fill(sum, sum + N, 0);
    std::fill(mts, mts + N, 0);
}

// Knuth 2Sum.

inline static double Knuth2Sum(double a, double b, double & s)
{
    double r = a + b;
    double z = r - a;
    double t1 = r - z;
    t1 = a - t1;
    z = b - z;
    s = t1 + z;
    return r;
}


template<int N, typename TRAITS> UNROLL_ATTRIBUTE
void FPExpansionVectM1<N,TRAITS>::XAccumulate(double x)
{
    double s1, s2;
    double y;
    double _sum = 0.;
    double _mts = 0.;
    y = x;

    for(unsigned int i = 0; i != N; ++i) {

        sum[i] = Knuth2Sum(sum[i], x, s1);
        _sum += sum[i];
        x = s1;

        mts[i] = Knuth2Sum(mts[i], y, s2);
        _mts += mts[i];
        y = s2;

        if(TRAITS::EarlyExit && i != 0 && x == 0 && y == 0) break;
    }
// This could be optimized a lot...
    if(_mts < 0. - y) {
        for(unsigned int i = 0; i != N; ++i) {
            mts[i] = 0.;
        }
    }

    if(x != 0) sum_superacc.Accumulate(x);

    if(y != 0) mts_superacc.Accumulate(y);
}

template<int N, typename TRAITS>
void FPExpansionVectM1<N,TRAITS>::Flush()
{
    for(unsigned int i = 0; i != N; ++i)
    {
        sum_superacc.Accumulate(sum[i]);
        sum[i] = 0.;
        mts_superacc.Accumulate(mts[i]);
//        mts[i] = 0.;
    }
}

#endif // EXSUM_FPE_HPP_
