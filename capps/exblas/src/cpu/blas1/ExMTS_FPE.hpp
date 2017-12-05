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

    /**
     * This function is meant to be used for printing the floating-point expansion
     */
    void Dump() const;
private:
    void Insert(double& x);
    void Insert(double & x1, double & x2);
    static void Swap(double & x1, double & x2);
    static double twosum(double a, double b, double & s);

    Superaccumulator & sum_superacc;
    Superaccumulator & mts_superacc;

    // Most significant digits first!
    // Store the sum
    double sum[N] __attribute__((aligned(32)));
    // Store the mts
    double mts[N] __attribute__((aligned(32)));
    double mtsbuf[N] __attribute__((aligned(32)));

    double victim;
};

template<int N, typename TRAITS>
FPExpansionVectM1<N,TRAITS>::FPExpansionVectM1(Superaccumulator & sa, Superaccumulator & sa2) :
    sum_superacc(sa),
    mts_superacc(sa2),
    victim(0.)
{
    std::fill(sum, sum + N, 0);
    std::fill(mts, mts + N, 0);
    std::fill(mtsbuf, mtsbuf + N, 0);
}

// Knuth 2Sum.

inline static double Knuth2Sum(double a, double b, double & s)
{
    double r = a + b;
    double z = r - a;
    s = (a - (r - z)) + (b - z);
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

        sum[i] = twosum(sum[i], x, s1);
        _sum += sum[i];
        x = s1;

        mtsbuf[i] = twosum(mtsbuf[i], y, s2);
        _mts += mtsbuf[i];
        y = s2;

        if(TRAITS::EarlyExit && i != 0 && x == 0 && y == 0) break;
    }
// This could be optimized a lot...
    if(_mts > 0.) {
        for(unsigned int i = 0; i != N; ++i) {
            mts[i] = mtsbuf[i];
        }
    } else {
        for(unsigned int i = 0; i != N; ++i) {
            mtsbuf[i] = mts[i];
        }
    }

    if(x != 0) sum_superacc.Accumulate(x);

    if(y != 0) mts_superacc.Accumulate(y);
}


template<int N, typename TRAITS>
double FPExpansionVectM1<N,TRAITS>::twosum(double a, double b, double & s)
{
	return Knuth2Sum(a, b, s);
}

inline static void swap_if_nonzero(double & a, double & b)
{
    if (a != 0) std::swap<double>(a, b);
}

template<int N, typename TRAITS>
void FPExpansionVectM1<N,TRAITS>::Swap(double & x1, double & x2)
{
    if(TRAITS::ConditionalSwap) {
        swap_if_nonzero(x1, x2);
    }
    else {
        std::swap<double>(x1, x2);
    }
}

template<int N, typename TRAITS> UNROLL_ATTRIBUTE
void FPExpansionVectM1<N, TRAITS>::Insert(double & x)
{
    if(TRAITS::Sort) {
        // Insert at tail. Unconditional version.
        // Rotate accumulators:
        // x <= a[0]
        // a[0] <= a[1]
        // a[1] <= a[2]
        // ...
        // a[N-2] <= a[N-1]
        // a[N-1] <= x
        //T xb = a[0];
        double xb = sum[0];
        for(int i = 0; i != N-1; ++i)
        {
            sum[i] = sum[i+1];
        }
        //a[N-1] = x;
        sum[N-1] = x;
        x = xb;
    }
    else {
        // Insert at head
        // Conditional or unconditional
        Swap(x, sum[0]);
    }
}

template<int N, typename TRAITS> UNROLL_ATTRIBUTE
void FPExpansionVectM1<N,TRAITS>::Insert(double & x1, double & x2)
{
    if(TRAITS::Sort) {
        // x1 <= a[0]
        // x2 <= a[1]
        // a[0] <= a[2]
        // a[1] <= a[3]
        // a[i] <= a[i+2]
        // a[N-3] <= a[N-1]
        // a[N-2] <= x1
        // a[N-1] <= x2
        double x1b = sum[0];
        double x2b = sum[1];
        for(int i = 0; i != N-2; ++i) {
            sum[i] = sum[i+2];
        }
        sum[N-2] = x1;
        sum[N-1] = x2;
        x1 = x1b;
        x2 = x2b;
    }
    else {
        Swap(x1, sum[0]);
        Swap(x2, sum[1]);
    }
}



#undef IACA
#undef IACA_START
#undef IACA_END

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


template<int N, typename TRAITS>
void FPExpansionVectM1<N,TRAITS>::Dump() const
{
    for(unsigned int i = 0; i != N; ++i)
    {
        std::cout << sum[i];
    }
    std::cout << std::endl;
    for(unsigned int i = 0; i != N; ++i)
    {
        std::cout << sum[i];
    }
    std::cout << std::endl;
}

#endif // EXSUM_FPE_HPP_
