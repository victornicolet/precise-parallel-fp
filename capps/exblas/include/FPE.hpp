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
#ifndef EXSUM_FPE_HPP_
#define EXSUM_FPE_HPP_

/**
 * \struct FPExpansionTraits
 * \ingroup ExSUM
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
 * \ingroup ExSUM
 * \brief This struct is meant to introduce functionality for working with
 *  floating-point expansions in conjuction with superaccumulators
 */
template<int N, typename TRAITS=FPExpansionTraits<false,false> >
struct FPExpansionVect
{
    /**
     * Constructor
     * \param sa superaccumulator
     */
    FPExpansionVect(Superaccumulator & sa);

    /** 
     * This function accumulates value x to the floating-point expansion
     * \param x input value
     */
    void Accumulate(double x);

    /**
     * This function is used to flush the floating-point expansion to the superaccumulator
     */
    void Flush();
    Superaccumulator & superacc;

private:


    // Most significant digits first!
    double a[N] __attribute__((aligned(32)));
};

template<int N, typename TRAITS>
FPExpansionVect<N,TRAITS>::FPExpansionVect(Superaccumulator & sa) :
    superacc(sa)
{
    std::fill(a, a + N, 0);
}

inline static double Knuth2Sum(double,double,double&);

// Knuth 2Sum.
inline static double Knuth2Sum(double a, double b, double & s)
{
    double r = a + b;
    double z = r - a;
    s = (a - (r - z)) + (b - z);
    return r;
}

template<int N, typename TRAITS> UNROLL_ATTRIBUTE
void FPExpansionVect<N,TRAITS>::Accumulate(double x)
{
    double s;
    for(unsigned int i = 0; i != N; ++i) {
        a[i] = Knuth2Sum(a[i], x, s);
        x = s;
    }
    if(x != 0) superacc.Accumulate(x);
}



template<int N, typename TRAITS>
void FPExpansionVect<N,TRAITS>::Flush()
{
    for(unsigned int i = 0; i != N; ++i)
    {
        superacc.Accumulate(a[i]);
        a[i] = 0;
    }
}

#endif // EXSUM_FPE_HPP_
