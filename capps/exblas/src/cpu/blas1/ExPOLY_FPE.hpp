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

#define f64factor (1 << 27) + 1.

#include "superaccumulator.hpp"

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
struct Poly_FPExp
{
    double factor;
    /**
     * Constructor
     * \param sa superaccumulator
     */
    Poly_FPExp(Superaccumulator & sa, double f);

    /** 
     * This function accumulates value x to the floating-point expansion
     * \param x input value
     */
    void Accumulate(double x);

    /** 
     * This function accumulates two values x to the floating-point expansion
     * \param x1 input value
     * \param x2 input value
     */
    void Accumulate(double x1, double x2);

    /**
     * This function is used to flush the floating-point expansion to the superaccumulator
     */
    void Flush();

    /**
     * This function is meant to be used for printing the floating-point expansion
     */
//    void Dump() const;
private:
    void FlushVector(double x) const;
    static double twosum(double a, double b, double & s);
    static double twomul(double a, double b, double &s);
    Superaccumulator & superacc;

    // Most significant digits first!
    double a[N] __attribute__((aligned(32)));

    double init_f;
    double victim;
};


template<int N, typename TRAITS>
Poly_FPExp<N,TRAITS>::Poly_FPExp(Superaccumulator & sa, double _f) :
    superacc(sa),
    init_f(_f),
    victim(0),
    factor(1)
{
    factor = 1;
    std::fill(a, a + N, 0);
}

// Knuth 2Sum.
template<typename T>
inline static double Knuth2Sum(double a, double b, double & s)
{
    double r = a + b;
    double z = r - a;
    s = (a - (r - z)) + (b - z);
    return r;
}


#if INSTRSET > 7                       // AVX2 and later


// Knuth 2Sum.
inline static double FMA2Sum(double a, double b, double & s)
{
//    double r = a + b;
//    double z = _fms(1., r, a);
//    s = _fma(1., a - _fms(1., r, z), b - z);
    double r = a + b;
    double z = r -a;
    s = a - (r - z) + (b - z);
    return r;
}
#endif


template<int N, typename TRAITS>
double Poly_FPExp<N,TRAITS>::twosum(double a, double b, double & s)
{
#if INSTRSET > 7                       // AVX2 and later
	// Assume Haswell-style architecture with parallel Add and FMA pipelines
	return FMA2Sum(a, b, s);
#else
    return Knuth2Sum(a, b, s);
#endif
}

inline double split(double a, double &s) {
    double c, x;
    c = f64factor * a;
    x = c - (c - a);
    s = a -x;
    return x;
}

inline double Dekker2Mul(double a,double b, double &s){
    double x, a1, a2, b1, b2;
    x = a * b;
    a1 = split(a, a2); // a = a1 + a2 and a1, a2 non overlapping
    b1 = split(b, b2); // a = a1 + a2 and a1, a2 non overlapping
    s = a2 * b2 - (((x - a1 * b1) - a2 * b1) - a1 * b2);
    return x;
}

#if INSTRSET > 7
inline double FMA2Mul(double a, double b, double &s){
    double x;
    x = a * b;
    s = fma(a, b, -x);
    return x;
}
#endif

template <int N, typename TRAITS>
double Poly_FPExp<N,TRAITS>::twomul(double a, double b, double &s) {
#if INSTRSEdouble > 7
    FMA2Mul(a, b, s);
#else
    Dekker2Mul(a, b, s);
#endif
}


template<int N, typename TRAITS> UNROLL_ATTRIBUTE
void Poly_FPExp<N,TRAITS>::Accumulate(double x)
{
    double s;
    factor *= init_f;
    x = factor * x;
    for(unsigned int i = 0; i != N; ++i) {
        a[i] = twosum(a[i], x, s);
        x = s;
        if(TRAITS::EarlyExit && i != 0 && x == 0) return;
    }
    if(TRAITS::EarlyExit || x != 0) {
        FlushVector(x);
    }
}

template<int N, typename TRAITS> UNROLL_ATTRIBUTE INLINE_ATTRIBUTE
void Poly_FPExp<N,TRAITS>::Accumulate(double x1, double x2)
{
    double s1, s2;
    factor *= init_f;
    x1 *= factor;
    factor *= init_f;
    x2 *= factor;
    for(unsigned int i = 0; i != N; ++i) {
        //double ai = a[i];
        a[i] = twosum(a[i], x1, s1);
        a[i] = twosum(a[i], x2, s2);
        //a[i] = ai;
        x1 = s1;
        x2 = s2;
        if(TRAITS::EarlyExit && i != 0 && !horizontal_or(Vec4d(x1)|Vec4d(x2))) return;
    }

    // Separate checks
    if(unlikely(x1 != 0)) {
            FlushVector(x1);
    }
    if(unlikely(x2 != 0)) {
            FlushVector(x2);
    }
}
#undef IACA
#undef IACA_START
#undef IACA_END

template<int N, typename TRAITS>
void Poly_FPExp<N,TRAITS>::Flush()
{
    for(unsigned int i = 0; i != N; ++i)
    {
        FlushVector(a[i]);
        a[i] = 0;
    }
    if(TRAITS::Victimcache) {
        FlushVector(victim);
    }
}

template<int N, typename TRAITS> inline
void Poly_FPExp<N,TRAITS>::FlushVector(double x) const
{
    superacc.Accumulate(x);
}

#endif // EXSUM_FPE_HPP_
