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
template<typename T, int N, typename TRAITS=FPExpansionTraits<false,false> >
struct Poly_FPExp
{
    /**
     * Constructor
     * \param sa superaccumulator
     */
    Poly_FPExp(Superaccumulator & sa, double f);

    /** 
     * This function accumulates value x to the floating-point expansion
     * \param x input value
     */
    void Accumulate(T x);

    /** 
     * This function accumulates two values x to the floating-point expansion
     * \param x1 input value
     * \param x2 input value
     */
    void Accumulate(T x1, T x2);

    /**
     * This function is used to flush the floating-point expansion to the superaccumulator
     */
    void Flush();

    /**
     * This function is meant to be used for printing the floating-point expansion
     */
//    void Dump() const;
private:
    void FlushVector(T x) const;
    static T twosum(T a, T b, T & s);
    static T twomul(T a, T b, T &s);
    static T split(T a, T &s);
    Superaccumulator & superacc;

    // Most significant digits first!
    T a[N] __attribute__((aligned(32)));
    T factor;
    double f;
    T victim;
};


template<typename T, int N, typename TRAITS>
Poly_FPExp<T,N,TRAITS>::Poly_FPExp(Superaccumulator & sa, double _f) :
    superacc(sa),
    f(_f),
    victim(0),
    factor(1)
{
    factor = 1;
    std::fill(a, a + N, 0);
}

// Knuth 2Sum.
template<typename T>
inline static T Knuth2Sum(T a, T b, T & s)
{
    T r = a + b;
    T z = r - a;
    s = (a - (r - z)) + (b - z);
    return r;
}

// Vector impl with test for fast path
template<typename T>
inline static T BiasedSIMD2Sum(T a, T b, T & s)
{
    T r = a + b;
    auto doswap = abs(b) > abs(a);
    //if(unlikely(!_mm256_testz_pd(doswap, doswap)))
    //asm("nop");
    if(/*unlikely*/(!_mm256_testz_si256(_mm256_castpd_si256(doswap), _mm256_castpd_si256(b))))  // any(doswap && b != +0)
    {
        // Slow path
        T a2 = select(doswap, b, a);
        T b2 = select(doswap, a, b);
        a = a2;
        b = b2;
    }
    s = (a - r) + b;
    return r;
}

#if INSTRSET > 7                       // AVX2 and later


// Knuth 2Sum.
template<typename T>
inline static T FMA2Sum(T a, T b, T & s)
{
//    T r = a + b;
//    T z = _fms(1., r, a);
//    s = _fma(1., a - _fms(1., r, z), b - z);
    T r = a + b;
    T z = r -a;
    s = a - (r - z) + (b - z);
    return r;
}
#endif

template<typename T, int N, typename TRAITS> UNROLL_ATTRIBUTE
void Poly_FPExp<T,N,TRAITS>::Accumulate(T x)
{
    T s;
    double elt[4];
    factor.store_a(elt);
    // elt = {x, x*x, x*x*x, x^4}
    for(unsigned int i = 0; i < 4; ++i) {
        elt[i] *= elt[3];
    }
    // elt = {x*x, x*x, x*x*x, x^4}
    factor.load_a(elt);
    x = factor * x;
    for(unsigned int i = 0; i != N; ++i) {
        a[i] = twosum(a[i], x, s);
        x = s;
        if(TRAITS::EarlyExit && i != 0 && !horizontal_or(x)) return;
    }
    if(TRAITS::EarlyExit || horizontal_or(x)) {
        FlushVector(x);
    }
}


template<typename T, int N, typename TRAITS>
T Poly_FPExp<T,N,TRAITS>::twosum(T a, T b, T & s)
{
#if INSTRSET > 7                       // AVX2 and later
	// Assume Haswell-style architecture with parallel Add and FMA pipelines
	return FMA2Sum(a, b, s);
#else
    if(TRAITS::Biased2Sum) {
        return BiasedSIMD2Sum(a, b, s);
    }
    else {
        return Knuth2Sum(a, b, s);
    }
#endif
}

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

template <typename T, int N, typename TRAITS>
T Poly_FPExp<T,N,TRAITS>::twomul(T a, T b, T &s) {
#if INSTRSET > 7
    FMA2Mul(a, b, s);
#else
    Dekker2Mul(a, b, s);
#endif
}


template<typename T, int N, typename TRAITS> UNROLL_ATTRIBUTE INLINE_ATTRIBUTE
void Poly_FPExp<T,N,TRAITS>::Accumulate(T x1, T x2)
{

    T s1, s2;
    for(unsigned int i = 0; i != N; ++i) {
        T ai = Vec4d().load_a((double*)(a+i));
        //T ai = a[i];
        ai = twosum(ai, x1, s1);
        ai = twosum(ai, x2, s2);
        ai.store_a((double*)(a+i));
        //a[i] = ai;
        x1 = s1;
        x2 = s2;
        if(TRAITS::EarlyExit && i != 0 && !horizontal_or(x1|x2)) return;
    }

    // Separate checks
    if(unlikely(horizontal_or(x1))) {
            FlushVector(x1);
    }
    if(unlikely(horizontal_or(x2))) {
            FlushVector(x2);
    }
}
#undef IACA
#undef IACA_START
#undef IACA_END

template<typename T, int N, typename TRAITS>
void Poly_FPExp<T,N,TRAITS>::Flush()
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

template<typename T, int N, typename TRAITS> inline
void Poly_FPExp<T,N,TRAITS>::FlushVector(T x) const
{
    // TODO: update status, handle Inf/Overflow/NaN cases
    // TODO: make it work for other values of 4
    double v[4];
    x.store(v);
    
//    _mm256_zeroupper();
    for(unsigned int j = 0; j != 4; ++j) {
        superacc.Accumulate(v[j]);
    }
}

#endif // EXSUM_FPE_HPP_
