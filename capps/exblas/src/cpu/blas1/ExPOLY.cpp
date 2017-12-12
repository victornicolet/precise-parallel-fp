/*
 *  Copyright (c) 2016 Inria and University Pierre and Marie Curie
 *  All rights reserved.
 */

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "ExPOLY.hpp"
#include "blas1.hpp"

#ifdef EXBLAS_TIMING
#define iterations 50
#endif


/*
 * Parallel summation using our algorithm
 * If fpe < 2, use superaccumulators only,
 * Otherwise, use floating-point expansions of size FPE with superaccumulators when needed
 * early_exit corresponds to the early-exit technique
 */
__mts expoly(int Ng, double *ag, double factor, int fpe, bool early_exit) {
#ifdef EXBLAS_MPI
    int np = 1, p, err;
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
#endif
    if (fpe < 0) {
        fprintf(stderr, "Size of floating-point expansion should be a positive number. Preferably, it should be in the interval [2, 8]\n");
        exit(1);
    }

    int N;
    double *a;
#ifdef EXBLAS_MPI
    Superaccumulator acc, acc_fin;
    N = Ng / np + Ng % np;

    a = (double *)_mm_malloc(N * sizeof(double), 32);
    if (!a)
        fprintf(stderr, "Cannot allocate memory for per process array\n");

    int i;
    if (p == 0) {
        //distribute
        a = ag;
        ag = ag + N;
        for (i = 1; i < np; i++) {
            err = MPI_Send(ag + (i - 1)  * (N - Ng % np), N - Ng % np, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            if (err != MPI_SUCCESS)
                fprintf(stderr, "MPI_Send does not word properly %d\n", err);
        }
    } else {
        MPI_Status status;
        err = MPI_Recv(a, N - Ng % np, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        if (err != MPI_SUCCESS)
            fprintf(stderr, "MPI_Recv does not word properly %d\n", err);
    }
#else
    N = Ng;
    a = ag;
#endif
    // thee

    // with superaccumulators only
    if (fpe < 2) {
        int nthread = tbb::task_scheduler_init::automatic;
        tbb::task_scheduler_init tbbinit(nthread);
        return ExPOLYSuperacc(N, a, factor, 0, 0);
    }

    if (early_exit) {
        if (fpe <= 4)
            return (ExPOLYFPE<Poly_FPExp<4, FPExpansionTraits<true> > >)(N, a, factor, 0, 0);
        if (fpe <= 6)
            return (ExPOLYFPE<Poly_FPExp<6, FPExpansionTraits<true> > >)(N, a, factor, 0, 0);
        if (fpe <= 8)
            return (ExPOLYFPE<Poly_FPExp<8, FPExpansionTraits<true> > >)(N, a, factor, 0, 0);
    } else { // ! early_exit
        if (fpe == 2)
            return (ExPOLYFPE<Poly_FPExp<2> >)(N, a, factor, 0, 0);
        if (fpe == 3)
            return (ExPOLYFPE<Poly_FPExp<3> >)(N, a, factor, 0, 0);
        if (fpe == 4)
            return (ExPOLYFPE<Poly_FPExp<4> >)(N, a, factor, 0, 0);
        if (fpe == 5)
            return (ExPOLYFPE<Poly_FPExp<5> >)(N, a, factor, 0, 0);
        if (fpe == 6)
            return (ExPOLYFPE<Poly_FPExp<6> >)(N, a, factor, 0, 0);
        if (fpe == 7)
            return (ExPOLYFPE<Poly_FPExp<7> >)(N, a, factor, 0, 0);
        if (fpe == 8)
            return (ExPOLYFPE<Poly_FPExp<8> >)(N, a, factor, 0, 0);
    }

    return {0.0, 0.0};
}

/*
 * Our alg with superaccumulators only
 */
__mts ExPOLYSuperacc(int N, double *a, double factor, int inca, int offset) {
    double dacc, final_factor;
#ifdef EXBLAS_TIMING
    double t, mint = 10000;
    uint64_t tstart, tend;
    for(int iter = 0; iter != iterations; ++iter) {
    	tstart = rdtsc();
#endif

    TBBlongPoly tbblongPoly(a, factor);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, N, inca), tbblongPoly);
#ifdef EXBLAS_MPI
    tbBlongPoly.acc.Normalize();
        std::vector<int64_t> result(tbBlongPoly.acc.get_f_words() + tbBlongPoly.acc.get_e_words(), 0);
        //MPI_Reduce((int64_t *) &tbBlongPoly.acc.accumulator[0], (int64_t *) &acc_fin.accumulator[0], get_f_words() + get_e_words(), MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&(tbBlongPoly.acc.get_accumulator()[0]), &(result[0]), tbBlongPoly.acc.get_f_words() + tbBlongPoly.acc.get_e_words(), MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        Superaccumulator acc_fin(result);
        dacc = acc_fin.Round();
#else
    dacc = tbblongPoly.acc.Round();
    final_factor = tbblongPoly.factor;
#endif

#ifdef EXBLAS_TIMING
    tend = rdtsc();
        t = double(tend - tstart) / N;
        mint = std::min(mint, t);
    }
    fprintf(stderr, "%f ", mint);
#endif

    return {dacc, final_factor};
}

/**
 * \brief Parallel reduction step
 *
 * \param step step among threads
 * \param tid1 id of the first thread
 * \param tid2 id of the second thread
 * \param acc1 superaccumulator of the first thread
 * \param acc2 superaccumulator of the second thread
 */
inline static void ReductionStep(int step, int tid1, int tid2,
                                 Superaccumulator * acc1, Superaccumulator * acc2,
                                 double * fac1, double * fac2,
                                 int volatile * ready1, int volatile * ready2)
{
    _mm_prefetch((char const*)ready2, _MM_HINT_T0);
    // Wait for thread 2
    while(*ready2 < step) {
        // wait
        _mm_pause();
    }
    double cr;
    cr = (*fac1) * (acc2->Round());
    acc1->Accumulate(cr);
    (*fac1) *= (*fac2);
}

/**
 * \brief Final step of summation -- Parallel reduction among threads
 *
 * \param tid thread ID
 * \param tnum number of threads
 * \param acc superaccumulator
 */
inline static void Reduction(unsigned int tid, unsigned int tnum, std::vector<int32_t>& ready,
                             std::vector<Superaccumulator>& acc,
                             std::vector<double>& factors, int const linesize)
{
    // Custom reduction
    for(unsigned int s = 1; (1 << (s-1)) < tnum; ++s)
    {
        int32_t volatile * c = &ready[tid * linesize];
        ++*c;
        if(tid % (1 << s) == 0) {
            unsigned int tid2 = tid | (1 << (s-1));
            if(tid2 < tnum) {
                //acc[tid2].Prefetch(); // No effect...
                ReductionStep(s, tid, tid2, &acc[tid], &acc[tid2],
                              &factors[tid], &factors[tid2],
                              &ready[tid * linesize], &ready[tid2 * linesize]);
            }
        }
    }
}

template<typename CACHE> __mts ExPOLYFPE(int N, double *a, double factor, int inca, int offset) {
    // OpenMP sum+reduction
    int const linesize = 16;    // * sizeof(int32_t)
    int maxthreads = omp_get_max_threads();
    double dacc;
#ifdef EXBLAS_TIMING
    double t, mint = 10000;
    uint64_t tstart, tend;
    for(int iter = 0; iter != iterations; ++iter) {
        tstart = rdtsc();
#endif
    std::vector<Superaccumulator> acc(maxthreads);
    std::vector<double> factors(maxthreads, factor);
    std::vector<int32_t> ready(maxthreads * linesize);

#pragma omp parallel
    {
        unsigned int tid = omp_get_thread_num();
        unsigned int tnum = omp_get_num_threads();

        CACHE cache(acc[tid], factors[tid]);
        *(int32_t volatile *)(&ready[tid * linesize]) = 0;  // Race here, who cares?

        int l = ((tid * int64_t(N)) / tnum) & ~7ul;
        int r = ((((tid+1) * int64_t(N)) / tnum) & ~7ul) - 1;

        for(int i = l; i < r; i+=2) {
            asm ("# myloop");

            cache.Accumulate(a[i], a[i+1]);
        }
        cache.Flush();
        acc[tid].Normalize();
        factors[tid] = cache.factor;

        Reduction(tid, tnum, ready, acc, factors, linesize);
    }
#ifdef EXBLAS_MPI
    acc[0].Normalize();
        std::vector<int64_t> result(acc[0].get_f_words() + acc[0].get_e_words(), 0);
        MPI_Reduce(&(acc[0].get_accumulator()[0]), &(result[0]), acc[0].get_f_words() + acc[0].get_e_words(), MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        //MPI_Reduce((int64_t *) &acc[0].accumulator[0], (int64_t *) &acc_fin.accumulator[0], get_f_words() + get_e_words(), MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        Superaccumulator acc_fin(result);
        dacc = acc_fin.Round();
#else
    dacc = acc[0].Round();
    factor = factors[0];
#endif

#ifdef EXBLAS_TIMING
    tend = rdtsc();
        t = double(tend - tstart) / N;
        mint = std::min(mint, t);
    }
    fprintf(stderr, "%f ", mint);
#endif

    return {dacc, factor};
}

