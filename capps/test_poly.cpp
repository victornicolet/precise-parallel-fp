//
// Created by nicolet on 05/12/17.
//

#include "test_poly.h"
#include <mm_malloc.h>
#include <string>
#include <fstream>
#include <iostream>


using namespace std;


__mts sequential_poly(int N, double* a, double f){
    double x = 1.;
    double sum = 0.;
    for(int i = 0; i  < N; i++){
        x *= f;
        sum += x * a[i];
    }
    return {sum, x};
}


__mts polyjoin(__mts r,__mts n) {
    // r is the already reduced value
    // n is the new value
    // sum is the sum, mts the factor
    r.sum = r.sum + r.mts * n.sum;
    r.mts = r.mts * n.mts;
    return r;
}

__mts inexact_parallel_poly(int N, double* a, double factor){

    __mts m = {0., 1.};

#pragma omp declare reduction \
      (poly_reduction:__mts:omp_out=polyjoin(omp_out,omp_in)) \
      initializer(omp_priv={0.,1.})

#pragma omp parallel for reduction(poly_reduction:m)
    for (int idata=0; idata<N; idata++){
        m.mts *= factor;
        m.sum += m.mts * a[idata];
    }

    return m;
}

void m_test_poly(int argc, char** argv) {
    string outputcsv, err_outputcsv;
    outputcsv = string(__FUNCTION__).append(".csv");
    err_outputcsv = string(__FUNCTION__).append("_errlog.csv");
    int NUM_RUNS = 20;

    double factor = 0.99;

    fstream  fp, fperr;
    fp.open(outputcsv, ios::app);
    fperr.open(err_outputcsv, ios::app);

    double eps = 1e-16;
    int dec = 20;


    if (argc > 2) {
        dec = atoi(argv[2]);
    }

    int N = 1 << dec;
    dec = ((1 << (dec + 1)) - N) / 4;

    for(int intermediate_size = 0; intermediate_size < 4 * dec; intermediate_size += dec ) {
        double *a;

        N += intermediate_size;

        printf("------------Size : %i --------------\n", N);

        a = (double *) _mm_malloc(N * sizeof(double), 32);

        if (!a)
            fprintf(stderr, "Cannot allocate memory for the main array\n");

        for (int initmode = 0; initmode < 3; initmode++) {
            switch (initmode) {
                case 0:
                    init_naive(N, a);
                    break;
                case 1:
                    init_fpuniform(N, a, 10, -3);
                    break;
                case 2:
                    init_ill_cond(N, a, 0.1);
                    break;
                default:
                    break;
            }

            fprintf(stderr, "%d ", N);

            bool is_pass = true;
            __mts inex_poly, expoly_acc, expoly_fpe2, expoly_fpe4, expoly_fpe4ee, expoly_fpe6ee, expoly_fpe8ee;
            double time_expoly[7] = {0.};
            double wtime_expoly[7] = {0.};
            double start, wstart;

            //        Time a sequential implementation for reference
            __mts seq_res;
            double seqtime, seqwtime;
            PFP_WTIME(seq_res = sequential_poly(N, a, factor), start, seqtime, wstart, seqwtime)

            for (int run_no = 0; run_no < NUM_RUNS; run_no++) {
                ACC_PFP_WTIME(inex_poly = inexact_parallel_poly(N, a, factor), start, time_expoly[0], wstart,
                              wtime_expoly[0])
                ACC_PFP_WTIME(expoly_acc = expoly(N, a, factor, 0, false), start, time_expoly[1], wstart,
                              wtime_expoly[1])
                ACC_PFP_WTIME(expoly_fpe2 = expoly(N, a, factor, 2, false), start, time_expoly[2], wstart,
                              wtime_expoly[2])
                ACC_PFP_WTIME(expoly_fpe4 = expoly(N, a, factor, 4, false), start, time_expoly[3], wstart,
                              wtime_expoly[3])
                ACC_PFP_WTIME(expoly_fpe4ee = expoly(N, a, factor, 4, true), start, time_expoly[4], wstart,
                              wtime_expoly[4])
                ACC_PFP_WTIME(expoly_fpe6ee = expoly(N, a, factor, 6, true), start, time_expoly[5], wstart,
                              wtime_expoly[5])
                ACC_PFP_WTIME(expoly_fpe8ee = expoly(N, a, factor, 8, true), start, time_expoly[6], wstart,
                              wtime_expoly[6])

                double expoly_sacc_esum, expoly_fpe2_esum, expoly_fpe4_esum, expoly_fpe4ee_esum, expoly_fpe6ee_esum, expoly_fpe8ee_esum;
                double inexpoly_esum;


                //            Compare the results to the sequential sum
                /// Error reporting --------------------------------------------------------------------------------------
                inexpoly_esum = fabs(seq_res.sum - inex_poly.sum) / fabs(seq_res.sum);
                expoly_sacc_esum = fabs(seq_res.sum - expoly_acc.sum) / fabs(seq_res.sum);
                expoly_fpe2_esum = fabs(seq_res.sum - expoly_fpe2.sum) / fabs(seq_res.sum);
                expoly_fpe4_esum = fabs(seq_res.sum - expoly_fpe4.sum) / fabs(seq_res.sum);
                expoly_fpe4ee_esum = fabs(seq_res.sum - expoly_fpe4ee.sum) / fabs(seq_res.sum);
                expoly_fpe6ee_esum = fabs(seq_res.sum - expoly_fpe6ee.sum) / fabs(seq_res.sum);
                expoly_fpe8ee_esum = fabs(seq_res.sum - expoly_fpe8ee.sum) / fabs(seq_res.sum);
                if ((expoly_fpe2_esum > eps) || (expoly_fpe4_esum > eps) || (expoly_fpe4ee_esum > eps) ||
                    (expoly_fpe6ee_esum > eps) || (expoly_fpe8ee_esum > eps)) {
                    is_pass = false;
                    printf("ERROR for inexact: %.16g\n", inexpoly_esum);
                    printf("FAILED:\t %.16g %.16g \t %.16g \n\t\t%.16g \t %.16g \t %.16g\n",
                           expoly_sacc_esum, expoly_fpe2_esum, expoly_fpe4_esum, expoly_fpe4ee_esum,
                           expoly_fpe6ee_esum, expoly_fpe8ee_esum);
                    fperr << N << "," << initmode << ","
                          << inexpoly_esum << ","
                          << expoly_sacc_esum << ","
                          << expoly_fpe2_esum << ","
                          << expoly_fpe4_esum << ","
                          << expoly_fpe4ee_esum << ","
                          << expoly_fpe6ee_esum << ","
                          << expoly_fpe8ee_esum << endl;
                }

                fprintf(stderr, "\n");

                if (is_pass)
                    printf("TestPassed; ALL OK!\n");
                else
                    printf("TestFailed!\n");
            }
            /// Speedup reporting -----------------------------------------------------------------------------------------
            fp << N << "," << initmode;
            for (int exno = 0; exno < 7; exno++) {
                fp << "," << seqwtime * (NUM_RUNS / wtime_expoly[exno]);
            }
            fp << endl;
        }
        free(a);
    }

    fp.flush();
    fp.close();

    fperr.flush();
    fperr.close();

}