//
// Created by nicolet on 05/12/17.
//

#include "test_mts.h"
#include <mm_malloc.h>
#include <string>
#include <fstream>
#include <iostream>


using namespace std;


__mts sequential_mts(int N, double* a){
    double mts = 0.;
    double sum = 0.;
    for(int i = 0; i  < N; i++){
        sum += a[i];
        mts = max(mts + a[i], 0.);
    }
    return {sum, mts};
}


__mts mtsjoin(__mts r,__mts n) {
    // r is the already reduced value
    // n is the new value
    double new_mts;
    new_mts = n.sum + r.mts;
    if (new_mts < n.mts){
        r.mts = n.mts;
    } else {
        r.mts = new_mts;
    }
    r.sum = r.sum + n.sum;
    return r;
}

__mts inexact_parallel_mts(int N, double* a){

    __mts m = {0., 0.};

#pragma omp declare reduction \
      (mts_reduction:__mts:omp_out=mtsjoin(omp_out,omp_in)) \
      initializer(omp_priv={0,0})

#pragma omp parallel for reduction(mts_reduction:m)
    for (int idata=0; idata<N; idata++){
        m.sum += a[idata];
        m.mts = (m.mts + a[idata] > 0) ? m.mts + a[idata] : 0;
    }

    return m;
}

void m_test_mts(int argc, char** argv) {
    string outputcsv, err_outputcsv;
    outputcsv = string(__FUNCTION__).append(".csv");
    err_outputcsv = string(__FUNCTION__).append("_errlog.csv");
    int NUM_RUNS = 10;

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

        printf("------------Size : %i --------------", N);

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
            __mts inex_mts, exmts_acc, exmts_fpe2, exmts_fpe4, exmts_fpe4ee, exmts_fpe6ee, exmts_fpe8ee;
            double time_exmts[7] = {0.};
            double wtime_exmts[7] = {0.};
            double start, wstart;

            //        Time a sequential implementation for reference
            __mts seq_res;
            double seqtime, seqwtime;
            PFP_WTIME(seq_res = sequential_mts(N, a), start, seqtime, wstart, seqwtime)

            //        ExSUM for reference on SUM
            double exsum_res = 0.;
            exsum_res = exsum(N, a, 0, 0, 4, false);

            for (int run_no = 0; run_no < NUM_RUNS; run_no++) {
                ACC_PFP_WTIME(inex_mts = inexact_parallel_mts(N, a), start, time_exmts[0], wstart, wtime_exmts[0])
                ACC_PFP_WTIME(exmts_acc = exmts(N, a, 0, false), start, time_exmts[1], wstart, wtime_exmts[1])
                ACC_PFP_WTIME(exmts_fpe2 = exmts(N, a, 2, false), start, time_exmts[2], wstart, wtime_exmts[2])
                ACC_PFP_WTIME(exmts_fpe4 = exmts(N, a, 4, false), start, time_exmts[3], wstart, wtime_exmts[3])
                ACC_PFP_WTIME(exmts_fpe4ee = exmts(N, a, 4, true), start, time_exmts[4], wstart, wtime_exmts[4])
                ACC_PFP_WTIME(exmts_fpe6ee = exmts(N, a, 6, true), start, time_exmts[5], wstart, wtime_exmts[5])
                ACC_PFP_WTIME(exmts_fpe8ee = exmts(N, a, 8, true), start, time_exmts[6], wstart, wtime_exmts[6])

                printf("  exsum whit FPE4                 = %.16g\n", exsum_res);
                printf("  exmts sequential                = %.16g  mts = %.16g\n", seq_res.sum, seq_res.mts);
                printf("  exmts naive inexact             = %.16g  mts = %.16g\n", inex_mts.sum, inex_mts.mts);
                printf("  exmts with superacc             = %.16g  mts = %.16g\n", exmts_acc.sum, exmts_acc.mts);
                printf("  exmts with FPE2 and superacc    = %.16g  mts = %.16g\n", exmts_fpe2.sum, exmts_fpe2.mts);
                printf("  exmts with FPE4 and superacc    = %.16g  mts = %.16g\n", exmts_fpe4.sum, exmts_fpe2.mts);
                printf("  exmts with FPE4ee and superacc  = %.16g  mts = %.16g\n", exmts_fpe4ee.sum,
                       exmts_fpe4ee.mts);
                printf("  exmts with FPE6ee and superacc  = %.16g  mts = %.16g\n", exmts_fpe6ee.sum,
                       exmts_fpe6ee.mts);
                printf("  exmts with FPE8ee and superacc  = %.16g  mts = %.16g\n", exmts_fpe8ee.sum,
                       exmts_fpe8ee.mts);

                double exmts_sacc_esum, exmts_fpe2_esum, exmts_fpe4_esum, exmts_fpe4ee_esum, exmts_fpe6ee_esum, exmts_fpe8ee_esum;
                double inexmts_esum;
                double exsum_res_esum;
                //            Compare the results to the sequential sum

                exsum_res_esum = fabs(seq_res.sum - exsum_res) / fabs(seq_res.sum);
                if (exsum_res_esum > eps) {
                    printf("FAILED for EXSUM: error of %.16g\n", exsum_res_esum);
                }

                inexmts_esum = fabs(seq_res.sum - inex_mts.sum) / fabs(seq_res.sum);
                exmts_sacc_esum = fabs(seq_res.sum - exmts_acc.sum) / fabs(seq_res.sum);
                exmts_fpe2_esum = fabs(seq_res.sum - exmts_fpe2.sum) / fabs(seq_res.sum);
                exmts_fpe4_esum = fabs(seq_res.sum - exmts_fpe4.sum) / fabs(seq_res.sum);
                exmts_fpe4ee_esum = fabs(seq_res.sum - exmts_fpe4ee.sum) / fabs(seq_res.sum);
                exmts_fpe6ee_esum = fabs(seq_res.sum - exmts_fpe6ee.sum) / fabs(seq_res.sum);
                exmts_fpe8ee_esum = fabs(seq_res.sum - exmts_fpe8ee.sum) / fabs(seq_res.sum);
                if ((exmts_fpe2_esum > eps) || (exmts_fpe4_esum > eps) || (exmts_fpe4ee_esum > eps) ||
                    (exmts_fpe6ee_esum > eps) || (exmts_fpe8ee_esum > eps)) {
                    is_pass = false;
                    printf("ERROR for inexact: %.16g\n", inexmts_esum);
                    printf("FAILED:\t %.16g %.16g \t %.16g \n\t\t%.16g \t %.16g \t %.16g\n",
                           exmts_sacc_esum, exmts_fpe2_esum, exmts_fpe4_esum, exmts_fpe4ee_esum,
                           exmts_fpe6ee_esum, exmts_fpe8ee_esum);
                    fperr << N << "," << initmode << ","
                          << inexmts_esum << ","
                          << exmts_sacc_esum << ","
                          << exmts_fpe2_esum << ","
                          << exmts_fpe4_esum << ","
                          << exmts_fpe4ee_esum << ","
                          << exmts_fpe6ee_esum << ","
                          << exmts_fpe8ee_esum << endl;
                }

                fprintf(stderr, "\n");

                if (is_pass)
                    printf("TestPassed; ALL OK!\n");
                else
                    printf("TestFailed!\n");
            }
            fp << N << "," << initmode;
            for (int exno = 0; exno < 7; exno++) {
                fp << "," << seqwtime * (NUM_RUNS / wtime_exmts[exno]);
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