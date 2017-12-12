//
// Created by nicolet on 05/12/17.
//

#include "test_poly_multicore.h"
#include "test_poly.h"
#include <mm_malloc.h>
#include <string>
#include <fstream>
#include <iostream>

#define NUM_IMPLEMS 2

using namespace std;

void m_test_poly_multicore(int argc, char** argv) {

    int maxcores = omp_get_max_threads();
    printf("MAXCORES = %i\n", maxcores);

    string outputcsv, err_outputcsv;
    outputcsv = string(__FUNCTION__).append(".csv");
    err_outputcsv = string(__FUNCTION__).append("_errlog.csv");
    int NUM_RUNS = 10;

    double factor = 0.99;

    fstream  fp, fperr;
    fp.open(outputcsv, ios::app);
    fperr.open(err_outputcsv, ios::app);

    double eps = 1e-16;



    int N = 10000; // 10e7


    if (argc > 2) {
        N = pow(10, atoi(argv[2]));
    }


    int testcores[11] = {1,2,4,6,8,12,16,24,32,48,64};

    for(int icore = 0; testcores[icore] <= maxcores; icore ++){
        double *a;
        int numcores = testcores[icore];
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

            fprintf(stderr, "%d ", numcores);

            bool is_pass = true;
            __mts inex_poly, expoly_fpe2;
            double time_expoly[NUM_IMPLEMS] = {0.};
            double wtime_expoly[NUM_IMPLEMS] = {0.};
            double start, wstart;

            //        Time a sequential implementation for reference
            __mts seq_res;
            double seqtime, seqwtime;
            PFP_WTIME(seq_res = sequential_poly(N, a, factor), start, seqtime, wstart, seqwtime)

            for (int run_no = 0; run_no < NUM_RUNS; run_no++) {
                omp_set_num_threads(numcores);
                ACC_PFP_WTIME(inex_poly = inexact_parallel_poly(N, a, factor), start, time_expoly[0], wstart,
                              wtime_expoly[0])
                ACC_PFP_WTIME(expoly_fpe2 = expoly(N, a, factor, 2, false), start, time_expoly[1], wstart,
                              wtime_expoly[1])

                double expoly_fpe2_esum;
                double inexpoly_esum;


                //            Compare the results to the sequential sum
                /// Error reporting --------------------------------------------------------------------------------------
                inexpoly_esum = fabs(seq_res.sum - inex_poly.sum) / fabs(seq_res.sum);
                expoly_fpe2_esum = fabs(seq_res.sum - expoly_fpe2.sum) / fabs(seq_res.sum);
                if (expoly_fpe2_esum > eps) {
                    is_pass = false;
                    fperr << numcores << "," << initmode << ","
                          << inexpoly_esum << ","
                          << expoly_fpe2_esum << endl;
                }
                if (is_pass)
                    printf("TestPassed; ALL OK!\n");
                else
                    printf("TestFailed!\n");
            }
            /// Speedup reporting -----------------------------------------------------------------------------------------
            fp << numcores << "," << initmode;
            for (int exno = 0; exno < NUM_IMPLEMS; exno++) {
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