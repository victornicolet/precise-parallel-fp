#include "pfpdefs.hpp"
#include <stdio.h>
#include <stdbool.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <mm_malloc.h>
#include "omp.h"
#include <getopt.h>
#include <time.h>
#include "par_precise_fp.hpp"
#include "pfpdefs.hpp"
#ifdef EXBLAS_MPI
#include <mpi.h>

#endif
// exblas
#include "blas1.hpp"
#include "common.hpp"

using namespace std;

#define TIME(start, call, mem)\
start = omp_get_wtime();\
call;\
mem = (double) (omp_get_wtime() - start);

#define VERBOSE(v,c)\
if(v){\
c;\
}

bool file_exists(const char * filename) {
    FILE* file;
    if (file = fopen(filename, "r")) {
        fclose(file);
        return true;
    }
    return false;
}

bool f_is_empty(std::fstream& pFile)
{
    return pFile.peek() == std::fstream::traits_type::eof();
}

void usage(char* arg){
    printf("Usage: precise_parallel_fp [-n size] [-o output-file] [-t tests-numbers] [-v]");
}

// _exblas_ usage: command [N=20] [stddev=1.] [mean=1.] [ill_cond_or_lognormal_or_uniform= (i | n)]
int _exblas_(int argc, char** argv) {
    double eps = 1e-16;
    int N = 1 << 20;
    bool lognormal = false;
    if(argc > 1) {
        N = 1 << atoi(argv[1]);
    }
    if(argc > 4) {
        if(argv[4][0] == 'n') {
            lognormal = true;
        }
    }

    int range = 1;
    int emax = 0;
    double mean = 1., stddev = 1.;
    if(lognormal) {
        stddev = strtod(argv[2], 0);
        mean = strtod(argv[3], 0);
    }
    else {
        if(argc > 2) {
            range = atoi(argv[2]);
        }
        if(argc > 3) {
            emax = atoi(argv[3]);
        }
    }

    double *a;

    a = (double*)_mm_malloc(N*sizeof(double), 32);
    if (!a)
        fprintf(stderr, "Cannot allocate memory for the main array\n");
    if(lognormal) {
        fprintf(stderr, "Data initialization: lognormal.\n");
        init_lognormal(N, a, mean, stddev);
    } else if ((argc > 4) && (argv[4][0] == 'i')) {
        fprintf(stderr, "Data initialization: ill conditioned.\n");
        init_ill_cond(N, a, range);
    } else {
        if(range == 1){
            fprintf(stderr, "Data intialization: naive.\n");
            init_naive(N, a);
        } else {
            fprintf(stderr, "Data initialization: fp uniform.\n");
            init_fpuniform(N, a, range, emax);
        }
    }

    fprintf(stderr, "%d ", N);

    if(lognormal) {
        fprintf(stderr, "%f ", stddev);
    } else {
        fprintf(stderr, "%d ", range);
    }

    bool is_pass = true;
    double exsum_acc, exsum_fpe2, exsum_fpe4, exsum_fpe4ee, exsum_fpe6ee, exsum_fpe8ee;
    exsum_acc = exsum(N, a, 1, 0, false);
    exsum_fpe2 = exsum(N, a, 1, 2, false);
    exsum_fpe4 = exsum(N, a, 1, 4, false);
    exsum_fpe4ee = exsum(N, a, 1, 4, true);
    exsum_fpe6ee = exsum(N, a, 1, 6, true);
    exsum_fpe8ee = exsum(N, a, 1, 8, true);

    printf("  exsum with superacc = %.16g\n", exsum_acc);
    printf("  exsum with FPE2 and superacc = %.16g\n", exsum_fpe2);
    printf("  exsum with FPE4 and superacc = %.16g\n", exsum_fpe4);
    printf("  exsum with FPE4 early-exit and superacc = %.16g\n", exsum_fpe4ee);
    printf("  exsum with FPE6 early-exit and superacc = %.16g\n", exsum_fpe6ee);
    printf("  exsum with FPE8 early-exit and superacc = %.16g\n", exsum_fpe8ee);


    exsum_fpe2 = fabs(exsum_acc - exsum_fpe2) / fabs(exsum_acc);
    exsum_fpe4 = fabs(exsum_acc - exsum_fpe4) / fabs(exsum_acc);
    exsum_fpe4ee = fabs(exsum_acc - exsum_fpe4ee) / fabs(exsum_acc);
    exsum_fpe6ee = fabs(exsum_acc - exsum_fpe6ee) / fabs(exsum_acc);
    exsum_fpe8ee = fabs(exsum_acc - exsum_fpe8ee) / fabs(exsum_acc);
    if ((exsum_fpe2 > eps) || (exsum_fpe4 > eps) || (exsum_fpe4ee > eps) || (exsum_fpe6ee > eps) || (exsum_fpe8ee > eps)) {
        is_pass = false;
        printf("FAILED: %.16g \t %.16g \t %.16g \t %.16g \t %.16g\n", exsum_fpe2, exsum_fpe4, exsum_fpe4ee, exsum_fpe6ee, exsum_fpe8ee);
    }

    fprintf(stderr, "\n");

    if (is_pass)
        printf("TestPassed; ALL OK!\n");
    else
        printf("TestFailed!\n");
#ifdef EXBLAS_MPI
    }
    MPI_Finalize();
#endif

    return 0;
}

void test_my_implems(int argc, char** argv){
    //    Parse options
    int verbose = 0;
    int number = 0;
    int opt;
    msize_t n = 2 << 10;
    char* filename = "output.csv";
    FILE *fp;

    while ((opt = getopt(argc, argv, "n:o:t:v")) != -1)    {
        switch (opt)
        {
            case 'o':
                filename = optarg;
                break;
            case 'n':
                number = atoi(optarg);
                n = number;
                break;
            case 'v':
                verbose = 1;
                break;
            default:
                usage(argv[0]);
        }
    }

    bool filext = file_exists(filename);

    fp = fopen(filename, "a");

    if(!filext){
        fprintf(fp, "n, seqtime, par_time, pfx_tot_time, pfx_cor_time, error\n");
    }

    srand((unsigned int)time(NULL));
    double precision = 1e-18;
    double start;

    double *a = (double*)malloc(sizeof(a)* n);
    double *p = (double*)calloc(n, sizeof(p));
    double *ep = (double*)calloc(n, sizeof(ep));
    double *epp = (double*)calloc(n, sizeof(epp));

    VERBOSE(verbose, printf("n = %i\n", n))
    frand_init(a, n);

// Sequential
    double seq_sum = 0.0;
    double seqtime = 0.0;

    TIME(start, seq_sum = dsum(a, n), seqtime)

//    Parallel
    double par_time = 0.0;
    double par_sum;
    TIME(start, par_sum = dsum_par(a, n), par_time)
    VERBOSE(verbose, printf("Error: %f \n", (seq_sum - par_sum) / seq_sum))

//    Prefix sum + error correction
    double pfx_par_time = 0.0;
    double pfx_cor_time = 0.0;
    double pfx_tot_time = 0.0;
    int reclevel = 0;
    TIME(start, dpxsum_par(a, p, n), pfx_par_time)
    TIME(start, dprecise_parallel(a, p, ep, epp, precision, n, &reclevel), pfx_cor_time)
    pfx_tot_time = pfx_par_time + pfx_cor_time;
    double pre_sum = p[n-1];
    double error = (seq_sum - pre_sum)/seq_sum;

    fprintf(fp, "%i,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f\n", n, seqtime, par_time, pfx_tot_time, pfx_cor_time, error);

    if(verbose) {
        printf("Error with correction: %f\n", (seq_sum - pre_sum) / seq_sum);
        printf("\nSequential time: %f seconds.\n", seqtime);
        printf("Parallel wtime (no correction): %f\n", par_time);
        printf("Parallel wtime (wt correction): %f\n", pfx_tot_time);
        printf("Prefix calculation time: %f\n", pfx_par_time);
        printf("Prefix correction time: %f\n", pfx_cor_time);
        printf("Number of recursions: %i\n", reclevel);
    }
    fclose(fp);
    free(p);
    free(ep);
    free(epp);
    free(a);
}

void small_tests(int argc, char** argv){
    int N = 1 << 10;
    int TRIALS = 10;
    int NUM_TESTS = 10;

    if(argc > 1){
        N = 1 << atoi(argv[1]);
    }

    double* a;
    _mts_ mtsres, refmtsres;

    fstream fp;
    fp.open("./output.csv", ios::app);
    a = (double*)_mm_malloc(N*sizeof(double), 32);


    for(int testno = 0; testno < NUM_TESTS; testno++) {
        double testline[6] = {0.};
        for (int initmode = 0; initmode < 3; initmode++) {
            switch (initmode) {
                case 0:
                    init_naive(N, a);
                    break;
                case 1:
                    init_fpuniform(N, a, 10, 3);
                    break;
                case 2:
                    init_ill_cond(N, a, 0.23);
                    break;
                default:
                    break;
            }
//            Execute once for reference and study variation
            refmtsres = myparallelmts(a, N);

            for (int i = 0; i < TRIALS; ++i) {
                mtsres = myparallelmts(a, N);
                testline[2 * initmode] += (refmtsres.sum - mtsres.sum) / TRIALS;
                testline[2 * initmode + 1] += (refmtsres.mts - mtsres.mts) / TRIALS;
            }
        }
        fp << N;
        for(int j = 0; j < 6; j++){
            fp << "," << testline[j];
        }
        fp << endl;
    } // End test number testno
    fp.flush();
    fp.close();
    free(a);
}

int main(int argc, char** argv) {
#ifdef TEST_CUSTOM
    test_my_implems(argc, argv);
#endif
//    Test exblas
//    test_my_exblas(argc, argv);
    small_tests(argc, argv);
    return 0;
}
