#include "pfpdefs.h"
#include "par_precise_fp.h"
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>
#include "omp.h"

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

void usage(){
    printf("Usage: precise_parallel_fp [-n size] [-o output-file] [-t tests-numbers] [-v]");
}

int main(int argc, char** argv) {
//    Parse options
    int verbose = 0;
    int number = 0;
    int opt;
    msize_t n = 2 << 26;
    char* filename = "output.csv";
    FILE *fp;

    while ((opt = getopt(argc, argv, "n:o:t:v")) != -1)
    {
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
    return 0;
}