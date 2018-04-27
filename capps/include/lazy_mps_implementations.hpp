/* Implementation of lasy parallel maximum prefix sum with interval arithmetic and superaccumulators.
 * Author: Raphael Dang-Nhu 
 * Date: 04/19/2018 */

#include <emmintrin.h>
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"

using namespace tbb;
using namespace std;

// Same structure with lazy computation, i.a. and superaccumulators. 
struct __mps{
    // pointer to the array
    double* array;
    // Both sum and maximum prefix sum are represented as intervals  
    __m128d sum_interval;
    __m128d mps_interval;
    // Position of the maximum prefix sum
    int position;
    // size
    int size;
    int left;
    int right;
    // Constructor
    __mps(double* a);
    // Splitting constructor
    __mps(__mps&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<int>&);
    // Join operation for the reduction
    void join(__mps& rightMps); 
    // Printing function
    void print_mps();
};

/* MPS computation with floats */
void parallel_mps_float(double*,int);

/* Function computing in parallel the mps with superaccumulators, same arguments as previous function */
void parallel_mps_superacc(double*,int);

/* Function computing mps with mpfr */
void parallel_mps_mpfr(double*,int);

/* Function computimg mps with Collange mixed FPE-superacc solution (not working yet)*/
void parallel_mps_Collange(double*,int);

/* Lazy computation of mps, with superaccs for precise computations */
void parallel_mps_superacc_lazy(double*,int);

/* Lazy computations of mps, with mpfr for precise computations */
void parallel_mps_mpfr_lazy(double*,int);

