/* Implementation of parallel maximum prefix sum with interval arithmetic.
 * Author: Raphael Dang-Nhu 
 * Date: 04/19/2018 */

#include <emmintrin.h>
#include "tbb/blocked_range.h"

using namespace tbb;

// Structure containing the result of a mps evaluation. 
struct __mps{
    // Both sum and maximum prefix sum are represented as intervals  
    __m128d sum_interval;
    __m128d mps_interval;
    // Position of the maximum prefix sum
    int position;
    // Int size
    int size;
    // Constructor
    __mps();
    // Splitting constructor
    __mps(__mps&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<double*>&);
    // Join operation for the reduction
    void join(__mps& rightMps); 
    // Printing function
    void print_mps();
};

/* Function computing in parallel the maximum prefix sum with interval arithmetic.
 * Arguments: input array, size of the array
 * Returns: __mps struct
 */
__mps parallel_ia_mps(double*,int);
