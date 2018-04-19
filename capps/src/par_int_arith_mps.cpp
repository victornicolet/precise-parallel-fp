/* Implementation of parallel maximum prefix sum with interval arithmetic.
 * Author: Raphael Dang-Nhu 
 * Date: 04/19/2018 */

#include "par_int_arith_mps.hpp"
#include "interval_arithmetic.hpp"

void __mps::mps_join(__mps& rightMps){
   // computing sum-l + mps-r
   __m128d mpsCandidate = in2_add(sum_interval,rightMps.mps_interval);
   // adding two sums 
   sum_interval = in2_add(sum_interval,rightMps.sum_interval);
   // comparison of mpsCandidate and mps-l

}

__mps parallel_ia_mps(double* array, int size){

}
