/* Sequential mps implementations 
 * Author: Raphaël Dang-Nhu
 * DAte: May 7th */


#include "interval_arithmetic.hpp"

/* This function computes the maximum prefix sum of an array with superaccumulators */
void sequential_mps_superacc(double*,int,double*,double*,int*);

/* This function computes the maximum prefix sum of an array with doubles */
void sequential_mps_double(double*,int,double*,double*,int*);

/* These function compute the maximum prefix sum/position in a lazy way
   Last argument is the memorization of decisions */
void sequential_mps_interval_memorized(double*,int,double*,double*,int*,memo**);

void sequential_mps_iterate_reverse_mps(double*,int,double*,double*,int*,memo**);

void sequential_mps_iterate_reverse_pos(double*,int,double*,double*,int*,memo**);

void sequential_mps_lazy_superacc(double*,int,double*,double*,int*,memo**);

void sequential_mps_lazy_mpfr(double*,int,double*,double*,int*,memo**);
