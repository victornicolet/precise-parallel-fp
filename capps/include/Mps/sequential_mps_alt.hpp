/* Sequential mps implementations 
 * Author: Raphaël Dang-Nhu
 * DAte: May 7th */


#include "interval_arithmetic.hpp"

/* This function computes the maximum prefix sum of an array with superaccumulators */
void sequential_mps_superacc_alt(double*,int,double*,double*,int*);

/* This function computes the maximum prefix sum of an array with mpfr */
void sequential_mps_mpfr_alt(double*,int,double*,double*,int*);

/* This function computes the maximum prefix sum of an array with doubles */
void sequential_mps_double_alt(double*,int,double*,double*,int*);

/* These function compute the maximum prefix sum/position in a lazy way
   Last argument is the memorization of decisions */
void sequential_mps_interval_memorized_alt(double*,int,double*,double*,int*,boolean*&);

void sequential_mps_iterate_reverse_mps_alt(double*,int,double*,double*,int*,boolean*);

void sequential_mps_iterate_reverse_pos_alt(double*,int,double*,double*,int*,boolean*);

void sequential_mps_lazy_superacc_alt(double*,int,double*,double*,int*,boolean*);

void sequential_mps_lazy_mpfr_alt(double*,int,double*,double*,int*,boolean*);
