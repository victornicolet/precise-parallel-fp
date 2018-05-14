/* Sequential mts implementations 
 * Author: Raphaël Dang-Nhu
 * DAte: May 7th */


/* This function computes the maximum prefix sum of an array with superaccumulators */
void sequential_mts_superacc(double*,long,double*,long*);

/* This function computes the maximum prefix sum of an array with doubles */
void sequential_mts_double(double*,long,double*,long*);

/* This function computes the maximum prefix sum in a lazy way
 * Last argument is for optional optimization method */
void sequential_mts_lazy(double*,long,double*,long*,long);
