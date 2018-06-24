/* Sequential mps implementations 
 * Author: Raphaël Dang-Nhu
 * DAte: May 7th */


/* This function computes the maximum prefix sum of an array with superaccumulators */
void sequential_mps_superacc(double*,int,double*,double*,int*);

/* This function computes the maximum prefix sum of an array with doubles */
void sequential_mps_double(double*,int,double*,double*,int*);

/* This function computes the maximum prefix sum in a lazy way
 * Last argument is for optional optimization method */
void sequential_mps_lazy(double*,int,double*,double*,int*,int);
