/* Optimized functions for interval arithmetic 
 * Source: https://locklessinc.com/articles/interval_arithmetic/ 
 * */

// Function to get the lower bound of an interval
double in2_min(__m128d);

// Function to get the upper bound of an interval
double in2_max(__m128d);

// Function to get the real value stored as the lower bound (opposite of the lower bound)
double in2_rmin(__m128d);

// Function to create an interval, first argument is lower bound, second argument is upper bound
__m128d in2_create(double,double);

// Same function, except that the first argument must be the opposite of the lower bound
__m128d in2_raw_create(double,double);

// Function to add to intervals. CAUTION: THE ROUNDING MODE MUST BE SET TOWARDS INFINITY, with  _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
__m128d in2_add(__m128d,__m128d);

// Functions below have been implemented by Raphael Dang-Nhu 

// Function to compare two intervals
