/* Parallel mss implementations.
 * Author: Raphaël Dang-Nhu
 * Date: 11/05/18 */

using namespace tbb;
using namespace std;

// Naive mss structure
struct __mss_naive{
    // pointer to the array
    double* array;
    // Superaccumulators for sum and mss
    double sum;
    double mps;
    double mts;
    double mss;
    // Position of the maximum prefix sum
    long posmps;
    long posmts;
    long posmssl;
    long posmssr;
    
    // Constructor
    __mss_naive(double* a);
    // Splitting constructor
    __mss_naive(__mss_naive&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<int>&);
    // Join operation for the reduction
    void join(__mss_naive& rightmss); 
    // Printing function
    void print_mss();
};
