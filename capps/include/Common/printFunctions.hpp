#include <vector>
#include <emmintrin.h>

// Function to print a vector (the template type t must be std::cout - able 
template<class t> void  printVector(std::vector<t> v){
    for(typename std::vector<t>::iterator it = v.begin(); it != v.end(); it++){
        std::cout << *it << ",";
    }
    std::cout << std::endl;
}


// Function to print a vector of intervals
void printVectorInt(std::vector<__m128d> v);
