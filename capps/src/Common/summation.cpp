
#include <iostream>

#include "summation.hpp"
#include "superaccumulator.hpp"
#include "debug.hpp"


using namespace std;

void sequential_summation_superacc(double* array, int size, double* sum){
    Superaccumulator sumA = Superaccumulator();
    for(int i = 0; i != size; i++){
        sumA.Accumulate(array[i]);
    }
    *sum = sumA.Round();

    if(PRINT){
        cout << endl << "Summation with superaccumulators" << endl;
        cout << "Sum: " << *sum << endl;
    }

}
