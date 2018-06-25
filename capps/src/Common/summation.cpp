
#include <iostream>

#include "summation.hpp"
#include "superaccumulator.hpp"
#include "ExSUM_FPE.hpp"
#include "debug.hpp"


using namespace std;

void sequential_summation_superacc(double* array, int size, double* sum){
    Superaccumulator sumA = Superaccumulator();
    FPExpansionVect<double,4> sumC(sumA);
    for(int i = 0; i != size; i++){
        sumC.Accumulate(array[i]);
    }
    sumC.Flush();
    *sum = sumA.Round();

    if(PRINT){
        cout << endl << "Summation with superaccumulators" << endl;
        cout << "Sum: " << *sum << endl;
    }

}
