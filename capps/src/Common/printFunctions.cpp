#include <iostream>

#include <vector>
#include <emmintrin.h>
#include "printFunctions.hpp"
#include "interval_arithmetic.hpp"

using namespace std;


void printVectorInt(vector<__m128d> v){
    for(vector<__m128d>::iterator it = v.begin(); it != v.end(); it++){
        print(*it);
        cout << endl;
    }
}
