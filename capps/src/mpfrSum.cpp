#include <algorithm>
#include <iostream>
#include <common.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>

using namespace std;

int debug = 0;
int print = 0;

struct mpsResult{
    mpsResult(){
        mpfr_init2(sum,10);
        mpfr_init2(mps,10);
    }
    ~mpsResult(){
        mpfr_clear(sum);
        mpfr_clear(mps);
    }
    mpfr_t sum;
    mpfr_t mps;
};

struct mpsResult recursiveMps(double array[], int l, int r){
    
   struct mpsResult result;
   if (r == l) {
       mpfr_set_d(result.sum,0.,MPFR_RNDD);
       mpfr_set_d(result.mps,0.,MPFR_RNDD);
   }
   else if (l+1 == r){
       mpfr_set_d(result.sum,array[l],MPFR_RNDD);
       mpfr_set_d(result.mps,(array[l] > 0) ? array[l] : 0,MPFR_RNDD);
   }
   else{
       int middle = (r+l)/2;
       result = recursiveMps(array, l, middle);
       struct mpsResult rightResult = recursiveMps(array, middle, r);
       mpfr_add(result.sum,rightResult.sum,MPFR_RNDD); 
       /*result.mps = (aux > result.mps) ? aux : result.mps;*/
   }

   if(debug){
       cout << "Trace: l = " << l << ", r = " << r << ", sum = " << result.sum << endl;
   }

   return result;
}

double mps(double array[], int size){
   struct mpsResult res = recursiveMps(array, 0, size);
   return res.sum;
}

void printArray(double array[], int size){
    cout << endl;
    for(int i = 0; i < size; i++){
        cout << array[i] << ", ";
    }
    cout << endl;
}

double simpleSum(double array[], int size){
    double sum = 0;
    for(int i = 0; i < size; i++){
        sum += array[i];
    }
    return sum;
}

void test(){
   int size = 100000000;
   cout << size << endl;
   double* a = new double[size];
   init_ill_cond(size, a, 1);
   /*a[0] = 1;
   a[1] = -1;
   a[2] = 1e-14;
    */
   
   srand(time(NULL));
   for(int i = 0; i < size ; i++){
        a[i] = (rand() % 2) ? a[i] : -a[i];
   }
   

   if( print) printArray(a,size);
   cout << endl << simpleSum(a,size) << endl;
   
   for(int i = 0; i < 20; i++){
       random_shuffle(&a[0],&a[size]);
       if( print) printArray(a,size);
       cout << endl << simpleSum(a,size) << endl;
   }

   if( print) printArray(a,size);
   cout << endl << mps(a,size) << endl;
   
   for(int i = 0; i < 20; i++){
       random_shuffle(&a[0],&a[size]);
       if( print) printArray(a,size);
       cout << endl << mps(a,size) << endl;
   }
}

int main(){
    test();
    return 0;
}
