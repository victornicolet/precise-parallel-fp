/* Implementation of lazy parallel maximum prefix sum with interval arithmetic and superaccumulators.
 * Author: Raphael Dang-Nhu 
 * Date: 04/19/2018 */

#include <iostream>
#include <math.h>

#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/task_group.h"
#include "tbb/task_scheduler_init.h"

#include "parallel_mps.hpp"
#include "sequential_mps.hpp"
#include "lazy_mps_implementations.hpp"
#include "2_lazy_mps_implementations.hpp"
#include "interval_arithmetic.hpp"
#include "debug.hpp"

using namespace std;
using namespace tbb;

void printA(double* array, int size){
    // printing the array
    cout << "{";
    if(size <= 10){
        for(int i = 0; i < size; i++){
            cout << array[i] << ",";
        }
        cout << "}" << endl;
    }
    else{
        for(int i = 0; i < 10; i++){
            cout << array[i] << ",";
        }
        cout << "... (size: " << size << ")" << endl;
    }

}

void parallel_mps_float(double* array, int size){
    task_scheduler_init init;
    __mps_naive result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    if(PRINT){
        cout << endl << "Parallel double" << endl;
        //printA(array,size);
        result.print_mps();
    }
}

void parallel_mps_Collange(double* array, int size){
    task_scheduler_init init;
    __mps_precise<4> result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    if(PRINT){
        cout << endl << "Parallel Collange" << endl;
        printA(array,size);
        result.print_mps();
    }
}

void parallel_mps_superacc(double* array, int size){
    task_scheduler_init init;
    __mps_acc result(array);
    parallel_reduce(blocked_range<int>(0,size),result);
    if(PRINT){
        cout << endl << "Parallel superacc" << endl;
        //printA(array,size);
        result.print_mps();
    }
}

void parallel_mps_superacc_lazy(double* array, int size){
    task_scheduler_init init;
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);

    // First step of computation with interval arithmetic
    __mps result(array);
    parallel_reduce(blocked_range<int>(0,size),result);

    // Print intermediate result
    if(PRINT){
        //printA(array,size);
        cout << endl << "Parallel superacc lazy, first results" << endl;
        result.print_mps();
    }
    init.terminate();

    // Second step of computation
    if(result.lposition != result.rposition){
        task_scheduler_init init;
        _MM_SET_ROUNDING_MODE(0);
        __mps_acc resultAcc(array);
        parallel_reduce(blocked_range<int>(result.lposition,result.rposition),resultAcc);
        if(PRINT){
            cout <<  "Parallel superacc lazy, precise results" << endl;
            resultAcc.print_mps();
        }
    }
}

void parallel_mps_mpfr(double* array, int size){
    task_scheduler_init init;
    __mps_mpfr result(array);
    parallel_reduce(blocked_range<long>(0,size),result);
    if(PRINT){
        cout << endl << "Parallel mpfr" << endl;
        //printA(array,size);
        result.print_mps();
    }
} 

void parallel_mps_mpfr_lazy(double* array, int size){
    task_scheduler_init init;
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);

    // First step of computation with interval arithmetic
    __mps result(array);
    parallel_reduce(blocked_range<int>(0,size),result);

    // Print intermediate result
    if(PRINT){
        //printA(array,size);
        cout << endl << "Parallel mpfr lazy, first results" << endl;
        result.print_mps();
    }

    // Second step of computation
    if(result.lposition != result.rposition){
        __mps_mpfr resultM(array);
        parallel_reduce(blocked_range<long>(result.lposition,result.rposition),resultM);
        if(PRINT){
            cout <<  "Parallel mpfr lazy, precise results" << endl;
            resultM.print_mps();
        }
    }
}

void parallel_mps_mpfr_lazy_2(double* array, int size, int grainsize){

    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    
    // Parameters
    int Cutoff = grainsize;
    
    // Create variables for result
    mps_struct result;
    
    // Create matrix to store the information
    double ratio = (double) (size)/(double)Cutoff;
    int ceilRatio = ceil(log2(ratio));
    if(ceilRatio < 1) ceilRatio = 1;
    const int maxDepth = ceilRatio + 1;
    
    int maxIndex = 1;
    Status** memo = new Status*[maxDepth]; 
    for(int i = 0; i < maxDepth;i++){
        memo[i] = new Status[maxIndex];
        maxIndex = 2*maxIndex;
    }

    int validity = 0;
    MpsTask1& root = *new(task::allocate_root()) MpsTask1(Cutoff,array,size,&result,memo,0,0,0,size);

    task::spawn_root_and_wait(root);

    // Get precise position
    Superaccumulator sum = Superaccumulator();
    Superaccumulator mps = Superaccumulator();
    int position2 = 0;

    MpsTask2& root2 = *new(task::allocate_root()) MpsTask2(Cutoff,array,size,&sum,&mps,&position2,memo);

    task::spawn_root_and_wait(root2);

    if(PRINT){
        // Printing result
        cout << endl << "Parallel mpfr lazy 2"<< endl;
        //printA(array,size);

        cout <<"First result" << endl;
        cout << "Sum: ";
        print(result.sum) ;
        cout << endl << "Mps: ";
        print(result.mps) ;
        cout << endl << "Position: " << result.pos << endl;
        if(result.val == 0){
            cout << "Valid position" << endl;
        } else{
            cout << "Unvalid position" << endl;
        }
        cout << "Precise result" << endl;
        cout << "Sum: " << sum.Round();
        cout << endl << "Mps: " << mps.Round();
        cout << endl << "Position: " << position2 << endl;
        

        // Printing memo
        cout << endl;
        maxIndex= 1;
        for(int i = 0; i!=maxDepth; i++){
            for(int j = 0; j!= maxIndex; j++){
                switch(memo[i][j]){
                    case undefinedComparison :
                        cout << "0";
                        break;
                    case cutoffPrecise : 
                        cout << "0";
                        break;
                    case leftChild :
                        cout << "|";
                        break;
                    case rightChild :
                        cout << "\\";
                        break;
                    case cutoff :
                        cout << "_";
                        break;
                
                }
            }
            maxIndex = 2*maxIndex;
            cout << endl;
        }
    }
}


__mps::__mps(double* a) :
    array(a),
    lposition(0),
    rposition(0)
{
    sum_interval = in2_create(0.,0.);
    mps_interval = in2_create(0.,0.);
}

__mps::__mps(__mps& x, split) :
    array(x.array),
    lposition(0),
    rposition(0)
{
    sum_interval = in2_create(0.,0.);
    mps_interval = in2_create(0.,0.);
}

void __mps::print_mps(){
    cout << "sum: ";
    print(sum_interval);
    cout << endl << "mps: ";
    print(mps_interval);
    cout << endl << "left position: " << lposition;
    cout << endl << "right position: " << rposition << endl;
}

void __mps::join(__mps& rightMps){
   // computing sum-l + mps-r
   __m128d mpsCandidate = in2_add(sum_interval,rightMps.mps_interval);
   // adding two sums
   sum_interval = in2_add(sum_interval,rightMps.sum_interval);
   // comparison of mpsCandidate and mps-l
   boolean b = inferior(mps_interval,mpsCandidate);
   if (b == True){
       mps_interval = mpsCandidate;
       lposition = rightMps.lposition;
       rposition = rightMps.rposition;
   }
   // In case of undefined comparison:
   else if(b == Undefined){
       mps_interval = in2_max(mps_interval,mpsCandidate);
       rposition = rightMps.rposition;
   } 
}

void __mps::operator()(const blocked_range<int>& r){

    if(lposition == 0 && rposition == 0){
        lposition = r.begin();
        rposition = r.begin();
    }
    // iterating over the subrange
    for(int i = r.begin(); i != r.end(); i++){
        sum_interval = in2_add_double(sum_interval,array[i]);
        boolean b = inferior(mps_interval,sum_interval);
        if (b == True){
            mps_interval = sum_interval;
            lposition = i+1; 
            rposition = i+1; 
        }
        else if(b == Undefined){
            mps_interval = in2_max(mps_interval,sum_interval);
            rposition = i+1; 
        }
    }
}
