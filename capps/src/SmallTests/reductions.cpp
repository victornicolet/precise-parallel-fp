/* Function comparing the different implementations of reduction operations
 * Author: Raphael Dang-Nhu
 * Date: 04/24/2018 */


#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

#include "tbb/task_group.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"

#include<common.hpp>
#include "pfpdefs.hpp"

using namespace tbb;
using namespace std;

#define DEBUG 0

/* Homemade reduction operation */
class MpsTask1: public task {
    public:
        MpsTask1(int C, double* a,int size_, double*s, double*m,int*p, int left_ = 0, int right_ = -1) : 
            Cutoff(C),
            array(a),
            size(size_),
            sum(s),
            mps(m),
            position(p),
            left(left_)
        {
            if(right_ == -1){
                right_ = size;
            }
            right = right_;

        }
        ~MpsTask1(){}
        
        task* execute(){
            if(size <= Cutoff){

                for(int i = left; i != right; i++){
                    *sum += array[i];
                    if(*sum >= *mps){
                       *mps = *sum;
                       *position = i+1;
                    }
                }

            }else{
                // Parameters for subtasks
                int middle = (right+left)/2;

                // Variables for results
                int rPos = middle;
                int sizel = left - middle;
                int sizer = middle - right;
                double rsum = 0, rmps = 0;

                // Create subtasks
                set_ref_count(3);

                MpsTask1& lTask = *new(allocate_child()) MpsTask1(Cutoff,array,sizel,sum,mps,position,left,middle);
                
                spawn(lTask);

                MpsTask1 &rTask = *new(allocate_child()) MpsTask1(Cutoff,array,sizer,&rsum, &rmps,&rPos,middle,right);
                
                spawn_and_wait_for_all(rTask);

                rmps = *sum + rmps;
                *sum = *sum + rsum;
                if(rmps >= *mps){
                    *mps = rmps;
                    *position = rPos;
                }
            }
            return NULL;
        }

    private:
        /* Below this size, the mps and sum are computed sequentially */
        int Cutoff;
        /* Input array and its size */
        double* array;
        int size;
        /* Intervals for sum and mps */
        double* sum;
        double* mps;
        /* Position of the mps */
        int* position;
        /* Bounds */
        int left;
        int right;
};

// Parallel reduce implementation 
struct __mps{
    // pointer to the array
    double* array;
    // Superaccumulators for sum and mps
    double sum;
    double mps;
    // Position of the maximum prefix sum
    int position;
    // Constructor
    __mps(double* a);
    // Splitting constructor
    __mps(__mps&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<int>&);
    // Join operation for the reduction
    void join(__mps& rightMps); 
    // Printing function
    void print_mps();
};

__mps::__mps(double* a):
    array(a),
    sum(0),
    mps(0),
    position(0)
{}

__mps::__mps(__mps& x, split) :
    array(x.array),
    sum(0),
    mps(0),
    position(0)
{}

void __mps::print_mps(){
    cout << "sum: " << sum;
    cout << endl << "mps: " << mps;
    cout << endl << "position: " << position << endl;
}

void __mps::operator()(const blocked_range<int>& r){
    if(position == 0){
        position = r.begin();
    }
    for(int i = r.begin(); i != r.end(); i++){
        sum += array[i];
        if(sum >= mps){
            mps = sum;
            position = i+1; 
        }
    }
}

void __mps::join(__mps& rightMps){
   // computing sum-l + mps-r
   rightMps.mps += sum;
   // adding two sums
   sum += rightMps.sum;
   // comparison of mpsCandidate and mps-l
   if(rightMps.mps >= mps){
       mps = rightMps.mps;
       position = rightMps.position;
   }
}

// Parallel reduce main function
void tbb_main(double* a, int size){
    __mps result = __mps(a);
    parallel_reduce(blocked_range<int>(0,size),result);
    if(DEBUG){
        cout << endl << "Parallel Reduce" << endl;
        result.print_mps();
    }
}

// Homemade main function
void homemade_main(double *a, int size,int grainsize){
    
    double sum = 0, mps = 0;
    int pos = 0;

    MpsTask1& root = *new(task::allocate_root()) MpsTask1(grainsize,a,size,&sum,&mps,&pos);

    task::spawn_root_and_wait(root);

    if(DEBUG){
        cout << endl <<"Homemade reduction" << endl;
        cout << "Sum: " << sum << endl << "Mps " << mps << endl << "Pos: " << pos << endl;
    }
}

// Function to compare the implementations
void runtime_comparison(){
    
    // Variables declaration and initialisation 
    double start;
    int size = pow(10,7);
    int N = 2;

    // for each dynamic range
    vector<int> grainsSizes  {100,300,1000,2000,3000,6000,10000,20000,30000};
    int s = grainsSizes.size();
    
    // Store results to plot
    fstream results;
    results.open("Plots/reductions.csv", ofstream::out | ofstream::trunc);
    vector<double> x(s),r0(s),r1(s);

    // Random seed
    srand(time(NULL));
    
    for(int r = 0; r < s; r++){

        // initialization of means
        double mean_tbb = 0., mean_homemade = 0.; 

        for(int i = 0; i < N; i++){
            
            // Generating array
            double* drray = new double[size];
            init_fpuniform(size, drray, 100, 50);

            // Randomly change signs
            for(int j = 0; j < size ; j++){
                 drray[j] = (rand() % 2) ? drray[j] : -drray[j];
            }
            
            // Declare result variables
            double time_tbb = 0.0;
            PFP_TIME(tbb_main(drray,size),start,time_tbb);
            double time_homemade = 0.0;
            PFP_TIME(homemade_main(drray,size,grainsSizes[r]),start,time_homemade);
       
            mean_tbb += time_tbb;
            mean_homemade += time_homemade;
            
            delete[] drray;
        }
        // Finalize mean computation
        mean_tbb += mean_tbb / N;
        mean_homemade += mean_homemade /N;
        
        x[r]= grainsSizes[r];
        r0[r]= mean_tbb/mean_tbb;
        r1[r]= mean_homemade/ mean_tbb;

        // Writing results to a file
        results << to_string(x[r]) << "," << to_string(r0[r]) << "," << to_string(r1[r]) << endl;
        
        
    }

    results.close();
}

int main(){
    runtime_comparison();
    return 0;
}
