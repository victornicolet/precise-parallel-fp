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
#include "tbb/partitioner.h"
#include "tbb/task_scheduler_init.h"

#include <common.hpp>
#include "pfpdefs.hpp"

using namespace tbb;
using namespace std;

#define DEBUG 1

struct mps_struct{
    mps_struct(){};
    double sum;
    double mps;
    int pos;
};

/* Continuation passing class */
class MpsContinuation: public task {
    public :
    mps_struct* lres;
    mps_struct rres;

    MpsContinuation(mps_struct* lres_):
        lres(lres_)
    {}
    
    task* execute(){
        double s = lres->sum;
        double rmps = s + rres.mps;
        lres->sum = s + rres.sum;
        if(rmps >= lres->mps){
            lres->mps = rmps;
            lres->pos = rres.pos;
        }
        return NULL;
    }
};

/* Homemade reduction operation */
class MpsTask1: public task {
    public:
        MpsTask1(int D, double* a, mps_struct* res_, unsigned int left_, unsigned int right_) : 
            depth(D),
            array(a),
            res(res_),
            left(left_),
            right(right_)
        {
        }
        ~MpsTask1(){}
        
        task* execute(){
            if(depth == 0){
                double s = 0.;
                double m = 0.;
                int p = left;
                for(int i = left; i != right; i++){
                    s += array[i];
                    if(s >= m){
                       m = s;
                       p = i+1;
                    }
                }
                res->sum = s;
                res->mps = m;
                res->pos = p;
                return NULL;

            }else{
                // Parameters for subtasks
                unsigned int middle = (right+left) >> 1;
                // Variables for results
                depth--;

                // Create subtasks
                MpsContinuation& c = *new(allocate_continuation()) MpsContinuation(res);

                //MpsTask1& lTask = *new(c.allocate_child()) MpsTask1(Cutoff,array,sizel,sum,mps,position,left,middle);
                MpsTask1 &rTask = *new(c.allocate_child()) MpsTask1(depth,array,&c.rres,middle,right);

                c.set_ref_count(2);
                spawn(rTask);

                recycle_as_child_of(c);
                right = middle;


                return this;
            }
        }

    private:
        /* Below this size, the mps and sum are computed sequentially */
        int depth;
        /* Input array and its size */
        double* array;
        /* Intervals for sum and mps */
        mps_struct* res;
        /* Bounds */
        unsigned int left;
        unsigned int right;
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
void tbb_main(double* a, int size, int grainsize){
    task_scheduler_init init;
    __mps result = __mps(a);
    parallel_reduce(blocked_range<int>(0,size),result);
    if(DEBUG){
        cout << endl << "Parallel Reduce" << endl;
        result.print_mps();
    }
    init.terminate();
}

// Parallel deterministic reduce main function
void tbb_deterministic_main(double* a, int size, int grainsize){
    task_scheduler_init init;
    __mps result = __mps(a);
    parallel_deterministic_reduce(blocked_range<int>(0,size,grainsize),result);
    if(DEBUG){
        cout << endl << "Parallel Deterministic Reduce" << endl;
        result.print_mps();
    }
    init.terminate();
}

// Homemade main function
void homemade_main(double *a, int size,int grainsize){
    
    mps_struct result;
    
    // Computing maxdepth
    double ratio = (double) (size)/(double)grainsize;
    int depth = ceil(log2(ratio));
    if(depth < 0) depth = 0;

    MpsTask1& root = *new(task::allocate_root()) MpsTask1(depth,a,&result,0,size);

    task::spawn_root_and_wait(root);

    if(DEBUG){
        cout << endl <<"Homemade reduction" << endl;
        cout << "Sum: " << result.sum << endl << "Mps " << result.mps << endl << "Pos: " << result.pos << endl;
    }
}

/*
// Reduction step
inline static void ReductionStep(int step, int tid1, int tid2, Superaccumulator * acc1, Superaccumulator * acc2,
    int volatile * ready1, int volatile * ready2)
{
    _mm_prefetch((char const*)ready2, _MM_HINT_T0);
    // Wait for thread 2
    while(*ready2 < step) {
        // wait
        _mm_pause();
    }
    acc1->Accumulate(*acc2);
}

// Reduction function
void Reduction(unsigned int tid, unsigned int tnum){
    // Custom reduction
    for(unsigned int s = 1; (1 << (s-1)) < tnum; ++s) 
    {
        int32_t volatile * c = &ready[tid * linesize];
        ++*c;
        if(tid % (1 << s) == 0) {
            unsigned int tid2 = tid | (1 << (s-1));
            if(tid2 < tnum) {
                ReductionStep(s, tid, tid2, &acc[tid], &acc[tid2],
                    &ready[tid * linesize], &ready[tid2 * linesize]);
            }
        }
    }
}


// Openmp reduction operation
void omp_main(double* a, int size){

    #pragma omp parallel
    {
        // Get thread index
        unsigned int tid = omp_get_thread_num();
        unsigned int tnum = omp_get_num_threads();

        int l = ((tid * int64_t(size)) / tnum) & ~7ul;
        int r = ((((tid+1) * int64_t(size)) / tnum) & ~7ul) - 1;

        // Variables
        double sum = 0.;
        double mps = 0.;
        int pos = l;

        for(int i=l; i < r; i++){
            sum += a[i];

            if(sum >= mps){
                mps = sum;
                pos = i+1;
            }
        }

        Reduction(tid, tnum)
    }
}*/

// Function to compare the implementations
void runtime_comparison(){
    
    // Variables declaration and initialisation 
    double start;
    int size = pow(10,9);
    int N = 1;

    // for each dynamic range
    vector<int> grainsSizes  {100,600,1000,3000,10000,20000,30000};
    //vector<int> grainsSizes  {10};
    int s = grainsSizes.size();
    
    // Store results to plot
    fstream results;
    results.open("Plots/reductions.csv", ofstream::out | ofstream::trunc);
    vector<double> x(s),r0(s),r1(s),r2(s);

    // Random seed
    srand(time(NULL));
    // Generating array
    double* drray = new double[size];
    init_fpuniform(size, drray, 100, 50);

    // Randomly change signs
    for(int j = 0; j < size ; j++){
         drray[j] = (rand() % 2) ? drray[j] : -drray[j];
    }
    
    for(int r = 0; r < s; r++){

        // initialization of means
        double mean_tbb = 0., mean_homemade = 0.;

        for(int i = 0; i < N; i++){
            
            // Declare result variables
            double time_homemade = 0.0;
            //PFP_TIME(homemade_main(drray,size,grainsSizes[r]),start,time_homemade);
            double time_tbb = 0.0;
            PFP_TIME(tbb_main(drray,size,grainsSizes[r]),start,time_tbb);
       
            mean_tbb += time_tbb;
            mean_homemade += time_homemade;
            
        }
        // Finalize mean computation
        mean_tbb += mean_tbb / N;
        mean_homemade += mean_homemade /N;
        
        x[r]= grainsSizes[r];
        r0[r]= mean_tbb;
        r1[r]= mean_homemade;

    }
    for(int r = 0; r < s; r++){

        // initialization of means
        double mean_deterministic_tbb = 0.;

        for(int i = 0; i < N; i++){
            
            // Declare result variables
            double time_deterministic_tbb = 0.0;
            PFP_TIME(tbb_deterministic_main(drray,size,grainsSizes[r]),start,time_deterministic_tbb);
       
            mean_deterministic_tbb += time_deterministic_tbb;
            
        }
        mean_deterministic_tbb += mean_deterministic_tbb / N;
        
        r2[r]= mean_deterministic_tbb;

    }
    for(int r = 0; r < s; r++){
        // Writing results to a file
        results << to_string(x[r]) << "," << to_string(r0[r]) << "," << to_string(r1[r]) << "," << to_string(r2[r]) << endl;
    }    

    results.close();
}

int main(){
    runtime_comparison();
    return 0;
}
