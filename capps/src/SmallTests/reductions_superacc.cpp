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

#include "common.hpp"
#include "superaccumulator.hpp"
#include "pfpdefs.hpp"

using namespace tbb;
using namespace std;

#define DEBUGSUPERACC 0

/*
// Continuation passing class
class MpsContinuationSuperacc: public task {
    public :
    Superaccumulator* sum;
    Superaccumulator* mps;
    Superaccumulator* rsum;
    Superaccumulator* rmps;
    int* position;
    int rposition;

    MpsContinuationSuperacc(Superaccumulator* sum_, Superaccumulator* mps_, int* pos_, int middle):
        sum(sum_),
        mps(mps_),
        rsum(),
        rmps(),
        position(pos_),
        rposition(middle)
    {}
    
    task* execute(){
        rmps.Accumulate(*sum);
        sum->Accumulate(*sum);
        if(!rmps.comp(*mps)){
            mps = rmps;
            *position = rposition;
        }
        return NULL;
    }
};

// Homemade reduction operation 
class MpsTaskSuperacc: public task {
    public:
        MpsTaskSuperacc(int C, double* a,int size_, double*s, double*m,int*p, int left_ = 0, int right_ = -1) : 
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
        ~MpsTaskSuperacc(){}
        
        task* execute(){
            if(size <= Cutoff){

                for(int i = left; i != right; i++){
                    *sum += array[i];
                    if(*sum >= *mps){
                       *mps = *sum;
                       *position = i+1;
                    }
                }
                return NULL;

            }else{
                // Parameters for subtasks
                int middle = (right+left)/2;

                // Variables for results
                int sizel = left - middle;
                int sizer = middle - right;

                // Create subtasks
                MpsContinuationSuperacc& c = *new(allocate_continuation()) MpsContinuationSuperacc(sum, mps, position, middle);

                //MpsTaskSuperacc& lTask = *new(c.allocate_child()) MpsTaskSuperacc(Cutoff,array,sizel,sum,mps,position,left,middle);
                MpsTaskSuperacc &rTask = *new(c.allocate_child()) MpsTaskSuperacc(Cutoff,array,sizer,&c.rsum, &c.rmps,&c.rposition,middle,right);

                recycle_as_child_of(c);
                right = middle;
                size = sizel;

                c.set_ref_count(2);
                spawn(rTask);

                return this;
            }
        }

    private:
        // Below this size, the mps and sum are computed sequentially
        int Cutoff;
        // Input array and its size
        double* array;
        int size;
        // Intervals for sum and mps
        Superaccumulator* sum;
        Superaccumulator* mps;
        // Position of the mps
        int* position;
        // Bounds 
        int left;
        int right;
};
*/
// Parallel reduce implementation 
struct __mps_acc{
    // pointer to the array
    double* array;
    // Superaccumulators for sum and mps
    Superaccumulator sum;
    Superaccumulator mps;
    // Position of the maximum prefix sum
    int position;
    // Constructor
    __mps_acc(double* a);
    // Splitting constructor
    __mps_acc(__mps_acc&,split);
    // Accumulate result for subrange
    void operator()(const blocked_range<int>&);
    // Join operation for the reduction
    void join(__mps_acc& rightMps); 
    // Printing function
    void print_mps();
};

__mps_acc::__mps_acc(double* a):
    array(a),
    position(0)
{
    sum = Superaccumulator();
    mps = Superaccumulator();
}

__mps_acc::__mps_acc(__mps_acc& x, split) :
    array(x.array),
    position(0)
{
    sum = Superaccumulator();
    mps = Superaccumulator();
}

void __mps_acc::print_mps(){
    cout << "sum: " << sum.Round();
    cout << endl << "mps: " << mps.Round();
    cout << endl << "position: " << position << endl;
}

void __mps_acc::operator()(const blocked_range<int>& r){
    if(position == 0){
        position = r.begin();
    }
    for(int i = r.begin(); i != r.end(); i++){
        sum.Accumulate(array[i]);
        if(!sum.comp(mps)){
            mps = Superaccumulator(sum.get_accumulator());
            position = i+1; 
        }
    }
}

void __mps_acc::join(__mps_acc& rightMps){
   //Set rounding mode
   _MM_SET_ROUNDING_MODE(0);
   // computing sum-l + mps-r
   rightMps.mps.Accumulate(sum);
   // adding two sums
   sum.Accumulate(rightMps.sum);
   // comparison of mpsCandidate and mps-l
   double candidate = rightMps.mps.Round();
   double mpsAux = mps.Round();
   if(candidate >= mpsAux){
       mps = rightMps.mps;
       position = rightMps.position;
   }
}

// Parallel reduce main function
void tbb_main_s(double* a, int size, int grainsize){
    __mps_acc result = __mps_acc(a);
    parallel_reduce(blocked_range<int>(0,size,grainsize),result);
    if(DEBUGSUPERACC){
        cout << endl << "Parallel Reduce" << endl;
        result.print_mps();
    }
}

// Parallel deterministic reduce main function
void tbb_deterministic_main_s(double* a, int size, int grainsize){
    __mps_acc result = __mps_acc(a);
    parallel_deterministic_reduce(blocked_range<int>(0,size,grainsize),result);
    if(DEBUGSUPERACC){
        cout << endl << "Parallel Deterministic Reduce" << endl;
        result.print_mps();
    }
}

/*
// Homemade main function
void homemade_main_s(double *a, int size,int grainsize){
    
    double sum = 0, mps = 0;
    int pos = 0;

    MpsTaskSuperacc& root = *new(task::allocate_root()) MpsTaskSuperacc(grainsize,a,size,&sum,&mps,&pos);

    task::spawn_root_and_wait(root);

    if(DEBUGSUPERACC){
        cout << endl <<"Homemade reduction" << endl;
        cout << "Sum: " << sum << endl << "Mps " << mps << endl << "Pos: " << pos << endl;
    }
}
*/

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
void runtime_comparison_s(){
    
    // Variables declaration and initialisation 
    double start;
    int size = pow(10,7);
    int N = 1;

    // for each dynamic range
    vector<int> grainsSizes  {100,300,600,1000,2000,3000,6000,10000,20000,30000};
    int s = grainsSizes.size();
    
    // Store results to plot
    fstream results;
    results.open("Plots/reductions_superacc.csv", ofstream::out | ofstream::trunc);
    vector<double> x(s),r0(s),r1(s),r2(s);

    // Random seed
    srand(time(NULL));
    
    // Warm-up
    for(int r = 0; r < s; r++){
        // Generating array
        double* drray = new double[size];
        init_fpuniform(size, drray, 100, 50);

        // Randomly change signs
        for(int j = 0; j < size ; j++){
             drray[j] = (rand() % 2) ? drray[j] : -drray[j];
        }
        
        // Declare result variables
        double time_tbb = 0.0;
        PFP_TIME(tbb_main_s(drray,size,grainsSizes[r]),start,time_tbb);
        double time_homemade = 0.0;
        //PFP_TIME(homemade_main(drray,size,grainsSizes[r]),start,time_homemade);
        double time_deterministic_tbb = 0.0;
        PFP_TIME(tbb_deterministic_main_s(drray,size,grainsSizes[r]),start,time_deterministic_tbb);
   
        
        delete[] drray;
    }
    
    for(int r = 0; r < s; r++){

        // initialization of means
        double mean_tbb = 0., mean_homemade = 0., mean_deterministic_tbb = 0.;

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
            PFP_TIME(tbb_main_s(drray,size,grainsSizes[r]),start,time_tbb);
            double time_homemade = 0.0;
            //PFP_TIME(homemade_main_s(drray,size,grainsSizes[r]),start,time_homemade);
            double time_deterministic_tbb = 0.0;
            PFP_TIME(tbb_deterministic_main_s(drray,size,grainsSizes[r]),start,time_deterministic_tbb);
       
            mean_tbb += time_tbb;
            mean_homemade += time_homemade;
            mean_deterministic_tbb += time_deterministic_tbb;
            
            delete[] drray;
        }
        // Finalize mean computation
        mean_tbb += mean_tbb / N;
        mean_homemade += mean_homemade /N;
        mean_deterministic_tbb += mean_deterministic_tbb / N;
        
        x[r]= grainsSizes[r];
        r0[r]= mean_tbb;
        r1[r]= mean_homemade/ mean_tbb;
        r2[r]= mean_deterministic_tbb;

        // Writing results to a file
        results << to_string(x[r]) << "," << to_string(r0[r]) << "," << to_string(r1[r]) << "," << to_string(r2[r]) << endl;
        
    }

    results.close();
}
/*
int main(){
    runtime_comparison_s();
    return 0;
}*/
