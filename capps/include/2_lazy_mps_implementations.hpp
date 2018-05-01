/* Lazy computation of mps.
 * Author: Raphael Dang-Nhu
 * Date: 04/24/2018 */

#ifndef PAR_LAZY_MPS_H
#define PAR_LAZY_MPS_H

#include <mpfr.h>
#include <emmintrin.h>

#include "interval_arithmetic.hpp"
#include "superaccumulator.hpp"

#include "tbb/task_group.h"

using namespace tbb;

// Struct to store computation preliminary results
enum Status {
    leftChild,
    rightChild,
    undefinedComparison,
    cutoff,
    cutoffPrecise
};

class MpsTask1: public task {
    public:
        MpsTask1(int C, double* a, int s, __m128d* sum_i, __m128d* mps_i, int* p,Status** m, int d = 0, int i = 0, int l = 0, int r = -1) : 
            Cutoff(C),
            array(a),
            size(s),
            depth(d),
            index(i),
            left(l),
            right(r),
            sum_interval(sum_i),
            mps_interval(mps_i),
            position(p),
            memo(m)
        {
            if(r == -1) right = size;
            //cout << endl << "Constructed MpsTask." << endl;
            //cout << "left: " << left << endl;
            //cout << "right: " << right << endl;
        }
        ~MpsTask1(){}
        
        task* execute(){
            if(size <= Cutoff){
                Status s = cutoff;
                __m128d delta_sum = in2_create(0.,0.);

                _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);

                for(int i = left; i != right; i++){
                    *sum_interval = in2_add_double(*sum_interval,array[i]);
                    delta_sum = in2_add_double(delta_sum,array[i]);
                    boolean b = inferior(*mps_interval,*sum_interval);
                    //boolean b = inferior_double(0,delta_sum);
                    if (b == True){
                        *mps_interval = in2_add(*mps_interval,delta_sum);
                        delta_sum = in2_create(0.,0.);
                        *position = i+1; 
                    }
                    else if(b == Undefined){
                        s = cutoffPrecise;
                    }
                }

                memo[depth][index] = s;
                
            }
            else{
                // Parameters for subtasks
                int middle = (right+left)/2;
                int sizel = middle - left;
                int sizer = right - middle;
                int newDepth = depth + 1;
                int lIndex = 2*index;
                int rIndex = lIndex + 1;

                // Variables for results
                int lPos = left;
                int rPos = middle;
                __m128d lsum = in2_create(0.,0.);
                __m128d rsum = in2_create(0.,0.);
                __m128d lmps = in2_create(0.,0.);
                __m128d rmps = in2_create(0.,0.);

                // Create subtasks
                set_ref_count(3);

                MpsTask1& lTask = *new(allocate_child()) MpsTask1(Cutoff,array,sizel,&lsum,&lmps,&lPos,memo,newDepth,lIndex,left,middle);
                
                spawn(lTask);

                MpsTask1 &rTask = *new(allocate_child()) MpsTask1(Cutoff,array,sizer,&rsum, &rmps,&rPos,memo,newDepth,rIndex,middle,right);
                
                spawn_and_wait_for_all(rTask);
                
                // Gather results
                _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
                rmps = in2_add(lsum,rmps);
                *sum_interval = in2_add(lsum,rsum);
                boolean b = inferior(*mps_interval,rmps);
                if(b == True){
                    *mps_interval = rmps;
                    *position = rPos;
                    memo[depth][index] = rightChild;
                }
                else if(b == False){
                    *mps_interval = lmps;
                    *position = lPos;
                    memo[depth][index] = leftChild;
                }
                else{
                    // Merge the two intervals
                    *mps_interval = in2_merge(lmps,rmps);
                    memo[depth][index] = undefinedComparison;
                }
                // if undefined comparison, call the precise method
                 
                // Gather results and continue
                return NULL;
            }
        }

    private:
        /* Below this size, the mps and sum are computed sequentially */
        int Cutoff;
        /* Input array and its size */
        double* array;
        int size;    
        /* Identification of the task */
        int depth;
        int index;
        int left;
        int right;
        /* Intervals for sum and mps */
        __m128d* sum_interval;
        __m128d* mps_interval;
        /* Position of the mps */
        int* position;
        Status** memo;
};


// Main function to perform the lazy reduction
void parallel_mps_mpfr_lazy_2(double*,int,int);         

#endif

