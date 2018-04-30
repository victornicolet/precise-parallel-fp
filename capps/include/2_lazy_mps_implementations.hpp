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

class MpsTask1: public task {
    public:
        MpsTask1(int Cutoff, double* array, int size, __m128d* sum_interval, __m128d* mps_interval, int* position, int depth = 0, int index = 0, int left = 0, int right = -1) : 
            Cutoff(Cutoff),
            array(array),
            size(size),
            depth(depth),
            index(index),
            left(left),
            sum_interval(sum_interval),
            mps_interval(mps_interval),
            position(position)
        {
            if(right == -1) right = size;
        }
        ~MpsTask1(){}
        
        // Method to process the chunk if it is smaller than the Cutoff size 
        void processChunk(){
            for(int i = left; i != right; i++){
                *mps_interval = in2_add_double(*sum_interval,array[i]);
                boolean b = inferior(*mps_interval,*sum_interval);
                if (b == True){
                    mps_interval = sum_interval;
                    *position = i; 
                }
                else if(b == Undefined){
                    // Call precise computation of chunk
                    // Gather results and continue
                }
            }
        }

        task* execute(){
            if(size <= Cutoff){
                processChunk();
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
                MpsTask1& lTask = *new(allocate_child()) MpsTask1(Cutoff,array,sizel,&lsum,&lmps,&lPos,newDepth,lIndex,left,middle);

                MpsTask1 &rTask = *new(allocate_child()) MpsTask1(Cutoff,array,sizer,&rsum, &rmps,&rPos,newDepth,rIndex,middle,right);
                
                set_ref_count(3);
                spawn(lTask);
                spawn_and_wait_for_all(rTask);
                
                // Gather results
                rmps = in2_add(lsum,rmps);
                *sum_interval = in2_add(lsum,rsum);
                boolean b = inferior(*mps_interval,rmps);
                if(b == True){
                    *mps_interval = rmps;
                    *position = rPos;
                }
                else{
                    *mps_interval = lmps;
                    *position = lPos;
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
};

// Main function to perform the lazy reduction
void parallel_mps_mpfr_lazy_2(double*,int);         

#endif

