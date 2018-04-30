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
        MpsTask1(int Cutoff, double* array, int size, __m128d* sum_interval, __m128d* mps_interval, int* p, int depth = 0, int index = 0, int l = 0, int r = -1) : 
            Cutoff(Cutoff),
            array(array),
            size(size),
            depth(depth),
            index(index),
            left(l),
            right(r),
            sum_interval(sum_interval),
            mps_interval(mps_interval),
            position(p)
        {
            if(r == -1) right = size;
            //cout << endl << "Constructed MpsTask." << endl;
            //cout << "left: " << left << endl;
            //cout << "right: " << right << endl;
        }
        ~MpsTask1(){}
        
        // Method to process the chunk if it is smaller than the Cutoff size 
        void processChunk(){
            _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
            for(int i = left; i != right; i++){
                *sum_interval = in2_add_double(*sum_interval,array[i]);
                boolean b = inferior(*mps_interval,*sum_interval);
                if (b == True){
                    *mps_interval = *sum_interval;
                    *position = i+1; 
                }
                else if(b == Undefined){
                    // Call precise computation of chunk
                    // Gather results and continue
                }
            }
            
        }

        task* execute(){
            if(size <= Cutoff){
                // This is really weird
                cout << "0" << endl;
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
                _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
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

class MpsTask2: public task {
    public:
        MpsTask2(int Cutoff, double* array, int size, int* sum_interval, int* mps_interval, int* p, int depth = 0, int index = 0, int l = 0, int r = -1) : 
            Cutoff(Cutoff),
            array(array),
            size(size),
            depth(depth),
            index(index),
            left(l),
            right(r),
            sum_interval(sum_interval),
            mps_interval(mps_interval),
            position(p)
        {
            if(r == -1) right = size;
            //cout << endl << "Constructed MpsTask." << endl;
            //cout << "left: " << left << endl;
            //cout << "right: " << right << endl;
        }
        ~MpsTask2(){}
        
        // Method to process the chunk if it is smaller than the Cutoff size 
        void processChunk(){
            _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
            for(int i = left; i != right; i++){
                *sum_interval = in2_add_double(*sum_interval,array[i]);
                boolean b = inferior(*mps_interval,*sum_interval);
                if (b == True){
                    *mps_interval = *sum_interval;
                    *position = i+1; 
                }
                else if(b == Undefined){
                    // Call precise computation of chunk
                    // Gather results and continue
                }
            }
            
        }

        task* execute(){
            if(size <= Cutoff){
                // This is really weird
                cout << "0" << endl;
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
                MpsTask2& lTask = *new(allocate_child()) MpsTask2(Cutoff,array,sizel,&lsum,&lmps,&lPos,newDepth,lIndex,left,middle);

                MpsTask2 &rTask = *new(allocate_child()) MpsTask2(Cutoff,array,sizer,&rsum, &rmps,&rPos,newDepth,rIndex,middle,right);
                
                set_ref_count(3);
                spawn(lTask);
                spawn_and_wait_for_all(rTask);
                
                // Gather results
                _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
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

// Main function to perform the lazy reduction
void parallel_mps_mpfr_lazy_2(double*,int);         

#endif

