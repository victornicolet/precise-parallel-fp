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

#define DEBUG 1

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
        MpsTask1(int C, double* a, int s,int* v, __m128d* sum_i, __m128d* mps_i, int* p,Status** m, int d = 0, int i = 0, int l = 0, int r = -1) : 
            Cutoff(C),
            array(a),
            size(s),
            validity(v),
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
            /*cout << endl << "Constructed MpsTask." << endl;
            cout << "left: " << left << endl;
            cout << "right: " << right << endl;
            */
        }
        ~MpsTask1(){}
        
        task* execute(){
            if(size <= Cutoff){
                Status stat = cutoff;
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
                        *mps_interval = in2_add(*mps_interval,delta_sum);
                        delta_sum = in2_create(0.,0.);
                        *position = i+1; 
                        stat = cutoffPrecise;
                        *validity = 1;
                    }
                }

                memo[depth][index] = stat;
                return NULL;
                
            }else{
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
                int vall = 0;
                int valr = 0;

                // Create subtasks
                set_ref_count(3);

                MpsTask1& lTask = *new(allocate_child()) MpsTask1(Cutoff,array,sizel,&vall,&lsum,&lmps,&lPos,memo,newDepth,lIndex,left,middle);
                
                spawn(lTask);

                MpsTask1 &rTask = *new(allocate_child()) MpsTask1(Cutoff,array,sizer,&valr,&rsum, &rmps,&rPos,memo,newDepth,rIndex,middle,right);
                
                spawn(rTask);
                wait_for_all();
                
                
                // Gather results
                _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
                rmps = in2_add(lsum,rmps);
                *sum_interval = in2_add(lsum,rsum);
                boolean b = inferior(lmps,rmps);
                if(b == True){
                    *mps_interval = rmps;
                    *position = rPos;
                    memo[depth][index] = rightChild;
                    *validity = valr;
                }
                else if(b == False){
                    *mps_interval = lmps;
                    *position = lPos;
                    memo[depth][index] = leftChild;
                    *validity = vall;
                }
                else{
                    // Merge the two intervals
                    *mps_interval = in2_merge(lmps,rmps);
                    *position = rPos;
                    *validity = 1;
                    memo[depth][index] = undefinedComparison;
                }
                return NULL;
            }
        }

    private:
        /* Below this size, the mps and sum are computed sequentially */
        int Cutoff;
        /* Input array and its size */
        double* array;
        int size;    
        /* Flag for validity of the result */
        int* validity;
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

// Task when asked for sum only
class MpsTask3: public task {
    public:
        MpsTask3(int C, double* a, int s, Superaccumulator* sum_,  int l = 0, int r = -1) : 
            Cutoff(C),
            array(a),
            size(s),
            left(l),
            right(r),
            sum(sum_)
        {
            if(r == -1) right = size;
            /*cout << endl << "Constructed MpsTask3." << endl;
            cout << "left: " << left << endl;
            cout << "right: " << right << endl;
            */
            
        }
        ~MpsTask3(){}
        
        task* execute(){
            if(size <= Cutoff){
                for(int i = left; i != right; i++){
                    sum->Accumulate(array[i]);
                }
              /*  cout << endl << "Chunk";
                cout << endl << "Left: " << left;
                cout << endl << "Right:" << right;
                cout << endl <<"Sum: ";
                mpfr_out_str(stdout,10,10,*sum,MPFR_RNDN);
                cout << endl << "Mps: ";
                mpfr_out_str(stdout,10,10,*mps,MPFR_RNDN);
                cout << endl;
                */

                
            }else{
                // Parameters for subtasks
                int middle = (right+left)/2;
                int sizel = middle - left;
                int sizer = right - middle;

                // Variables for results
                Superaccumulator lsum = Superaccumulator();
                Superaccumulator rsum = Superaccumulator();

                set_ref_count(3);
                MpsTask3& lTask = *new(allocate_child()) MpsTask3(Cutoff,array,sizel,&lsum,left,middle);
            
                spawn(lTask);
                
                // To be replaced by sum task
                MpsTask3 &rTask = *new(allocate_child()) MpsTask3(Cutoff,array,sizer,&rsum,middle,right);
            
                spawn_and_wait_for_all(rTask);
                
                // Gather results
                sum->Accumulate(lsum);
                sum->Accumulate(rsum);
                                
            }
            return NULL;
        }

    private:
        /* Below this size, the mps and sum are computed sequentially */
        int Cutoff;
        /* Input array and its size */
        double* array;
        int size;    
        /* Identification of the task */
        int left;
        int right;
        /* Intervals for sum and mps */
        Superaccumulator* sum;
};
#endif


// Task when asked for mps
class MpsTask2: public task {
    public:
        MpsTask2(int C, double* a, int s, Superaccumulator* sum_, Superaccumulator* mps_, int* p,Status** m, int d = 0, int i = 0, int l = 0, int r = -1) : 
            Cutoff(C),
            array(a),
            size(s),
            depth(d),
            index(i),
            left(l),
            right(r),
            sum(sum_),
            mps(mps_),
            position(p),
            memo(m)
        {
            if(r == -1) right = size;
            /*cout << endl << "Constructed MpsTask2." << endl;
            cout << "left: " << left << endl;
            cout << "right: " << right << endl;
            */
            
        }
        ~MpsTask2(){}
        
        task* execute(){
            if(size <= Cutoff){
                for(int i = left; i != right; i++){
                    sum->Accumulate(array[i]);
                    if (!sum->comp(*mps)){
                        *mps = Superaccumulator(sum->get_accumulator());
                        *position = i+1; 
                    }
                }
              /*  cout << endl << "Chunk";
                cout << endl << "Left: " << left;
                cout << endl << "Right:" << right;
                cout << endl <<"Sum: ";
                mpfr_out_str(stdout,10,10,*sum,MPFR_RNDN);
                cout << endl << "Mps: ";
                mpfr_out_str(stdout,10,10,*mps,MPFR_RNDN);
                cout << endl;
                */

                
            }else{
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
                // Check status
                Status stat = memo[depth][index];

                if(stat ==leftChild){
                    
                    Superaccumulator rsum = Superaccumulator();
                    Superaccumulator lsum = Superaccumulator();
                    Superaccumulator lmps = Superaccumulator();
                    
                    set_ref_count(3);
                    MpsTask2& lTask = *new(allocate_child()) MpsTask2(Cutoff,array,sizel,&lsum,&lmps,&lPos,memo,newDepth,lIndex,left,middle);
                
                    spawn(lTask);
                    
                    MpsTask3 &rTask = *new(allocate_child()) MpsTask3(Cutoff,array,sizer,&rsum, middle,right);
                
                    spawn_and_wait_for_all(rTask);
                    
                    // Gather results
                    sum->Accumulate(lsum);
                    sum->Accumulate(rsum);
                    *mps = Superaccumulator(lmps.get_accumulator());
                    *position = lPos;
                } 
                else if(stat == rightChild){
                    Superaccumulator rsum = Superaccumulator();
                    Superaccumulator lsum = Superaccumulator();
                    Superaccumulator rmps = Superaccumulator();

                    set_ref_count(3);
                    MpsTask3& lTask = *new(allocate_child()) MpsTask3(Cutoff,array,sizel,&lsum,left,middle);
                
                    spawn(lTask);

                    MpsTask2 &rTask = *new(allocate_child()) MpsTask2(Cutoff,array,sizer,&rsum, &rmps,&rPos,memo,newDepth,rIndex,middle,right);
                
                    spawn_and_wait_for_all(rTask);
                    
                    // Gather results
                    sum->Accumulate(lsum);
                    sum->Accumulate(rsum);
                    *mps = Superaccumulator(rmps.get_accumulator());
                    *position = rPos;
                }
                else if(stat == undefinedComparison){
                    Superaccumulator rsum = Superaccumulator();
                    Superaccumulator lsum = Superaccumulator();
                    Superaccumulator rmps = Superaccumulator();
                    Superaccumulator lmps = Superaccumulator();

                    set_ref_count(3);
                    MpsTask2& lTask = *new(allocate_child()) MpsTask2(Cutoff,array,sizel,&lsum,&lmps,&lPos,memo,newDepth,lIndex,left,middle);
                
                    spawn(lTask);

                    MpsTask2 &rTask = *new(allocate_child()) MpsTask2(Cutoff,array,sizer,&rsum, &rmps,&rPos,memo,newDepth,rIndex,middle,right);
                
                    spawn_and_wait_for_all(rTask);
                    
                    // Gather results
                    sum->Accumulate(lsum);
                    sum->Accumulate(rsum);
                    rmps.Accumulate(lsum);

                    if(!rmps.comp(*mps)){
                        *mps = Superaccumulator(rmps.get_accumulator());
                        *position = rPos;
                    }
                    else {
                        *mps = Superaccumulator(lmps.get_accumulator());
                        *position = lPos;
                    }
                }
                /*cout << endl << "Join";
                cout << endl << "Left: " << left;
                cout << endl << "Right:" << right;
                cout << endl <<"Sum: ";
                mpfr_out_str(stdout,10,10,*sum,MPFR_RNDN);
                cout << endl << "Mps: ";
                mpfr_out_str(stdout,10,10,*mps,MPFR_RNDN);
                cout << endl;
                */
                                
            }
            return NULL;
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
        Superaccumulator* sum;
        Superaccumulator* mps;
        /* Position of the mps */
        int* position;
        Status** memo;
};
