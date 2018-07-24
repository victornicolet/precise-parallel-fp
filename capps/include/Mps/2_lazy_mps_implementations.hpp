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

// Mps result
struct mps_struct{
    mps_struct(){};
    __m128d sum;
    __m128d mps;
    int pos;
    bool val;
};

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
        MpsTask1(int C, double* a, int s, mps_struct* result_,Status** m, int d = 0, int i = 0, int l = 0, int r = -1) : 
            Cutoff(C),
            array(a),
            size(s),
            depth(d),
            index(i),
            left(l),
            right(r),
            result(result_),
            memo(m)
        {
            //cout << endl << "Constructed MpsTask." << endl;
            //cout << "left: " << left << endl;
            //cout << "right: " << right << endl;
            
        }
        ~MpsTask1(){}
        
        task* execute(){
            if(size <= Cutoff){
                Status stat = cutoff;
                _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);

                __m128d s = in2_create(0.,0.);  
                __m128d m = in2_create(0.,0.);  
                int v = 0;
                int p = left;

                for(int i = left; i != right; i++){
                    s = in2_add(s,array[i]);
                    boolean b = in2_le(m,s);
                    if (b == True){
                        m = s;
                        p = i+1; 
                        v = 0;
                    }
                    else if(b == Undefined){
                        m = in2_merge(m,s);
                        p = i+1;
                        v = 1;
                    }
                }

                result->sum = s;
                result->mps = m;
                result->pos = p;
                result->val = v;

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
                mps_struct rresult;

                // Create subtasks
                set_ref_count(3);

                MpsTask1& lTask = *new(allocate_child()) MpsTask1(Cutoff,array,sizel,result,memo,newDepth,lIndex,left,middle);
                
                spawn(lTask);

                MpsTask1 &rTask = *new(allocate_child()) MpsTask1(Cutoff,array,sizer,&rresult,memo,newDepth,rIndex,middle,right);
                
                spawn(rTask);
                wait_for_all();
                
                
                // Gather results
                _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);

                rresult.mps = in2_add(result->sum,rresult.mps);
                result->sum = in2_add(result->sum,rresult.sum);
                boolean b = in2_le(result->mps,rresult.mps);
                if(b == True){
                    result->mps = rresult.mps;
                    result->pos = rresult.pos;
                    memo[depth][index] = rightChild;
                    result->val = rresult.val;
                }
                else if(b == Undefined){
                    // Merge the two intervals
                    result->mps = in2_merge(result->mps,rresult.mps);
                    result->pos = rresult.pos;
                    result->val = 1;
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
        /* Identification of the task */
        int depth;
        int index;
        int left;
        int right;
        mps_struct* result;
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
