/* Lazy steep implementation with hybrid reduction operation */

#include <iostream>

#include "tbb/task_group.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_reduce.h"

#include "pfpdefs.hpp"
#include "debug.hpp"
#include "parallel_steep.hpp"
#include "interval_arithmetic.hpp"

using namespace tbb;
using namespace std;


struct steep_struct {
    double sum;
    double capacity;
    bool b;
};

class HybridSteepReduction : public task {
    public :
        HybridSteepReduction(int depth_,int index_, double* a_, steep_struct* res_, unsigned long left_, unsigned long right_) :
            depth(depth_),
            index(index_),
            a(a_),
            res(res_),
            left(left_),
            right(right_)
    {}
    ~HybridSteepReduction(){}

    task* execute(){
        if(depth == 0){
            // Call parallel_reduce
            __steep_naive result = __steep_naive(a);
            parallel_reduce(blocked_range<long>(left,right),result);
            res -> sum = result.sum;
            res -> capacity = result.capacity;
            res -> b = result.b;

        }
        else{
            int newDepth = depth - 1;
            steep_struct resLeft;
            steep_struct resRight;
            long middle = (left + right) / 2;
            int lIndex = 2*index;
            int rIndex = lIndex + 1;

            // Call subtasks and return result;
                HybridSteepReduction& lTask = *new(allocate_child()) HybridSteepReduction(newDepth,lIndex,a,&resLeft,left,middle);
                HybridSteepReduction &rTask = *new(allocate_child()) HybridSteepReduction(newDepth,rIndex,a,&resRight,middle,right);

                set_ref_count(3);
                spawn(lTask);
                spawn_and_wait_for_all(rTask);

                // Join
                double aux = resRight.capacity - resLeft.sum;
                res -> capacity = min(resLeft.capacity,aux);
                res -> sum = resLeft.sum + resRight.sum;
                res -> b &= aux >=0.;

        }
        return NULL;
    }
    private:
        int depth;
        int index;
        double* a;
        steep_struct* res;
        unsigned long left;
        unsigned long right;
        boolean** memo;
};



void parallel_steep_hybrid(double* a, long size, int maxDepth){
    task_scheduler_init init;
    
    steep_struct res; 

    HybridSteepReduction& root = *new(task::allocate_root()) HybridSteepReduction(maxDepth-1,0,a,&res,0,size);

    task::spawn_root_and_wait(root);
    if(PRINT){
        cout << endl << "Hybrid steep" << endl;
        cout << "Sum: " << res.sum << endl;
        cout << "Capacity: " << res.capacity << endl;
        cout << "Boolean: " << res.b << endl;
    }
    init.terminate();

}

struct steep_struct_bool {
    steep_struct_bool() :
        sum(false),
        capacity(false),
        b(false)
    {}
    bool sum;
    bool capacity;
    bool b;

    void print(){
        cout << "Sum: " << sum << endl;
        cout << "Capacity: " << capacity << endl;
        cout << "Boolean: " << b << endl;
    }

};

struct steep_struct_interval {
    __m128d sum;
    __m128d capacity;
    boolean b;
};

struct steep_struct_mpfr {
    mpreal sum;
    mpreal capacity;
    bool b;
};

class HybridSteepReductionInterval : public task {
    public :
        HybridSteepReductionInterval(int depth_,int index_, double* a_, steep_struct_interval* res_, unsigned long left_, unsigned long right_,boolean*** memo_) :
            depth(depth_),
            index(index_),
            a(a_),
            res(res_),
            left(left_),
            right(right_),
            memo(memo_)
    {}
    ~HybridSteepReductionInterval(){}

    task* execute(){
        if(depth == 0){
            // Call parallel_reduce
            __steep_interval result = __steep_interval(this->a);

            parallel_reduce(blocked_range<long>(left,right),result);
            res -> sum = result.sum;
            res -> capacity = result.capacity;
            res -> b = in2_ge(result.capacity,0.);

        }
        else{

            int newDepth = depth - 1;
            steep_struct_interval resLeft;
            steep_struct_interval resRight;
            long middle = (left + right) / 2;
            int lIndex = 2*index;
            int rIndex = lIndex + 1;

            // Call subtasks and return result;
                HybridSteepReductionInterval& lTask = *new(allocate_child()) HybridSteepReductionInterval(newDepth,lIndex,a,&resLeft,left,middle,memo);
                HybridSteepReductionInterval &rTask = *new(allocate_child()) HybridSteepReductionInterval(newDepth,rIndex,a,&resRight,middle,right,memo);

                set_ref_count(3);
                spawn(lTask);
                spawn_and_wait_for_all(rTask);

                // Join
                boolean dec;
                
                __m128d aux = in2_sub(resRight.capacity,resLeft.sum);

                // First decision
                dec = in2_le(resLeft.capacity,aux);
                switch(dec){
                    case True:
                        res -> capacity = resLeft.capacity;
                        break;
                    case False:
                        res -> capacity = aux;
                        break;
                    case Undefined:
                        res -> capacity = in2_merge(resLeft.capacity,aux); 
                }
                memo[depth][index][0] = dec;
    
                // Second decision
                dec = resLeft.b;
                memo[depth][index][1] = dec;
                switch(dec){
                    case False:
                        res -> b = False;
                        break;
                    case Undefined:
                        // If b false not undefined
                        dec = in2_ge(aux,0.);
                        memo[depth][index][2] = dec;
                        switch(dec){
                            case False:
                                res -> b = False;
                            break;
                            default:
                                res -> b = Undefined;
                        }
                        break;
                    case True:
                        res -> b = in2_ge(aux,0.);
                }

                res -> sum = in2_add(resLeft.sum,resRight.sum);

        }
        return NULL;
    }
    private:
        int depth;
        int index;
        double* a;
        steep_struct_interval* res;
        unsigned long left;
        unsigned long right;
        boolean*** memo;
};

class HybridSteepReductionMpfr : public task {
    public :
        HybridSteepReductionMpfr(int depth_,int index_, double* a_, steep_struct_mpfr* res_,steep_struct_bool* res_bool_, unsigned long left_, unsigned long right_,boolean*** memo_) :
            depth(depth_),
            index(index_),
            a(a_),
            res(res_),
            res_bool(res_bool_),
            left(left_),
            right(right_),
            memo(memo_)
    {}
    ~HybridSteepReductionMpfr(){}

    task* execute(){

        // Debug
        /*cout << endl;
        cout << "Depth: " << depth << endl;
        cout << "Index: " << index << endl;
        res_bool->print();
        */
        

        if(depth == 0){
            // Call parallel_reduce
            bool aux = res_bool->capacity || res_bool->sum || res_bool -> b;
            if(aux){

                // Debug
                /*cout << endl;
                cout << "Computation started " << depth << endl;
                cout << "Index: " << index << endl;
                cout << "capacity status: " << res_bool->capacity << endl;
                cout << "sum status: " << res_bool->sum << endl;
                cout << "bool status: " << res_bool->b << endl;
                */

                __steep_mpfr result = __steep_mpfr(a);
                parallel_reduce(blocked_range<long>(left,right),result);
                res -> sum = result.sum;
                res -> capacity = result.capacity;
                res -> b = result.b;
            }
        }
        else{
            int newDepth = depth - 1;
            steep_struct_mpfr resLeft;
            steep_struct_mpfr resRight;
            unsigned long middle = (left + right) / 2;
            int lIndex = 2*index;
            int rIndex = lIndex + 1;
            
            steep_struct_bool resBoolLeft, resBoolRight;
            bool t;
            boolean dec;
            bool aux;

            // Sum
            if(res_bool -> sum){
                //cout << "test" << endl;
                res_bool->sum = false;
                resBoolRight.sum = true;
                resBoolLeft.sum = true;
            }

            // Second decision
            // dec is b, if false we don't need to know aux, else we need to
            if(res_bool -> b){
                dec = memo[depth][index][1];
                switch(dec){
                    case False:
                        res_bool -> b = false;
                        break;
                    case Undefined:
                        // If b false not undefined
                        // Third decision
                        dec = memo[depth][index][2];
                        switch(dec){
                            case False:
                                res_bool -> b = false;
                            break;
                            default:
                                aux = true;
                        }
                        break;
                    case True:
                        res_bool -> b = false;
                        aux = true;
                }
            }

            // First decision 
            dec = memo[depth][index][0];
            switch(dec){
                case True:
                    if(res_bool -> capacity){
                        res_bool -> capacity = false;
                        resLeft.capacity = true;
                    }
                    break;
                case False:
                    if(res_bool -> capacity){
                        res_bool -> capacity = false;
                        aux = true;
                    }
                    break;
                case Undefined:
                    if(res_bool -> capacity){
                        res_bool -> capacity = false;
                        aux = true;
                        resLeft.capacity = true;
                    }
            }

            if(aux){
                resRight.capacity = true;
                resLeft.sum = true;
            }

            // Call subtasks and return result;
            HybridSteepReductionMpfr& lTask = *new(allocate_child()) HybridSteepReductionMpfr(newDepth,lIndex,a,&resLeft,&resBoolLeft,left,middle,memo);
            HybridSteepReductionMpfr &rTask = *new(allocate_child()) HybridSteepReductionMpfr(newDepth,rIndex,a,&resRight,&resBoolRight,middle,right,memo);

            set_ref_count(3);
            spawn(lTask);
            spawn_and_wait_for_all(rTask);

            // Join
            mpreal aux2 = resRight.capacity - resLeft.sum;
            res -> capacity = min(resLeft.capacity,aux2);
            res -> sum = resLeft.sum + resRight.sum;
            res -> b &= aux2 >=0.;

        }
        return NULL;
    }
    private:
        int depth;
        int index;
        double* a;
        steep_struct_mpfr* res;
        steep_struct_bool* res_bool;
        unsigned long left;
        unsigned long right;
        boolean*** memo;
};

void parallel_steep_hybrid_lazy(double* a, long size,int maxDepth,double& time_hybrid_interval, double& time_hybrid_total){

    double start;
    PFP_TIME(boolean*** memo = parallel_steep_hybrid_interval(a, size, maxDepth),start,time_hybrid_interval);

    double time_exact;
    PFP_TIME(parallel_steep_hybrid_exact(a, size, maxDepth,memo),start,time_exact);

    time_hybrid_total = time_hybrid_interval + time_exact;
}

boolean*** parallel_steep_hybrid_interval(double* a, long size,int maxDepth){
    
    // Set rounding mode
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    task_scheduler_init init;
    
    int maxIndex = 1;
    boolean*** memo = new boolean**[maxDepth]; 
    for(int i = maxDepth-1; i >= 0;i--){
        memo[i] = new boolean*[maxIndex];
        for(int j = 0; j != maxIndex; j++){
            memo[i][j] = new boolean[3];
        }
        maxIndex = 2*maxIndex;
    }

    steep_struct_interval res; 

    HybridSteepReductionInterval& root = *new(task::allocate_root()) HybridSteepReductionInterval(maxDepth-1,0,a,&res,0,size,memo);

    task::spawn_root_and_wait(root);
    if(PRINT){
        cout << endl << "Hybrid steep interval" << endl;
        cout << "Sum: " ;
        print(res.sum);
        cout << endl;
        cout << "Capacity: ";
        print(res.capacity);
        cout << endl;
        cout << "Boolean: ";
        print(res.b);
        cout << endl;
    }
    init.terminate();
    _MM_SET_ROUNDING_MODE(0);
    return memo;

}

void parallel_steep_hybrid_exact(double* a, long size,int maxDepth,boolean*** memo){
    
    // Second step
    task_scheduler_init init2(1);

    steep_struct_mpfr res_mpfr; 
    steep_struct_bool res_bool; 
    res_bool.b = true;

    HybridSteepReductionMpfr& rootMpfr = *new(task::allocate_root()) HybridSteepReductionMpfr(maxDepth-1,0,a,&res_mpfr,&res_bool,0,size,memo);

    task::spawn_root_and_wait(rootMpfr);
    if(PRINT){
        cout << endl << "Hybrid steep Mpfr" << endl;
        cout << "Sum: " << res_mpfr.sum << endl;
        cout << "Capacity: " << res_mpfr.capacity << endl;
        cout << "Boolean: " << res_mpfr.b << endl;
    }
}


