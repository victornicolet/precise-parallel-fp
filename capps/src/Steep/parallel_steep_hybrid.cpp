/* Lazy steep implementation with hybrid reduction operation */

#include <iostream>

#include "tbb/task_group.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_reduce.h"

#include "debug.hpp"
#include "parallel_steep.hpp"
#include "interval_arithmetic.hpp"

using namespace tbb;
using namespace std;


struct steep_struct {
    double sum;
    double capacity;
    boolean b;
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
                double aux = resRight.capacity + resLeft.sum;
                capacity = min(resRight.capacity,aux);
                res -> sum = resLeft.sum + resRight.sum;
                b &= aux >=0.;

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
        b(false),
    {}
    bool sum;
    bool capacity;
    bool sum;

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
        HybridSteep(int depth_,int index_, double* a_, steep_struct_interval* res_, unsigned long left_, unsigned long right_,boolean*** memo_) :
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
            __steep_interval_without_pos result = __steep_interval_without_pos(a);
            parallel_reduce(blocked_range<long>(left,right),result);
            res -> sum = result.sum;
            res -> capacity = result.capacity;
            res -> b = result.b;

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
                
                __m128d aux = in2_sub(resRight.capacity,sum);

                // First decision
                dec = in2_le(resLeft.capacity,aux);
                switch(dec){
                    case True:
                        capacity = resLeft.capacity;
                        break;
                    case False:
                        capacity = aux;
                        break;
                    case Undefined:
                        capacity = in2_merge(resLeft.capacity,aux); 
                }
                memo[depth][index][0] = dec;

                // Second decision
                dec = resLeft.b;
                memo[depth][index][1] = dec;
                switch(dec){
                    case False:
                        return False;
                    case Undefined:
                        // If b false not undefined
                        // Third decision
                        dec = in2_ge(aux,0.);
                        memo[depth][index][2] = dec;
                        switch(dec){
                            case False:
                                b = False;
                            break;
                            default:
                                b = Undefined;
                        }
                        break;
                    case True:
                        b = True;
                }

                sum = in2_add(resLeft.sum,resRight.sum);

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
                cout << endl;
                cout << "Computation started " << depth << endl;
                cout << "Index: " << index << endl;

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

            boolean b = memo[depth][index][3];
            switch(b){
                case True:
                
                if(res_bool->steep){
                    res_bool->steep = false;
                    resBoolRight.mps = true;
                    resBoolLeft.mts = true;
                }
                if(res_bool->posl){
                    res_bool->posl = false;
                    resBoolLeft.posMps = true;
                }
                if(res_bool->posr){
                    res_bool->posr = false;
                    resBoolRight.posMts = true;
                }

                break;
                case Undefined:
                // Variables used in the conditionnal
                t = res_bool->steep || res_bool->posl || res_bool->posr;

                // Assignments
                if(res_bool->steep){
                    resBoolRight.mps = true;
                    resBoolLeft.mts = true;
                }
                if(res_bool->posl){
                    resBoolLeft.posMps = true;
                }
                if(res_bool->posr){
                    resBoolRight.posMts = true;
                }
                if(t){
                    res_bool->steep = true;
                    resBoolLeft.mts = true;
                    resBoolRight.mps = true;
                }
                default:
                break;
            }

            b = memo[depth][index][2];
            switch(b){
                case True:
                if(res_bool->steep){
                    res_bool->steep = false;
                    resBoolLeft.steep = true;
                }
                if(res_bool->posl){
                    res_bool->posl = false;
                    resBoolLeft.posl = true;
                }
                if(res_bool->posr){
                    res_bool->posr = false;
                    resBoolLeft.posr = true;
                }
                break;
                case False:
                if(res_bool->steep){
                    res_bool->steep = false;
                    resBoolRight.steep = true;
                }
                if(res_bool->posl){
                    res_bool->posl = false;
                    resBoolRight.posl = true;
                }
                if(res_bool->posr){
                    res_bool->posr = false;
                    resBoolRight.posr = true;
                }
                break;
                default:
                // Variables used in the conditionnal
                t = res_bool->steep || res_bool->posl || res_bool->posr;

                if(res_bool->steep){
                    res_bool->steep = false;
                    resBoolRight.steep = true;
                    resBoolLeft.steep = true;
                }
                if(res_bool->posl){
                    res_bool->posl = false;
                    resBoolRight.posl = true;
                    resBoolLeft.posl = true;
                }
                if(res_bool->posr){
                    res_bool->posr = false;
                    resBoolRight.posr = true;
                    resBoolLeft.posr = true;
                }

                if(t){
                    resBoolLeft.steep = true;
                    resBoolRight.steep = true;
                }
            }

            // Mts
            b = memo[depth][index][1];
            switch(b){
                case True:
                if(res_bool->mts){
                    res_bool->mts = false;
                    resBoolRight.sum = true;
                    resBoolLeft.mts = true;
                }
                if(res_bool->posMts){
                    res_bool->posMts = false;
                    resBoolLeft.posMts = true;
                }
                break;
                case False:
                if(res_bool->mts){
                    res_bool->mts = false;
                    resBoolRight.mts = true;
                }
                if(res_bool->posMts){
                    res_bool->posMts = false;
                    resBoolRight.posMts = true;
                }
                break;
                default:
                t = res_bool->mts || res_bool->posMts;

                if(res_bool->mts){
                    res_bool->mts = false;
                    resBoolRight.sum = true;
                    resBoolLeft.mts = true;
                    resBoolRight.mts = true;
                }
                if(res_bool->posMts){
                    res_bool->posMts = false;
                    resBoolRight.posMts = true;
                    resBoolLeft.posMts = true;
                }

                if(t){
                    resBoolRight.sum = true;
                    resBoolLeft.mts = true;
                    resBoolRight.mts = true;
                }
            }

            // Mps
            b = memo[depth][index][0];
            switch(b){
                case True:
                if(res_bool->mps){
                    res_bool->mps = false;
                    resBoolLeft.sum = true;
                    resBoolRight.mps = true;
                }
                if(res_bool->posMps){
                    res_bool->posMps = false;
                    resBoolRight.posMps = true;
                }
                break;
                case False:
                if(res_bool->mps){
                    res_bool->mps = false;
                    resBoolLeft.mps = true;
                }
                if(res_bool->posMps){
                    res_bool->posMps = false;
                    resBoolLeft.posMps = true;
                }
                break;
                default:
                t = res_bool->mps || res_bool->posMps;

                if(res_bool->mps){
                    res_bool->mps = false;
                    resBoolLeft.mps = true;
                    resBoolLeft.sum = true;
                    resBoolRight.mps = true;
                }
                if(res_bool->posMps){
                    res_bool->posMps = false;
                    resBoolLeft.posMps = true;
                    resBoolRight.posMps = true;
                }

                if(t){
                    resBoolLeft.mps = true;
                    resBoolLeft.sum = true;
                    resBoolRight.mps = true;
                }

            }

            // Sum
            if(res->sum){
                res->sum = false;
                resBoolLeft.sum = true;
                resBoolRight.sum = true;
            }

            // Call subtasks and return result;
            HybridSteepReductionMpfr& lTask = *new(allocate_child()) HybridSteepReductionMpfr(newDepth,lIndex,a,&resLeft,&resBoolLeft,left,middle,memo);
            HybridSteepReductionMpfr &rTask = *new(allocate_child()) HybridSteepReductionMpfr(newDepth,rIndex,a,&resRight,&resBoolRight,middle,right,memo);

            set_ref_count(3);
            spawn(lTask);
            spawn_and_wait_for_all(rTask);

            // Join
            // Sum
            res->sum = resLeft.sum + resRight.sum;

            // Mps
            mpreal mpsAux = resLeft.sum+resRight.mps;
            if(mpsAux >= resLeft.mps){
                res -> mps = mpsAux;
                res -> posMps = resRight.posMps;
            }else{
                res -> mps = resLeft.mps;
                res -> posMps = resLeft.posMps;
            }

            // Mts
            mpreal mtsAux = resRight.sum+resLeft.mts;
            if(mtsAux>=resRight.mts){
                res->mts = mtsAux;
                res->posMts = resLeft.posMts;
            }else{
                res->mts = resRight.mts;
                res->posMts = resRight.posMts;
            }

            // Mss
            if(resLeft.steep>=resRight.steep){
                res->steep = resLeft.steep; 
                res->posl = resLeft.posl;
                res->posr = resLeft.posr;
            }else{
                res->steep = resRight.steep; 
                res->posl = resRight.posl;
                res->posr = resRight.posr;
            }

            mpreal steepAux = resRight.mps+resLeft.mts;
            if(steepAux>=res->steep){
                res->steep = steepAux;
                res->posl = resLeft.posMts;
                res->posr = resRight.posMps;
            }

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

void parallel_steep_hybrid_interval(double* a, long size,int maxDepth){
    
    // Set rounding mode
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    task_scheduler_init init;
    
    int maxIndex = 1;
    boolean*** memo = new boolean**[maxDepth]; 
    for(int i = maxDepth-1; i >= 0;i--){
        memo[i] = new boolean*[maxIndex];
        for(int j = 0; j != maxIndex; j++){
            memo[i][j] = new boolean[4];
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
        cout << "Steep: ";
        print(res.steep);
        cout << endl;
        cout << "Mps: ";
        print(res.mps);
        cout << endl;
        cout << "Mts: ";
        print(res.mts);
        cout << endl;
    }
    init.terminate();
    _MM_SET_ROUNDING_MODE(0);

    // Second step
    task_scheduler_init init2;

    steep_struct_mpfr res_mpfr; 
    steep_struct_bool res_bool; 
    res_bool.posl = true;
    res_bool.posr = true;

    HybridSteepReductionMpfr& rootMpfr = *new(task::allocate_root()) HybridSteepReductionMpfr(maxDepth-1,0,a,&res_mpfr,&res_bool,0,size,memo);

    task::spawn_root_and_wait(rootMpfr);
    if(PRINT){
        cout << endl << "Hybrid steep Mpfr" << endl;
        cout << "Mss Left pos: " << res_mpfr.posl << endl;
        cout << "Mss Right pos: " << res_mpfr.posr << endl;
        cout << "Mps Pos: " << res_mpfr.posMps << endl;
        cout << "Mts Pos: " << res_mpfr.posMts << endl;
    }


}
