/* Lazy mss implementation with hybrid reduction operation */

#include <iostream>

#include "tbb/task_group.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_reduce.h"

#include "debug.hpp"
#include "parallel_mss.hpp"
#include "interval_arithmetic.hpp"

using namespace tbb;
using namespace std;

struct mss_struct {
    double sum;
    double mss;
    double mts;
    double mps;
    long posl;
    long posr;
    long posMps;
    long posMts;
};

class HybridMssReduction : public task {
    public :
        HybridMssReduction(int depth_,int index_, double* a_, mss_struct* res_, unsigned long left_, unsigned long right_) :
            depth(depth_),
            index(index_),
            a(a_),
            res(res_),
            left(left_),
            right(right_)
    {}
    ~HybridMssReduction(){}

    task* execute(){
        if(depth == 0){
            // Call parallel_reduce
            __mss_naive result = __mss_naive(a);
            parallel_reduce(blocked_range<long>(left,right),result);
            res -> sum = result.sum;
            res -> mss = result.mss;
            res -> mts = result.mts;
            res -> mps = result.mps;
            res -> posl = result.posmssl;
            res -> posr = result.posmssr;
            res -> posMps = result.posmps;
            res -> posMts = result.posmts;

        }
        else{
            int newDepth = depth - 1;
            mss_struct resLeft;
            mss_struct resRight;
            long middle = (left + right) / 2;
            int lIndex = 2*index;
            int rIndex = lIndex + 1;

            // Call subtasks and return result;
                HybridMssReduction& lTask = *new(allocate_child()) HybridMssReduction(newDepth,lIndex,a,&resLeft,left,middle);
                HybridMssReduction &rTask = *new(allocate_child()) HybridMssReduction(newDepth,rIndex,a,&resRight,middle,right);

                set_ref_count(3);
                spawn(lTask);
                spawn_and_wait_for_all(rTask);

                // Join
                // Sum
                res -> sum = resLeft . sum + resRight . sum;

                // Mps
                double mpsAux = resLeft . sum + resRight . mps;
                if(mpsAux >= resLeft.mps){
                    res -> mps = mpsAux;
                    res -> posMps = resRight.posMps;
                }else{
                    res -> mps = resLeft.mps;
                    res -> posMps = resLeft.posMps;
                }

                // Mts
                double mtsAux = resRight.sum + resLeft.mts;
                if(mtsAux >= resRight.mts){
                    res -> mts = mtsAux;
                    res -> posMts = resLeft.posMts;
                }else{
                    res->mts = resRight.mts;
                    res->posMts = resRight.posMts;
                }

                // Mss
                if(resLeft.mss >= resRight.mss){
                    res->mss = resLeft.mss; 
                    res->posl = resLeft.posl;
                    res->posr = resLeft.posr;
                }else{
                    res->mss = resRight.mss; 
                    res->posl = resRight.posl;
                    res->posr = resRight.posr;
                }
                double mssAux = resRight.mps + resLeft.mts;
                if(mssAux >= res->mss){
                    res->mss = mssAux;
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
        mss_struct* res;
        unsigned long left;
        unsigned long right;
        boolean** memo;
};


void parallel_mss_hybrid(double* a, long size){
    task_scheduler_init init;
    const int maxDepth = 5;
    
    mss_struct res; 

    HybridMssReduction& root = *new(task::allocate_root()) HybridMssReduction(maxDepth-1,0,a,&res,0,size);

    task::spawn_root_and_wait(root);
    if(PRINT){
        cout << endl << "Hybrid mss" << endl;
        cout << "Sum: " << res.sum << endl;
        cout << "Mss: " << res.mss << endl;
        cout << "Left pos: " << res.posl << endl;
        cout << "Right pos: " << res.posr << endl;
        cout << "Mps: " << res.mps << endl;
        cout << "Pos: " << res.posMps << endl;
        cout << "Mts: " << res.mts << endl;
        cout << "Pos: " << res.posMts << endl;
    }
    init.terminate();

}

struct mss_struct_bool {
    mss_struct_bool() :
        sum(false),
        mss(false),
        mts(false),
        mps(false),
        posl(false),
        posr(false),
        posMps(false),
        posMts(false)
    {}
    bool sum;
    bool mss;
    bool mts;
    bool mps;
    bool posl;
    bool posr;
    bool posMps;
    bool posMts;
};

struct mss_struct_interval {
    __m128d sum;
    __m128d mss;
    __m128d mts;
    __m128d mps;
    long posl;
    long posr;
    long posMps;
    long posMts;
};

struct mss_struct_mpfr {
    mpreal sum;
    mpreal mss;
    mpreal mts;
    mpreal mps;
    long posl;
    long posr;
    long posMps;
    long posMts;
};

class HybridMssReductionInterval : public task {
    public :
        HybridMssReductionInterval(int depth_,int index_, double* a_, mss_struct_interval* res_, unsigned long left_, unsigned long right_,boolean*** memo_) :
            depth(depth_),
            index(index_),
            a(a_),
            res(res_),
            left(left_),
            right(right_),
            memo(memo_)
    {}
    ~HybridMssReductionInterval(){}

    task* execute(){
        if(depth == 0){
            // Call parallel_reduce
            __mss_interval_without_pos result = __mss_interval_without_pos(a);
            parallel_reduce(blocked_range<long>(left,right),result);
            res -> sum = result.sum;
            res -> mss = result.mss;
            res -> mts = result.mts;
            res -> mps = result.mps;

        }
        else{
            int newDepth = depth - 1;
            mss_struct_interval resLeft;
            mss_struct_interval resRight;
            long middle = (left + right) / 2;
            int lIndex = 2*index;
            int rIndex = lIndex + 1;

            // Call subtasks and return result;
                HybridMssReductionInterval& lTask = *new(allocate_child()) HybridMssReductionInterval(newDepth,lIndex,a,&resLeft,left,middle,memo);
                HybridMssReductionInterval &rTask = *new(allocate_child()) HybridMssReductionInterval(newDepth,rIndex,a,&resRight,middle,right,memo);

                set_ref_count(3);
                spawn(lTask);
                spawn_and_wait_for_all(rTask);

                // Join
                boolean b;
                // Sum
                res->sum = in2_add(resLeft.sum,resRight.sum);

                // Mps
                __m128d mpsAux = in2_add(resLeft.sum,resRight.mps);
                b = in2_ge(mpsAux,resLeft.mps);
                switch(b){
                    case True:
                    res -> mps = mpsAux;
                    res -> posMps = resRight.posMps;
                    break;
                    case False:
                    res -> mps = resLeft.mps;
                    res -> posMps = resLeft.posMps;
                    break;
                    default:
                    res -> mps = in2_merge(mpsAux,resLeft.mps);
                }
                memo[depth][index][0] = b;

                // Mts
                __m128d mtsAux = in2_add(resRight.sum,resLeft.mts);
                b = in2_ge(mtsAux,resRight.mts);
                switch(b){
                    case True:
                    res -> mts = mtsAux;
                    res -> posMts = resLeft.posMts;
                    break;
                    case False:
                    res->mts = resRight.mts;
                    res->posMts = resRight.posMts;
                    break;
                    default:
                    res->mts = in2_merge(mtsAux,resRight.mts);
                }
                memo[depth][index][1] = b;

                // Mss
                b = in2_ge(resLeft.mss,resRight.mss);
                switch(b){
                    case True:
                    res->mss = resLeft.mss; 
                    res->posl = resLeft.posl;
                    res->posr = resLeft.posr;
                    break;
                    case False:
                    res->mss = resRight.mss; 
                    res->posl = resRight.posl;
                    res->posr = resRight.posr;
                    break;
                    default:
                    res->mss = in2_merge(resLeft.mss,resRight.mss);
                }
                memo[depth][index][2] = b;

                __m128d mssAux = in2_add(resRight.mps,resLeft.mts);
                b = in2_ge(mssAux,res->mss);
                switch(b){
                    case True:
                    res->mss = mssAux;
                    res->posl = resLeft.posMts;
                    res->posr = resRight.posMps;
                    break;
                    case Undefined:
                    res->mss = in2_merge(res->mss,mssAux);
                    break;
                    default:
                    break;
                }
                memo[depth][index][3] = b;

        }
        return NULL;
    }
    private:
        int depth;
        int index;
        double* a;
        mss_struct_interval* res;
        unsigned long left;
        unsigned long right;
        boolean*** memo;
};

class HybridMssReductionMpfr : public task {
    public :
        HybridMssReductionMpfr(int depth_,int index_, double* a_, mss_struct_mpfr* res_,mss_struct_bool* res_bool_, unsigned long left_, unsigned long right_,boolean*** memo_) :
            depth(depth_),
            index(index_),
            a(a_),
            res(res_),
            res_bool(res_bool_),
            left(left_),
            right(right_),
            memo(memo_)
    {}
    ~HybridMssReductionMpfr(){}

    task* execute(){
        if(depth == 0){
            // Call parallel_reduce
            bool auxMss = res_bool->mss || res_bool->posl || res_bool -> posr;
            bool auxMps = res_bool->mps || res_bool -> posMps;
            bool auxMts = res_bool->mts || res_bool -> posMts;
            bool auxSum = res_bool->sum;
            if(auxMss || auxMps || auxMts || auxSum){
                __mss_mpfr result = __mss_mpfr(a);
                parallel_reduce(blocked_range<long>(left,right),result);
                res -> sum = result.sum;
                res -> mss = result.mss;
                res -> mts = result.mts;
                res -> mps = result.mps;
                res -> posMps = result.posmps;
                res -> posMts = result.posmts;
                res -> posl = result.posmssl;
                res -> posr = result.posmssr;
            }
        }
        else{
            int newDepth = depth - 1;
            mss_struct_mpfr resLeft;
            mss_struct_mpfr resRight;
            unsigned long middle = (left + right) / 2;
            int lIndex = 2*index;
            int rIndex = lIndex + 1;
            
            mss_struct_bool resBoolLeft, resBoolRight;
            bool t;

            boolean b = memo[depth][index][3];
            switch(b){
                case True:
                if(res_bool->mss){
                    res_bool->mss = false;
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
                t = res_bool->mss || res_bool->posl || res_bool->posr;

                // Assignments
                if(res_bool->mss){
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
                    res_bool->mss = true;
                    resBoolLeft.mts = true;
                    resBoolRight.mps = true;
                }
                default:
                break;
            }

            b = memo[depth][index][2];
            switch(b){
                case True:
                if(res_bool->mss){
                    res_bool->mss = false;
                    resBoolLeft.mss = true;
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
                if(res_bool->mss){
                    res_bool->mss = false;
                    resBoolRight.mss = true;
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
                t = res_bool->mss || res_bool->posl || res_bool->posr;

                if(res_bool->mss){
                    res_bool->mss = false;
                    resBoolRight.mss = true;
                    resBoolLeft.mss = true;
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
                    resBoolLeft.mss = true;
                    resBoolRight.mss = true;
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
            HybridMssReductionMpfr& lTask = *new(allocate_child()) HybridMssReductionMpfr(newDepth,lIndex,a,&resLeft,&resBoolLeft,left,middle,memo);
            HybridMssReductionMpfr &rTask = *new(allocate_child()) HybridMssReductionMpfr(newDepth,rIndex,a,&resRight,&resBoolRight,middle,right,memo);

            set_ref_count(3);
            spawn(lTask);
            spawn_and_wait_for_all(rTask);

            // Join





        }
        return NULL;
    }
    private:
        int depth;
        int index;
        double* a;
        mss_struct_mpfr* res;
        mss_struct_bool* res_bool;
        unsigned long left;
        unsigned long right;
        boolean*** memo;
};

void parallel_mss_hybrid_interval(double* a, long size){
    const int maxDepth = 5;
    
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

    mss_struct_interval res; 

    HybridMssReductionInterval& root = *new(task::allocate_root()) HybridMssReductionInterval(maxDepth-1,0,a,&res,0,size,memo);

    task::spawn_root_and_wait(root);
    if(PRINT){
        cout << endl << "Hybrid mss interval" << endl;
        cout << "Sum: " ;
        print(res.sum);
        cout << endl;
        cout << "Mss: ";
        print(res.mss);
        cout << endl;
        cout << "Left pos: " << res.posl << endl;
        cout << "Right pos: " << res.posr << endl;
        cout << "Mps: ";
        print(res.mps);
        cout << endl;
        cout << "Pos: " << res.posMps << endl;
        cout << "Mts: ";
        print(res.mts);
        cout << endl;
        cout << "Pos: " << res.posMts << endl;
    }
    init.terminate();
    _MM_SET_ROUNDING_MODE(0);

}
