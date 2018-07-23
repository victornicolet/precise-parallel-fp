/* Lazy mss implementation with hybrid reduction operation */

#include <iostream>

#include "tbb/task_group.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_reduce.h"

#include "debug.hpp"
#include "parallel_mss.hpp"

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
    const int maxDepth = 10;
    
    int maxIndex = 1;
    boolean** memo = new boolean*[maxDepth]; 
    for(int i = maxDepth-1; i >= 0;i--){
        memo[i] = new boolean[maxIndex];
        maxIndex = 2*maxIndex;
    }

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
