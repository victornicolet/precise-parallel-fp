/* Lazy mss implementation with hybrid reduction operation */

#include <iostream>

#include "tbb/task_group.h"
#include "tbb/parallel_reduce.h"

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
        HybridMssReduction(int depth_, double* a_, mss_struct* res_, unsigned long left_, unsigned long right_) :
            depth(depth_),
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

            // Call subtasks and return result;
                HybridMssReduction& lTask = *new(allocate_child()) HybridMssReduction(newDepth,a,&resLeft,left,middle);
                HybridMssReduction &rTask = *new(allocate_child()) HybridMssReduction(newDepth,a,&resRight,middle,right);

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
                    res -> mps = mtsAux;
                    res -> posMts = resLeft.posMts;
                }else{
                    res->mts = resRight.mts;
                    res->posMts = resRight.posMts;
                }

                // Mss

        }
        return NULL;
    }
    private:
        int depth;
        double* a;
        mss_struct* res;
        unsigned long left;
        unsigned long right;
        boolean** memo;
};


void parallel_mps_interval_hybrid(){
}
