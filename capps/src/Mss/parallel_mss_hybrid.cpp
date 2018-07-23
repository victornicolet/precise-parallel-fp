/* Lazy mss implementation with hybrid reduction operation */

#include <iostream>

using namespace std;

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
            res->mss = mss;
            res->mts = mts;
            res->pos = pos;
            return NULL;
        }
        else{
            // Call subtasks and return result;
                MpsTask1& lTask = *new(c.allocate_child()) MpsTask1(Cutoff,array,sizel,sum,mps,position,left,middle);
                MpsTask1 &rTask = *new(c.allocate_child()) MpsTask1(depth,array,&c.rres,middle,right);

                c.set_ref_count(3);
        }
    }
}


void parallel_mps_interval_hybrid(){
}
