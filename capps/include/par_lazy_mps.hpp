/* Lazy computation of mps.
 * Author: Raphael Dang-Nhu
 * Date: 04/24/2018 */

#ifndef PAR_LAZY_MPS_H
#define PAR_LAZY_MPS_H

#include "interval_arithmetic.hpp"
#include "superaccumulator.hpp"

#include "tbb/task_group.h"

using namespace tbb;

class MpsTask: public task {
    public:
        MpsTask(int Cutoff, double* array, __m128d** mps_values, __m128d** sum_values, Superaccumulator** mps_values2, Superaccumulator** sum_values2, int depth, int index) : 
            Cutoff(Cutoff),
            array(array),
            depth(depth),
            index(index),
            mps_values(mps_values),
            sum_values(sum_values),
            mps_values2(mps_values2),
            sum_values2(sum_values2)
        {}
        ~MpsTask();
        task* execute();

    private:
        /* Below this size, the mps and sum are computed sequentially */
        int Cutoff;
        /* Input array and its size */
        double* array;
        int size;    
        /* Identification of the task */
        int depth;
        int index;
        /* Matrix containing double precision, interval arithmetic results */
        __m128d** mps_values;
        __m128d** sum_values;
        /* Matrix containing superaccumulators */
        Superaccumulator** mps_values2;
        Superaccumulator** sum_values2;
        /* Size of the matrixs */
        int maxDepth;
        int maxIndex;
};
        

#endif

