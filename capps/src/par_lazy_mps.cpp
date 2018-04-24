#include "par_lazy_mps.hpp"

MpsTask::~MpsTask(){
    for(int i = 0; i < maxIndex; i++){
        delete[] mps_values[i];
        delete[] mps_values2[i];
        delete[] sum_values[i];
        delete[] sum_values2[i];
    } 
    delete[] mps_values;
    delete[] sum_values;
    delete[] mps_values2;
    delete[] sum_values2;
}
        
task* MpsTask::execute(){
}
