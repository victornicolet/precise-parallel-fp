
#include <iostream>
#include <algorithm>

#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/task_group.h"
#include "tbb/task_scheduler_init.h"

#include "parallel_mss.hpp"
#include "parallel_mps.hpp"
#include "sequential_mts.hpp"

#include "debug.hpp"

using namespace tbb;
using namespace std;


void printA(double* array, int size){
    // printing the array
    cout << "{";
    if(size <= 10){
        for(int i = 0; i < size; i++){
            cout << array[i] << ",";
        }
        cout << "}" << endl;
    }
    else{
        for(int i = 0; i < 10; i++){
            cout << array[i] << ",";
        }
        cout << "... (size: " << size << ")" << endl;
    }

}

void parallel_mss_double(double* array, long size){
    task_scheduler_init init;
    __mss_naive result(array);
    parallel_reduce(blocked_range<long>(0,size),result);
    if(PRINT){
        cout << endl << "Parallel double" << endl;
        printA(array,size);
        result.print_mss();
    }
    init.terminate();
}

void parallel_mss_interval(double* array, long size){
    _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
    task_scheduler_init init;
    __mss_interval result(array);
    parallel_reduce(blocked_range<long>(0,size),result);
    init.terminate();
    _MM_SET_ROUNDING_MODE(0);

    // If non precise result redo computation
    long posl = -1;
    long posr = -1;
    if(result.posmssl1 != result.posmssl2 && result.posmssr1 != result.posmssr2){

       // First case : posmssl2 <= posmssr1
       if(result.posmssl2 <= result.posmssr1){
           task_scheduler_init init2(2);
           __mps_mpfr resultmps(array);
                       
           parallel_reduce(blocked_range<long>(result.posmssr1,result.posmssr2),resultmps);

           posr = resultmps.position;
           
           // Creation of parameters
           long length = result.posmssl2 - result.posmssl1;
           double* a = new double[length];
           for(long i = 0; i != length; i++){
               a[i] = array[i+result.posmssl1];
           }
           double mts;
           sequential_mts_superacc(a,length,&mts,&posl);
           posl += result.posmssl1;
        }
       // Second case : posmssl2 > posmssr1
       else {
           cout << endl <<  "Recomputation not yet implemented, because super rare" << endl;
        }
    }
    else if(result.posmssr1 != result.posmssr2){
           task_scheduler_init init2(2);
           __mps_mpfr resultmps(array);
                       
           parallel_reduce(blocked_range<long>(result.posmssr1,result.posmssr2),resultmps);

           posr = resultmps.position;
    }
    else if(result.posmssl1 != result.posmssl2){
           // Creation of parameters
           long length = result.posmssl2 - result.posmssl1;
           double* a = new double[length];
           for(long i = 0; i != length; i++){
               a[i] = array[i+result.posmssl1];
           }
           double mts;
           sequential_mts_superacc(a,length,&mts,&posr);
           posl += result.posmssl1;
    }
    if(PRINT){
        cout << endl << "Dynamic lazy computation first results" << endl;
        printA(array,size);
        result.print_mss();
        cout << endl << "Precise results" << endl;
        cout << "Left pos: " << posl << endl;
        cout << "Right pos: " << posr << endl;
    }
}

__mss_naive::__mss_naive(double* a) :
    array(a),
    sum(0.),
    mps(0.),
    mts(0.),
    mss(0.),
    posmps(0),
    posmts(0),
    posmssl(0),
    posmssr(0)
{}

__mss_naive::__mss_naive(__mss_naive& x, split) :
    array(x.array),
    sum(0.),
    mps(0.),
    mts(0.),
    mss(0.),
    posmps(0),
    posmts(0),
    posmssl(0),
    posmssr(0)
{}

void __mss_naive::operator()(const blocked_range<long>& r){
    if(posmps == 0){
        posmps = r.begin();
    }
    if(posmts == 0){
        posmts = r.begin();
    }
    if(posmssl == 0 && posmssr == 0){
        posmssl = r.begin();
        posmssr = r.begin();
    }

    for(int i = r.begin(); i != r.end(); i++){

        // Update sum
        sum += array[i];
        
        // Update mts
        mts += array[i];
        if(mts < 0){
            mts = 0.;
            posmts = i+1;
        }
        
        // Update mps
        if(sum >= mps){
            mps = sum;
            posmps = i+1;
        }

        // Update mss
        if(mts >= mss){
            mss = mts;
            posmssl = posmts;
            posmssr = i+1;
        }
    }
}

void __mss_naive::join(__mss_naive& right){

    // Join mss
    double aux = mts + right.mps;
    bool test1 = aux >= mss;
    bool test2 = right.mss >= aux;
    bool test3 = right.mss >= mss;

    if(test2 && test3){
        mss = right.mss;
        posmssl = right.posmssl;
        posmssr = right.posmssr;
    }
    else if (test1 && !test2){
        mss = aux;
        posmssl = posmts;
        posmssr = right.posmps;
    }
    
    // Join mps
    right.mps += sum;
    if(right.mps >= mps){
        mps = right.mps;
        posmps = right.posmps;
    }

    // Join mts
    mts += right.sum;
    if(right.mts >= mts){
        mts = right.mts;
        posmts = right.posmts;
    }

    // Join sum
    sum += right.sum;
}

void __mss_naive::print_mss(){
    cout << "Sum: " << sum << endl;
    cout << "Mss: " << mss << endl;
    cout << "Left pos: " << posmssl << endl;
    cout << "Right pos: " << posmssr << endl;
    cout << "Mps: " << mps << endl;
    cout << "Pos: " << posmps << endl;
    cout << "Mts: " << mts << endl;
    cout << "Pos: " << posmts << endl;

}

__mss_mpfr::__mss_mpfr(double* a) :
    array(a),
    sum(0.),
    mps(0.),
    mts(0.),
    mss(0.),
    posmps(0),
    posmts(0),
    posmssl(0),
    posmssr(0)
{}

__mss_mpfr::__mss_mpfr(__mss_mpfr& x, split) :
    array(x.array),
    sum(0.),
    mps(0.),
    mts(0.),
    mss(0.),
    posmps(0),
    posmts(0),
    posmssl(0),
    posmssr(0)
{}

void __mss_mpfr::operator()(const blocked_range<long>& r){
    if(posmps == 0){
        posmps = r.begin();
    }
    if(posmts == 0){
        posmts = r.begin();
    }
    if(posmssl == 0 && posmssr == 0){
        posmssl = r.begin();
        posmssr = r.begin();
    }

    for(int i = r.begin(); i != r.end(); i++){

        // Update sum
        sum += array[i];
        
        // Update mts
        mts += array[i];
        if(mts < 0){
            mts = 0.;
            posmts = i+1;
        }
        
        // Update mps
        if(sum >= mps){
            mps = sum;
            posmps = i+1;
        }

        // Update mss
        if(mts >= mss){
            mss = mts;
            posmssl = posmts;
            posmssr = i+1;
        }
    }
}

void __mss_mpfr::join(__mss_mpfr& right){

    // Join mss
    mpreal aux = mts + right.mps;
    bool test1 = aux >= mss;
    bool test2 = right.mss >= aux;
    bool test3 = right.mss >= mss;

    if(test2 && test3){
        mss = right.mss;
        posmssl = right.posmssl;
        posmssr = right.posmssr;
    }
    else if (test1 && !test2){
        mss = aux;
        posmssl = posmts;
        posmssr = right.posmps;
    }
    
    // Join mps
    right.mps += sum;
    if(right.mps >= mps){
        mps = right.mps;
        posmps = right.posmps;
    }

    // Join mts
    mts += right.sum;
    if(right.mts >= mts){
        mts = right.mts;
        posmts = right.posmts;
    }

    // Join sum
    sum += right.sum;
}

void __mss_mpfr::print_mss(){
    cout << "Sum: " << sum << endl;
    cout << "Mss: " << mss << endl;
    cout << "Left pos: " << posmssl << endl;
    cout << "Right pos: " << posmssr << endl;
    cout << "Mps: " << mps << endl;
    cout << "Pos: " << posmps << endl;
    cout << "Mts: " << mts << endl;
    cout << "Pos: " << posmts << endl;

}

__mss_interval::__mss_interval(double* array_) :
    array(array_),
    posmssl1(0),
    posmssl2(0),
    posmssr1(0),
    posmssr2(0),
    posmts1(0),
    posmts2(0),
    posmps1(0),
    posmps2(0)
{
    sum = in2_create(0,0);
    mss = in2_create(0,0);
    mps = in2_create(0,0);
    mts = in2_create(0,0);
}

__mss_interval::__mss_interval(__mss_interval& x ,split s) :
    array(x.array),
    posmssl1(0),
    posmssl2(0),
    posmssr1(0),
    posmssr2(0),
    posmts1(0),
    posmts2(0),
    posmps1(0),
    posmps2(0)
{
    sum = in2_create(0,0);
    mss = in2_create(0,0);
    mps = in2_create(0,0);
    mts = in2_create(0,0);
}

void __mss_interval::operator()(const blocked_range<long>& r){
    if(posmts1 == 0 && posmts2 == 0){
        posmts1 = r.begin();
        posmts2 = r.begin();
    }
    if(posmps1 == 0 && posmps2 == 0){
        posmps1 = r.begin();
        posmps2 = r.begin();
    }
    if(posmssl1 == 0 && posmssl2 == 0 && posmssr1 == 0 && posmssr2 == 0){
        posmssl1 = r.begin();
        posmssl2 = r.begin();
        posmssr1 = r.begin();
        posmssr2 = r.begin();
    }
    for(int i = r.begin(); i != r.end(); i++){

        // Update sum
        sum = in2_add(sum,array[i]);
        
        // Update mts
        mts = in2_add(mts,array[i]);
        boolean b = in2_le(0,mts); 
        if(b == False){
            mts = in2_create(0.,0.);
            posmts1 = i+1;
            posmts2 = i+1;
        } else if (b == Undefined){
            posmts2 = i+1;
        }
        
        // Update mps
        boolean b0 = in2_le(mps,sum);
        if(b0 == True){
            mps = sum;
            posmps1 = i+1;
            posmps2 = i+1;
        } else if (b0 == Undefined){
            mps = in2_merge(mps,sum);
            posmps2 = i+1;
        }

        // Update mss
        boolean b1 = in2_le(mss,mts);
        if(b1 == True){
            mss = mts;
            posmssl1 = posmts1;
            posmssl2 = posmts2;
            posmssr1 = i+1;
            posmssr2 = i+1;
        } else if(b1 == Undefined){
            mss = in2_merge(mss,mts);
            if(posmssl1 > posmts1) posmssl1 = posmts1; 
            if(posmssl2 < posmts2) posmssl2 = posmts2;
            posmssr2 = i+1;
        }
    }
}

void __mss_interval::join(__mss_interval& right){

    // Join mss
    __m128d aux = in2_add(mts,right.mps);
    boolean test1 = in2_le(mss,aux);
    boolean test2 = in2_le(aux,right.mss);
    boolean test3 = in2_le(mss,right.mss);

    if(test2 == True && test3 == True){
        // Then right.mss
        mss = right.mss;
        posmssl1 = right.posmssl1;
        posmssl2 = right.posmssl2;
        posmssr1 = right.posmssr1;
        posmssr2 = right.posmssr2;
    }
    else if (test1 == True && test2 == False){
        // Then aux
        mss = aux;
        posmssl1 = posmts1;
        posmssl2 = posmts2;
        posmssr1 = right.posmps1;
        posmssr2 = right.posmps2;
    }
    else if (test1 == False && test3 == False){}
    else if (test1 == True || test3 == True){
        // Then either aux or right.mss
        mss = in2_merge(aux,right.mss);
        posmssl1 = posmts1;
        posmssl2 = right.posmssl2;
        posmssr1 = min(right.posmssr1,right.posmps1);
        posmssr2 = max(right.posmssr2,right.posmps2);
    }
    else if (test1 == False || test2 == True){
        // Then either mss or right.mss
        mss = in2_merge(mss,right.mss);
        posmssl2 = right.posmssl2;
        posmssr2 = right.posmssr2;
    }
    else if (test2 == False || test3 == False){
        // Then either aux or mss
        mss = in2_merge(mss,aux);
        posmssl1 = min(posmssl1,posmts1);
        posmssl2 = max(posmssl2,posmts2);
        posmssr2 = right.posmps2;
    }
    else {
        // Then merge all three intervals
        mss = in2_merge(mss,aux);
        mss = in2_merge(mss,right.mps);
        // Merge all three positions
        posmssl1 = min(posmssl1,posmts1);
        posmssl2 = right.posmssl2;
        posmssr1 = min(right.posmssr1,right.posmps1);
        posmssr2 = max(right.posmssr2,right.posmps2);
    }
    
    // Join mps
    right.mps = in2_add(right.mps,sum);
    boolean b = in2_le(mps,right.mps); 
    if(b == True){
        mps = right.mps;
        posmps1 = right.posmps1;
        posmps2 = right.posmps2;
    }
    else if(b == Undefined){
        mps = in2_merge(mps,right.mps);
        posmps2 = right.posmps2;
    }

    // Join mts
    mts = in2_add(mts,right.sum);
    boolean b0 = in2_le(mts,right.mts);
    if(b0 == True){
        mts = right.mts;
        posmts1 = right.posmts1;
        posmts2 = right.posmts2;
    }
    else if (b0 == Undefined){
        mts = in2_merge(mts,right.mps);
        posmts2 = right.posmts2;
    }


    // Join sum
    sum = in2_add(sum,right.sum);
}

void __mss_interval::print_mss(){
    cout << "Sum: ";
    print(sum);
    cout << endl << "Mss: ";
    print(mss);
    cout << endl << "Left pos: " << posmssl1 <<", " << posmssl2;
    cout << endl << "Right pos: " << posmssr1 << ", " << posmssr2;
    cout << endl << "Mps: ";
    print(mps);
    cout << endl << "Pos: " << posmps1 << ", " << posmps2;
    cout << endl << "Mts: ";
    print(mts);
    cout << endl << "Pos: " << posmts1 << ", " << posmts2;
    cout << endl;
}

// Without position
__mss_interval_without_pos::__mss_interval_without_pos(double* array_) :
    array(array_)
{
    sum = in2_create(0,0);
    mss = in2_create(0,0);
    mps = in2_create(0,0);
    mts = in2_create(0,0);
}

__mss_interval_without_pos::__mss_interval_without_pos(__mss_interval_without_pos& x ,split s) :
    array(x.array)
{
    sum = in2_create(0,0);
    mss = in2_create(0,0);
    mps = in2_create(0,0);
    mts = in2_create(0,0);
}

void __mss_interval_without_pos::operator()(const blocked_range<long>& r){
    for(int i = r.begin(); i != r.end(); i++){

        // Update sum
        sum = in2_add(sum,array[i]);
        
        // Update mts
        mts = in2_add(mts,array[i]);
        boolean b = in2_le(0,mts); 
        if(b == False){
            mts = in2_create(0.,0.);
        } else if (b == Undefined){
            mts = in2_merge(0.,mts);
        }
        
        // Update mps
        boolean b0 = in2_le(mps,sum);
        if(b0 == True){
            mps = sum;
        } else if (b0 == Undefined){
            mps = in2_merge(mps,sum);
        }

        // Update mss
        boolean b1 = in2_le(mss,mts);
        if(b1 == True){
            mss = mts;
        } else if(b1 == Undefined){
            mss = in2_merge(mss,mts);
        }
    }
}

void __mss_interval_without_pos::join(__mss_interval_without_pos& right){

    // Join mss
    __m128d aux = in2_add(mts,right.mps);
    boolean test1 = in2_le(mss,aux);
    boolean test2 = in2_le(aux,right.mss);
    boolean test3 = in2_le(mss,right.mss);

    if(test2 == True && test3 == True){
        // Then right.mss
        mss = right.mss;
    }
    else if (test1 == True && test2 == False){
        // Then aux
        mss = aux;
    }
    else if (test1 == False && test3 == False){}
    else if (test1 == True || test3 == True){
        // Then either aux or right.mss
        mss = in2_merge(aux,right.mss);
    }
    else if (test1 == False || test2 == True){
        // Then either mss or right.mss
        mss = in2_merge(mss,right.mss);
    }
    else if (test2 == False || test3 == False){
        // Then either aux or mss
        mss = in2_merge(mss,aux);
    }
    else {
        // Then merge all three intervals
        mss = in2_merge(mss,aux);
        mss = in2_merge(mss,right.mps);
    }
    
    // Join mps
    right.mps = in2_add(right.mps,sum);
    boolean b = in2_le(mps,right.mps); 
    if(b == True){
        mps = right.mps;
    }
    else if(b == Undefined){
        mps = in2_merge(mps,right.mps);
    }

    // Join mts
    mts = in2_add(mts,right.sum);
    boolean b0 = in2_le(mts,right.mts);
    if(b0 == True){
        mts = right.mts;
    }
    else if (b0 == Undefined){
        mts = in2_merge(mts,right.mps);
    }


    // Join sum
    sum = in2_add(sum,right.sum);
}

void __mss_interval_without_pos::print_mss(){
    cout << "Sum: ";
    print(sum);
    cout << endl << "Mss: ";
    print(mss);
    cout << endl << "Mps: ";
    print(mps);
    cout << endl << "Mts: ";
    print(mts);
    cout << endl;
}


