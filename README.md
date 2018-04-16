# Precise parallel floating point applications

## How to modify

This project is meant to act as a rough template to define other reduction
that will still guarantee reproducibility. 

### How to add your own reduction

To add you own example using floating point expansions and superaccumulators,
you can follow these steps. 

Let's say we want to add an experiment for a reduction called 'mps'. First, you
can add the source files in `capps/exblas/src/cpu/blas1`. We will take the 
'mts' implementation as an example. 

#### Floating point expansion
You will need one class defining the floating point expansion, such as 
`ExMPS_FPE.hpp`. For that, just copy one of the files already present 
(like `ExMTS.hpp`) and adapt the `FPExpansionVectM1` class to your needs.

There are two main methods you should look into:
```c++
      void XAccumulate(double x);
      void Flush();
```
The first one defines the accumulation operation of your reduction operator. 
The second one defines how the values of the expansion are flushed into the
superaccumulator. It should be called at least once per thread, since only the
superaccumulators will be passed to the join.

#### Superaccumulator only version

The file `ExMTS.hpp / ExMTS.cpp` contains a TBB reduction implementation using
only the superaccumulator.
You should also copy and adapt this file. Look at the last two methods it 
contains: one of them defines the superaccumulator only function, the
other one is the floating point expansion function:
```C++
template<typename CACHE> __mts ExMTSFPE(int N, double *a);
```
You do not need to worry about the CACHE type. The first parameter defines the
size of the floating point expansion to be used, the second one the input
array. The implementation if in `ExMTS.cpp`


#### Adding your reduction to ExBLAS

Once you have implemented your components, go in `capps/exblas/blas1.hpp` and
add your main function there. For example. this main function for MTS has been
defined in `ExMTS.cpp`.

#### Use your reduction

You can look at `capps/test_mts.cpp` for an example on how to use the reduction
you just defined. You should just need to import `blas1.hpp` in order to 
use the main function of your reduction.


*To be continued...*