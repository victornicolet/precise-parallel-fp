cmake_minimum_required(VERSION 2.8.8)
project(precise_parallel_fp)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -O2 -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5")

add_subdirectory(exblas)

set(SOURCE_FILES main.c par_precise_fp.c pfpdefs.h par_precise_fp.h pfpdefs.c)
add_executable(precise_parallel_fp ${SOURCE_FILES})