cmake_minimum_required(VERSION 3.1)
project(precise_parallel_fp)
enable_language(CXX)
set (EXBLAS_HOME "./exblas/")
add_subdirectory(exblas)
# add exblas
include_directories ("exblas/include/")
link_directories ("exblas/lib/")
link_libraries(gomp)
link_libraries(exblas)

set (CMAKE_CXX_STANDARD_REQUIRED 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu++11 -march=native -fabi-version=0 -ltbb
-O3 -Wall -std=c++11 -fopenmp -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5")

set(SOURCE_FILES main.cpp par_precise_fp.cpp pfpdefs.hpp par_precise_fp.hpp pfpdefs.cpp)
add_executable(precise_parallel_fp ${SOURCE_FILES})

