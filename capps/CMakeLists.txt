cmake_minimum_required(VERSION 3.1)
project(precise_parallel_fp)
enable_language(CXX)
set (CMAKE_CXX_STANDARD_REQUIRED 11)

set (EXBLAS_HOME "./exblas/")

# Add compilation database for ycm/libclang
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

# Compilation flags
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11 -g -msse3 -march=native -fabi-version=0 -O3 -Wall -fopenmp -masm=intel")

# This line seems to be useless, I comment it
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")

# Define output directory for binaries
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# source directories without includes
add_subdirectory(exblas)

include_directories (exblas/include/)
include_directories (exblas/src/common/)
include_directories(/usr/lib/gcc/x86_64-linux-gnu/6/include/)
include_directories(include/)
include_directories(include/Victor/)
link_directories (exblas/lib/)
link_libraries(gomp)
link_libraries(exblas)
link_libraries(mpfr)
link_libraries(gmp)
# Add OpenMP
find_package(OpenMP REQUIRED)

# source subdirectories with includes
add_subdirectory(src)


