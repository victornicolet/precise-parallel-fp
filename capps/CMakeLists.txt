cmake_minimum_required(VERSION 3.1)
project(precise_parallel_fp)
enable_language(CXX)
set (CMAKE_CXX_STANDARD_REQUIRED 11)

####################
# General parameters
####################

# exblas setting
set (EXBLAS_HOME "./exblas/")

# Add compilation database for ycm/libclang
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )


# Compilation flags
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11 -g -msse3 -march=native -fabi-version=0 -O3 -Wall -fopenmp -masm=intel")

# Define output directory for binaries
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

###################
# Subdirectories ##
###################

add_subdirectory(exblas)
add_subdirectory(src)



