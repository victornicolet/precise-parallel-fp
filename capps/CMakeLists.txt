cmake_minimum_required(VERSION 3.1)
project(precise_parallel_fp)
enable_language(CXX)
set (CMAKE_CXX_STANDARD_REQUIRED 11)


set (EXBLAS_HOME "./exblas/")
add_subdirectory(exblas)
# add exblas
include_directories (exblas/include/)
include_directories (./)
link_directories (exblas/lib/)
link_directories (./)
link_libraries(gomp)
link_libraries(exblas)
# Add OpenMP
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES pfpdefs.hpp main.cpp par_precise_fp.cpp par_precise_fp.hpp pfpdefs.cpp)
add_executable(precise_parallel_fp ${SOURCE_FILES})

