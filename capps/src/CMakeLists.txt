
##################################
# Settings for precise parallel fp
##################################

# Include directories
include_directories (../exblas/include/)
include_directories (../exblas/src/common/)
include_directories(/usr/lib/gcc/x86_64-linux-gnu/6/include/)
include_directories(Victor/)
include_directories(../include/)
include_directories(../include/Victor/)
include_directories(../include/Common/)
include_directories(../include/Mps/)
include_directories(../include/Mss/)
include_directories(../include/Mts/)
include_directories(../include/Graphs/)

# Link directories
link_directories (../exblas/lib/)

# Link libraries
link_libraries(gomp)
link_libraries(exblas)
link_libraries(mpfr)
link_libraries(gmp)

# Add OpenMP
find_package(OpenMP REQUIRED)

# add matplotlibcpp
find_package(PythonLibs 2.7)
include_directories(PRIVATE ${PYTHON_INCLUDE_DIRS})
link_libraries(${PYTHON_LIBRARIES})



file(GLOB Common_sources "Common/*.cpp")
file(GLOB Mps_sources "Mps/*.cpp" "Mts/sequential_mts.cpp")
file(GLOB Mss_sources "Mss/*.cpp" "Mps/parallel_mps.cpp" "Mts/sequential_mts.cpp")
file(GLOB Graph_sources "Graphs/bellman_ford.cpp" "Graphs/test_graph.cpp")
file(GLOB Reductions_sources "*SmallTests/*.cpp" )

add_executable(test_mss ${Common_sources} ${Mss_sources} )
add_executable(test_mps ${Common_sources} ${Mps_sources})
add_executable(test_graph ${Common_sources} ${Graph_sources})
add_executable(test_reductions ${Reductions_sources})

