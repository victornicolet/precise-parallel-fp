# Graph subproject
add_subdirectory(Graphs/)

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

# Define targets
file(GLOB SOURCE_FILES "*.cpp" "*.hpp" "*.h")
add_executable(precise_parallel_fp ${SOURCE_FILES})
