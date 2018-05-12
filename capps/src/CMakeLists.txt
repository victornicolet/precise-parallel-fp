# Graph subproject (single directory project without includes)
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
include_directories(../include/Common/)
include_directories(../include/Mps/)
include_directories(../include/Mss/)

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

# Single directory project with includes
add_subdirectory(SmallTests/)


file(GLOB Common_sources "Common/*.cpp")
file(GLOB Mps_sources "Mps/*.cpp")
file(GLOB Mss_sources "Mss/*.cpp")

add_executable(test_mss ${Common_sources} ${Mss_sources})
add_executable(test_mps ${Common_sources} ${Mps_sources})

