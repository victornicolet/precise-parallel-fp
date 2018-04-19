# Define targets
file(GLOB SOURCE_FILES "*.cpp" "*.hpp" "*.h")
add_executable(precise_parallel_fp ${SOURCE_FILES})
