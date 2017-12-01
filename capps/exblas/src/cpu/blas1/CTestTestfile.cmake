# CMake generated Testfile for 
# Source directory: /home/nicolet/projects/precise_parallel_fp/capps/exblas/src/cpu/blas1
# Build directory: /home/nicolet/projects/precise_parallel_fp/capps/exblas/src/cpu/blas1
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(TestSumNaiveNumbers "test.exsum" "24")
set_tests_properties(TestSumNaiveNumbers PROPERTIES  PASS_REGULAR_EXPRESSION "TestPassed; ALL OK!")
add_test(TestSumStdDynRange "test.exsum" "24" "2" "0" "n")
set_tests_properties(TestSumStdDynRange PROPERTIES  PASS_REGULAR_EXPRESSION "TestPassed; ALL OK!")
add_test(TestSumLargeDynRange "test.exsum" "24" "50" "0" "n")
set_tests_properties(TestSumLargeDynRange PROPERTIES  PASS_REGULAR_EXPRESSION "TestPassed; ALL OK!")
add_test(TestSumIllConditioned "test.exsum" "24" "1e+50" "0" "i")
set_tests_properties(TestSumIllConditioned PROPERTIES  PASS_REGULAR_EXPRESSION "TestPassed; ALL OK!")
