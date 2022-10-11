cmake_minimum_required(VERSION 3.5)
project (lab)

include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/include
)

add_executable(lab src/lab3.cpp src/matrices.cpp src/direct_methods/direct_sqrt.cpp src/iteration_methods/simple_iteration.cpp src/iteration_methods/calc_params.cpp src/iteration_methods/method_log/method_info.cpp)