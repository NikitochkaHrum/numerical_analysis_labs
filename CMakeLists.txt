cmake_minimum_required(VERSION 3.5)
project (lab)

include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/include
)

add_executable(lab src/lab2.cpp src/matrices.cpp src/methods.cpp)