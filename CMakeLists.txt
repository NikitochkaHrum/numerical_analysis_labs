cmake_minimum_required(VERSION 3.5)
project (lab)

include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/include
)

add_executable(lab  src/lab6.cpp
                    src/matrices.cpp
                    src/direct_methods/direct_sqrt.cpp
                    src/iteration_methods/simple_iteration.cpp
                    src/iteration_methods/gradient_descent.cpp
                    src/iteration_methods/simple_relaxation.cpp
                    src/iteration_methods/conj_grad.cpp
                    src/iteration_methods/calc_params.cpp
                    src/iteration_methods/method_log/method_info.cpp
                    src/non_lynear_systems/calc_params.cpp
                    src/non_lynear_systems/input_data.cpp
                    src/non_lynear_systems/newton_method.cpp
                    src/non_lynear_systems/simple_iteration.cpp
                    src/non_lynear_systems/gradient_descent.cpp
                    src/non_lynear_systems/method_log/method_info.cpp
                    src/integration/function_tool.cpp
                    src/integration/gauss.cpp
                    src/integration/simpson.cpp
                    src/integration/trapez.cpp
                    src/integration/trapez_spline.cpp
                    src/integration/simpson.cpp
                    src/integration/gauss.cpp
                    src/integration/method_log/method_info.cpp
                    )