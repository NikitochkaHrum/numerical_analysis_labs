#pragma once
#include "calc_params.h"
#include "method_log/method_info.h"

const ld GAMMA = 0.9;

Vec simple_iteration(Mat a, Vec b, Vec init, Vec solution, ld eps);

int simple_iteration_estimate_n_iterations(ld eps, ld cond);