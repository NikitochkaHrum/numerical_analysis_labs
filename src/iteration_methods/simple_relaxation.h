#pragma once
#include "../matrices.h"
#include "method_log/method_info.h"
#include "calc_params.h"

Vec simple_relaxation(Mat a, Vec b, Vec init, Vec solution, ld EPS);

Vec get_solution_simple_relaxation(Mat a, Vec b, Vec init, Vec solution, ld EPS, ld omega, int& iter_count, bool log);

int simple_relaxation_estimate_n_iterations(ld EPS, ld cond);
