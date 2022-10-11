#pragma once
#include "../matrices.h"
#include "method_log/method_info.h"
#include "calc_params.h"

Vec gradient_descent(Mat a, Vec b, Vec init, Vec solution, ld EPS);

int gradient_descent_estimate_n_iterations(ld EPS, ld cond);
