#pragma once
#include "../matrices.h"
#include "method_log/method_info.h"
#include "calc_params.h"

Vec conj_grad(Mat a, Vec b, Vec init, Vec solution, ld EPS);

int conj_grad_estimate_n_iterations(ld EPS, ld cond);
