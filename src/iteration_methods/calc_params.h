#pragma once
#include "../matrices.h"

std::tuple<ld, ld, ld, ld> calc_params(Vec* x1, Vec* x2, Vec* x3, Mat& a, Vec& b, Vec& solution, ld solution_norm);
