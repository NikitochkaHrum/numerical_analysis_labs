#pragma once

#include "../matrices.h"

void cubic_spline_interpolation(function<ld(ld)> f, function<ld(ld)> der_f, Mat& values, ld m, ld m1);
