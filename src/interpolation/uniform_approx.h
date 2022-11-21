#pragma once

#include "../matrices.h"

void uniform_approximation_first_order_polynomial(function<ld(ld)> f, function<ld(ld)> der_f, Mat& values);
