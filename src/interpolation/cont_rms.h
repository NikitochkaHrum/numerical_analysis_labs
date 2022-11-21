#pragma once

#include "../matrices.h"

void continuous_rms(function<ld(ld)> f, vector<function<ld(ld)>>& basis, Mat& values);
ld integral_value(function<ld(ld)> f, ld a, ld b, const int n);
