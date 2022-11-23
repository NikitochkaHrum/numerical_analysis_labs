#pragma once

#include "../matrices.h"

const double A = 1;
const double B = 2;

const int N = 5;

ld F_var9(ld x);

ld Der_var9(ld x);

ld Der4_var9(ld x);

ld Der5_var9(ld x);

ld Der6_var9(ld x);

ld F_var25(ld x);

ld Der_var25(ld x);

ld Der4_var25(ld x);

ld Der5_var25(ld x);

ld Der6_var25(ld x);

ld omega(ld x, vector<ld> &points);

Mat table_of_values(function<ld(ld)> f, int n, ld a, ld b);