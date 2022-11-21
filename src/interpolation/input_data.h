#pragma once

#include "../matrices.h"

const double A = 1;
const double B = 2;

const int N = 5;

function<ld(ld)> F_var9 = [](ld x) -> ld {
    return exp(x) + 6 * x + 3;
};

function<ld(ld)> Der_var9 = [](ld x) -> ld {
    return exp(x) + 6;
};

function<ld(ld)> Der4_var9 = [](ld x) -> ld {
    return exp(x);
};

function<ld(ld)> Der5_var9 = [](ld x) -> ld {
    return exp(x);
};

function<ld(ld)> Der6_var9 = [](ld x) -> ld {
    return exp(x);
};

function<ld(ld)> F_var25 = [](ld x) -> ld {
    return pow(3, x) - x + 2;
};

function<ld(ld)> Der_var25 = [](ld x) -> ld {
    return pow(3, x) * log(3) - 1;
};

function<ld(ld)> Der4_var25 = [](ld x) -> ld {
    return pow(3, x) * pow(log(3), 4);
};

function<ld(ld)> Der5_var25 = [](ld x) -> ld {
    return pow(3, x) * pow(log(3), 5);
};

function<ld(ld)> Der6_var25 = [](ld x) -> ld {
    return pow(3, x) * pow(log(3), 6);
};

ld omega(ld x, vector<ld> &points);

Mat table_of_values(function<ld(ld)> f, int n, ld a, ld b);