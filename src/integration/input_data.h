#include "../matrices.h"

function<ld(ld)> F_var9 = [](ld x) -> ld {
    return exp(x) + 6 * x + 3;
};

function<ld(ld)> Der_var9 = [](ld x) -> ld {
    return exp(x) + 6;
};

function<ld(ld)> F_var25 = [](ld x) -> ld {
    return pow(3, x) - x + 2;
};

function<ld(ld)> Der_var25 = [](ld x) -> ld {
    return pow(3, x) * log(3) - 1;
};

const ld J_EXACT_var9 = exp(2) - exp(1) + 12.;
const ld J_EXACT_var25 = 6./log(3) + 0.5;
const ld A = 1;
const ld B = 2;