#pragma once

#include "../matrices.h"

const ld PI = acos(-1);
const ld HI = 0.1;

function<ld(ld)> A_9 = [](ld x) -> ld {
	return 50 * (x + 1);
};

function<ld(ld)> B_9 = [](ld x) -> ld {
	return -x*x+2;
};

function<ld(ld)> C_9 = [](ld x) -> ld {
	return x + 1;
};

function<ld(ld)> A_25 = [](ld x) -> ld {
	return 50 * (x + 1);
};

function<ld(ld)> B_25 = [](ld x) -> ld {
	return x*x+1;
};

function<ld(ld)> C_25 = [](ld x) -> ld {
	return x + 1;
};

function<ld(ld)> Y_9 = [](ld x) -> ld {
	return 1 + x + 10 * log(9/*вариант*/ + 1) * x * x * x * (1 - x) * (1 - x) * (1 - x);
};

function<ld(ld)> Y_der1_9 = [](ld x) -> ld {
	return 1 - 30 * (1 - x) * (1 - x) * x * x * (2 * x - 1) * log(9 + 1);
};

function<ld(ld)> Y_der2_9 = [](ld x) -> ld {
	return -60 * (5 * x * x * x * x- 10 * x * x * x + 6 * x * x - x) * log(9 + 1);
};

function<ld(ld)> Y_25 = [](ld x) -> ld {
	return 1 + x + 10 * log(25/*вариант*/ + 1) * x * x * x * (1 - x) * (1 - x) * (1 - x);
};

function<ld(ld)> Y_der1_25 = [](ld x) -> ld {
	return 1 - 30 * (1 - x) * (1 - x) * x * x * (2 * x - 1) * log(25 + 1);
};

function<ld(ld)> Y_der2_25 = [](ld x) -> ld {
	return -60 * (5 * x * x * x * x- 10 * x * x * x + 6 * x * x - x) * log(25 + 1);
};

function<ld(ld)> F_9 = [](ld x) -> ld {
	return Y_der2_9(x) + A_9(x) * Y_der1_9(x) - B_9(x) * Y_9(x) + C_9(x) * sin(Y_9(x));
};

function<ld(ld)> F_25 = [](ld x) -> ld {
	return Y_der2_25(x) + A_25(x) * Y_der1_25(x) - B_25(x) * Y_25(x) + C_25(x) * sin(Y_25(x));
};

function<ld(ld)> PHI = [](ld x) -> ld {
    return x;
};

function<ld(ld, ld, ld)> Z_9 = [](ld x, ld y, ld z) -> ld {
	return F_9(x) - A_9(x) * z + B_9(x) * y - C_9(x) * sin(y);
};

function<ld(ld, ld, ld)> Z_25 = [](ld x, ld y, ld z) -> ld {
	return F_25(x) - A_25(x) * z + B_25(x) * y - C_25(x) * sin(y);
};

function<ld(ld, ld)> U_9 = [](ld t, ld x) -> ld {
    return x + 0.1 * sin(M_PI * x) * t * 9/*вариант*/;
};

function<ld(ld, ld)> U_25 = [](ld t, ld x) -> ld {
    return x + 0.1 * sin(M_PI * x) * t * 25/*вариант*/;
};

function<ld(ld, ld)> FU_9 = [](ld t, ld x) -> ld {
    return 0.1 * 9/*вариант*/ * sin(M_PI * x) * (1 + HI * pow(M_PI, 2) * t);
};

function<ld(ld, ld)> FU_25 = [](ld t, ld x) -> ld {
    return 0.1 * 25/*вариант*/ * sin(M_PI * x) * (1 + HI * pow(M_PI, 2) * t);
};
