#include "input_data.h"

ld f1_v9(ld x, ld y) {
    return cos(x + 0.5) - y - 2.;
}

ld f2_v9(ld x, ld y) {
    return sin(y) - 2 * x - 1;
}

ld f1_v25(ld x, ld y) {
    return cos(y + 0.5) + x - 0.8;
}

ld f2_v25(ld x, ld y) {
    return sin(x) - 2 * y - 1.6;
}

ld f1_v20(ld x, ld y) {
    return sin(y + 2) - x - 1.5;
}

ld f2_v20(ld x, ld y) {
    return y + cos(x - 2) - 0.5;
}

ld fi1_v9(ld x, ld y) {
    return (sin(y) - 1.) / 2.;
}

ld fi2_v9(ld x, ld y) {
    return cos(x + 0.5) - 2.;
}

ld fi1_v25(ld x, ld y) {
    return 0.8 - cos(y + 0.5);
}

ld fi2_v25(ld x, ld y) {
    return (sin(x) - 1.6) / 2.;
}

ld fi1_v20(ld x, ld y) {
    return sin(y + 2) - 1.5;
}

ld fi2_v20(ld x, ld y) {
    return 0.5 - cos(x - 2);
}

ld F_v9(ld x, ld y) {
    return pow(f1_v9(x, y), 2) + pow(f2_v9(x, y), 2);
}

ld F_v25(ld x, ld y) {
    return pow(f1_v25(x, y), 2) + pow(f2_v25(x, y), 2);
}

ld F_v20(ld x, ld y) {
    return pow(f1_v20(x, y), 2) + pow(f2_v20(x, y), 2);
}

ld der_fi1_x_v9(ld x, ld y) {
    return 0;
}

ld der_fi1_y_v9(ld x, ld y) {
    return cos(y) / 2.;
}

ld der_fi2_x_v9(ld x, ld y) {
    return -sin(x + 0.5);
}

ld der_fi2_y_v9(ld x, ld y) {
    return 0;
}

ld der_fi1_x_v25(ld x, ld y) {
    return 0;
}

ld der_fi1_y_v25(ld x, ld y) {
    return sin(y + 0.5);
}

ld der_fi2_x_v25(ld x, ld y) {
    return cos(x) / 2.;
}

ld der_fi2_y_v25(ld x, ld y) {
    return 0;
}

ld der_fi1_x_v20(ld x, ld y) {
    return 0;
}

ld der_fi1_y_v20(ld x, ld y) {
    return cos(y + 2);
}

ld der_fi2_x_v20(ld x, ld y) {
    return sin(x - 2);
}

ld der_fi2_y_v20(ld x, ld y) {
    return 0;
}

ld der_f1_x_v9(ld x, ld y) {
    return -sin(x + 0.5);
}

ld der_f1_y_v9(ld x, ld y) {
    return -1;
}

ld der_f2_x_v9(ld x, ld y) {
    return -2;
}

ld der_f2_y_v9(ld x, ld y) {
    return cos(y);
}

ld der_f1_x_v25(ld x, ld y) {
    return 1;
}

ld der_f1_y_v25(ld x, ld y) {
    return -sin(y + 0.5);
}

ld der_f2_x_v25(ld x, ld y) {
    return cos(x);
}

ld der_f2_y_v25(ld x, ld y) {
    return -2;
}

ld der_f1_x_v20(ld x, ld y) {
    return -1;
}

ld der_f1_y_v20(ld x, ld y) {
    return cos(y + 2);
}

ld der_f2_x_v20(ld x, ld y) {
    return -sin(x - 2);
}

ld der_f2_y_v20(ld x, ld y) {
    return 1;
}

ld der_F_x_v9(ld x, ld y) {
    return 2 * f1_v9(x, y) * der_f1_x_v9(x, y) + 2 * f2_v9(x, y) * der_f2_x_v9(x, y);
}

ld der_F_y_v9(ld x, ld y) {
    return 2 * f1_v9(x, y) * der_f1_y_v9(x, y) + 2 * f2_v9(x, y) * der_f2_y_v9(x, y);
}

ld der_F_x_v25(ld x, ld y) {
    return 2 * f1_v25(x, y) * der_f1_x_v25(x, y) + 2 * f2_v25(x, y) * der_f2_x_v25(x, y);
}

ld der_F_y_v25(ld x, ld y) {
    return 2 * f1_v25(x, y) * der_f1_y_v25(x, y) + 2 * f2_v25(x, y) * der_f2_y_v25(x, y);
}

ld der_F_x_v20(ld x, ld y) {
    return 2 * f1_v20(x, y) * der_f1_x_v20(x, y) + 2 * f2_v20(x, y) * der_f2_x_v20(x, y);
}

ld der_F_y_v20(ld x, ld y) {
    return 2 * f1_v20(x, y) * der_f1_y_v20(x, y) + 2 * f2_v20(x, y) * der_f2_y_v20(x, y);
}