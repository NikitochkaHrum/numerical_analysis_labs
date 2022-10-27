#include "input_data.h"

ld F(function<ld(ld, ld)> a1, function<ld(ld, ld)> a2, ld x, ld y){
    ld val1 = a1(x, y), val2= a2(x, y);
    return val1*val1 + val2*val2;
}

ld f1_v9(ld x, ld y){
    return cos(x + 0.5) - y - 2.;
}

ld f2_v9(ld x, ld y){
    return sin(y) - 2*x - 1;
}

ld f1_v25(ld x, ld y){
    return cos(y + 0.5) + x - 0.8;
}

ld f2_v25(ld x, ld y){
    return sin(x) - 2*y - 1.6;
}

ld fi1_v9(ld x, ld y){
    return (sin(y) - 1.)/2.;
}

ld fi2_v9(ld x, ld y){
    return cos(x + 0.5) - 2.;
}

ld fi1_v25(ld x, ld y){
    return 0.8 - cos(y + 0.5);
}

ld fi2_v25(ld x, ld y){
    return (sin(x) - 1.6)/2.;
}

ld der_fi1_x_v9(ld x, ld y){
    return 0;
}

ld der_fi1_y_v9(ld x, ld y){
    return cos(y)/2.;
}

ld der_fi2_x_v9(ld x, ld y){
    return -sin(x + 0.5);
}

ld der_fi2_y_v9(ld x, ld y){
    return 0;
}

ld der_fi1_x_v25(ld x, ld y){
    return 0;
}

ld der_fi1_y_v25(ld x, ld y){
    return sin(y + 0.5);
}

ld der_fi2_x_v25(ld x, ld y){
    return cos(x)/2.;
}

ld der_fi2_y_v25(ld x, ld y){
    return 0;
}

ld der_f1_x_v9(ld x, ld y){
    return -sin(x + 0.5);
}

ld der_f1_y_v9(ld x, ld y){
    return -1;
}

ld der_f2_x_v9(ld x, ld y){
    return -2;
}

ld der_f2_y_v9(ld x, ld y){
    return cos(y);
}

ld der_f1_x_v25(ld x, ld y){
    return 1;
}

ld der_f1_y_v25(ld x, ld y){
    return -sin(y + 0.5);
}

ld der_f2_x_v25(ld x, ld y){
    return cos(x);
}

ld der_f2_y_v25(ld x, ld y){
    return -2;
}