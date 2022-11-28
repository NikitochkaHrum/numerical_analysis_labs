#pragma once

#include "../matrices.h"

typedef tuple<ld, ld, ld> triplet;

// Метод Рунге-Кутты 4-го порядка
pair<ld, ld> runge_kutta_step(function<ld(ld, ld, ld)> &Z, ld h, ld x, ld y, ld z);
vector<triplet> runge_kutta_method(function<ld(ld, ld, ld)> &Z, function<ld(ld)> &Y_r, ld x_first, ld x_last, ld y, ld z, ld h, ld EPS, bool log=true);