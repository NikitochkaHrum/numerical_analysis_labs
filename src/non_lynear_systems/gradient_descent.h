#pragma once

#include "../matrices.h"

// Метод градиентного спуска (поиск alpha методом наискорейшего спуска)
Vec gradient_descent(Vec init, Vec solution, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2,
                        function<ld(ld, ld)> der_f1_x, function<ld(ld, ld)> der_f1_y,
                        function<ld(ld, ld)> der_f2_x, function<ld(ld, ld)> der_f2_y,
                        function<ld(ld, ld)> F,
                        function<ld(ld, ld)> der_F_x, function<ld(ld, ld)> der_F_y,
                        ld eps);

pair<ld, int> golden_ratio(const function<ld(ld)> f, ld a, ld b, ld eps); // Метод золотого сечения
static const ld sqrt_5 = std::sqrt(5);
ld gr_left(ld a, ld b); // Метод для вычисления c'
ld gr_right(ld a, ld b); // Метод для вычисления c
