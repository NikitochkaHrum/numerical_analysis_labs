#pragma once

#include "method_log/method_info.h"

// Таблица разделенных разностей
Mat table_separated_diff(const Mat& table_of_values);

// Интерполяция по формуле Ньютона
ld newton_method(const Mat& table, ld x, ld m, ld& error_est);

// Решение уравнения f(x) = c обратным интерполированием
ld root_of_equation(const Mat& table_of_values, ld c);

// Вывод результатов
void newton_method(function<ld(ld)> f, ld m, const Mat& values);
