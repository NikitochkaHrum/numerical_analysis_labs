#pragma once

#include "../../matrices.h"

void heat_eq_header(int n, ld h, ld hi, Vec &s, bool expl);
void heat_eq_iter(Vec &u, ld t, ld err);
void heat_eq_max_err(ld err);

void log_shooting_method_header();
void log_shooting_method_iteration(int iter, ld alpha, ld y, ld B);
void log_runge_kutta_header();
void log_runge_kutta_iteration(ld x, ld y, ld ex_y, ld z);
