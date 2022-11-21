#pragma once

#include "../../matrices.h"

void print_error(function<ld(ld)> g, Mat& values);
void print_polynom(Vec &coefficients);

void print_interpolation_results_header();
void print_interpolation_results(Vec& x, function<ld(ld)> f, function<ld(ld)> g, ld error_est);
void print_interpolation_results(Vec& x, function<ld(ld)> f, Vec& g_values, ld error_est);
void print_interpolation_results(Vec& x, function<ld(ld)> f, function<ld(ld)> g, Vec& errors);
