#pragma once
#include "../../matrices.h"
#include "../calc_params.h"

void print_method_header(const string& method_name, const int iteration_cnt);
void print_method_iteration_info(const int iteration, const ld tau, const ld q, const ld norm_a, const ld norm_error, const ld estimate_error, Vec* x);
void print_omega_search_header();
void print_omega_search(ld omega, int iterations_cnt);
void print_found_omega(ld omega, int iterations_cnt);