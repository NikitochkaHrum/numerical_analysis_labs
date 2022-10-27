#pragma once

#include "../../matrices.h"

void print_method_header(const string& method_name, bool end = true);
void print_method_header_with_error(const string& method_name, bool end = true);
void print_method_header_with_jacobi(const string& method_name);
void print_method_header_with_additional_iterations(const string& method_name);

void print_method_iteration_info(const int iteration, const ld x, const ld y, const ld residual_error_norm, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2 , bool switch_to_new_line);
void print_method_iteration_info(const int iteration, const ld x, const ld y, const ld residual_error_norm, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2 , const ld error, bool switch_to_new_line);
void print_method_iteration_info(const int iteration, const ld x, const ld y, const ld residual_error_norm, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2 , const ld error, const ld jacobi_matrix_norm);
void print_method_iteration_info(const int iteration, const ld x, const ld y, const ld residual_error_norm, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2 , const ld error, const ld alpha, const int additional_iterations_count);
