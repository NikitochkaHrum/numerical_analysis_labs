#include "method_info.h"

void print_method_header(const std::string& method_name, bool end)
{
	cout << method_name << '\n';
    cout << setw(6) << "Itr";
    cout << "|" << setw(12) << "X";
    cout << "|" << setw(12) << "Y";
	cout << "|   Норма невязки";
    cout << "|" << setw(16) << "f1";
    cout << "|" << setw(16) << "f2";
    if (end) {
        cout << "\n-----------------------------------------------------------------------------------------------------------------------------------------\n";
    }
}

void print_method_header_with_error(const string& method_name, bool end){
    print_method_header(method_name, false);
    cout << "|" << setw(19) << "Погрешность решения";
    if (end) {
        cout << "\n-----------------------------------------------------------------------------------------------------------------------------------------\n";
    }
}

void print_method_header_with_jacobi(const string& method_name)
{
    print_method_header_with_error(method_name, false);
    cout << "|" << setw(19) << "Норма матрицы Якоби\n";
    cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";
}

void print_method_header_with_additional_iterations(const std::string& method_name)
{
    print_method_header_with_error(method_name, false);
	cout << "|" << setw(16) << "Альфа"
		<< "|" << setw(6)  << "k" << '\n';
	cout << "-----------------------------------------------------------------------------------------------------------------------------------------\n";
}

void print_method_iteration_info(const int iteration, const ld x, const ld y, const ld residual_error_norm, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2 , bool switch_to_new_line){
    cout << setw(6) << iteration;
    cout << "|" << setw(12) << setprecision(8) << x;
    cout << "|" << setw(12) << setprecision(8) << y;
    cout << "|" << setw(16) << setprecision(10) << residual_error_norm;
    cout << "|" << setw(16) << setprecision(10) << f1(x, y);
    cout << "|" << setw(16) << setprecision(10) << f2(x, y);
    if (switch_to_new_line) {
        cout << '\n';
    }
}

void print_method_iteration_info(const int iteration, const ld x, const ld y, const ld residual_error_norm, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2, const ld error, bool switch_to_new_line){
    print_method_iteration_info(iteration, x, y, residual_error_norm, f1, f2, false);
    cout << "|" << setw(19) << error;
    if (switch_to_new_line) {
        cout << '\n';
    }
}

void print_method_iteration_info(const int iteration, const ld x, const ld y, const ld residual_error_norm, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2 , const ld error, const ld jacobi_matrix_norm){
    print_method_iteration_info(iteration, x, y, residual_error_norm, f1, f2, error, false);
    cout << "|" << setw(19) << jacobi_matrix_norm << '\n';
}

void print_method_iteration_info(const int iteration, const ld x, const ld y, const ld residual_error_norm, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2 , const ld error, const ld alpha, const int additional_iterations_count){
    print_method_iteration_info(iteration, x, y, residual_error_norm, f1, f2, error, false);
    cout << "|" << setw(16) << alpha
        << "|" << setw(6) << additional_iterations_count << '\n';
}
