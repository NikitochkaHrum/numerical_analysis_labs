#include "method_info.h"

void print_method_header(const string& method_name, const int iteration_cnt){
    cout << method_name << " | " << "Теоретическая оценка кол-ва итераций: " << iteration_cnt << '\n';
	cout << "|" << setw(6) << "Itr"
		 << "|" << setw(6) << "Tau"
		 << "|" << setw(6) << "q"
		 << "|" << setw(16) << "Норма невязки"
		 << "|" << setw(20) << "Норма погрешности"
		 << "|" << setw(22) << "Оценка погрешности"
		 << "|" << setw(10) << "x[1]"
		 << " " << setw(10) << "x[2]"
		 << " " << setw(10) << "x[3]"
		 << " " << setw(10) << "x[4]" << endl;
	cout << "------------------------------------------------------------------------------------------------------------------------------\n";
}

void print_method_iteration_info(const int iteration, const ld tau, const ld q, const ld norm_a, const ld norm_error, const ld estimate_error, Vec* x){
    cout << "|" << setw(6) << iteration << "|"
         << setw(6) << fixed << setprecision(4) << tau << "|"
		 << setw(6) << fixed << setprecision(4) << q << "|"
		 << setw(16) << fixed << setprecision(6) << norm_a << "|"
		 << setw(20) << fixed << setprecision(6) << norm_error << "|"
		 << setw(22) << fixed << setprecision(6) << estimate_error << "|"
		 << setw(10) << fixed << setprecision(6) << (*x)[0] << " "
		 << setw(10) << fixed << setprecision(6) << (*x)[1] << " "
		 << setw(10) << fixed << setprecision(6) << (*x)[2] << " "
		 << setw(10) << fixed << setprecision(6) << (*x)[3] << " \n";
}

void print_omega_search_header(){
    cout << "Выбор оптимального w для метода ПВР:\n";
}

void print_omega_search(ld omega, int iterations_cnt){
    cout << "w = " << omega << ", количество итераций = "<< iterations_cnt << '\n';
}

void print_found_omega(ld omega, int iterations_cnt){
    cout << "Оптимальное значение w* = " << omega << ", количество итераций = "<< iterations_cnt << '\n';
}