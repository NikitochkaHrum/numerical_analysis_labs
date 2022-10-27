#include "calc_params.h"

tuple<ld, ld, Vec> calc_params(Vec* x, Vec* x_next, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2){

    // разница между ответами
	ld error_norm = ((*x_next) - (*x)).euclidian_norm();
	Vec residual_error(2);

    // Невязка
	residual_error.data[0] = f1(x_next->data[0], x_next->data[1]);
	residual_error.data[1] = f2(x_next->data[0], x_next->data[1]);
	
	ld residual_error_norm = residual_error.euclidian_norm();

	return { error_norm, residual_error_norm, residual_error };
}

tuple<ld, ld, ld, Vec> calc_params(Vec* x, Vec* x_next, Vec* solution, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2){

    // разница между ответами
	ld error_norm = ((*x_next) - (*x)).euclidian_norm();

	// Погрешность ответа
	ld error_from_solution_norm = ((*x_next) - (*solution)).euclidian_norm();

    // Невязка
	Vec residual_error(2);
	residual_error.data[0] = f1(x_next->data[0], x_next->data[1]);
	residual_error.data[1] = f2(x_next->data[0], x_next->data[1]);
	
	ld residual_error_norm = residual_error.euclidian_norm();

	return { error_norm, residual_error_norm, error_from_solution_norm, residual_error };
}