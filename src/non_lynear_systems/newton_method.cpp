#include "newton_method.h"
#include "method_log/method_info.h"
#include "calc_params.h"

Vec newton_method(Vec init, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2,
                        function<ld(ld, ld)> der_f1_x, function<ld(ld, ld)> der_f1_y, function<ld(ld, ld)> der_f2_x, function<ld(ld, ld)> der_f2_y,
                        ld eps){
	const string method_name = " Метод Ньютона"; // Название метода

	auto* x = new Vec(init.data); // Вектор решения 

	Mat values_inv_jacobian = get_matrix_of_inversed_jacobian(der_f1_x, der_f1_y, der_f2_x, der_f2_y, *x);// Матрица значений обрат. Якобиана при опред. X
	Vec values_function = Vec(2);
	values_function.data[0] = f1(x->data[0], x->data[1]);
	values_function.data[1] = f2(x->data[0], x->data[1]);


    Mat values_jacobian(2, 2); // Матрица Якоби
    values_jacobian.data[0][0] = der_f1_x(x->data[0], x->data[1]);
    values_jacobian.data[0][1] = der_f1_y(x->data[0], x->data[1]);
    values_jacobian.data[1][0] = der_f2_x(x->data[0], x->data[1]);
    values_jacobian.data[1][1] = der_f2_y(x->data[0], x->data[1]);

	//double norm_jacobian = values_jacobian.euclideanNorm(); // Якобиан

	Vec residual_error(2); // Вектор невязки
	ld residual_error_norm; // Норма вектора невязки
	ld error_norm; // Норма погрешности на k-ом шаге

	print_method_header(method_name, true);

	int iter = 0; // Номер итерации
	do {
		// Очередное приближение
		Vec* x_next = new Vec(2);
		*x_next = *x - values_inv_jacobian * values_function;
		
		tie(error_norm, residual_error_norm, residual_error) = calc_params(x, x_next, f1, f2);

		delete x;
		x = x_next;
		
		values_inv_jacobian = get_matrix_of_inversed_jacobian(der_f1_x, der_f1_y, der_f2_x, der_f2_y, *x);
        values_function.data[0] = f1(x->data[0], x->data[1]);
        values_function.data[1] = f2(x->data[0], x->data[1]);
        
		// values_jacobian.data[0][0] = der_f1_x(x->data[0], x->data[1]);
        // values_jacobian.data[0][1] = der_f1_y(x->data[0], x->data[1]);
        // values_jacobian.data[1][0] = der_f2_x(x->data[0], x->data[1]);
        // values_jacobian.data[1][1] = der_f2_y(x->data[0], x->data[1]);
		//norm_jacobian = values_jacobian.euclideanNorm();

		iter++;

		print_method_iteration_info(iter, x->data[0], x->data[1], residual_error_norm, f1, f2, true);
	} while (residual_error_norm >= eps || iter < 3);

	return *x;
}

Mat get_matrix_of_inversed_jacobian(function<ld(ld, ld)> der_f1_x, function<ld(ld, ld)> der_f1_y, function<ld(ld, ld)> der_f2_x, function<ld(ld, ld)> der_f2_y,
                                        Vec vals){
    Mat r(2, 2);
    r.data[0][0] = der_f1_x(vals[0], vals[1]);
    r.data[0][1] = der_f1_y(vals[0], vals[1]);
    r.data[1][0] = der_f2_x(vals[0], vals[1]);
    r.data[1][1] = der_f2_y(vals[0], vals[1]);
    
	ld det = 1. / (r[0][0] * r[1][1] - r[0][1] * r[1][0]);

	swap(r.data[0][0], r.data[1][1]);
	r.data[0][1] = -r[0][1];
	r.data[1][0] = -r[1][0];
	
    r.data[0][0] *=det;
    r.data[0][1] *=det;
    r.data[1][0] *=det;
    r.data[1][1] *=det;
	

	return r;
}
