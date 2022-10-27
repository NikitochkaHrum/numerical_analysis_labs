#include "simple_iteration.h"
#include "method_log/method_info.h"
#include "calc_params.h"

Vec simple_iteration(Vec init, Vec solution, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2,
                        function<ld(ld, ld)> fi1, function<ld(ld, ld)> fi2,
                        function<ld(ld, ld)> der_fi1_x, function<ld(ld, ld)> der_fi1_y,
                        function<ld(ld, ld)> der_fi2_x, function<ld(ld, ld)> der_fi2_y,
                        ld eps){
    const string method_name = " Метод простой итерации"; // Название метода
    auto* x = new Vec(init.data); // Вектор решения 
    auto* sol = new Vec(solution.data); // Вектор правильногоо ответа
    Mat values_jacobian(2, 2); // Матрица Якоби
    values_jacobian.data[0][0] = der_fi1_x(x->data[0], x->data[1]);
    values_jacobian.data[0][1] = der_fi1_y(x->data[0], x->data[1]);
    values_jacobian.data[1][0] = der_fi2_x(x->data[0], x->data[1]);
    values_jacobian.data[1][1] = der_fi2_y(x->data[0], x->data[1]);

    ld j = values_jacobian.matrixnormeuk(); // Якобиан

    Vec residual_error(2); // Вектор невязки
    ld residual_error_norm; // Норма вектора невязки
    ld error_norm; // Норма погрешности на k-ом шаге
    ld error_from_solution_norm; // Норма погрешности по сравнению с ответом
    int iter = 0; // Номер итерации

    // Вывод шапки метода
    print_method_header_with_jacobi(method_name);

    do {
        // Очередное приближение
		Vec* x_next = new Vec(2);
		x_next->data[0] = fi1(x->data[0], x->data[1]);
		x_next->data[1] = fi2(x->data[0], x->data[1]);

		tie(error_norm, residual_error_norm, error_from_solution_norm, residual_error) = calc_params(x, x_next, sol, f1, f2);

        delete x;
        x = x_next;

        values_jacobian.data[0][0] = der_fi1_x(x->data[0], x->data[1]);
        values_jacobian.data[0][1] = der_fi1_y(x->data[0], x->data[1]);
        values_jacobian.data[1][0] = der_fi2_x(x->data[0], x->data[1]);
        values_jacobian.data[1][1] = der_fi2_y(x->data[0], x->data[1]);

        j = values_jacobian.matrixnormeuk();

        iter++;

        print_method_iteration_info(iter, x->data[0], x->data[1], residual_error_norm, f1, f2, error_from_solution_norm, j);
    } while (residual_error_norm >= eps);

    return *x;
}