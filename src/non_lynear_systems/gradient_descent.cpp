#include "gradient_descent.h"
#include "method_log/method_info.h"
#include "calc_params.h"

Vec gradient_descent(Vec init, Vec solution, function<ld(ld, ld)> f1, function<ld(ld, ld)> f2,
                        function<ld(ld, ld)> der_f1_x, function<ld(ld, ld)> der_f1_y,
                        function<ld(ld, ld)> der_f2_x, function<ld(ld, ld)> der_f2_y,
                        function<ld(ld, ld)> F,
                        function<ld(ld, ld)> der_F_x, function<ld(ld, ld)> der_F_y,
                        ld eps){
    const string method_name = " Метод наискорейшего спуска";
    ld a = eps, b = 1; // Отрезок, на котором будем искать alpha
    auto* x = new Vec(init.data); // Вектор решения
    auto* sol = new Vec(solution.data); // Вектор правильногоо ответа
    Vec residual_error(2); // Вектор невязки
    ld residual_error_norm; // Норма вектора невязки
    ld error_norm; // Норма погрешности на k-ом шаге
    ld error_from_solution_norm; // Норма погрешности по сравнению с ответом
    
    int additional_iterations_count = 1;

    print_method_header_with_additional_iterations(method_name);

    int iter = 0; // Номер итерации
    int sum_iters = 0;
    do {
        auto* x_next = new Vec(2); // Значения вектора решения на итерации iter + 1

        ld alpha;
        // Функция от alpha (Функция сервтки F = f^2 + g^2)
        function<ld(ld)> fi_alpha = [=, &f1, &f2, &der_F_x, &der_F_y](ld _alpha) {
            Vec _x(2);
            _x.data[0] = x->data[0] - _alpha * der_F_x(x->data[0], x->data[1]);
            _x.data[1] = x->data[1] - _alpha * der_F_y(x->data[0], x->data[1]);
            ld res = pow(f1(_x[0], _x[1]), 2) + pow(f2(_x[0], _x[1]), 2);
            return res;
        };

        // Находим минимум функции от alpha
        tie(alpha, additional_iterations_count) = golden_ratio(fi_alpha, a, b, eps * 100.);

        // Очередное приближение
        x_next->data[0] = x->data[0] - alpha * der_F_x(x->data[0], x->data[1]);
        x_next->data[1] = x->data[1] - alpha * der_F_y(x->data[0], x->data[1]);

		tie(error_norm, residual_error_norm, error_from_solution_norm, residual_error) = calc_params(x, x_next, sol, f1, f2);

        iter++;
        delete x;
        x = x_next;
        sum_iters += additional_iterations_count;
        print_method_iteration_info(iter, x->data[0], x->data[1], residual_error_norm, f1, f2, error_from_solution_norm, alpha, additional_iterations_count);
    } while(residual_error_norm >= eps || error_norm >= eps);
    cout << "Всего итераций: " << sum_iters << '\n';
    return *x;
}

pair<ld, int> golden_ratio(const std::function<ld(ld)> f, ld a, ld b, ld eps) {
    // Находим c' и f(c')

    int iterations_count = 0;

    ld c_left = gr_left(a, b);
    ld f_left = f(c_left);

    // Находим c и f(c)
    ld c_right = gr_right(a, b);
    ld f_right = f(c_right);

    while (abs(c_right - a) > eps) {
        iterations_count++;
        if (f_left < f_right) {
            b = c_right;
            c_right = c_left;
            f_right = f_left;
            c_left = gr_left(a, b);
            f_left = f(c_left);
        } else {
            a = c_left;
            c_left = c_right;
            f_left = f_right;
            c_right = gr_right(a, b);
            f_right = f(c_right);
        }
    }

    return { c_right, iterations_count };
}

inline ld gr_left(ld a, ld b) {
    return a + (b - a) * (3 - sqrt_5) / 2;
}

inline ld gr_right(ld a, ld b) {
    return a + (b - a) * (sqrt_5 - 1) / 2;
}