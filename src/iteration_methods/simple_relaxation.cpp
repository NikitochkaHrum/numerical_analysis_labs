#include "simple_relaxation.h"

Vec simple_relaxation(Mat a, Vec b, Vec init, Vec solution, ld EPS){
    // Найдем оптимальное значение омега
    ld omega = 0.1, bst_omega = 0, m_eps = EPS * 100.; // Искать будем для точности эпсилон * 100
    int cur_itr, min_itr = INT_MAX;
    print_omega_search_header();
    while (omega < 2.)
    {
        Vec tmp = get_solution_simple_relaxation(a, b, init, solution, m_eps, omega, cur_itr, false);
        if (cur_itr < min_itr) {
            min_itr = cur_itr;
            bst_omega = omega;
        }
        print_omega_search(omega, cur_itr);
        omega += 0.1;
    }

    print_found_omega(bst_omega, min_itr);

    return get_solution_simple_relaxation(a, b, init, solution, EPS, bst_omega, cur_itr, true);
}

Vec get_solution_simple_relaxation(Mat a, Vec b, Vec init, Vec solution, ld EPS, ld omega, int& iter_count, bool log){
    if (log)
		print_method_header("Метод ПВР", simple_relaxation_estimate_n_iterations(EPS, a.cond(3)));
	

    Vec * x3 = new Vec(init.data); // Задаем решению начальное приблежение
    Vec * x2 = nullptr, * x1 = nullptr; // Вектор x на предыдущей и две итерации назад

    ld q, residual_norm, relative_error, error_estimate;
    ld solution_norm = solution.euclidian_norm();
    iter_count = 0;
    do {
        iter_count++;

        delete x1;
        x1 = x2;
        x2 = x3;

        // Очередное приближение
        x3 = new Vec(x2->n);
        for (int i = 0; i < b.n; i++) {
            (*x3).data[i] = b[i];
            for (int j = i - 1; j >= 0; j--)
                (*x3).data[i] -= a[i][j] * (*x3)[j];
            for (int j = i + 1; j < b.n; j++)
                (*x3).data[i] -= a[i][j] * (*x2)[j];
            (*x3).data[i] /= a[i][i];
            (*x3).data[i] = (*x2)[i] + omega * ((*x3)[i] - (*x2)[i]);
        }

        tie(q, residual_norm, relative_error, error_estimate) = calc_params(x1, x2, x3, a, b, solution, solution_norm);
        if (log && iter_count >= 2)
			print_method_iteration_info(iter_count-1, omega, q, residual_norm, relative_error, error_estimate, x3);
    } while (relative_error >= EPS); // Критерий остановки — относительная погрешность решения

    delete x2; delete x1;
    return *x3;
}

int simple_relaxation_estimate_n_iterations(ld EPS, ld cond){
    return log(1 / EPS) / 4 * sqrt(cond);
}