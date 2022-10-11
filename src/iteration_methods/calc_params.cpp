#include "calc_params.h"

tuple<ld, ld, ld, ld> calc_params(Vec* x1, Vec* x2, Vec* x3, Mat& a, Vec& b, Vec& solution, ld solution_norm){
    // || x(k+1) - x(k) ||
    ld iter_rel_error = ((*x3) - (*x2)).euclidian_norm();
    // || x(k) - x(k-1) ||
    ld prev_iter_rel_error = x1 && x2 ? ((*x2) - (*x1)).euclidian_norm() : 1;
    // Оценка нормы матрицы переходов
    ld transition_matrix_estimate = iter_rel_error / prev_iter_rel_error;
    // Оценка погрешность решения
    ld error_estimate = transition_matrix_estimate * prev_iter_rel_error / (1 - transition_matrix_estimate);

    // Относительная погрешность решения
    ld relative_error = (solution - (*x3)).euclidian_norm() / solution_norm;
    // Норма вектора невязки
    double residual_norm = (a * (*x3) - b).euclidian_norm();

    return { transition_matrix_estimate, residual_norm, relative_error, error_estimate };
}