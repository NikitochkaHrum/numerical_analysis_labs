#include "heat_equation.h"
#include "method_log/method_info.h"
#include "../direct_methods/sweep.h"

void heat_eq_explicit_schema(function<ld(ld, ld)> &f, ld a, ld b, function<ld(ld)> &phi, ld hi, int n, ld end_t, function<ld(ld, ld)> &exact_solution){
    // Считаем нужные параметры
    ld h = (b - a) / n;
    ld h2 = pow(h, 2); // Запомним, пригодиться
    ld tau = h2 / (4 * hi);

    // Инициализируем сетку при t(0)
    Vec prev(n + 1), cur(n + 1);
    prev.data[0] = a;
    cur.data[0] = a;
    prev.data[n] = b;
    cur.data[n] = b;
    for (int i = 1; i < n; i++)
        prev.data[i] = phi(i * h);

    heat_eq_header(n, h, hi, prev);

    ld t; // Текущее время t(n), у нас t(i)
    ld max_err = 0; // Максимальная погрешность
    for (int i = 1; (t = i * tau) <= end_t; i++) {
        for (int j = 1; j < n; j++) {
            ld tmp = (prev[j + 1] - 2 * prev[j] + prev[j - 1]) / h2;
            ld x = j * h;
            cur.data[j] = prev[j] + tau * (hi * tmp + f(t - tau, x));
        }

        ld cur_err = 0; // Смотрим на максимальную погрешность для t(n)
        for (int j = 0; j < n; j++) {
            ld x = j * h;
            ld err = abs(cur[j] - exact_solution(t, x));
            if (err > cur_err)
                cur_err = err;
        }

        heat_eq_iter(cur, t, cur_err);

        if (cur_err > max_err)
            max_err = cur_err;

        prev = cur;
    }

    heat_eq_max_err(max_err);

}

void heat_eq_implicit_schema(function<ld(ld, ld)> &f, ld a, ld b, function<ld(ld)> &phi, ld hi, int n, ld end_t, function<ld(ld, ld)> &exact_solution){
    // Считаем нужные параметры
    ld h = (b - a) / n, tau = h;
    ld d = tau * hi / pow(h, 2);

    int mn = n - 1;
    Mat left(mn); // матрица коэф-в левой части
    Vec right(mn); // вектор свободных членов правой части

    for (int i = 0; i < mn; i++) {
        if (i != 0)
            left.data[i][i - 1] = -d;
        left.data[i][i] = 1 + 2 * d;
        if (i != mn - 1)
            left.data[i][i + 1] = -d;
    }

    // Инициализируем сетку при t(0)
    Vec prev(n + 1), cur(n + 1);
    prev.data[0] = a;
    cur.data[0] = a;
    prev.data[n] = b;
    cur.data[n] = b;
    for (int i = 1; i < n; i++)
        prev.data[i] = phi(i * h);

    heat_eq_header(n, h, hi, prev);

    ld t; // Текущее время t(n), у нас t(i)
    ld max_err = 0; // Максимальная погрешность
    for (int i = 1; (t = i * tau) <= end_t; i++) {
        for (int j = 0; j < mn; j++)
            right.data[j] = prev[j + 1] + tau * f(t, h * (j + 1));

        right.data[0] += d * a;
        right.data[mn - 1] += d * b;

        prev = sweep_method(left, right); // Решаем СЛАУ

        for (int j = 0; j < mn; j++)
            cur.data[j + 1] = prev[j];

        ld cur_err = 0; // Смотрим на максимальную погрешность для t(n)
        for (int j = 0; j < n; j++) {
            ld x = j * h;
            ld err = abs(cur[j] - exact_solution(t, x));
            if (err > cur_err)
                cur_err = err;
        }

        heat_eq_iter(cur, t, cur_err);

        if (cur_err > max_err)
            max_err = cur_err;

        prev = cur;
    }

    heat_eq_max_err(max_err);
}