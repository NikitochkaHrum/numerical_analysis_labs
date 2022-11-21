#include "newton.h"


Mat table_separated_diff(const Mat& table_of_values) {
    int n = table_of_values.n - 1;
    Mat table(table_of_values.n, table_of_values.m + n);
    for (int i = 0; i < table_of_values.n; i++)
        for (int j = 0; j < table_of_values.m; j++)
            table.data[i][j] = table_of_values.data[i][j];

    // Индекс первого столбца разделенных разностей и их количество в каждом столбце
    int offset = table_of_values.m, rows_in_col = n;
    for (int j = 0; j < n; j++) {
        int col_ind = offset + j;
        for (int i = 1; i <= rows_in_col; i++)
            table.data[i - 1][col_ind] = (table.data[i][col_ind - 1] - table[i - 1][col_ind - 1]) / (table[i + j][0] - table[i - 1][0]);
        rows_in_col--;
    }
    return table;
}

ld newton_method(const Mat& table, ld x, ld m, ld& error_est) {
    int n = table.n - 1;
    ld omega = x - table.data[n][0]; // омега для оценки погрешности
    // Схема Горнера
    ld result = table.data[0][table.m- 1]; // P(k = 0)
    for (int k = 1; k <= n; k++) {
        ld tmp = (x - table.data[n - k][0]); // x - x(n - k)
        result = table.data[0][table.m - k - 1] + tmp * result;
        // Сразу считаем омега
        omega *= tmp;
    }
    // Оценка погрешности
    error_est = m * abs(omega);
    ld fact = 1;
    for (int i = n + 1; i > 1; i--)
        fact *= i;
    error_est = error_est / fact;
    return result;
}

ld root_of_equation(const Mat& table_of_values, ld c) {
    int n = table_of_values.n - 1;
    Mat rev_table_of_values(table_of_values.n, table_of_values.m);
    for (int i = 0; i <= n; i++) {
        ld tmp = table_of_values.data[i][1] - c;
        rev_table_of_values.data[i][1] = table_of_values.data[i][0];
        rev_table_of_values.data[i][0] = tmp;
    }

    Mat table = table_separated_diff(rev_table_of_values);
    table.print();

    ld err;
    return newton_method(table, 0, 0, err);
}

void newton_method(function<ld(ld)> f, ld m, const Mat& values)
{
    const int N = values.n - 1;
    const ld A = values.data[0][0];
    const ld B = values.data[N][0];

    ld h = (B - A) / N, x_ = A + h / 2;

    Mat table_diffs = table_separated_diff(values);

    Vec errors(N);
    Vec xs(N);

    for (int i = 0; i < N; i++) {
        xs.data[i] = x_ + h * i;
        ld x = xs[i];
        ld f_newton_ = newton_method(table_diffs, x, m, errors.data[i]);
        ld f_ = values.data[i][1];
    }
    ld err;// KocTblJlU
    auto g = [&](ld x) -> ld { return newton_method(table_diffs, x, m, err); };

    cout << "Интерполяционная формула Ньютона" << endl;
    cout << "Таблица разделенных разностей" << endl;

    table_diffs.print();
    cout << '\n';

    cout << "M6 = " << fixed << m << '\n';
    print_interpolation_results(xs, f, g, errors);

    cout << '\n';

    // Пример вывода
    /*
    Интерполяционная формула Ньютона
 Таблица разделенных разностей
 1.00  6.000000  8.685964  2.264389  0.927384  0.284859  0.069999
 1.20  7.737193  9.591720  2.820819  1.155271  0.354857
 1.40  9.655537 10.720047  3.513981  1.439157
 1.60 11.799546 12.125640  4.377475
 1.80 14.224674 13.876630
 2.00 17.000000

 M6=15.8237431
   x    f(x)       Pn(x)           Delta             Оценка
 1.10  6.84837   6.84838  1.13639971166535E-05  2.07686628835692E-05
 1.30  8.67117   8.67116  3.90437515918052E-06  6.92288762785639E-06
 1.50 10.69615  10.69616  2.87656676434267E-06  4.94491973418315E-06
 1.70 12.97301  12.97300  4.15693947886098E-06  6.9228876278564E-06
 1.90 15.56363  15.56364  1.2882429841099E-05   2.07686628835692E-05
    */
}