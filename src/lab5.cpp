#include "interpolation/input_data.h"
#include "interpolation/newton.h"
#include "interpolation/spline.h"
#include "interpolation/discrete_rms.h"
#include "interpolation/cont_rms.h"
#include "interpolation/uniform_approx.h"


void run_var(function<ld(ld)> f, function<ld(ld)> Der, function<ld(ld)> Der4, function<ld(ld)> Der5, function<ld(ld)> Der6){
    Mat val_table = table_of_values(f, N, A, B);
    Mat table = table_separated_diff(val_table);

    ld m_der6_f = max(Der6(A), Der6(B)), m_der4_f = max(Der4(A), Der4(B)), m_der5_f = max(Der5(A), Der5(B));

    newton_method(f, m_der6_f, val_table);

    cubic_spline_interpolation(f, Der, val_table, m_der4_f, m_der5_f);

    vector<function<ld(ld)>> basis;
    basis.push_back([](ld x) -> ld { return 1; });
    basis.push_back([](ld x) -> ld { return x; });
    basis.push_back([](ld x) -> ld { return x*x; });

    discrete_rms(basis, val_table);
    continuous_rms(f, basis, val_table);

    uniform_approximation_first_order_polynomial(f, Der, val_table);

    cout << "Решение уравнения методом обратной интерполяции\n";
    cout << "Таблица разделенных разностей\n";

    ld c = (val_table[N][1] + val_table[0][1]) / 2;
    ld root = root_of_equation(val_table, c);

    cout << '\n';

    cout << "c = " << c << '\n';
    cout << "Корень: " << root << '\n';
    cout << "Невязка: " << abs(f(root) - c) << '\n';
}

int main(){
    freopen("/home/pna/Documents/study/numerical_analysis_labs/output.txt", "w", stdout);

    cout << "Вариант: 9\n";
    run_var(F_var9, Der_var9, Der4_var9, Der5_var9, Der6_var9);
    cout << "Вариант: 25\n";
    run_var(F_var25, Der_var25, Der4_var25, Der5_var25, Der6_var25);

}