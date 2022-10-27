#include "non_lynear_systems/newton_method.h"
#include "non_lynear_systems/simple_iteration.h"
#include "non_lynear_systems/gradient_descent.h"
#include "non_lynear_systems/input_data.h"

int main(){
    freopen("/home/pna/Documents/study/numerical_analysis_labs/output.txt", "w", stdout);
    
    const ld EPS1 = 1e-12;
    const ld EPS2 = 1e-4;
    // Vec init({-0.95, -1.1}); // v9
    //Vec init({-0.15, -0.85}); // v25
    Vec init({-1.7, 1.3}); // vovavar
    cout << "Начальное приближение:\nx = " << init[0] << "\ny = " << init[1] << '\n';

    auto x = newton_method(init, f1_v20, f2_v20, der_f1_x_v20, der_f1_y_v20, der_f2_x_v20, der_f2_y_v20, EPS1);
    cout << fixed << setprecision(10) << "x = " << x[0] << "\ny = " << x[1] << '\n';

    auto x_it = simple_iteration(init, x, f1_v20, f2_v20, fi1_v20, fi2_v20, der_fi1_x_v20, der_fi1_y_v20, der_fi2_x_v20, der_fi2_y_v20, EPS2);
    cout << fixed << setprecision(10) << "x = " << x_it[0] << "\ny = " << x_it[1] << '\n';

    x_it = gradient_descent(init, x, f1_v20, f2_v20, der_f1_x_v20, der_f1_y_v20, der_f2_x_v20, der_f2_y_v20, F_v20, der_F_x_v20, der_F_y_v20, EPS2);
    cout << fixed << setprecision(10) << "x = " << x_it[0] << "\ny = " << x_it[1] << '\n';
}