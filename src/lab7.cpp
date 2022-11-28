#include "ode/input_data.h"
#include "ode/runge_kutt.h"
#include "ode/shooting.h"
#include "ode/heat_equation.h"

const ld _A = 1., _B = 2.;
const ld x_left = 0., x_right = 1.;
const ld h = 0.1, h_eps = 1e-5;
const ld y_eps = 1e-4;


void run_var(function<ld(ld)> &Y_r, function<ld(ld, ld)> &U, function<ld(ld, ld)> &FU, function<ld(ld, ld, ld)> &Z){
    ld alpha = shooting_method(Z, Y_r, _A, _B, x_left, x_right, h, h_eps, y_eps);
    runge_kutta_method(Z, Y_r, x_left, x_right, _A, alpha, h, h_eps);

    heat_eq_explicit_schema(FU, 0, 1, PHI, HI, 8, 1, U);
    heat_eq_explicit_schema(FU, 0, 1, PHI, HI, 16, 1, U);
    heat_eq_explicit_schema(FU, 0, 1, PHI, HI, 32, 1, U);
    cout << '\n';
    heat_eq_implicit_schema(FU, 0, 1, PHI, HI, 8, 1, U);
    heat_eq_implicit_schema(FU, 0, 1, PHI, HI, 16, 1, U);
    heat_eq_implicit_schema(FU, 0, 1, PHI, HI, 32, 1, U);
    cout << "\n\n";
}

int main(){
    freopen("/home/pna/Documents/study/numerical_analysis_labs/output.txt", "w", stdout);
    cout << "Вариант: 9\n";
    run_var(Y_9, U_9, FU_9, Z_9);
    cout << "Вариант: 25\n";
    run_var(Y_25, U_25, FU_25, Z_25);
}