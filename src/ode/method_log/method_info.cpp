#include "method_info.h"

using namespace std;

void heat_eq_header(int n, ld h, ld hi, Vec &s) {
    cout << "Метод конечных разностей (явная схема)\n"
         << "N = " << n
         << "\nH = " << h
         << "\nHi = " << hi << '\n';

    //    t    delta       x:
    cout << setw(16) << "              t|"
         << setw(22) << "            delta|";

    cout << " x: ";
    s.print();
}

void heat_eq_iter(Vec &u, ld t, ld err) {
    cout << setw(10) << setprecision(6) << t << "   |"
        << setw(20) << setprecision(12) << err  << " | ";

    cout << "x: ";
    u.print();
}

void heat_eq_max_err(ld err) {
    cout << "Max error = " << err << "\n\n";
}

void log_shooting_method_iteration(int iter, ld alpha, ld y, ld B) {
    cout << setw(8) << iter << " | "
        << setw(12) << setprecision(6) << alpha << " | "
        << setw(12) << setprecision(6) << y << " | "
        << setw(20) << setprecision(12) << abs(y-B) << " | \n";
}

void log_runge_kutta_header() {
    cout << "\nМетод Рунге-Кутты\n";
    cout << setw(8) << "x" << " | "
        << setw(12) << "y(x)" << " | "
        << setw(12) << "y_exact(x)" << " | "
        << setw(14) << "z(x)" << " | "
        << setw(20) << "Delta" << " | \n";
}

void log_runge_kutta_iteration(ld x, ld y, ld ex_y, ld z) {
    cout << setw(8) << setprecision(6) << x << " | "
         << setw(12) << setprecision(8) << y << " | "
         << setw(12) << setprecision(8) << ex_y << " | "
         << setw(14) << setprecision(8) << z << " | "
         << setw(20) << setprecision(12) << abs(y-ex_y) << " | \n";
}

void log_shooting_method_header() {
    cout << "Метод стрельб\n";
    cout << setw(8) << "Iter" << " | "
         << setw(12) << "Alpha" << " | "
         << setw(12) << "y(1;alpha)" << " | "
         << setw(20) << "Delta" << " | \n";
}