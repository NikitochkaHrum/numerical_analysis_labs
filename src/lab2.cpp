#include <bits/stdc++.h>
#include "matrices.h"

int main() {
    freopen("/home/pna/Documents/study/numerical_analysis_labs/src/input.txt", "r", stdin);
    freopen("/home/pna/Documents/study/numerical_analysis_labs/src/output.txt", "w", stdout);
    int n;
    cin >> n;

    vector<vector<ld>> a(n, vector<ld>(n));
    vector<ld> X_correct(n);

    for (int i = 0; i < n; i++) {
        X_correct[i] = i + 1;
    }
    cout << "X:\n";
    for (int i = 0; i < n; i++) cout << X_correct[i] << ' ';
    cout << "\n\n";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) cin >> a[i][j];
    Mat aa(a);

    aa.calc_LUP();
    cout << "Final P:\n";
    for (int i = 0; i < n; i++) cout << aa.P[i] << ' ';
    cout << '\n';
    cout << "Final L:\n";
    Mat(aa.L).print();
    cout << "Final U:\n";;
    Mat(aa.U).print();
    cout << "\nLU - PA:\n";
    (Mat(aa.L) * Mat(aa.U) - aa.permutted_data()).print();

    cout << "\nDeterminant = " << aa.determinant() << '\n';
    Vec b = aa * Vec(X_correct);
    Vec x = solve_equation(aa, Vec(b));
    cout << "\nB:\n";
    Vec(b).print();
    cout << "\nX:\n";
    x.print();
    cout << "\nAX - B:\n";
    (aa * x - Vec(b)).print();
    cout << "\nA^(-1):\n";
    aa.inverse().print();
    auto a_ainv = aa * aa.inverse();
    for (int i = 0; i < n; i++) a_ainv.data[i][i]--;
    cout << "\nA*A^(-1)-E:\n";
    a_ainv.print();

    // Числа обусловленности
    cout << "\nCond numbers:\n\tFirst norm:\n\t";
    cout << aa.cond(1) << '\n';
    cout << "\tSecond norm:\n\t";
    cout << aa.cond(2) << '\n';
    cout << "\tEuclidian norm:\n\t";
    cout << aa.cond(3) << '\n';

    cout << "Uncertainity:\n"; //погрешность
    (Vec(X_correct) - x).print();
}