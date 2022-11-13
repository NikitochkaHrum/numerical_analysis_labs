#include "method_info.h"

// Оценка погрешности по правилу Рунге
ld error_estimation(ld s1, ld s2, ld h, int k) {
    return (s2 - s1) * pow(alpha, k) / (1 - pow(alpha, k));
    return ((s2 - s1) / (pow(h / alpha, k) - pow(h, k))) * pow(h, 2);
}

// Оценка порядка погрешности
ld degree_estimation(ld s1, ld s2, ld s3) {
    return d * log((s3 - s1) / (s2 - s1) - 1);
}

void print_header() {
    //    N   h   I   Погрешность   Оценка погрешности   k
    cout << setw(6)  << "N" << " |"
        << setw(12) << "h" << " |"
        << setw(16) << "I" << " |"
        << setw(24) << "             Погрешность" << " |"
        << setw(24) << "      Оценка погрешности" << " |"
        << setw(12) << "    k\n";
}

// Вывод информации об очередной итерации
void print_iteration(int n, ld h, ld s, ld s1, ld s2, ld s3, int k) {
    const int integral_precision = 10, precision = 8;
    auto old_precision = cout.precision();

    cout << setw(6) << n << " |";

    cout << setprecision(precision) << fixed
        << setw(12) << h << " |";

    cout << setprecision(integral_precision) << fixed
        << setw(16) << s3 << " |";

    cout << resetiosflags(ios_base::fixed);
    cout.precision(old_precision);

    cout << scientific;
    cout << setw(24) << abs(s3 - s) << " |"
        << setw(24);

    if (s2)
        cout << error_estimation(s2, s3, h, k);
    else
        cout << " ";

    cout << " |";
    cout << resetiosflags(ios_base::scientific);

    cout << setprecision(precision) << fixed
        << setw(12) << (s1 ? to_string(degree_estimation(s1, s2, s3)) : " ") << '\n';
    cout << resetiosflags(ios_base::fixed);
    cout.precision(old_precision);
}

