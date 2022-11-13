#include "method_log/method_info.h"
#include "trapez_spline.h"

ld trapezoidal_rule_spline(FunctionTool & f, const function<ld(ld)> & der_f, ld a, ld b, ld eps, ld exact_value)
{
    // sum = sum(j = 0..n-1) (f(x(j))*1/2 + f(x(j+1))*1/2 + h*m(j)*1/12 - h * m(j+1) * 1/12)
    // sum = 1/2 * (f(x(j)) + f(x(j+1)) + h/12 * (m(j) - m(j+1))
    // sum = sum(j = 0..n-1) ( 1/2 * (f(x(j) + f(x(j+1)))    + h/12 * m(0) - h/12 * m(n)
    // sum = sum(j = 1..n-1) (f(x(j)))  +    1/2 * (f(0) + f(n))    + h/12 * (der_f(a) - der_f(b)) - получили такую формулу для оптимальных вычислений
    int n = 1; // количество отрезков
    ld h = (b - a) / n; // длина отрезка
    ld sk = (f(a) + f(b)) / 2; // накапливаемая сумма
    ld der_sk = (der_f(a) - der_f(b)) / 12; // сумма с производными
    ld s3 = 0, s2 = 0, s1 = 0; // предыдущие значения сумм

    print_header();

    do {
        for (int i = n - 1; i > 0; i -= 2) //  Каждый раз n := 2n => лишняя проверка на четность не нужна
            sk += f(a + i * h);

        s1 = s2;
        s2 = s3;
        s3 = h * sk + h * h * der_sk;

        print_iteration(n, h, exact_value, s1, s2, s3, 4);

        n *= 2;
        h /= 2;
    } while (abs(s3 - exact_value) >= eps * abs(s3));

    return s3;
}
