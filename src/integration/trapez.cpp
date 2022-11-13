#include "method_log/method_info.h"
#include "trapez.h"

ld trapezoidal_rule(FunctionTool& f, ld a, ld b, ld eps, ld exact_value){
    int n = 1; // количество отрезков
    ld h = (b - a) / n; // длина отрезка
    ld sk = (f(a) + f(b)) / 2; // накапливаемая сумма
    ld s3 = 0, s2 = 0, s1 = 0; // предыдущие значения сумм

    print_header();
    do
    {
        for (int i = n - 1; i > 0; i -= 2) // Каждый раз n := 2n => Лишняя проверка на четность не нужна
            sk += f(a + i * h);

        s1 = s2;
        s2 = s3;
        s3 = h * sk;

        print_iteration(n, h, exact_value, s1, s2, s3, 2);

        n *= 2;
        h /= 2;
    } while (abs(s3 - exact_value) >= eps * abs(s3));

    return s3;
}