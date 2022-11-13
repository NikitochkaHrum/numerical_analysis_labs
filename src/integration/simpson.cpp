#include "simpson.h"
#include "method_log/method_info.h"

ld simpson_rule(FunctionTool& f, ld a, ld b, ld eps, ld exact_value) {
    int n = 2;
    ld h = (b - a) / n;
    ld sk = f(a) + f(b);
    ld s_even = 0;
    ld s_odd = 0;
    ld s3 = 0, s2 = 0, s1 = 0;

    print_header();
    do
    {
        s_even += s_odd;
        s_odd = 0;
        for (int i = n / 2; i > 0; i--)
            s_odd += f(a + (2 * i - 1) * h);

        s1 = s2;
        s2 = s3;
        s3 = h / 3 * (sk + 4 * s_odd + 2 * s_even);

        print_iteration(n / 2, h * 2, exact_value, s1, s2, s3, 4);

        n *= 2;
        h /= 2;
    } while (abs(s3 - exact_value) >= eps * abs(s3));

    return s3;
}