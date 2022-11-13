#include "gauss.h"
#include "method_log/method_info.h"

// Трехточечная квадратурная формула Гаусса
ld gauss(FunctionTool& f, ld a, ld b, ld eps, ld exact_value){
    // Константа корень 3 / 5
    const ld q = sqrt(3. / 5.);

    int n = 2;
    ld h = (b - a) / n;
    ld s3 = 0, s2 = 0, s1 = 0;

    print_header();

    do {
        ld tmp_s0 = 0, tmp_s1 = 0, tmp_s2 = 0;
        for (int i = n / 2; i > 0; i--) {
            tmp_s0 += f(a + (2 * i - 1) * h - q * h);
            tmp_s1 += f(a + (2 * i - 1) * h);
            tmp_s2 += f(a + (2 * i - 1) * h + q * h);
        }

        s1 = s2;
        s2 = s3;
        s3 = (h / 9) * (5 * (tmp_s0 + tmp_s2) + 8 * tmp_s1);

        print_iteration(n / 2, h * 2, exact_value, s1, s2, s3, 6);
        n *= 2;
        h /= 2;
    } while (abs(s3 - exact_value) >= eps * abs(s3) || n<=8);
    return s3;
}
