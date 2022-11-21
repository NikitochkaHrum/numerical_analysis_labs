#include "solve_nle.h"

ld solve_nle(const std::function<ld(ld)>& f, ld a, ld b, ld eps){
    ld c;
    while ((b - a) / 2 > eps) {
        c = (a + b) / 2;
        if ((f(a) * f(c)) > 0) a = c;
        else b = c;
    }
    return c;
}
