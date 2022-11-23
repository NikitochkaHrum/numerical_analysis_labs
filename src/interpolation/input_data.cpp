#include "input_data.h"


ld F_var9(ld x){
    return exp(x) + 6 * x + 3;
}

ld Der_var9 (ld x) {
    return exp(x) + 6;
}

ld Der4_var9(ld x) {
    return exp(x);
}

ld Der5_var9(ld x) {
    return exp(x);
}

ld Der6_var9(ld x) {
    return exp(x);
}

ld F_var25(ld x) {
    return pow(3, x) - x + 2;
}

ld Der_var25(ld x) {
    return pow(3, x) * log(3) - 1;
}

ld Der4_var25(ld x) {
    return pow(3, x) * pow(log(3), 4);
}

ld Der5_var25(ld x) {
    return pow(3, x) * pow(log(3), 5);
}

ld Der6_var25(ld x){
    return pow(3, x) * pow(log(3), 6);
}

ld omega(ld x, vector<ld>& points) {
	ld result = 1;
	for (ld xi : points)
		result *= (x - xi);

	return result;
}

Mat table_of_values(function<ld(ld)> f, int n, ld a, ld b) {
    Mat table(n + 1, 2);

    ld h = (b - a) / n;
    for (int i = 0; i <= n; i++) {
        table.data[i][0] = a + i * h;
        table.data[i][1] = f(table[i][0]);
    }

    return table;
}
