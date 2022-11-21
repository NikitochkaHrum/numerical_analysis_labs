#include "input_data.h"

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
