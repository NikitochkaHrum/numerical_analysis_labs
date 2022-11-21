#include "sweep.h"

Vec sweep_method(Mat& A, Vec& B) {
	int n = A.n;

	Vec x(n);

	Vec a(n), c(n), b(n);

	for (int i = 0; i < n; i++) {
		if (i - 1 >= 0)
			a.data[i] = -A[i][i - 1];
		else
			a.data[i] = 0;

		c.data[i] = A.data[i][i];

		if (i + 1 < n)
			b.data[i] = -A.data[i][i + 1];
		else
			b.data[i] = 0;
	}

	Vec alphas(n), betas(n);

	alphas.data[1] = b[0] / c[0];
	betas.data[1] = B[0] / c[0];

	for (int i = 1; i < n-1; i++) {
		alphas.data[i + 1] = b[i] / (c[i] - a[i] * alphas[i]);
		betas.data[i + 1]  = (B[i] + a[i] * betas[i]) / (c[i] - a[i] * alphas[i]);
	}

	x.data[n - 1] = (B[n - 1] + a[n - 1] * betas[n - 1]) / (c[n - 1] - a[n - 1] * alphas[n - 1]);

	for (int i = n - 2; i >= 0; i--) {
		x.data[i] = alphas[i + 1] * x[i + 1] + betas[i + 1];
	}

	return x;
}
