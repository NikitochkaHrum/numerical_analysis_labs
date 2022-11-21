#include "cont_rms.h"
#include "method_log/method_info.h"

// вычисление значения интеграла по формуле трапеций, разбивая отрезок на N равных частей
ld integral_value(function<ld(ld)> f, ld a, ld b, const int n) {
	ld h = (b - a) / n;

	ld result = 0;

	for (int i = 0; i < n; i++)
	{
		ld l = a + i * h;
		ld r = l + h;

		result += h * (f(l) + f(r)) / 2;
	}

	return result;
}

void continuous_rms(function<ld(ld)> f, vector<function<ld(ld)>>& basis, Mat& values) {
	const int N = values.n - 1;
	const ld A = values[0][0];
	const ld B = values[N][0];

	const int M = basis.size();
	const int COUNT_SEGMENT_SPLIT = 10000;

	Mat basis_values(M);

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			auto product_basis_function = [&](ld x) -> ld { return basis[i](x) * basis[j](x); };
			basis_values.data[i][j] = integral_value(product_basis_function, A, B, COUNT_SEGMENT_SPLIT);
		}
	}

	Vec scalar_with_f(M);
	for (int i = 0; i < M; i++) {
		for (int k = 0; k <= N; k++) {
			auto product_basis_function = [&](ld x) -> ld { return basis[i](x) * f(x); };
			scalar_with_f.data[i] = integral_value(product_basis_function, A, B, COUNT_SEGMENT_SPLIT);
		}
	}


	cout << integral_value([](ld x) -> ld { return x * x; }, 0, 5, 100) << endl;

	Vec c = solve_equation(basis_values, scalar_with_f);

	cout << "Непрерывный вариант\n";

	cout << "Матрица коэффициентов в правой части\n";
	basis_values.print();
	cout << '\n';

	cout << "Вектор правой части\n";
	scalar_with_f.print();
	cout << '\n';

	print_polynom(c);

	auto g = [&](ld x) -> ld { 
		return c[0] + c[1] * x + c[2] * x * x;
	};

	auto product_f = [&](ld x) -> ld {
		return f(x) * f(x);
	};

	auto product_g = [&](ld x) -> ld {
		return g(x) * g(x);
	};

	ld norm_f = integral_value(product_f, A, B, COUNT_SEGMENT_SPLIT),
		   norm_g = integral_value(product_g, A, B, COUNT_SEGMENT_SPLIT);

    cout << "Норма погрешности: " << sqrt(norm_f - norm_g) << "\n\n";

	print_error(g, values);

	cout << '\n';
	// Пример вывода
	/*
	 Непрерывный вариант
 Матрица 
  1.00000  1.50000  2.33333
  1.50000  2.33333  3.75000
  2.33333  3.75000  6.20000

 Вектор правых частей 

  10.9614353597607  17.3490423679259  28.3151499248069

 P2(x)=(1.57921) + (1.27098)*x + (3.20390)*x^2
 Норма погрешности: 0.0221201698386038

     x         Погрешность
    1.00   0.0540824053781916
    1.20   -0.0191997202979595
    1.40   -0.0173210826705557
    1.60   0.0152038921636972
    1.80   0.0229222063568955
    2.00   -0.0632456552045042
	*/
}
