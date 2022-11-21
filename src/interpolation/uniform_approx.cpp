#include "uniform_approx.h"
#include "method_log/method_info.h"
#include "../non_lynear_systems/solve_nle.h"

void uniform_approximation_first_order_polynomial(function<ld(ld)> f, function<ld(ld)> der_f, Mat& values){
    const int N = values.n - 1;
	const ld A = values[0][0];
	const ld B = values[N][0];

	const ld PRECISION = 1e-7; // точность для алгоритма золотого сечения
	Vec coefficients(2);

	// равномерное приближение с помощью полинома первого порядка
	// g(x) = a0 + a1 * x
	coefficients.data[1] = (f(B) - f(A)) / (B - A);

	auto r = [&](ld x) -> ld { return der_f(x) - coefficients[1]; }; // f'(d) = a1 -> f'(d) - a1 = 0 -> r(x) = f'(x) - a1

	ld d  = solve_nle(r, A, B, PRECISION); // находим исходя из f'(d) = a
	coefficients.data[0] = (f(A) + f(d) - coefficients[1] * (A + d)) / 2;

	auto g = [&](ld x) -> ld { return coefficients[0] + coefficients[1] * x; };

	Vec l_values(3);
	l_values.data[0] = f(A) - g(A);
	l_values.data[1] = f(d) - g(d);
	l_values.data[2] = f(B) - g(B);

	// Вывод результатов
	cout << "Равномерное приближение" << endl;
	cout << endl;
	print_polynom(coefficients);
	cout.precision(8);
	cout << fixed;

	cout << "d = " << d << endl;

	cout << "L(a) = " << l_values[0] << ", ";
	cout << "L(d) = " << l_values[1] << ", ";
	cout << "L(b) = " << l_values[2] << "\n\n";

	print_error(g, values);
	// 
	// Пример вывода
	/*
	 Равномерное приближение

 P1(x)=-5.4053+11.0000*x,   d= 1.5453
 L(a)=0.40525,   L(d)=-0.40525,   L(b)=0.40525

     x         Погрешность
    1.00   -0.405253514144292
    1.20   0.0575536670091559
    1.40   0.339209764109629
    1.60   0.395200351060421
    1.80   0.170072430013629
    2.00   -0.405253514144292*/
}
