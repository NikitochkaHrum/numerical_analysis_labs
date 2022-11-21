#include "spline.h"
#include "../direct_methods/sweep.h"
#include "method_log/method_info.h"

void cubic_spline_interpolation(function<ld(ld)> f, function<ld(ld)> der_f, Mat& values, ld m, ld m1)
{
	const int N = values.n - 1;
	const ld A = values[0][0];
	const ld B = values[N][0];

	Vec x(N + 1);

	for (int i = 0; i <= N; i++)
		x.data[i] = values[i][0];

	// значения производной исходной функции в крайней левой и в крайней правой точках отрезка
	ld m0 = der_f(A), mN = der_f(B);

	ld h = (B - A) / N; // шаг

	// m(0), m(N) - знаем
	// 1/2 * m(i-1) + 2 * m(i) + 1/2 * m(i+1) = 3 / 2h * (f(i+1)-f(i-1))
	// m(i-1) + 4 * m(i) + m(i+1) = 3 / h * (f(i+1)-f(i-1))

	// N = 5
	// i = 1: m0   + 4 * m(1) + m(2) = ...
	// i = 2: m(1) + 4 * m(2) + m(3) = ...
	// i = 3: m(2) + 4 * m(3) + m(4) = ...
	// i = 4: m(3) + 4 * m(4) + mN   = ...

	// m0, mN - знаем, найдем m(1), m(2), m(3), m(4)

	const int mCount = N - 1;

	Vec right(mCount); // матрица коэф-в левой части
	Mat left(mCount); // матрица свободных членов правой части

	for (int i = 1; i <= mCount; i++) {
		right.data[i - 1] = 3 * (values.data[i + 1][1] - values.data[i - 1][1]) / h;
	}

	for (int i = 0; i < mCount; i++) {
		if (i - 1 >= 0)
			left.data[i][i - 1] = 1;
		left.data[i][i] = 4;
		if (i + 1 < mCount)
			left.data[i][i + 1] = 1;
	}

	right.data[0] -= m0;
	right.data[N - 2] -= mN;

	// Решаем систему left * (m1, m2, m3, m4) = right с помощью метода прогонки, т.к. left - матрица трехдиаг.
	Vec found_m_values = sweep_method(left, right); // m1, ..., mN-1
	Vec all_m_values(N + 1);
	for (int i = 1; i < N; i++)
		all_m_values.data[i] = found_m_values[i - 1];

	all_m_values.data[0] = m0;
	all_m_values.data[N] = mN;

	auto fi0 = [](ld tau) -> ld
	{
		return (1 + 2 * tau) * (1 - tau) * (1 - tau);
	};

	auto fi1 = [](ld tau) -> ld
	{
		return tau * (1 - tau) * (1 - tau);
	};

	// функция, полученная при интерполяции сплайнами
	auto g = [&](ld x) -> ld
	{
		if (x < A || x > B) return 0;

		ld h = (B - A) / N;

		int i = min(int((x - A) / h), N - 1);

		ld x_left = values.data[i][0];
		ld x_right = values.data[i + 1][0];

		ld f_left = values.data[i][1];
		ld f_right = values.data[i + 1][1];

		ld tau = (x - x_left) / h;

		return fi0(tau) * f_left + fi0(1 - tau) * f_right + h * (fi1(tau) * all_m_values[i] - fi1(1 - tau) * all_m_values[i + 1]);
	};

	cout << "Интерполяция кубическим сплайном\n";
	cout << "M5 = " << fixed << m1 << '\n';

	print_interpolation_results(x, der_f, all_m_values, m1 * pow(h, 4) / 60);
	cout << '\n';

	Vec x2(N); int i = 0;
	for (ld t = 1.1; t <= B; t += h) {
		x2.data[i] = t;
		i++;
	}

	cout << "M4 = " << fixed << m << '\n';
	print_interpolation_results(x2, f, g, (m / 384 + m1 * h / 240) * pow(h, 4));

	cout << '\n';
	// Пример вывода
	/*
	   Интерполяция кубическим сплайном

 M5=14.4033917
   x[i]   df/dx(x[i])    m[i]          Delta          Оценка
    1.00  8.2958369    8.2958369           0            0.000384090446043703
    1.20  9.1057260    9.1056623  6.36193151724029E-05  0.000384090446043703
    1.40  10.1146299  10.1145646  6.52390365036837E-05  0.000384090446043703
    1.60  11.3714527  11.3713789  7.37050220855906E-05  0.000384090446043703
    1.80  12.9371157  12.9369796  0.000136091195779287  0.000384090446043703
    2.00  14.8875106  14.8875106           0            0.000384090446043703

 M4=13.1105321

     x      f(x)     S31(f;x)   Abs(f(x)-S31(f;x))       Оценка
    1.10   6.84837   6.84835  1.87494431251878E-05  7.38317395774316E-05
    1.30   8.67117   8.67114  2.52975821020129E-05  7.38317395774316E-05
    1.50  10.69615  10.69612  3.1352774506388E-05   7.38317395774316E-05
    1.70  12.97301  12.97297  3.77611254052113E-05  7.38317395774316E-05
    1.90  15.56363  15.56357  5.23853917950845E-05  7.38317395774316E-05
	*/
}