#include "method_info.h"

void print_error(function<ld(ld)> g, Mat& values)
{
	const int N = values.n;
	cout << setw(5) << "x" << setw(3) << "" << " | " << setw(14) << "Погрешность\n";
	for (int i = 0; i < N; i++)
	{
		ld x = values[i][0];
		ld f = values[i][1];
		ld error = f - g(x);

		cout.precision(4);
		cout << fixed << setw(8) << x;
		cout << "   ";
		cout.precision(10);
		cout << scientific << setw(15) << error << '\n';
	}
}

void print_polynom(Vec& coefficients)
{
	const int N = coefficients.n;

	cout << "P" << N-1 << "(x) = ";
	cout.precision(5);
	cout << fixed;
	for (int i = N - 1; i >= 0; i--)
	{
		ld c = coefficients[i];
		cout << " ";
		if (c > 0 && i != N - 1) {
			cout << "+";
		}
		cout << c;
		if (i != 0) {
			cout << "*x";
			if (i > 1)
				cout << "^" << i;
		}
	}

	cout << '\n';
}

void print_interpolation_results_header()
{
	//    x    f(x)       Pn(x)           Delta             Оценка
	cout << setw(6)  << "x"					<< " |"
		<< setw(12) << "f(x)"				<< " |"
		<< setw(12) << "g(x)"				<< " |"
		<< setw(24) << "Погрешность"		<< " |"
		<< setw(24) << "Оценка погрешности\n";
}

void print_interpolation_results(Vec& x, function<ld(ld)> f, function<ld(ld)> g, ld error_est)
{
	print_interpolation_results_header();
	int n = x.n;

	for (int i = 0; i < n; i++) {
		ld xv = x[i];
		cout.precision(2);
		cout << fixed;
		cout << setw(6) << xv << " |";
		cout.precision(8);
		cout << setw(12) << f(xv) << " |";
		cout << setw(12) << g(xv) << " |";
		cout << scientific;
		cout.precision(14);
		cout << setw(24) << abs(f(xv) - g(xv)) << " |";
		cout << setw(24) << error_est << '\n';
	}
}

void print_interpolation_results(Vec& x, function<ld(ld)> f, Vec& g_values, ld error_est)
{
	print_interpolation_results_header();
	int n = x.n;

	for (int i = 0; i < n; i++) {
		ld xv = x[i];
		cout.precision(2);
		cout << fixed;
		cout << setw(6) << xv << " |";
		cout.precision(8);
		cout << setw(12) << f(xv) << " |";
		cout << setw(12) << g_values[i] << " |";
		cout << scientific;
		cout.precision(14);
		cout << setw(24) << abs(f(xv) - g_values[i]) << " |";
		cout << setw(24) << error_est << '\n';
	}
}

void print_interpolation_results(Vec& x, function<ld(ld)> f, function<ld(ld)> g, Vec& errors)
{
	print_interpolation_results_header();
	int n = x.n;

	for (int i = 0; i < n; i++) {
		ld xv = x[i];
		cout.precision(2);
		cout << fixed;
		cout << setw(6) << xv << " |";
		cout.precision(8);
		cout << setw(12) << f(xv) << " |";
		cout << setw(12) << g(xv) << " |";
		cout << scientific;
		cout.precision(14);
		cout << setw(24) << abs(f(xv) - g(xv)) << " |";
		cout << setw(24) << errors[i] << '\n';
	}
}
