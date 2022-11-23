#include "discrete_rms.h"
#include "method_log/method_info.h"


void discrete_rms(vector<function<ld(ld)>>& basis, Mat& values){
    const int N = values.n - 1;
	const ld A = values[0][0];
	const ld B = values[N][0];

	const int M = basis.size();

	Mat basis_values(M);

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			for (int k = 0; k <= N; k++) {
				ld x = values[k][0];
				basis_values.data[i][j] += basis[i](x) * basis[j](x);
			}
		}
	}

	Vec scalar_with_f(M);
	for (int i = 0; i < M; i++) {
		for (int k = 0; k <= N; k++) {
			ld x = values[k][0];
			ld f = values[k][1];
			scalar_with_f.data[i] += basis[i](x) * f;
		}
	}


	Vec c = solve_equation(basis_values, scalar_with_f);
	// dodelat'
	cout << c[0] << " + " << c[1] << " * x + " << c[2] << " * x^2\n";

	cout << c[0] << " + " << c[1] << " * x + " << c[2] << " * x^2\n";

	cout << "Дискретный вариант:\n";

	cout << "Матрица коэффициентов в правой части:\n";
	basis_values.print();
	cout << '\n';

	cout << "Вектор правой части" << endl;
	
	scalar_with_f.print();
	cout << '\n';

	print_polynom(c);

    auto g = [&](ld x) -> ld
    {
        return c[0] + c[1] * x + c[2] * x * x;
    };

	ld norm_f = 0, norm_g = 0;

	for (int k = 0; k <= N; k++)
	{
		ld x = values[k][0];
		ld f = values[k][1];
		norm_f += f * f;
		norm_g += g(x) * g(x);
	}

    cout << "Норма погрешности: " << sqrt(norm_f - norm_g) << "\n\n";

	print_error(g, values);

	cout << '\n';

	// Пример вывода
	/*
	 Дискретный вариант 
 Матрица 
  6.00000  9.00000 14.20000
  9.00000 14.20000 23.40000
 14.20000 23.40000 39.96640

 Вектор правых частей 
  66.41694973123  107.286069909249  180.361191679766

 P2(x)=(1.52792) + (1.27708)*x + (3.22223)*x^2
 Норма погрешности 0.0755427116996264

     x           Погрешность
    1.00   0.0272331665325982
    1.20   -0.0367611129948235
    1.40   -0.0241278444628481
    1.60   0.0206185460318196
    1.80   0.0420250606412758
    2.00   -0.0289878157480246
	*/
}
