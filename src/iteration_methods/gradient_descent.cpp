#include "gradient_descent.h"

Vec gradient_descent(Mat a, Vec b, Vec init, Vec solution, ld EPS){
    print_method_header("Метод наискорейшего спуска", gradient_descent_estimate_n_iterations(EPS, a.cond(3, true)));
    ld tau, r_err;

	Vec* x3 = new Vec(init.data);
	Vec* x2 = nullptr;
    Vec* x1 = nullptr;

	ld q, residual_norm, relative_error, error_estimate;
	ld solution_norm = solution.euclidian_norm();

	Vec rk; 
	Vec ark;
	int itr = -1;

	do
	{
		rk = a * (*x3) - b;
		ark = a * rk;
		tau = ark.scalar_prod(rk) / ark.scalar_prod(ark);

		delete x1;
		x1 = x2;
		x2 = x3;
		x3 = new Vec(((*x3) - rk * tau).data);

		std::tie(q, residual_norm, relative_error, error_estimate) = calc_params(x1, x2, x3, a, b, solution, solution_norm);

		itr++;
		if (itr >= 1)
			print_method_iteration_info(itr, tau, q, residual_norm, relative_error, error_estimate, x3);

	} while (relative_error >= EPS);

    delete x2; delete x1;

	return *x3;
}

int gradient_descent_estimate_n_iterations(ld EPS, ld cond){
    return log(1 / EPS) / 2 * cond;
}