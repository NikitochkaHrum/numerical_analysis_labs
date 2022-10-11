#include "conj_grad.h"

Vec conj_grad(Mat a, Vec b, Vec init, Vec solution, ld EPS){
    print_method_header("Метод сопряженных градиентов", conj_grad_estimate_n_iterations(EPS, a.cond(3)));

	ld tau = 0, tauprev = 0, alpha = 1, alphaprev = 1;

	Vec* x3 = new Vec(init.data);
	Vec* x2 = new Vec(4), * x1 = nullptr; 

	ld q, residual_norm, relative_error, error_estimate;
	ld solution_norm = solution.euclidian_norm();

	Vec rk, rkprev;
	Vec ark;
	int itr = 0;

	ld tmp;

	do
	{
		rk = a * (*x3) - b;
		rkprev = a * (*x2) - b;
		ark = a * rk;

		tauprev = tau;
		tau = rk.scalar_prod(rk) / ark.scalar_prod(rk);
		if (itr >= 1)
		{
			alpha = 1 - (tau / tauprev) * (1 / alphaprev) * (rk.scalar_prod(rk) / rkprev.scalar_prod(rkprev));
			alpha = 1 / alpha;

			alphaprev = alpha;
		}

		delete x1;
		x1 = x2;
		x2 = x3;
		x3 = new Vec(((*x2) * alpha + (*x1) * (1 - alpha) - rk * (tau * alpha)).data);

		tie(q, residual_norm, relative_error, error_estimate) = calc_params(x1, x2, x3, a, b, solution, solution_norm);
		itr++;
		if (itr >= 1)
			print_method_iteration_info(itr, tau, q, residual_norm, relative_error, error_estimate, x3);
	} while (relative_error >= EPS);

    delete x2; delete x1;
	return (*x3);
}

int conj_grad_estimate_n_iterations(ld EPS, ld cond){
    return log(2 / EPS) / 2 * sqrt(cond);
}
