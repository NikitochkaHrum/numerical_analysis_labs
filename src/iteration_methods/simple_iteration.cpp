#include "simple_iteration.h"

Vec simple_iteration(Mat a, Vec b, Vec init, Vec solution, ld EPS){
   	print_method_header("Метод простой итерации", simple_iteration_estimate_n_iterations(EPS, a.cond(3, true)));
    const ld tau = GAMMA * 2 / a.matrixnormeuk(true);
    Vec* x3 = new Vec(init.data);
    Vec* x2 = nullptr;
    Vec* x1 = nullptr;
    ld q, residual_norm, relative_error=1e9, error_estimate;
    ld solution_norm = solution.euclidian_norm();
    int itr = -1;

    do{
        delete x1;
        x1 = x2;
        x2 = x3;

        x3 = new Vec(x2->n);
        for(int i=0; i<b.n; i++){
            ld sum = 0;
            for (int j = 0; j < b.n; j++)
                sum += a[i][j] * (*x2)[j];
            (*x3).data[i] = (*x2)[i] + tau * (b[i] - sum);
        }
        
        tie(q, residual_norm, relative_error, error_estimate) = calc_params(x1, x2, x3, a, b, solution, solution_norm);

        itr++;
		if (itr >= 1)
			print_method_iteration_info(itr, tau, q, residual_norm, relative_error, error_estimate, x3);
    } while(relative_error>=EPS);
    delete x2;
    delete x1;
    return *x3;
}

int simple_iteration_estimate_n_iterations(ld EPS, ld cond){
    return (int)(log(1/EPS)/2*cond);
}