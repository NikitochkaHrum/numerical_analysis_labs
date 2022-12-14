#include "direct_methods/direct_sqrt.h"
#include "iteration_methods/simple_iteration.h"
#include "iteration_methods/gradient_descent.h"
#include "iteration_methods/simple_relaxation.h"
#include "iteration_methods/conj_grad.h"

int main(){
    const ld EPS = 1e-4;


    freopen("/home/pna/Documents/study/numerical_analysis_labs/input.txt", "r", stdin);
    freopen("/home/pna/Documents/study/numerical_analysis_labs/output.txt", "w", stdout);
    int n;
    cin >> n;
    Mat a (vector<vector<ld>>(n, vector<ld>(n)));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) cin >> a.data[i][j];
    vector<ld> bb(n);
    for(int i=0; i<n; i++) bb[i]=i+1;
    Vec b(bb);
                    cout << "Вариант 9А.\n"; // ВАРИАНТ 
    cout << "A:\n";
    a.print();
    cout << "\n";
    cout << "b:\n";
    b.print();
    cout << '\n';
    Vec solution = direct_sqrt(a, b);
    cout << "Метод прямого квадратного корня:\nX:\n";
    solution.print();
    cout << '\n';
    cout << "Норма матрицы:\n";
    cout << a.matrixnormeuk(true) << '\n';
    // cout << "\nМетод простой итерации:\n";
    Vec x = simple_iteration(a, b, b, solution, EPS);
    cout << '\n';
    x = gradient_descent(a, b, b, solution, EPS);
    cout << '\n';
    x = simple_relaxation(a, b, b, solution, EPS);
    cout << '\n';
    x = conj_grad(a, b, b, solution, EPS);
    ld cond = a.cond(3, true);
    cout << "\nЧисло обусловленности = " << cond << "\nТеоретическая оценка числа итераций:\n";
    cout << "\tМетоды простой итерации и градиентного спуска: " << simple_iteration_estimate_n_iterations(EPS, cond) << '\n';
    cout << "\tМетод простой релаксации: " << simple_relaxation_estimate_n_iterations(EPS, cond) << '\n';
    cout << "\tМетод сопряжённых градиентов: " << conj_grad_estimate_n_iterations(EPS, cond) << '\n';
    
}