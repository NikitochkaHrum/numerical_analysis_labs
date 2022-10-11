#include "direct_methods/direct_sqrt.h"
#include "iteration_methods/simple_iteration.h"

int main(){
    const double EPS = 1e-4;


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
    cout << "b:\n";
    b.print();
    cout << '\n';
    Vec solution = direct_sqrt(a, b);
    cout << "Метод прямого квадратного корня:\nX:\n";
    solution.print();
    cout << '\n';

    cout << "Метод простой итерации:\n";
    Vec x = simple_iteration(a, b, b, solution, EPS);
    cout << "X:\n";
    x.print();
    cout << '\n';
}