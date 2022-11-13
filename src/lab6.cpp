#include "integration/input_data.h"
#include "integration/function_tool.h"
#include "integration/trapez.h"
#include "integration/trapez_spline.h"
#include "integration/simpson.h"
#include "integration/gauss.h"

const ld eps = 1e-8;
const int RESULT_PRECISION = 16;

void run_var(FunctionTool f, function<ld(ld)> Der, const ld J_EXACT){
    ld result = 0;

    auto old_precision = cout.precision();
    cout << "J = " << std::setprecision(RESULT_PRECISION) << J_EXACT << "\n\n";
    cout.precision(old_precision);

    old_precision = cout.precision();
    cout << "Формула трапеций" << '\n';
    result = trapezoidal_rule(f, A, B, eps, J_EXACT);
    cout << std::setprecision(RESULT_PRECISION) << "Результат: " << J_EXACT - result << '\n';
    cout << "Количество обращений к f(x): " << f.get_calls_count() << '\n' << '\n';
    cout.precision(old_precision);

    f.reset();

    old_precision = cout.precision();
    cout << "Формула трапеций(модифицированной с помощью сплайна)" << '\n';
    result = trapezoidal_rule_spline(f, Der, A, B, eps, J_EXACT);
    cout << std::setprecision(RESULT_PRECISION) << "Результат: " << J_EXACT - result << '\n';
    cout << "Количество обращений к f(x): " << f.get_calls_count() << '\n';
    cout << "Количество обращений к f'(x): 2" << '\n' << '\n';
    cout.precision(old_precision);

    f.reset();

    old_precision = cout.precision();
    cout << "Формула Симпсона" << '\n';
    result = simpson_rule(f, A, B, eps, J_EXACT);
    cout << std::setprecision(RESULT_PRECISION) << "Результат: " << J_EXACT-result << '\n';
    cout << "Количество обращений к f(x): " << f.get_calls_count() << '\n' << '\n';
    cout.precision(old_precision);

    f.reset();

    old_precision = cout.precision();
    cout << "Формула Гаусса(трехточечной)" << '\n';
    result = gauss(f, A, B, eps, J_EXACT);
    cout << setprecision(RESULT_PRECISION) << "Результат: " << J_EXACT - result << '\n';
    cout << "Количество обращений к f(x): " << f.get_calls_count() << '\n' << '\n';
    cout.precision(old_precision);

}

int main() {
    freopen("/home/pna/Documents/study/numerical_analysis_labs/output.txt", "w", stdout);
    FunctionTool f9(F_var9);
    FunctionTool f25(F_var25);
    cout << "Вариант: 9\n";
    run_var(f9, Der_var9, J_EXACT_var9);
    cout << "Вариант: 25\n";
    run_var(f25, Der_var25, J_EXACT_var25);
    
}