using namespace std;

typedef long double ld;

class Vec
{
public:
    vector<ld> data;
    Vec();
    Vec(vector<ld> data);
    void print(); // вывод вектора
    ld cubic_norm(); // максимальный элемент
    ld octahedral_norm(); // сумма элементов
    ld euclidian_norm(); // корень из суммы квадратов
    Vec operator+(Vec b); // сложение векторов
    Vec operator -(Vec b); // вычитание векторов
};

class Mat;
Vec solve_equation(Mat A, Vec b);

class Mat
{
public:
    vector<vector<ld>> data;
    vector<vector<ld>> L, U;
    vector<int> P;
    int n_swaps = 0, rang = 0;
    Mat();
    Mat(vector<vector<ld>> data);
    void print(); // вывод матрицы
    vector<vector<ld>> permutted_data(); // P*A
    Mat transpose(); // транспонирование матрицы
    ld matrixnormone(); //считает сумму по строкам
    ld matrixnormtwo(); //считает сумму по столбцам
    ld matrixnormeuk(); //корень суммы квадратов всех элементов
    ld cond(int a); // число обусловленности
    pair<int, int> maxNonDiagonal();
    ld maxDiagonal();
    ld sumOfNonDiagonalSquares();
    void calc_LUP(bool log = true);
    ld determinant(); // определитель
    Mat inverse(); // обратная матрица
    Mat operator *(Mat r); // умножение матриц
    Mat operator +(Mat b); // сложение матриц
    Mat operator -(Mat b); // вычитание матрицы из матрицы
    Vec operator *(Vec v); // умножение матрицы на вектор
};

Vec solve_equation(Mat A, Vec b) { // решение СЛУ
    if (A.data.empty() || A.data.size() != A.data[0].size() || A.data.size() != b.data.size())
        throw invalid_argument("incomparable matrix and vector sizes");
    if (A.L.empty() || A.U.empty())
        A.calc_LUP(false);
    int n = b.data.size();
    vector<ld> bb(n);
    for (int i = 0; i < n; ++i)
        bb[i] = b.data[A.P[i]];
    b.data = bb;
    vector<ld> y(n);

    int idx = 0;
    for (int i = 0; i < n; i++) {
        ld s = 0;
        for (int j = 0; j < i; j++) {
            s += A.L[i][j] * y[j];
        }
        y[i] = (b.data[i] - s) / A.L[i][i];
    }
    vector<ld> x(n);

    for (int i = n - 1; i > -1; i--) {
        ld s = 0;
        for (int j = n - 1; j > i; --j) {
            s += A.U[i][j] * x[j];
        }
        x[i] = (y[i] - s);
    }
    return Vec(x);
}

int main() {
    // freopen_s("input.txt", "r", stdin);
    // freopen_s("output.txt", "w", stdout);
    int n;
    cin >> n;

    vector<vector<ld>> a(n, vector<ld>(n));
    vector<ld> X_correct(n);

    for (int i = 0; i < n; i++) {
        X_correct[i] = i + 1;
    }
    cout << "X:\n";
    for (int i = 0; i < n; i++) cout << X_correct[i] << ' ';
    cout << "\n\n";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) cin >> a[i][j];
    Mat aa(a);

    aa.calc_LUP();
    cout << "Final P:\n";
    for (int i = 0; i < n; i++) cout << aa.P[i] << ' ';
    cout << '\n';
    cout << "Final L:\n";
    Mat(aa.L).print();
    cout << "Final U:\n";;
    Mat(aa.U).print();
    cout << "\nLU - PA:\n";
    (Mat(aa.L) * Mat(aa.U) - aa.permutted_data()).print();

    cout << "\nDeterminant = " << aa.determinant() << '\n';
    Vec b = aa * Vec(X_correct);
    Vec x = solve_equation(aa, Vec(b));
    cout << "\nB:\n";
    Vec(b).print();
    cout << "\nX:\n";
    x.print();
    cout << "\nAX - B:\n";
    (aa * x - Vec(b)).print();
    cout << "\nA^(-1):\n";
    aa.inverse().print();
    auto a_ainv = aa * aa.inverse();
    for (int i = 0; i < n; i++) a_ainv.data[i][i]--;
    cout << "\nA*A^(-1)-E:\n";
    a_ainv.print();

    // Числа обусловленности
    cout << "\nCond numbers:\n\tFirst norm:\n\t";
    cout << aa.cond(1) << '\n';
    cout << "\tSecond norm:\n\t";
    cout << aa.cond(2) << '\n';
    cout << "\tEuclidian norm:\n\t";
    cout << aa.cond(3) << '\n';

    cout << "Uncertainity:\n"; //погрешность
    (Vec(X_correct) - x).print();
}