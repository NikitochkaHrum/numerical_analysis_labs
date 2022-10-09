#include <bits/stdc++.h>

using namespace std;

typedef long double ld;

class Vec
{
public:
    int n;
    vector<ld> data;
    Vec();
    Vec(int n);
    Vec(vector<ld> data);
    void print(); // вывод вектора
    ld cubic_norm(); // максимальный элемент
    ld octahedral_norm(); // сумма элементов
    ld euclidian_norm(); // корень из суммы квадратов
    Vec operator+(Vec b); // сложение векторов
    Vec operator-(Vec b); // вычитание векторов
    ld operator[](int idx); // индексация
};

class Mat
{
public:
    int n, m;
    vector<vector<ld>> data;
    vector<vector<ld>> L, U;
    vector<int> P;
    int n_swaps = 0, rang = 0;
    Mat();
    Mat(int n);
    Mat(int n, int m);
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
    void calc_LUP(bool log=true);
    ld determinant(); // определитель
    Mat inverse(); // обратная матрица
    Mat operator *(Mat r); // умножение матриц
    Mat operator +(Mat b); // сложение матриц
    Mat operator -(Mat b); // вычитание матрицы из матрицы
    Vec operator *(Vec v); // умножение матрицы на вектор
    vector<ld> operator [](int idx); // индексация
};

Vec solve_equation(Mat A, Vec b);