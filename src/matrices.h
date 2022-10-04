#include <bits/stdc++.h>

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
    Vec operator-(Vec b); // вычитание векторов
};

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
    void calc_LUP(bool log=true);
    ld determinant(); // определитель
    Mat inverse(); // обратная матрица
    Mat operator *(Mat r); // умножение матриц
    Mat operator +(Mat b); // сложение матриц
    Mat operator -(Mat b); // вычитание матрицы из матрицы
    Vec operator *(Vec v); // умножение матрицы на вектор
};


Vec solve_equation(Mat A, Vec b);