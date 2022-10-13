#include <bits/stdc++.h>
#include "matrices.h"

using namespace std;

typedef long double ld;
const ld pi = 3.14159265358979323846, eps = 1e-12;

Vec::Vec(){
    data.resize(0);
    n=0;
}

Vec::Vec(int n){
    this->n=n;
    data.resize(n, 0);
}

Vec::Vec(vector<ld> data) {
    this->data = data;
    n=data.size();
}

void Vec::print(){ // вывод вектора
    cout << fixed;
    for (auto el : this->data){
        //if (abs(el)<eps) cout << setw(12)  << 0 << ' ';
        //else 
            cout << setw(12)  << scientific << el << ' ';
    }
    cout << '\n';
}

ld Vec::cubic_norm() { // максимальный элемент
    if (data.empty()) {
        throw invalid_argument("empty vector");
    }
    ld mx = abs(this->data[0]);
    for (auto x : this->data) mx = max(mx, abs(x));
    return mx;
}

ld Vec::octahedral_norm() { // сумма элементов
    if (data.empty()) {
        throw invalid_argument("empty vector");
    }
    ld s = 0;
    for (auto x : data) s += abs(x);
    return s;
}

ld Vec::euclidian_norm() { // корень из суммы квадратов
    if (data.empty()) {
        throw invalid_argument("empty vector");
    }
    ld s = 0;
    for (auto x : data) s += x * x;
    return sqrt(s);
}

ld Vec::scalar_prod(Vec b){
    if (n != b.n) {
        throw invalid_argument("incomparable vector sizes");
    }
    ld res=0;
    for(int i=0; i<n; i++){
        res+=data[i]*b[i];
    }
    return res;
}

Vec Vec::operator +(const Vec b) { // сложение векторов
    if (data.size() != b.data.size()) {
        throw invalid_argument("incomparable vector sizes");
    }
    vector<ld> res(data.size());
    for (int i = 0; i < data.size(); i++)
        res[i] = data[i] + b.data[i];
    return Vec(res);
}

Vec Vec::operator -(const Vec b) { // вычитание векторов
    if (data.size() != b.data.size()) {
        throw invalid_argument("incomparable vector sizes");
    }
    vector<ld> res(data.size());
    for (int i = 0; i < data.size(); i++)
        res[i] = data[i] - b.data[i];
    return Vec(res);
}

Vec Vec::operator*(const ld x){
    vector<ld> res(n);
    for(int i=0; i<n; i++)
        res[i]=data[i]*x;
    return Vec(res);
}

ld Vec::operator[](int idx) { // индексация
    return data[idx];
}

Vec& Vec::operator=(const Vec& b){
    data=b.data;
    n=b.n;
    return *this;
}

Mat::Mat(){
    data.resize(0);
    n=0;
    m=0;
}

Mat::Mat(int n){
    this->n=n;
    this->m=n;
    data.resize(n, vector<ld>(n, 0));
}

Mat::Mat(int n, int m){
    this->n=n;
    this->m=m;
    data.resize(n, vector<ld>(m, 0));
}

Mat::Mat(vector<vector<ld>> data) {
    this->data = data;
    n=data.size();
    if(!n)
        m=0;
    else
        m=data[0].size();
}

void Mat::print(){ // вывод матрицы
    cout << fixed;
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            //if (abs(data[i][j])<eps ) cout << setw(20)  << 0 << ' ';
            //else 
                cout << setw(20)  << scientific << data[i][j] << ' ';
        }
        cout << '\n';
    }
}

vector<vector<ld>> Mat::permutted_data() { // P*A
    int n = data.size();
    vector<vector<ld>> res(n, vector<ld>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            res[i][j] = data[P[i]][j];
    }
    return res;
}

Mat Mat::transpose() { // транспонирование матрицы
    vector<vector<ld>> a = data;
    for (int i = 0; i < a.size(); i++) {
        for (int j = i + 1; j < a[0].size(); j++) {
            swap(a[i][j], a[j][i]);
        }
    }
    return Mat(a);
}

ld Mat::matrixnormone(){ //считает сумму по строкам
    if (data.empty() || data.size() != data[0].size()) {
        throw invalid_argument("not square matrix");
    }
    ld mx = 0;
    for (int i = 0; i < data.size(); i++) {
        ld sum = 0;
        for (int j = 0; j < data[0].size(); j++)
            sum += fabs(data[i][j]);
        if (sum > mx)
            mx = sum;
    }
    return mx;
}

ld Mat::matrixnormtwo(){ //считает сумму по столбцам
    if (data.empty() || data.size() != data[0].size()) {
        throw invalid_argument("not square matrix");
    }
    ld max = 0;
    for (int i = 0; i < data.size(); i++)
    {
        ld sum = 0;
        for (int j = 0; j < data[0].size(); j++)
            sum += abs(data[j][i]);
        if (sum > max)
            max = sum;
    }
    return max;
}

ld Mat::matrixnormeuk(bool symmetric){ //корень суммы квадратов всех элементов
    Mat A;
    if(symmetric)
        A=*this;
    else
        A = (this->transpose()) * (*this);
    while (sqrt(A.sumOfNonDiagonalSquares()) > eps) {
        ld theta, C, S, CS;
        auto p = A.maxNonDiagonal();
        int I = p.first, J = p.second;
        if (A.data[I][I] == A.data[J][J]) theta = pi / 4;
        else
            theta = atan(2 * A.data[I][J] / (A.data[I][I] - A.data[J][J])) / 2;
        C = cos(theta);
        S = sin(theta);
        Mat B = A;

        for (int k = 0; k < A.data.size(); k++) {
            B.data[k][I] = C * A.data[k][I] + S * A.data[k][J];
            B.data[k][J] = -S * A.data[k][I] + C * A.data[k][J];
        }

        for (int k = 0; k < A.data.size(); k++) {
            A.data[I][k] = A.data[k][I] = C * B.data[I][k] + S * B.data[J][k];
            A.data[J][k] = A.data[k][J] = -S * B.data[I][k] + C * B.data[J][k];
        }

        A.data[I][J] = A.data[J][I] = 0;
    }
    if(symmetric)
        return A.maxDiagonal();
    return sqrt(A.maxDiagonal());
}

ld Mat::cond(int a, bool symmetric){ // число обусловленности
    switch (a) {
    case 1:
        return (this->inverse()).matrixnormone() * this->matrixnormone();
    case 2:
        return (this->inverse()).matrixnormtwo() * this->matrixnormtwo();
    case 3:
        return symmetric ? (this->inverse()).matrixnormeuk(true) * this->matrixnormeuk(true) : (this->inverse()).matrixnormeuk() * this->matrixnormeuk();
    default:
        return -1;
    }
}

pair<int, int> Mat::maxNonDiagonal(){
    ld max = abs(data[0][1]);
    int I = 0;
    int J = 1;
    for (int i = 0; i < data.size(); i++)
        for (int j = i + 1; j < data.size(); j++) {
            if (abs(data[i][j]) > max) {
                max = abs(data[i][j]);
                I = i;
                J = j;
            }
        }
    return make_pair(I, J);
}

ld Mat::maxDiagonal() {
    ld mx = abs(data[0][0]);
    for (int i = 1; i < data.size(); i++) {
        if (abs(data[i][i]) > mx)
            mx = abs(data[i][i]);
    }
    return mx;
}

ld Mat::sumOfNonDiagonalSquares() {
    ld ret = 0;
    for (int i = 0; i < data.size(); i++)
        for (int j = i + 1; j < data.size(); j++)
            ret += data[i][j] * data[i][j];
    return 2 * ret;
}

void Mat::calc_LUP(bool log) {
    if (data.empty() || data.size() != data[0].size()) {
        throw invalid_argument("not square matrix");
    }
    n_swaps = 0;
    int n = data.size();
    rang = 0;
    vector<vector<ld>> E(n, vector<ld>(n)); // единичная матрица
    for (int i = 0; i < n; i++) E[i][i] = 1;
    P.assign(n, 0);
    L.assign(n, vector<ld>(n, 0));
    for (int i = 0; i < n; i++) P[i] = i;
    U = data;
    for (int i = 0; i < n; i++) {
        ld mx = abs(U[i][i]);
        int mx_idx = i;
        for (int row = i + 1; row < n; row++) { // ищем максимальный элемент
            if (fabs(U[row][i]) > mx) {
                mx = fabs(U[row][i]);
                mx_idx = row;
            }
        }
        if (mx > eps)
            rang++;

        ld norm = U[mx_idx][i];

        if (log)
            cout << "k = " << i + 1 << " m = " << mx_idx + 1 << '\n';

        for (int j = i; j < n; ++j) { // меняем местами текущую строку и строку с макс элементом, и делим элементы на макс элемент
            swap(U[i][j], U[mx_idx][j]);
            U[i][j] /= norm;
        }
        if (i != mx_idx) n_swaps++; // число перестановок
        swap(P[i], P[mx_idx]);
        ld kof;

        for (int j = i + 1; j < n; ++j) {
            kof = U[j][i] / U[i][i];
            for (int z = i; z < n; ++z)
                U[j][z] -= kof * U[i][z];
        }


        for (int j = 0; j <= i; ++j) { // считаем матрицу L
            L[i][j] = data[P[i]][j];
            for (int z = 0; z <= j - 1; ++z)
                L[i][j] -= L[i][z] * U[z][j];
        }

        if (log) {
            cout << "Current L:\n";
            Mat(L).print();
            cout << "Current U:\n";
            Mat(U).print();
            cout << "Current P:\n";
            for (int ii = 0; ii < n; ii++) cout << P[ii] << ' ';
            cout << "\n\n\n";
        }

    }
    if (log)
        cout << "rang = " << rang << "\n\n";
}

ld Mat::determinant() { // определитель
    if (L.empty()) {
        this->calc_LUP(false);
    }
    ld ans = 1;
    for (int i = 0; i < L.size(); i++)
        ans *= L[i][i];
    if (n_swaps % 2)
        return -ans;
    return ans;
}

Mat Mat::inverse() { // обратная матрица
    if (data.empty() || data.size() != data[0].size()) {
        throw invalid_argument("not square matrix");
    }
    int n = data.size();

    vector<vector<ld>> inv(n, vector<ld>(n));
    for (int i = 0; i < n; i++)
    {
        vector<ld> y(n);
        y[i] = 1;
        Vec x = solve_equation(*this, Vec(y));
        for (int j = 0; j < n; j++) {
            inv[j][i] = x.data[j];
        }
    }
    return Mat(inv);
}

Mat Mat::operator +(const Mat b) { // сложение матриц
    if (data.size() != b.data.size() || data.empty() || data[0].size() != b.data[0].size()) {
        throw invalid_argument("incomparable matrix sizes");
    }

    vector<vector<ld>> res(data.size(), vector<ld>(data[0].size()));
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            res[i][j] = data[i][j] + b.data[i][j];
        }
    }
    return Mat(res);
}

Mat Mat::operator -(const Mat b) { // сложение матриц
    if (data.size() != b.data.size() || data.empty() || data[0].size() != b.data[0].size()) {
        throw invalid_argument("incomparable matrix sizes");
    }

    vector<vector<ld>> res(data.size(), vector<ld>(data[0].size()));
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            res[i][j] = data[i][j] - b.data[i][j];
        }
    }
    return Mat(res);
}

Mat Mat::operator *(const Mat r) { // умножение матриц
    if (data.empty() || r.data.empty() || data[0].size() != r.data.size()) {
        throw invalid_argument("incomparable matrix sizes");
    }
    vector<vector<ld>> res(data.size(), vector<ld>(r.data[0].size()));
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < r.data[0].size(); j++) {
            for (int k = 0; k < r.data.size(); k++) {
                res[i][j] += data[i][k] * r.data[k][j];
            }
        }
    }
    return Mat(res);
}

Vec Mat::operator *(const Vec v) { // умножение матрицы на вектор
    if (data.empty() || data[0].size() != v.data.size()) {
        throw invalid_argument("incomparable matrix and vector sizes");
    }
    vector<ld> res(data.size());
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < v.data.size(); j++)
            res[i] += data[i][j] * v.data[j];
    }
    return Vec(res);
}

vector<ld> Mat::operator[](int idx){ // индексация
    return data[idx];
}

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
