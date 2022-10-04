#include "methods.h"

int main(){
    freopen("/home/pna/Documents/study/numerical_analysis_labs/src/input.txt", "r", stdin);
    freopen("/home/pna/Documents/study/numerical_analysis_labs/src/output.txt", "w", stdout);
    int n;
    cin >> n;
    vector<vector<ld>> a(n, vector<ld>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) cin >> a[i][j];
    vector<ld> bb(n);
    for(int i=0; i<n; i++) bb[i]=i+1;
    Vec b(bb);
    Vec ans = direct_sqrt(Mat(a), b);
    ans.print();
}