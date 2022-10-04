#include "methods.h"

Vec direct_sqrt(Mat a, Vec b){
    int n = b.n;
    vector<vector<ld>> g(n, vector<ld>(n, 0));
    g[0][0]=sqrt(a[0][0]);
    for(int i=1; i<n; i++)
        g[0][i]=a[0][i]/g[0][0];
    for(int i=1; i<n; i++){
        ld s=0;
        for(int j=0; j<i; j++)
            s+=g[j][i]*g[j][i];
        g[i][i]=sqrt(a[i][i]-s);
    }
    for(int i=1; i<n; i++){
        for(int j=i+1; j<n; j++){
            ld s=0;
            for(int k=0; k<i; k++){
                s+=g[k][i]*g[k][j];
            }
            g[i][j]=(a[i][j]-s)/g[i][i];
        }
    }
    // Mat(g).print();
    // (Mat(g)*(Mat(g).transpose())).print();
    vector<ld> y(n);
    for(int i=0; i<n; i++){
        ld s=0;
        for(int j=0; j<i; j++){
            s+=g[j][i]*y[j];
        }
        y[i]=(b[i]-s)/g[i][i];
    }
    vector<ld> x(n);
    for(int i=n-1; i>-1; i--){
        ld s=0;
        for(int j=i+1; j<n; j++){
            s+=g[i][j]*x[j];
        }
        x[i]=(y[i]-s)/g[i][i];
    }
    return Vec(x);
}