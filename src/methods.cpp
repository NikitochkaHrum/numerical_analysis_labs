#include "methods.h"

ld sign(ld x){
    return (x > 0) - (x < 0);
}

Vec direct_sqrt(Mat a, Vec b){
    int n = b.n;

    Mat s(n, n);
    Mat d(n, n);

    ld sum = 0;

    for(int i=0; i<n; i++){
        sum=0;
        for(int k=0; k<i; k++){
            sum+=s[k][i]*s[k][i]*d[k][k];
        }
        d.data[i][i] = sign(a[i][i]-sum);
        s.data[i][i] = sqrt(abs(a[i][i]-sum));

        for(int j=i+1; j<n; j++){
            sum=0;
            for(int k=0; k<i; k++){
                sum+=s[k][i]*s[k][j]*d[k][k];
            }
            s.data[i][j]=(a[i][j]-sum)/(s[i][i]*d[i][i]);
        }
    }
    auto st=s.transpose()*d;
    Vec y(n);
    for (int i = 0; i < n; i++) {
        sum = b[i];
        for (int j = 0; j < i; j++){
            sum -= st[i][j] * y[j];
        }
        y.data[i] = sum / st[i][i];
    }

    Vec x(n);
    for (int i = n - 1; i >= 0; i--){
        sum = y[i];
        for (int j = i + 1; j < n; j++){
            sum -= s[i][j] * x[j];
        }
        x.data[i] = sum / s[i][i];
    }

    return x;
}