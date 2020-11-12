#include <iostream>
#include "integrate.h"

void lorenz(double *x_0, int m, double h, double ch, double *k){
    double sigma = 10.0, b = 8.0/3.0, r = 28.0;
    k[0] = -sigma*x_0[0] + sigma*x_0[1];
    k[1] = r*x_0[0] - x_0[1] - x_0[0]*x_0[2];
    k[2] = x_0[1]*x_0[0] -b*x_0[2];
};


int main(){
    int m = 3;
    int nt = 100000;
    double dt = 0.001;
    double x[m][nt];
    double x0[m], x_n[m];

    runge_kutta rk;

    x[0][0] = 1.0; x[1][0] = 1.0; x[2][0] = 0.75;

    for (int j = 0; j < nt; j++){
        x0[0] = x[0][j]; x0[1] = x[1][j]; x0[2] = x[2][j];
        rk.RK5(x0, x_n, dt, m, (DERIV_TYP)lorenz);
        x[0][j+1] = x_n[0]; x[1][j+1] = x_n[1]; x[2][j+1] = x_n[2];
        std::cout << x[0][j] << "   " << x[1][j] << "   " << x[2][j] << std::endl;
    };

    return 0;
};
