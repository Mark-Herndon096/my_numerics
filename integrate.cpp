#include <iostream>
#include "integrate.h"

runge_kutta::runge_kutta(){};

runge_kutta::~runge_kutta(){};

void runge_kutta::RK5(double *x_0, double *x_new, double h, int m, DERIV_TYP DERIV){
    double  *y1, *y2, *y3, *y4, *y5, *y6;
    double  *k1, *k2, *k3, *k4, *k5, *k6, *k7;
    y1 = new double[m]; 
    y2 = new double[m]; 
    y3 = new double[m]; 
    y4 = new double[m]; 
    y5 = new double[m]; 
    y6 = new double[m];
 
    k1 = new double[m]; 
    k2 = new double[m]; 
    k3 = new double[m]; 
    k4 = new double[m]; 
    k5 = new double[m]; 
    k6 = new double[m]; 
    k7 = new double[m];

    ch = c1*h;

    DERIV(x_0,m,h,ch,k1);

    for (int j = 0; j < m; j++){
         y1[j] = x_0[j] + h*a21*k1[j];
    };

    ch = c2*h;

    DERIV(y1,m,h,ch,k2);

    for (int j = 0; j < m; j++){
        y2[j] = x_0[j] + h*(a31*k1[j] + a32*k2[j]);
    };

    ch = c3*h;

    DERIV(y2,m,h,ch,k3);

    for (int j = 0; j < m; j++){
        y3[j] = x_0[j] + h*(a41*k1[j] + a42*k2[j] + a43*k3[j]);
    };

    ch = c4*h;

    DERIV(y3,m,h,ch,k4);

    for (int j = 0; j < m; j++){
        y4[j] = x_0[j] + h*(a51*k1[j] + a52*k2[j] + a53*k3[j] + a54*k4[j]);
    };

    ch = c5*h;

    DERIV(y4,m,h,ch,k5);

    for (int j = 0; j < m; j++){
        y5[j] = x_0[j] + h*(a61*k1[j] + a62*k2[j] + a63*k3[j] + a64*k4[j] + a65*k5[j]);
    };

    ch = c6*h;

    DERIV(y5,m,h,ch,k6);

    for (int j = 0; j < m; j++){
        y6[j] = x_0[j] + h*(a71*k1[j] + a72*k2[j] + a73*k3[j] + a74*k4[j] + a75*k5[j] +a76*k6[j]);
    };

    ch = c7*h;

    DERIV(y6,m,h,ch,k7);
    
    for (int j = 0; j < m; j++){
        x_new[j]  = x_0[j] + h*(b1_2*k1[j] + b2_2*k2[j] + b3_2*k3[j] + b4_2*k4[j] + b5_2*k5[j] + b6_2*k6[j] + b7_2*k7[j]);
    };

    delete [] y1; 
    delete [] y2; 
    delete [] y3; 
    delete [] y4; 
    delete [] y5; 
    delete [] y6; 

    delete [] k1; 
    delete [] k2; 
    delete [] k3; 
    delete [] k4; 
    delete [] k5; 
    delete [] k6; 
    delete [] k7; 
};


const double  runge_kutta::c1   = 0.0;
const double  runge_kutta::c2   = 1.0/5.0;
const double  runge_kutta::c3   = 3.0/10.0;
const double  runge_kutta::c4   = 4.0/5.0;
const double  runge_kutta::c5   = 8.0/9.0;
const double  runge_kutta::c6   = 1.0;
const double  runge_kutta::c7   = 1.0;
const double  runge_kutta::a21  = 1.0/5.0;
const double  runge_kutta::a31  = 3.0/40.0;
const double  runge_kutta::a32  = 9.0/40.0;
const double  runge_kutta::a41  = 44.0/45.0;
const double  runge_kutta::a42  = -56.0/15.0;
const double  runge_kutta::a43  = 32.0/9.0;
const double  runge_kutta::a51  = 19372.0/6561.0;
const double  runge_kutta::a52  = -25360.0/2187.0;
const double  runge_kutta::a53  = 64448.0/6561.0;
const double  runge_kutta::a54  = -212.0/729.0;
const double  runge_kutta::a61  = 9017.0/3168.0;
const double  runge_kutta::a62  = -355.0/33.0;
const double  runge_kutta::a63  = 46732.0/5247.0;
const double  runge_kutta::a64  = 49.0/176.0;
const double  runge_kutta::a65  = -5103.0/18656.0;
const double  runge_kutta::a71  = 35.0/384.0;
const double  runge_kutta::a72  = 0.0;
const double  runge_kutta::a73  = 500.0/1113.0;
const double  runge_kutta::a74  = 125.0/192.0;
const double  runge_kutta::a75  = -2187/6784.0;
const double  runge_kutta::a76  = 11.0/84.0;
const double  runge_kutta::b1   = runge_kutta::a71;
const double  runge_kutta::b2   = runge_kutta::a72;
const double  runge_kutta::b3   = runge_kutta::a73;
const double  runge_kutta::b4   = runge_kutta::a74;
const double  runge_kutta::b5   = runge_kutta::a75;
const double  runge_kutta::b6   = runge_kutta::a76;
const double  runge_kutta::b7   = 0.0;
const double  runge_kutta::b1_2 = 5179.0/57600.0;
const double  runge_kutta::b2_2 = 0.0;
const double  runge_kutta::b3_2 = 7571.0/16695.0;
const double  runge_kutta::b4_2 = 393.0/640.0;
const double  runge_kutta::b5_2 = -92097.0/339200.0;
const double  runge_kutta::b6_2 = 187.0/2100.0;
const double  runge_kutta::b7_2 = 1.0/40.0;     
