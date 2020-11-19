#ifndef __INTEGRATE_H__
#define __INTEGRATE_H__

#include <iostream>

typedef void (* DERIV_TYP)(double *x_0, int m, double h, double ch, double *k_arr); 

/* RK5 aims to be adaptive step size 4th order runge-kutta method based on Dorman
prince method. Need to add step size control and interpolation at grid points. Need to specify 
TSPAN and RK5 should return solution interpolated for the N points requested. */
class runge_kutta{
    public:
        runge_kutta();
        ~runge_kutta();
//        typedef void (* DERIV_TYP)(double *x_0, int m, double h, double ch); 
        void RK5(double *x_0, double *x_new, double h, int m, DERIV_TYP DERIV);
    private:
    double ch;
    static const double *y1, *y2, *y3, *y4, *y5, *y6;
    static const double *k1, *k2, *k3, *k4, *k5, *k6, *k7; 
    static const double a21                                       ;   
    static const double a31,  a32                                 ;   
    static const double a41,  a42,  a43                           ;   
    static const double a51,  a52,  a53,  a54                     ;
    static const double a61,  a62,  a63,  a64,  a65               ;
    static const double a71,  a72,  a73,  a74,  a75,  a76         ;
    static const double b1,   b2,   b3,   b4,   b5,   b6,   b7    ;
    static const double b1_2, b2_2, b3_2, b4_2, b5_2, b6_2, b7_2  ;
    static const double c1,   c2,   c3,   c4,   c5,   c6,   c7    ;
};

#endif
