#ifndef RATE_RDP_COST
#define RATE_RDP_COST

double fun_R(double z,double zf, double R0,double a1,double a2,double a3,double a4);
double fun_U(double z,double zf, double u0,double b1,double b2,double b3,double b4);
double fun_dRdz(double z,double zf, double R0,double a1,double a2,double a3,double a4);
double fun_dudz(double z,double zf, double u0,double b1,double b2,double b3,double b4);
double fun_d2Rdz2(double z,double zf, double R0,double a1,double a2,double a3,double a4); 
double  CostFunction(double *X); 
#endif
