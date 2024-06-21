#ifndef RDP_HISTORY
#define RDP_HISTORY

double repcr(double taud1, double taud2,  double atube);
double ret(double l, double taus, double fe, double atube, int q);
double ccr(double l1,double l2,  double taus, double fe, double atube, int q);
int RDP(int N_filed, double *Z, double edot[][3], double *Vz);

#endif
