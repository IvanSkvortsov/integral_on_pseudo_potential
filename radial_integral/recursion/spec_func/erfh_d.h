#ifndef	D_ERFH_H
#define D_ERFH_H

#ifndef Inf
#define Inf (1./0.)
#endif
#ifndef NaN
#define NaN (0./0.)
#endif

static double erfh_y100(double y100, double x);

double erfh(double x);

#endif//D_ERFH_H
