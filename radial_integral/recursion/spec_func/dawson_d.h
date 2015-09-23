#ifndef	D_DAWSON_H
#define D_DAWSON_H

#ifndef Inf
#define Inf (1./0.)
#endif
#ifndef NaN
#define NaN (0./0.)
#endif

static double dawson_y100(double y100, double x);

double dawson(double x);

#endif//D_DAWSON_H
