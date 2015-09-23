#ifndef	LD_DAWSON_H
#define LD_DAWSON_H

#ifndef Inf
#define Inf (1./0.)
#endif
#ifndef NaN
#define NaN (0./0.)
#endif

static long double dawson_y100(long double y100, long double x);

long double dawson(long double x);

#endif//LD_DAWSON_H
