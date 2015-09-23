#ifndef	LD_ERFH_H
#define LD_ERFH_H

#ifndef Inf
#define Inf (1./0.)
#endif
#ifndef NaN
#define NaN (0./0.)
#endif

static long double erfh_y100(long double y100, long double x);

long double erfh(long double x);

#endif//LD_ERFH_H
