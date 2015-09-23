#ifndef VALUE_RECURS_H
#define VALUE_RECURS_H

//---------- Values of recursive functions --------//

// plus_minus
//
// These are two variables (values of function):
// p := ( + ) -> function...
// m := ( - ) -> function...
//
// where:
// +   := (ka + kb) / ( 2*alp )
// -   := |ka - kb| / ( 2*alp )
// ka  := -2 * alp_GBFa * CA
// kb  := -2 * alp_GBFb * CB
// alp := alp_GBFa + alp_ECP + alp_GBFb
//
// CA  := [ (C - A)_x, (C - A)_y, (C - A)_z ] -- 3d vector
// CB  := [ (C - B)_x, (C - B)_y, (C - B)_z ] -- 3d vector
//
template<class T>
struct plus_minus
{
	T p, m;
	plus_minus(): p(0), m(0) {}
};

// recursion_bc
//
// These are two recursively calculated functions
// b := (i, +/-) -> (i-1)/(2*alp) * b(i-2, +/-) + (+/-) * c(i-1, +/-)
// c := (i, +/-) -> (i-1)/(2*alp) * c(i-2, +/-) + (+/-) * b(i-1, +/-)
//
// where:
// (+/-) := see above
template<class T>
struct recursion_bc
{
	plus_minus<T> b, c;
};

#endif//VALUE_RECURS_H
