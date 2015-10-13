#ifndef __HANKEL_HPP__
#define __HANKEL_HPP__
//#include"../../../lib_math/math.functions.h"// Cnpk
#include"../../lib_math/math.functions.h"// Cnpk

//#define __aabb_values//

#ifdef  __aabb_values
#define __test_hankel_alpha_beta
#ifdef  __test_hankel_alpha_beta
#include<iostream>
#include<iomanip>
#endif//__test_hankel_alpha_beta
#endif//__aabb_values

namespace hankel
{
	template<class float_type> float_type coef  (int i, int l);
	template<class float_type> float_type alpha (int i, int l);
	template<class float_type> float_type beta  (int i, int l);
	template<class float_type> float_type bessel(int l, float_type const & x);
	template<class float_type> float_type bessel_e(int l, float_type const & x);
};

template<class T>
T hankel::bessel_e(int l, T const & x)
{
	if( x == T(0) ) return 0;
	T _2x = 2 * x;
	T exp_2x = (l%2?-exp(-_2x): exp(-_2x));
	T _i2x = T(1) / _2x, pow_i2x = _i2x;
	T v = T(0);
	for(int i = 0; i <= l; ++i)
	{
		v += Cnpk<T>(l,i) * ((i%2?-1:1) - exp_2x) * pow_i2x;
		pow_i2x *= _i2x;
	}
	return v;
}

template<class T>
T hankel::bessel(int l, T const & x)
{
	T v = T(0);
	T sinh_x = sinh(x), cosh_x = cosh(x);
	T _ix = T(1) / x, pow_ix = _ix;
	for(int i = 0; i <= l; ++i)
	{
		v += (hankel::beta<T>(i,l) * sinh_x - hankel::alpha<T>(i,l) * cosh_x) * pow_ix;
		pow_ix *= _ix;
	}
	return (l%2? -v : v);
}

template<class T>
T hankel::coef(int i, int l)
{
	return Cnpk<T>(l, i) / T(uint32_t(1)<<i);// (l+i)! / (i! * (l - i)! * (2^i))
	//return fact<T>(l + i) / (fact<T>(i) * fact<T>(l-i) * T(1<<i) );// (l+i)! / (i! * (l - i)! * (2^i))
}

template<class T>
T hankel::alpha(int i, int l)
{
	if( (l + i) % 2 )
		return hankel::coef<T>(i, l);
	return 0;
}

template<class T>
T hankel::beta(int i, int l)
{
	if( (l + i) % 2 )
		return 0;
	return hankel::coef<T>(i, l);
}

#ifdef  __aabb_values
// alpha * alpha + beta * beta
template<class T>
T aabb_p(int s, int l, int t, int l_)
{
	return hankel::alpha<T>(s, l) * hankel::alpha<T>(t, l_) + hankel::beta<T>(s, l) * hankel::beta<T>(t, l_);
}

// alpha * alpha - beta * beta
template<class T>
T aabb_m(int s, int l, int t, int l_)
{
	return hankel::alpha<T>(s, l) * hankel::alpha<T>(t, l_) - hankel::beta<T>(s, l) * hankel::beta<T>(t, l_);
}

// alpha * beta + beta * alpha
template<class T>
T abba_p(int s, int l, int t, int l_)
{
	return hankel::alpha<T>(s, l) * hankel::beta<T>(t, l_) + hankel::beta<T>(s, l) * hankel::alpha<T>(t, l_);
}

// alpha * beta - beta * alpha
template<class T>
T abba_m(int s, int l, int t, int l_)
{
	return hankel::alpha<T>(s, l) * hankel::beta<T>(t, l_) - hankel::beta<T>(s, l) * hankel::alpha<T>(t, l_);
}
#endif//__aabb_values

#ifdef  __aabb_values
#ifdef  __test_hankel_alpha_beta
template<class T>
void test_alpha(std::ostream & out, int l, int l_)
{
	int l1 = l+1, l1_ = l_+1;
	int w = 30, p = 20;
	out << std::setw(4) << l << std::setw(4) << l_ << std::endl << std::endl;
	out.setf( std::ios::scientific );
	for(int i = 1; i <= l1; ++i)
	{
		for(int j = 1; j <= l1_; ++j)
		{
			out << std::setw(4) << i << std::setw(4) << j << std::endl;
			out << std::setw(8) << "aabb_p" << std::setw( w ) << std::setprecision( p ) << aabb_p<T>(i, l, j, l_) << std::endl;
			out << std::setw(8) << "aabb_m" << std::setw( w ) << std::setprecision( p ) << aabb_m<T>(i, l, j, l_) << std::endl;
			out << std::setw(8) << "abba_p" << std::setw( w ) << std::setprecision( p ) << abba_p<T>(i, l, j, l_) << std::endl;
			out << std::setw(8) << "abba_m" << std::setw( w ) << std::setprecision( p ) << abba_m<T>(i, l, j, l_) << std::endl;
		}
		out << std::endl;
	}
	out << std::endl;
	out.unsetf( std::ios::scientific );
}
#endif//__test_hankel_alpha_beta
#endif//__aabb_values


#endif//__HANKEL_HPP__
