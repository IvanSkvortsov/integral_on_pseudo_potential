#ifndef __HANKEL_HPP__
#define __HANKEL_HPP__
#include"../../../lib_math/math.functions.h"// Cnpk

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
	template<class float_type> float_type bessel(int i, int l);
	template<class float_type> float_type alpha (int i, int l);
	template<class float_type> float_type beta  (int i, int l);
};

template<class T>
T hankel::bessel(int i, int l)
{
	int i1 = i - 1;
	return Cnpk<T>(l, i1) / T(1<<i1);// (l+i-1)! / ((i-1)! * (l - (i-1))!) / (2^(i-1))
}

template<class T>
T hankel::alpha(int i, int l)
{
	if( (l - i) % 2 == 0 )
		return hankel::bessel<T>(i, l);
	return 0;
}

template<class T>
T hankel::beta(int i, int l)
{
	if( (l - i) % 2 == 0 )
		return 0;
	return hankel::bessel<T>(i, l);
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
