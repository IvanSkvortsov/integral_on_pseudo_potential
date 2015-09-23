#ifndef __MATH_FUNCTIONS_H__
#define __MATH_FUNCTIONS_H__
#include"math.constants.h"

// power [n - positive]
template<class T>
T pow_int_pos(T const & v, unsigned int const & n)
{
	T value = v;
	for(int i = 1; i < n; ++i)
		value *= v;
	return value;
}
// power [main]
template<class T>
T pow_int(T const & v, int const & n)
{
	if( n < 0 ) return 1 / pow_int_pos<T>(v, -n);
	if( n == 0 ) return 1;
	return pow_int_pos<T>(v, n);
}
// 2^t fast version, t >= 0;
unsigned int _2pow(int const & t)
{
	if( t > 31 ) return -1;
	return 1u<<t;
}

template<class T>
T par_fact(int const & x, int const & n)// x * (x + 1) * ... * (x + n - 1)
{
	if( n > 11 )
	{
		T val = T(x);
		for(int i = 1; i < n; ++i)
			val *= (x+i);
		return val;
	}
	int val = x;
	for(int i = 1; i < n; ++i)
		val *= (x+i);
	return val;
}

template<class T>
T fact(int const & n)
{
	if( n > 11 )
	{
		T val = T(fact_11);
		for(int i = 12; i <= n; ++i)
			val *= i;
		return val;
	}
	int val_i = 1;
	for(int i = 2; i <= n; ++i)
		val_i *= i;
	return val_i;
}

template<class T>
T factt(int const & n)// (2 * n + 1)!!
{
	if( !( n%2 ) || n < 0)
		return 1;
	if( n == 1 ) return 1;
	if( n == 3 ) return 3;
	if( n > 19 )
	{
		T tmp_f19 = T(654729075);// TODO: include this constant to 'math.constants.h'
		for(int i = 21; i <= n; i += 2)
			tmp_f19 *= i;
		return tmp_f19;
	}
	int ival = 3;
	for(int i = 5; i <= n; i+=2)
		ival *= i;
	return T(ival);
}

int logi(int const & n, int const & base)
{
	if( n < 0 ) return -1;
	int a = n;
	int pow_base = 0;
	while( a >= base )
	{
		a /= base;
		++pow_base;
	}
	return pow_base;
}

template<class T>
T NewtonC_n18(int const & n, int const & max_k, int const & min_k)// n < 18
{
	int ival = 1;
	for(int i = n; i > max_k; --i)
		ival *= i;
	return ival / fact<int>(min_k);
}

template<class T>
T NewtonC_n18(int const & n, int const & k)// n < 18
{
	const int min_k = ( (n - k) < k ? n - k : k ), max_k = n - min_k;
	return NewtonC_n18<T>(n, max_k, min_k);
}

template<class T>
T NewtonC_min9(int const & n, int const & max_k, int const & min_k)// min_k <= 9 && n >= 18
{
	int int_maxpow = logi(int32_max, n);
	if( int_maxpow < min_k )
	{
		int ival = n;
		for(int i = 1; i < int_maxpow; ++i)
			ival *= (n - i);
		T fact_val = T(ival);
		for(int i = int_maxpow; i < min_k; ++i)
			fact_val *= (n - i);
		return fact_val / fact<int>(min_k);
	}
	int ival = 1;
	for(int i = n; i > max_k; --i)
		ival *= i;
	return ival / fact<int>(min_k);
}

template<class T>
T NewtonC_ming9(int const & n, int const & max_k, int const & min_k)// n >= 18 && min_k > 9
{
	T min_fact_val = T(fact_10);
	for(int i = 11; i <= min_k; ++i)
		min_fact_val *= i;
	int int_maxpow = logi(int32_max, n);// ~ 7
	int ival = n;
	for(int i = 1; i < int_maxpow; ++i)
		ival *= (n - i);
	T fact_val = T(ival);// n * (n - 1) * ... * (n - int_maxpow + 1);
	for(int i = int_maxpow; i < min_k; ++i)
		fact_val *= (n - i);
	return fact_val / min_fact_val;
}

template<class T>
T NewtonC(int const & n, int const & k)// n!/k!/(n-k)!
{
	const int min_k = ( (n - k) < k ? n - k : k ), max_k = n - min_k;
	if( min_k == 0 ) return 1;
	else if( min_k == 1 ) return n;
	//
	if( n < 18 )
	{
		return NewtonC_n18<T>(n, max_k, min_k);
	}
	else if( min_k <= 9 ) // n >= 18
	{
		return NewtonC_min9<T>(n, max_k, min_k);
	}
	else// min_k > 9 && n >= 20
	{
		return NewtonC_ming9<T>(n, max_k, min_k);
	}
	return 0;
}

template<class T>
T Cn2k(int const & n, int const & k)// n!/k!/(n - 2 * k)!
{
	T NewtC = NewtonC<T>(n, k);
	const int nmk = n - k, nm2k = nmk - k;
	if( nmk < 18 )
	{
		if( k <= 9 )
		{
			int ival = 1;
			for(int i = nmk; i > nm2k; --i)
				ival *= i;
			return T(ival) * NewtC;
		}
		else
		{
			int int_maxpow = logi(int32_max, nmk);// ~ 7
			int ival = nmk;
			for(int i = 1; i < int_maxpow; ++i)
				ival *= (nmk - i);
			T fact_val = T(ival);// n * (n - 1) * ... * (n - int_maxpow + 1);
			for(int i = int_maxpow; i < k; ++i)
				fact_val *= (nmk - i);
			return fact_val * NewtC;
		}
	}
	else if( nm2k <= 9 ) // n >= 18
	{
		int int_maxpow = logi(int32_max, nmk);
		if( int_maxpow < k )
		{
			int ival = nmk;
			for(int i = 1; i < int_maxpow; ++i)
				ival *= (nmk - i);
			T fact_val = T(ival);
			for(int i = int_maxpow; i < k; ++i)
				fact_val *= (nmk - i);
			return fact_val * NewtC;
		}
		int ival = 1;
		for(int i = nmk; i > nm2k; --i)
			ival *= i;
		return ival * NewtC;
	}
	else// nm2k > 9 && n >= 20
	{
		int int_maxpow = logi(int32_max, nmk);// ~ 7
		int ival = nmk;
		for(int i = 1; i < int_maxpow; ++i)
			ival *= (nmk - i);
		T fact_val = T(ival);// n * (n - 1) * ... * (n - int_maxpow + 1);
		for(int i = int_maxpow; i < k; ++i)
			fact_val *= (nmk - i);
		return fact_val * NewtC;
	}
	return 0;
}

template<class T>
T Cnpk(int const & n, int const & k)// (n+k)!/(k! * (n-k)!)
{
	return Cn2k<T>( n+k, k );// n!/(k! * (n-2*k)!), here n = n+k, k = k;
}

//#define tgamma_05_CHECK_IF_N_POSITIVE

#ifdef  tgamma_05_CHECK_IF_N_POSITIVE
#include<cstdlib>// exit
#include<iostream>// cerr
#endif

template<class T>
T tgamma_05(int const & n)// tgamma( n + 0.5 ) / sqrt( Pi )
{
	if( n < 0 )
	{
#ifdef  tgamma_05_CHECK_IF_N_POSITIVE
		std::cerr << "Error: [tgamma_05] in tgamma_05( const int & n ) n must be positive " << std::endl;
		std::cerr << "n : " << n << std::endl;
		std::exit(1);
#endif
		return 1;
	}
	if( n == 0 ) return T(1);
	//if( n == 1 ) return 1 / T(2);
	T tmp_v = T(1), t_05 = 1 / T(2);
	for(int i = n; i > 0; i--)
		tmp_v *= ( T(i) - t_05 );
	return tmp_v;
}

#endif//__MATH_FUNCTIONS_H__
