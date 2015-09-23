#ifndef __QU_H__
#define __Qu_H__
#include"hankel.h"
#include<vector>
#define __sqrt_pi 1.7724538509055160272981674833411

namespace qu_integral
{
//--- calculation ---//
// Qu(N, lmb_a, lmb_b) / Qu(N, lmb)
template<class T>
T Qu(int const & N, int const & lmb_a, int const & lmb_b, T const & ka, T const & kb, T const * tau, T const * sigma, T const * _sigma, T const * rho);
template<class T>
T Qu(int const & N, int const & lmb, T const & k, T const * kappa, T const * etta);
//----- resize ------//
// semi local
template<class T>
void Qu_resize_1(std::vector<T> & v, int & N_size, int & lmb_a_size, int & lmb_b_size, int const & la, int const & lb, int const& l);
template<class T>
void Qu_resize_2(std::vector<T> & v, int & N_size, int & lmb_size, int const & la, int const & lb, int const & l);
// local
template<class T>
void Qu_resize_1(std::vector<T> & v, int & N_size, int & lmb_size, int const & la, int const & lb);
template<class T>
void Qu_resize_2(std::vector<T> & v, int & N_size, int & lmb_size, int const & la, int const & lb);
//------- run -------//
template<class T>
void Qu_run(std::vector<T> & v, int const & N_size, int const & lmb_a_size, int const & lmb_b_size, int const& n, T const & ka, T const & kb,
	T const * tau, T const * sigma, T const * _sigma, T const * rho);
template<class T>
void Qu_run(std::vector<T> & v, int const & N_size, int const & lmb_size, int const& n, T const & k, T const * kappa, T const * etta);
};

//---------------------- resize -----------------------//
                //- ( A != C != B ) -//
// semi local
template<class T>
void qu_integral::Qu_resize_1(std::vector<T> & v, int & N_size, int & lmb_a_size, int & lmb_b_size, int const & la, int const & lb, int const& l)
{
	N_size = la + lb + 1;
	lmb_a_size = l + la + 1;
	lmb_b_size = l + lb + 1;
	v.resize( N_size * lmb_a_size * lmb_b_size );
}
// local
template<class T>
void qu_integral::Qu_resize_1(std::vector<T> & v, int & N_size, int & lmb_size, int const & la, int const & lb)
{
	N_size = la + lb + 1;
	lmb_size = la + lb + 1;
	v.resize( N_size * lmb_size );
}
                //- ( A == C != B ) -//
// semi local
template<class T>
void qu_integral::Qu_resize_2(std::vector<T> & v, int & N_size, int & lmb_size, int const & la, int const & lb, int const & l)
{
	N_size = lb + 1;
	lmb_size = l + lb + 1;
	v.resize( N_size * lmb_size );
}
// local
template<class T>
void qu_integral::Qu_resize_2(std::vector<T> & v, int & N_size, int & lmb_size, int const & la, int const & lb)
{
	N_size = lb + 1;
	lmb_size = la + lb + 1;
	v.resize( N_size * lmb_size );
}

//------------------------ run ------------------------//
                //- ( A != C != B ) -//
template<class T>
void qu_integral::Qu_run(std::vector<T> & v, int const & N_size, int const & lmb_a_size, int const & lmb_b_size, int const& n,
		T const & ka, T const & kb, T const * tau, T const * sigma, T const * _sigma, T const * rho)
{
	T * p = v.data();
	int N = n;
	for(int N_n = 0; N_n < N_size; ++N_n)
	{
		for(int lmb_a = 0; lmb_a < lmb_a_size; ++lmb_a)
		{
			for(int lmb_b = 0; lmb_b < lmb_b_size; ++lmb_b)
			{
				*p++ = qu_integral::Qu<T>(N, lmb_a, lmb_b, ka, kb, tau, sigma, _sigma, rho);
			}
		}
		++N;
	}
}
template<class T>
void qu_integral::Qu_run(std::vector<T> & v, int const & N_size, int const & lmb_size, int const& n, T const & k,
		T const * kappa, T const * etta)
{
	T * p = v.data();
	int N = n;
	for(int N_n = 0; N_n < N_size; ++N_n)
	{
		for(int lmb = 0; lmb < lmb_size; ++lmb)
		{
			*p++ = qu_integral::Qu<T>(N, lmb, k, kappa, etta);
		}
		++N;
	}
}
//-------------------- calculation --------------------//
template<class T>
T qu_integral::Qu(int const & N, int const & lmb_a, int const & lmb_b, T const & ka, T const & kb, T const * tau, T const * sigma, 
		T const * _sigma, T const * rho)
{
	//typedef long long int llz;
	int N_i, N_ij;
	int la_p = lmb_a + 1, lb_p = lmb_b + 1;
	T pow_ka_ = T(1), pow_kab_;
	T aa, ab, bb, ba, a_, b_, _a, _b;
	T val = 0;
	for(int i = 1; i <= la_p; ++i)
	{
		N_i = N - i;
		pow_ka_ *= ka;
		a_ = hankel::alpha<T>(i, lmb_a);
		b_ = hankel::beta <T>(i, lmb_a);
		pow_kab_ = pow_ka_;
		for(int j = 1; j <= lb_p; ++j)
		{
			N_ij = N_i - j;
			_a = hankel::alpha<T>(j, lmb_b);
			_b = hankel::beta <T>(j, lmb_b);

			aa = a_ * _a;
			ab = a_ * _b;
			bb = b_ * _b;
			ba = b_ * _a;

			pow_kab_ *= kb;
			val += ( aa * tau[N_ij] + bb * rho[N_ij] - ba * sigma[N_ij] - ab * _sigma[N_ij] ) / pow_kab_;
		}
	}
	return val * ((lmb_a - lmb_b)%2 ? -1 : 1);
}

template<class T> 
T qu_integral::Qu(int const & N, int const & lmb, T const & k, T const * kappa, T const * etta)
{
	int N_i, lp = lmb + 1;
	T pow_k_ = T(1);
	T a, b;
	T val = 0;
	for(int i = 1; i <= lp; ++i)
	{
		N_i = N - i;
		pow_k_ *= k;
		a = hankel::alpha<T>(i, lmb);
		b = hankel::beta <T>(i, lmb);
		val += (a * kappa[N_i] - b * etta[N_i]) / pow_k_;
	}
	return val * (lp%2 ? -1 : 1);
}


#endif//__QU_H__
