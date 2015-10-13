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
T Qu(int const & N, int const & lmb, T const & k, T const * kappa, T const * eta);
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
void Qu_run(std::vector<T> & v, int const & N_size, int const & lmb_size, int const& n, T const & k, T const * kappa, T const * eta);
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
	int N = n + 2;
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
		T const * kappa, T const * eta)
{
	T * p = v.data();
	int N = n + 2;
	for(int N_n = 0; N_n < N_size; ++N_n)
	{
		for(int lmb = 0; lmb < lmb_size; ++lmb)
		{
			*p++ = qu_integral::Qu<T>(N, lmb, k, kappa, eta);
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
	//int la_p = lmb_a + 1, lb_p = lmb_b + 1;
	T pow_ika_ = T(1), pow_ikab_, _ika = T(1)/ka, _ikb = T(1)/kb;
	T aa, ab, bb, ba, a_, b_, _a, _b;
	T val = 0;
	N_i = N - 2;
	for(int i = 0; i <= lmb_a; ++i)
	{
		pow_ika_ *= _ika;
		a_ = hankel::alpha<T>(i, lmb_a);
		b_ = hankel::beta <T>(i, lmb_a);
		pow_ikab_ = pow_ika_;
		N_ij = N_i;// N - i - j - 2
		for(int j = 0; j <= lmb_b; ++j)
		{
			_a = hankel::alpha<T>(j, lmb_b);
			_b = hankel::beta <T>(j, lmb_b);

			aa = a_ * _a;
			ab = a_ * _b;
			bb = b_ * _b;
			ba = b_ * _a;

			pow_ikab_ *= _ikb;
			val += ( aa * tau[N_ij] + bb * rho[N_ij] - ba * sigma[N_ij] - ab * _sigma[N_ij] ) * pow_ikab_;
			--N_ij;
		}
		--N_i;
	}
	return val * ((lmb_a - lmb_b)%2 ? -1 : 1);
}

#include<stdlib.h>

template<class T> 
T qu_integral::Qu(int const & N, int const & lmb, T const & k, T const * kappa, T const * eta)
{
	if( lmb == 0 ) return eta[N-1]/k;
	if( lmb == 1 )
	{
		T _ik = T(1) / k;
		return (kappa[N-1] - eta[N-2] * _ik) * _ik;
	}
	if( lmb == 2 )
	{
		T _ik = T(1) / k;
		//return (eta[N-1] + 3 * (-kappa[N-2] + eta[N-3] * _ik) * _ik ) * _ik;
		return ((eta[N-3] * _ik - kappa[N-2]) * T(3) * _ik + eta[N-1] ) * _ik;
	}
	if( lmb == 3 )
	{
		T _ik = T(1) / k;
		//return (((kappa[N-3] - eta[N-4] * _ik) * 15 * _ik - eta[N-2]) * 6 * _ik + kappa[N-1]) * _ik;
		return (kappa[N-1] + (-6 * eta[N-2] + 15 * (kappa[N-3] - eta[N-4] * _ik) * _ik) * _ik) * _ik;
		//return (kappa[N-1] - (6 * eta[N-2] - (15 * (kappa[N-3] - eta[N-4] * _ik) * _ik) * _ik) * _ik;
	}
	int N_i = N - 1;
	T pow_ik_ = T(1), _ik = T(1)/k;
	T a = 1, b = 1;
	T val = T(0);
	for(int i = 0; i <= lmb; ++i)
	{
		pow_ik_ *= _ik;
		a = hankel::alpha<T>(i, lmb);// i:=0 if( lmb%2 == 1 )
		b = hankel::beta <T>(i, lmb);// i:=0 if( lmb%2 == 0 )
		val += (b * eta[N_i] - a * kappa[N_i]) * pow_ik_;
		//val += (eta[N_i] + kappa[N_i]);
		//val += (a * kappa[N_i] - b * eta[N_i]) * pow_ik_;
		--N_i;
	}
	//return 1;
	//return val;
	//return val * (lmb%2 ? 1 : -1);
	return val * (lmb%2 ? -1 : 1);
}


#endif//__QU_H__
