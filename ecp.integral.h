#ifndef __ECP_INTEGRAL_H__
#define __ECP_INTEGRAL_H__

#include"angular_integral/omega.integral.h"
#include"radial_integral/radial.h"
#include"../lib_electron/primitive.h"
#include"ecp.primitive.h"

#ifndef __ECP_PRIMITIVE_H__
#define __ECP_PRIMITIVE_H__
template<class T>
struct ecp_primitive
{
	T alp, r[3];
	int n, l;
	ecp_primitive():alp(0), n(0), l(0)
	{
		for(int i = 0; i < 3; ++i) r[i] = T(0);
	}
};
#endif//__ECP_PRIMITIVE_H__

#include<iostream>
#include<iomanip>
#include<fstream>

template<class T>
void read_task(Primitive<T> & a, Primitive<T> & b, vector<ecp_primitive<T> > & v_ecp, char const * file)
{
	std::ifstream inp( file );
	for(int i = 0; i < 3; ++i) inp >> a.i[i];
	for(int i = 0; i < 3; ++i) inp >> a.r[i];
	inp >> a.alp;
	for(int i = 0; i < 3; ++i) inp >> b.i[i];
	for(int i = 0; i < 3; ++i) inp >> b.r[i];
	inp >> b.alp;
	int size = 1;
	inp >> size;
	v_ecp.resize( size );
	for(int k = 0; k < size; ++k)
	{
		inp >> v_ecp[k].l >> v_ecp[k].n;
		for(int i = 0; i < 3; ++i) inp >> v_ecp[k].r[i];
		inp >> v_ecp[k].alp;
	}
	inp.close();
	return;
}

template<class T>
bool is_v3_equal(T const * A, T const * B)
{
	return A[0] == B[0] && A[1] == B[1] && A[2] == B[2];
}

template<class T>
int integral_num(T const * A, T const * B, T const * C)
{
	if( !is_v3_equal<T>( A, C ) && !is_v3_equal<T>( B, C ) )// A != C && C != B
		return 1;
	if( is_v3_equal<T>( A, C ) && !is_v3_equal<T>( B, C ) )//  A == C && C != B
		return 2;
	if( !is_v3_equal<T>( A, C ) && is_v3_equal<T>( B, C ) )//  A != C && C == B
		return -2;
	if( is_v3_equal<T>( A, C ) && is_v3_equal<T>( B, C ) )//   A == C && C == B
		return 3;
}

template<class T>
void kx_calc(T * kx, T * CX, T const * X, T const * C, T const & alp_x)
{
	for(int i = 0; i < 3; ++i) CX[i] = C[i] - X[i];
	for(int i = 0; i < 3; ++i) kx[i] = CX[i] * (-2 * alp_x );
}

template<class T>
T ecp_local_integral_1(Primitive<T> const & a, Primitive<T> const & b, ecp_primitive<T> const & ecp)
{
	// angular
	std::vector<T> ang;
	int la = a.i[0] + a.i[1] + a.i[2];
	int lb = b.i[0] + b.i[1] + b.i[2];
	int l = ecp.l;
	int la_size, lb_size, lmb_size;
	omega_integral::omega_resize_1<T>( ang, la_size, lb_size, lmb_size, la, lb);
	T ka[3], CA[3];
	kx_calc<T>( ka, CA, a.r, ecp.r, a.alp );
	T kb[3], CB[3];
	kx_calc<T>( kb, CB, b.r, ecp.r, b.alp );
	T k[3];
	for(int i = 0; i < 3; ++i) k[i] = ka[i] + kb[i];
	omega_integral::omega_run_1<T>( ang, la_size, lb_size, lmb_size, a.i, b.i, k, CA, CB);
	// radial
	std::vector<T> rad;
	int N_size;//, lmb_size;
	int n = ecp.n;
	qu_run_1<T>( rad, N_size, lmb_size, a.r, b.r, ecp.r, la, lb, n, a.alp, b.alp, ecp.alp );
	// ecp integral
	T value = T(0);
	for(int na = 0; na < la_size; ++na)
	{
		for(int nb = 0; nb < lb_size; ++nb )
		{
			for(int lmb = 0; lmb <= na + nb; ++lmb )
			{
				value += ang[(na * lb_size + nb) * lmb_size + lmb] * rad[(na + nb) * lmb_size + lmb];
			}
		}
	}
	//return value;
	T _4Pi = 4 * T(Pi);
	T sqr_CA = T(0); for(int i = 0; i < 3; ++i) sqr_CA += CA[i] * CA[i];
	T sqr_CB = T(0); for(int i = 0; i < 3; ++i) sqr_CB += CB[i] * CB[i];
	T exp_ca = exp( -a.alp * sqr_CA ), exp_cb = exp( -b.alp * sqr_CB );
	return _4Pi * exp_ca * exp_cb * value;
}

template<class T>
T ecp_semi_local_integral_1(Primitive<T> const & a, Primitive<T> const & b, ecp_primitive<T> const & ecp)
{
	// angular
	std::vector<T> ang;
	int la = a.i[0] + a.i[1] + a.i[2];
	int lb = b.i[0] + b.i[1] + b.i[2];
	int l = ecp.l;
	int la_size, lb_size, lmb_a_size, lmb_b_size;
	omega_integral::omega_resize_1<T>( ang, la_size, lb_size, lmb_a_size, lmb_b_size, la, lb, l);
	T ka[3], CA[3];
	kx_calc<T>( ka, CA, a.r, ecp.r, a.alp );
	T kb[3], CB[3];
	kx_calc<T>( kb, CB, b.r, ecp.r, b.alp );
	omega_integral::omega_run_1<T>( ang, la_size, lb_size, lmb_a_size, lmb_b_size, a.i, b.i, l, ka, kb, CA, CB);
	// radial
	std::vector<T> rad;
	int N_size;//, lmb_a_size, lmb_b_size;
	int n = ecp.n;
	qu_run_1<T>( rad, N_size, lmb_a_size, lmb_b_size, a.r, b.r, ecp.r, la, lb, l, n, a.alp, b.alp, ecp.alp );
	// ecp integral
	T value = T(0);
	for(int na = 0; na < la_size; ++na)
	{
		for(int nb = 0; nb < lb_size; ++nb )
		{
			for(int lmb_a = 0; lmb_a <= l + na; ++lmb_a )
			{
				for(int lmb_b = 0; lmb_b <= l + nb; ++lmb_b )
					value += ang[( (na * lb_size + nb) * lmb_a_size + lmb_a ) * lmb_b_size + lmb_b] *
						rad[((na + nb) * lmb_a_size + lmb_a) * lmb_b_size + lmb_b];
			}
		}
	}
	//return value;
	T _4_pi = 4 * T(Pi);
	T sqr_CA = T(0); for(int i = 0; i < 3; ++i) sqr_CA += CA[i] * CA[i];
	T sqr_CB = T(0); for(int i = 0; i < 3; ++i) sqr_CB += CB[i] * CB[i];
	T exp_v = exp( -a.alp * sqr_CA - b.alp * sqr_CB );
	return _4_pi * exp_v * value;
}

template<class T>
T ecp_semi_local_integral_2(Primitive<T> const & a, Primitive<T> const & b, ecp_primitive<T> const & ecp)
{
	// angular
	std::vector<T> ang;
	int la = a.i[0] + a.i[1] + a.i[2];
	int lb = b.i[0] + b.i[1] + b.i[2];
	int l = ecp.l;
	int lb_size, lmb_size;
	omega_integral::omega_resize_2<T>( ang, lb_size, lmb_size, la, lb, l);
	T ka[3], CA[3];
	kx_calc<T>( ka, CA, a.r, ecp.r, a.alp );
	T kb[3], CB[3];
	kx_calc<T>( kb, CB, b.r, ecp.r, b.alp );
	omega_integral::omega_run_2<T>( ang, lb_size, lmb_size, a.i, b.i, l, kb, CB);
	// radial
	std::vector<T> rad;
	int N_size;//, lmb_a_size, lmb_b_size;
	int n = ecp.n;
	qu_run_2<T>( rad, N_size, lmb_size, b.r, ecp.r, la, lb, l, n, a.alp, b.alp, ecp.alp );
	// ecp integral
	T value = T(0);
	for(int nb = 0; nb < lb_size; ++nb )
	{
		for(int lmb_b = 0; lmb_b <= l + nb; ++lmb_b )
			value += ang[nb * lmb_size + lmb_b] * rad[nb * lmb_size + lmb_b];
	}
	//return value;
	T _4Pi = 4 * T(Pi);
	T sqr_CB = T(0); for(int i = 0; i < 3; ++i) sqr_CB += CB[i] * CB[i];
	T exp_cb = exp( -b.alp * sqr_CB );
	return _4Pi * exp_cb * value;
}

template<class T>
T ecp_semi_local_integral_3(Primitive<T> const & a, Primitive<T> const & b, ecp_primitive<T> const & ecp)
{
	// orbital momenta
	int la = a.i[0] + a.i[1] + a.i[2];
	int lb = b.i[0] + b.i[1] + b.i[2];
	int l = ecp.l;
	// angular integral
	vector<spherical<T> > v_slm;
	v_slm.resize( 2 * l + 1 );
	for(int i = 0; i < v_slm.size(); ++i)
	{
		v_slm[i].run(l, i-l);
		v_slm[i].optimize_ez();
	}
	angular_omega_xyz<T> omg_xyz_;
	int lx_max = la > lb ? la : lb;
	int i_sz = l + lx_max;
	omg_xyz_.run(i_sz, i_sz, i_sz);
	spherical<T> sph_buf;
	sph_buf.reserve( 40 );
	T sa, sb;
	typename spherical<T>::polynomial_type pol_a(T(1), a.i[0], a.i[1], a.i[2]), pol_b(T(1), b.i[0], b.i[1], b.i[2]), * p;
	T ang = 0;
	for(int i = 0; i < v_slm.size(); ++i)
	{
		sa = 0;
		sph_buf = v_slm[i];// m = i-l;
		sph_buf *= pol_a;
		p = &sph_buf[0];
		for(int j = 0; j < sph_buf.size(); ++j)
		{
			if( !p->is_even_x() ) continue;
			sa += p->d * omg_xyz_(p->x, p->y, p->z);
			++p;
		}
		sb = 0;
		sph_buf = v_slm[i];// m = i-l;
		sph_buf *= pol_b;
		p = &sph_buf[0];
		for(int j = 0; j < sph_buf.size(); ++j)
		{
			if( !p->is_even_x() ) continue;
			sb += p->d * omg_xyz_(p->x, p->y, p->z);
			++p;
		}
		ang += sa * sb;
	}
	T _4Pi = 4 * T(Pi);
	ang *= 1;//_4Pi * _4Pi;
	// radial integral
	int N = la + lb + ecp.n;
	T alp = a.alp + b.alp + ecp.alp;
	T rad = (N%2 ? fact<T>((N-1)/2)*T(0.5)/pow_int<T>(alp, (N+1)/2) : tgamma_05<T>(N/2)*T(sqrtPi) * T(0.5) / (pow_int<T>(alp, N/2) * sqrt(alp))); 
	return ang * rad;
}

template<class T>
T ecp_local_integral(Primitive<T> const & a, Primitive<T> const & b, ecp_primitive<T> const & ecp)
{
	int i_num = integral_num<T>( a.r, b.r, ecp.r );
	//std::cout << "integral_num : " << std::setw(3) << i_num << std::endl;
	if( i_num == 1 ) return ecp_local_integral_1<T>( a, b, ecp );
	if( i_num == 2 ) return -1;//ecp_local_integral_2<T>( a, b, ecp );
	if( i_num == -2 ) return -2;//ecp_local_integral_2<T>( b, a, ecp );
	if( i_num == 3 ) return -3;//ecp_local_integral_3<T>( b, a, ecp );
	return T(-1);
}

template<class T>
T ecp_semi_local_integral(Primitive<T> const & a, Primitive<T> const & b, ecp_primitive<T> const & ecp)
{
	int i_num = integral_num<T>( a.r, b.r, ecp.r );
	//std::cout << "integral_num : " << std::setw(3) << i_num << std::endl;
	if( i_num == 1 ) return ecp_semi_local_integral_1<T>( a, b, ecp );
	if( i_num == 2 ) return ecp_semi_local_integral_2<T>( a, b, ecp );
	if( i_num == -2 ) return ecp_semi_local_integral_2<T>( b, a, ecp );
	if( i_num == 3 ) return ecp_semi_local_integral_3<T>( b, a, ecp );
	return T(-1);
}

#endif//__ECP_INTEGRAL_H__
