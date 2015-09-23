#ifndef __RADIAL_H__
#define __RADIAL_H__

#include"lib_qu/qu.h"
#include"recursion/taus.h"

// semi local : A != C != B
template<class T>
void qu_run_1(std::vector<T> & v, int & N_size, int & lmb_a_size, int & lmb_b_size, T const * A, T const * B, T const * C,
		int const & la, int const & lb, int const & l, int const & n, T const & alp_i, T const & alp_j, T const & alp_k)
{
	qu_integral::Qu_resize_1(v, N_size, lmb_a_size, lmb_b_size, la, lb, l);
	recursion<T> recurs;
	recurs_resize_1<T>( recurs, la, lb, l, n );
	recursion_elem<T> elc;
	elc.init_a( alp_i + alp_j + alp_k );
	T ka[3];
	// ka
	for(int i = 0; i < 3; ++i) ka[i] = C[i] - A[i];
	T kb[3];
	// kb
	for(int i = 0; i < 3; ++i) kb[i] = C[i] - B[i];
	// ka_len
	T ka_len = T(0);
	for(int i = 0; i < 3; ++i) ka_len += (ka[i] * ka[i]);
	ka_len = sqrt( ka_len ) * 2 * alp_i;
	// kb_len
	T kb_len = T(0);
	for(int i = 0; i < 3; ++i) kb_len += (kb[i] * kb[i]);
	kb_len = sqrt( kb_len ) * 2 * alp_j;
	//
	for(int i = 0; i < 3; ++i) ka[i] *= (-2 * alp_i);
	for(int i = 0; i < 3; ++i) kb[i] *= (-2 * alp_j);
	// recursion_elem
	elc.run_k( ka_len, kb_len );
	elc.run_1();
	// recursion calculation
	recurs.run_1( elc );
	std::vector<T> taus;
	int begin = 0, end = 0;
	T * tau, * sigma, * _sigma, * rho;
	taus_run_1<T>( taus, begin, end, tau, sigma, _sigma, rho, recurs, elc.sgn__ka_m_kb );
	qu_integral::Qu_run(v, N_size, lmb_a_size, lmb_b_size, n, ka_len, kb_len, tau, sigma, _sigma, rho);
}
// semi local : A == C != B
template<class T>
void qu_run_2(std::vector<T> & v, int & N_size, int & lmb_size, T const * B, T const * C,
		int const & la, int const & lb, int const & l, int const & n, T const & alp_i, T const & alp_j, T const & alp_k)
{
	qu_integral::Qu_resize_2(v, N_size, lmb_size, la, lb, l);
	recursion<T> recurs;
	recurs_resize_2<T>( recurs, la, lb, l, n );
	recursion_elem<T> elc;
	elc.init_a( alp_i + alp_j + alp_k );
	T kb[3];
	// kb
	for(int i = 0; i < 3; ++i) kb[i] = C[i] - B[i];
	// kb_len
	T kb_len = T(0);
	for(int i = 0; i < 3; ++i) kb_len += (kb[i] * kb[i]);
	kb_len = sqrt( kb_len ) * 2 * alp_j;
	//
	for(int i = 0; i < 3; ++i) kb[i] *= (-2 * alp_j);
	// recursion_elem
	elc.run_k( 0, kb_len );
	elc.run_2();
	// recursion calculation
	recurs.run_2( elc );
	std::vector<T> taus;
	int begin = 0, end = 0;
	T * kappa, * etta;
	taus_run_2<T>( taus, begin, end, kappa, etta, recurs );
	qu_integral::Qu_run(v, N_size, lmb_size, n+la, kb_len, kappa, etta);
}

// local : A != C != B
template<class T>
void qu_run_1(std::vector<T> & v, int & N_size, int & lmb_size, T const * A, T const * B, T const * C,
		int const & la, int const & lb, int const & n, T const & alp_i, T const & alp_j, T const & alp_k)
{
	recursion<T> recurs;
	recurs_resize_1<T>( recurs, la, lb, n );
	recursion_elem<T> elc;
	elc.init_a( alp_i + alp_j + alp_k );
	T ka[3];
	// ka
	for(int i = 0; i < 3; ++i) ka[i] = C[i] - A[i];
	for(int i = 0; i < 3; ++i) ka[i] *= (-2 * alp_i);
	T kb[3];
	// kb
	for(int i = 0; i < 3; ++i) kb[i] = C[i] - B[i];
	for(int i = 0; i < 3; ++i) kb[i] *= (-2 * alp_j);
	// k
	for(int i = 0; i < 3; ++i) ka[i] += kb[i];
	// k_len
	T k_len = T(0);
	for(int i = 0; i < 3; ++i) k_len += (ka[i] * ka[i]);
	k_len = sqrt( k_len );
	// recursion_elem
	elc.run_k( 0, k_len );
	elc.run_2();
	// recursion calculation
	recurs.run_2( elc );
	std::vector<T> taus;
	int begin = 0, end = 0;
	T * kappa, * etta;
	taus_run_2<T>( taus, begin, end, kappa, etta, recurs );
	qu_integral::Qu_resize_1(v, N_size, lmb_size, la, lb);
	qu_integral::Qu_run(v, N_size, lmb_size, n, k_len, kappa, etta);
}
// local : A == C != B
template<class T>
void qu_run_2(std::vector<T> & v, int & N_size, int & lmb_size, T const * B, T const * C,
		int const & la, int const & lb, int const & n, T const & alp_i, T const & alp_j, T const & alp_k)
{
	recursion<T> recurs;
	recurs_resize_2<T>( recurs, la, lb, n );
	recursion_elem<T> elc;
	elc.init_a( alp_i + alp_j + alp_k );
	T kb[3];
	// kb
	for(int i = 0; i < 3; ++i) kb[i] = C[i] - B[i];
	// kb_len
	T kb_len = T(0);
	for(int i = 0; i < 3; ++i) kb_len += (kb[i] * kb[i]);
	kb_len = sqrt( kb_len ) * 2 * alp_j;
	//
	for(int i = 0; i < 3; ++i) kb[i] *= (-2 * alp_j);
	// recursion_elem
	elc.run_k( 0, kb_len );
	elc.run_2();
	// recursion calculation
	recurs.run_2( elc );
	std::vector<T> taus;
	int begin = 0, end = 0;
	T * kappa, * etta;
	taus_run_2<T>( taus, begin, end, kappa, etta, recurs );
	qu_integral::Qu_resize_2(v, N_size, lmb_size, la, lb);
	qu_integral::Qu_run(v, N_size, lmb_size, n+la, kb_len, kappa, etta);
}

#endif//__RADIAL_H__
