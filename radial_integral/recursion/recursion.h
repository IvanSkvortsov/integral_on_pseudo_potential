#ifndef RECURSION_H
#define RECURSION_H
#include"recursion.elem.h"// erfh.h dawson.h l_mymath.h defined_consts.h
#include"recursion.bc.h"// recursion_bc, plus_minus, (p_recursion_bc)
#include<vector>

// recursion - elements of Radial Integral reccursive realization (in Maple semantics) :
//
// B[+] := i -> ( i - 1 ) / ( 2 * alp ) * B[+](i-2) + a[+] * C[+](i-1)
// C[+] := i -> ( i - 1 ) / ( 2 * alp ) * C[+](i-2) + a[+] * B[+](i-1)
// G[+] := 0.25 * exp( a[+]^2 * alp )
//
// B[-] := i -> ( i - 1 ) / ( 2 * alp ) * B[-](i-2) + a[-] * C[-](i-1)
// C[-] := i -> ( i - 1 ) / ( 2 * alp ) * C[-](i-2) + a[-] * B[-](i-1)
// G[-] := 0.25 * exp( a[-]^2 * alp )
//
// a[+] := p
// a[-] := m
//
// p    := k_A + k_B
// m    := |k_A - k_B|
// k_A  := 2 * alp_i * |CA|
// k_B  := 2 * alp_j * |CB|
//
// P.S.
// alp  := alp_i + alp_j + alp_k
//
// alp_i -> Basis_A_alp(i)
// alp_j -> Basis_B_alp(j)
// alp_k -> ECP_alp(k)
// 
// v_bc   - the functions B and C (see above);
// sh_bc  - shifted pointer to &v_bc[0] - min_n;
// _g     - the function G (see above);
// _min_n - relative address of sh_bc[0] with reference address is &v_bc[0]
// _max_n   - sh_bc[_max_n] -eq v_bc[size], where size := _max_n - _min_n;
//
// recursion_bc<T> - structure that contains 2 variables of type plus_minus<T> : b, c;
// plus_minus<T> - structure that also contains 2 variables but type of T: p, m;

//#define __recursion_log__
//#define __recursion_warn_msg__
#define __recursion_err_msg__
#define __recursion_print__

#if defined(__recursion_print__) || defined(__recursion_err_msg__) || defined(__recursion_log__)
#include<iostream>
#include<iomanip>
#endif

template<class T>
class recursion
{
protected:
	std::vector<recursion_bc<T> > v_bc;
	recursion_bc<T> * sh_bc;
	plus_minus<T> _g;
	int _min_n, _max_n;
	void set_zero();
#ifdef  __recursion_log__
	void log(std::string const & s)const;
#endif
public:
	recursion();
	recursion(recursion<T> const & rcsn);
	~recursion();
	// operator=
	recursion<T> & operator=(recursion<T> const & rcsn);
	void reverse(){}
	// std::vector : resize, reserve, size, capacity
	void resize(std::size_t _sz ){v_bc.resize( _sz );}
	void reserve(std::size_t _sz){v_bc.reserve(_sz );}
	const std::size_t size()const{return v_bc.size();}
	const std::size_t capacity()const{return v_bc.capacity();}
	// min_n, max_n
	int min_n(int __min_n){return _min_n = __min_n;}
	int max_n(int __max_n ){return _max_n = __max_n;}
	int & min_n(){return _min_n;}
	int & max_n(){return _max_n;}
	const int & min_n()const{return _min_n;}
	const int & max_n()const{return _max_n;}
	// operator[]
	recursion_bc<T> & operator[](int i){return v_bc[i];}
	recursion_bc<T> const & operator[](int i)const{return v_bc[i];}
	// operator()
	recursion_bc<T> & operator()(int i){return sh_bc[i];}
	recursion_bc<T> const & operator()(int i)const{return sh_bc[i];}
	// p_bc -- pointer to bc
	recursion_bc<T> * p_bc(int i){return &sh_bc[i];}
	recursion_bc<T> const * p_bc(int i)const{return &sh_bc[i];}
	// 
	T & gp(){return _g.p;}
	T & gm(){return _g.m;}
	T const & gp()const{return _g.p;}
	T const & gm()const{return _g.m;}
	// resize
	int resize();
	int resize(int __min_n, int __max_n );
	// run_3
	void run_0_eq_3( recursion_elem<T> const & );
	void run_p1_geq_3( recursion_elem<T> const & );
	void run_m1_eq_3( recursion_elem<T> const &);
	void run_m2_eq_3( recursion_elem<T> const &);
	void run_m3_leq_3( recursion_elem<T> const &);
	void run_3(recursion_elem<T> const & elc );
	void run_2(recursion_elem<T> const & elc ){this->run_3( elc );}
	// run_1
	void run_0_eq( recursion_elem<T> const & );
	void run_p1_geq( recursion_elem<T> const & );
	void run_m1_eq( recursion_elem<T> const &);
	void run_m2_eq( recursion_elem<T> const &);
	void run_m3_leq( recursion_elem<T> const &);
	void run(recursion_elem<T> const & elc );
	void run_1(recursion_elem<T> const & elc){this->run( elc );}
	// print
#ifdef  __recursion_print__
	void print(std::ostream & out = std::cout)const;
#endif
};

template<class T>
void recurs_resize_acb_sl(recursion<T> & recurs, int const & la, int const & lb, int const & l, int const & n)
{
	recurs.resize( n - 2*(l + 1) + 2, la + lb + n - 2 + 2);
}
template<class T>
void recurs_resize_acb_l(recursion<T> & recurs, int const & la, int const & lb, int const & n)
{
	recurs.resize( n - 1 + 2, la + lb + n - 1 + 2);
}

template<class T>
void recurs_resize_ccb_sl(recursion<T> & recurs, int const & la, int const & lb, int const & l, int const & n)
{
	recurs.resize( la + n - l - 1 + 2, la + lb + n - 1 + 2);
}
template<class T>
void recurs_resize_ccb_l(recursion<T> & recurs, int const & la, int const & lb, int const & n)
{
	recurs.resize( n - 1 + 2, la + lb + n - 1 + 2);
}

template<class T>
void recurs_run_1(recursion<T> & recurs, T const * A, T const * B, T const * C, T const & alp_i, T const & alp_j, T const & alp_k)
{
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
}

template<class T>
void recurs_run_2(recursion<T> & recurs, T const * B, T const * C, T const & alp_i, T const & alp_j, T const & alp_k)
{
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
}

// set_zero()
template<class T>
void recursion<T>::set_zero()
{
	sh_bc = 0;
	_min_n = 0;
	_max_n = 0;
}
// recursion()
template<class T>
recursion<T>::recursion():v_bc(), sh_bc(0), _g(), _min_n(0), _max_n(0)
{
#ifdef  __recursion_log__
	log("recursion()");
#endif
}
// recursion( recursion<T> const & )
template<class T>
recursion<T>::recursion(recursion<T> const & rcsn):v_bc(rcsn.v_bc), sh_bc(v_bc.data() - rcsn._min_n), _g(rcsn._g), _min_n(rcsn._min_n), _max_n(rcsn._max_n)
{
#ifdef  __recursion_log__
	log("recursion(recursion<T> const &)");
#endif
}
// ~recursion()
template<class T>
recursion<T>::~recursion()
{
	this->set_zero();
#ifdef  __recursion_log__
	log("~recursion()");
#endif
}
// operator=
template<class T>
recursion<T> & recursion<T>::operator=(recursion<T> const & rcsn)
{
#ifdef  __recursion_log__
	log("operator=(recursion<T> const & )");
#endif
	if( this == &rcsn ) return *this;
	recursion<T> tmp(rcsn);
	swap(v_bc, tmp.v_bc);
	swap(_g, tmp._g);
	swap(_min_n, tmp._min_n);
	swap(_max_n, tmp._max_n);
	sh_bc = &v_bc[0] - _min_n;
	return *this;
}
// resize
template<class T>
int recursion<T>::resize()
{
	//if( _max_n <= _min_n )
	if( _max_n < _min_n )
	{
#ifdef  __recursion_err_msg__
		std::cerr << "Error: recursion<T>::resize()" << std::endl;
		std::cerr << "min_n: " << _min_n << std::endl;
		std::cerr << "max_n: " << _max_n << std::endl;
#endif
		return -1;
	}
	int sz = 0;
	//if( _max_n <= 0 )
	if( _max_n < 0 )
		sz = 1 - _min_n;// include zero's element
	else
		if( _min_n > 0 )
			sz = 1 + _max_n;// include zero's element
		else
			sz = 1 + _max_n - _min_n;// already contains zero's element
	v_bc.resize( sz );
	//sh_bc = &v_bc[0] - _min_n;
	sh_bc = &v_bc[0] - (_min_n > 0 ? 0 : _min_n );
	return 1;
}
// min_n and end values setting with corresponding memory allocation for vector<p_recursion_bc<T> > v_bc;
template<class T>
int recursion<T>::resize(int __min_n, int __max_n )
{
#ifdef  __recursion_log__
	log("resize(int ,int )");
#endif
	if( __min_n > __max_n )
	{
#ifdef  __recursion_err_msg__
		std::cerr << "Error: recursion<T>::resize(int, int)" << std::endl;
		std::cerr << "min_n: " << __min_n << std::endl;
		std::cerr << "max_n: " << __max_n << std::endl;
#endif
		return -1;
	}
	int sz = 0;
	_min_n = __min_n; _max_n = __max_n;
	if( _max_n < 0 )
		sz = 1 - _min_n;
	else
		if( _min_n > 0 )
			sz = 1 + _max_n;
		else
			sz = 1 + _max_n - _min_n;
	v_bc.resize( sz );
	sh_bc = &v_bc[0] - (_min_n > 0 ? 0 : _min_n );
	int __w = 16;
#ifdef  __recursion_log__
	std::cout << "------- recursion<T>::resize(int, int) ------" << std::endl;
	std::cout << std::setw( __w ) << "size : " << v_bc.size() << std::endl;
	std::cout << std::setw( __w ) << "min_n : " << _min_n << std::endl;
	std::cout << std::setw( __w ) << "max_n : " << _max_n << std::endl;
#endif
	if( v_bc.size() == 0 )
	{
#ifdef  __recursion_err_msg__
		std::cerr << "Error: recursion<T>::resize(int, int) -> size : " << v_bc.size() << std::endl;
#endif
		return -2;
	}
	//
	if( _max_n - _min_n != v_bc.size() - 1 )
	{
#ifdef  __recursion_warn_msg__
		std::cerr << "Warning: recursion<T>::resize(int, int)" << std::endl;
		std::cerr << std::setw(__w) << "max_n - min_n : " << _max_n-_min_n << std::endl;
		std::cerr << std::setw(__w) << "size - 1 : " << v_bc.size()-1 << std::endl;
		std::cerr << std::setw( __w ) << "min_n : " << _min_n << std::endl;
		std::cerr << std::setw( __w ) << "max_n : " << _max_n << std::endl;
#endif
		return -3;
	}
#ifdef  __recursion_log__
	std::cout << "-----------------------------------" << std::endl;
#endif
	return 1;
}
/*
 * Constants of elem_recurs structure from "elem_recurs.h"
 *
	T alp;
	T alp2;// 2 * alp;
	T sqrt_alp;// sqrt( alp );
	T ap_sqrtalp, am_sqrtalp;// a * sqrt( alp );
	T ap;// (ka + kb) / ( 2 * alp ); 
	T am;// |ka - kb| / ( 2 * alp ); 
	T sqrt_pialp;// sqrt( Pi / alp );
	T erf_ap, erf_am;//    erf( a * sqrt( alp ) );
	T erfh_ap, erfh_am;// erfh( a * sqrt( alp ) );
	T daw_ap, daw_am;// dawson( a * sqrt( alp ) );
	T exp_ap, exp_am;//    exp( -alp * a**2 );
*/


// calculate elements of Radial Integral reccursive implementation
//
// P.S.
// ka_p_kb :=  ka + kb
// ka_m_kb := |ka - kb|

template<class T>
struct recursion_lit
{
	recursion_bc<T> * ptr, *ptr_1, *ptr_2;
	plus_minus<T> * b, *c, *b1, *c1, *b2, *c2;
};

// sh_bc[0]
template<class T>
void recursion<T>::run_0_eq_3( recursion_elem<T> const & elc )
{
	recursion_bc<T> * ptr = &sh_bc[0];
	plus_minus<T> * b = &(ptr -> b), *c = &(ptr -> c);
	//_g.m = T(1)/(elc.exp_am * T(4));
	_g.p = T(1)/(elc.exp_ap * T(4));

	b = &(ptr -> b);
	c = &(ptr -> c);
	//b -> m = elc.sqrt_pialp;
	b -> p = elc.sqrt_pialp;

	//c -> m = elc.sqrt_pialp * elc.erf_am;
	c -> p = elc.sqrt_pialp * elc.erf_ap;
}

// _max_n >= 1
template<class T>
void recursion<T>::run_p1_geq_3( recursion_elem<T> const & elc )
{
	recursion_bc<T> * ptr = &sh_bc[1], *ptr_1 = ptr - 1, * ptr_2;
	plus_minus<T> * b = &(ptr -> b), *c = &(ptr -> c), *c1 = &( ptr_1 -> c ), *b1, *b2, *c2;

	b = &(ptr -> b);
	c = &(ptr -> c);

	//b -> m = elc.exp_am / elc.alp + elc.am * c1 -> m;
	//c -> m = elc.am * elc.sqrt_pialp;

	b -> p = elc.exp_ap / elc.alp + elc.ap * c1 -> p;
	c -> p = elc.ap * elc.sqrt_pialp;

	// ptr ---- &sh_bc[1]
	int im1 = 1;
	for( int i = 2; i <= _max_n; ++i )
	{
		ptr_1 = ptr;// &sh_bc[1]
		ptr_2 = ptr_1 - 1;// &sh_bc[0]

		++ptr;// &sh_bc[2]
		b = &(ptr -> b);
		c = &(ptr -> c);
		b1 = &(ptr_1 -> b);
		c1 = &(ptr_1 -> c);
		b2 = &(ptr_2 -> b);
		c2 = &(ptr_2 -> c);

		//b -> m = im1 / elc.alp2 * b2 -> m + elc.am * c1 -> m;
		//c -> m = im1 / elc.alp2 * c2 -> m + elc.am * b1 -> m;

		b -> p = im1 / elc.alp2 * b2 -> p + elc.ap * c1 -> p;
		c -> p = im1 / elc.alp2 * c2 -> p + elc.ap * b1 -> p;
		im1 = i;
	}
}

// ptr_1 - &sh_bc[-1]
template<class T>
void recursion<T>::run_m1_eq_3( recursion_elem<T> const & elc )
{
	recursion_bc<T> * ptr_1 = &sh_bc[-1];
	plus_minus<T> * b1 = &(ptr_1 -> b), *c1 = &(ptr_1 -> c);
	//
	//b1 -> m = T(_2_sqrtPi) * elc.erfh_am;
	//c1 -> m = T(_2_sqrtPi) * elc.daw_am;
	b1 -> p = T(_2_sqrtPi) * elc.erfh_ap;
	c1 -> p = T(_2_sqrtPi) * elc.daw_ap;
}

// ptr   - &sh_bc[0]
// ptr_1 - &sh_bc[-1]
// ptr_2 - &sh_bc[-2]
template<class T>
void recursion<T>::run_m2_eq_3( recursion_elem<T> const & elc )
{
	recursion_bc<T> * ptr_1 = &sh_bc[-1], * ptr = ptr_1 + 1, * ptr_2 = ptr_1 - 1;
	plus_minus<T> *b = &(ptr -> b), *c = &(ptr -> c), *b1 = &(ptr_1 -> b), *c1 = &(ptr_1 -> c), *b2 = &(ptr_2 -> b), *c2 = &(ptr_2 -> c);

	//b1 -> m = T(_2_sqrtPi) * elc.erfh_am;
	//c1 -> m = T(_2_sqrtPi) * elc.daw_am;
	b1 -> p = T(_2_sqrtPi) * elc.erfh_ap;
	c1 -> p = T(_2_sqrtPi) * elc.daw_ap;

	//b2 -> m = elc.alp2 * ( elc.am * c1 -> m - b -> m );
	//c2 -> m = elc.alp2 * ( elc.am * ( 2 * elc.exp_am + b1 -> m ) - c -> m );
	b2 -> p = elc.alp2 * ( elc.ap * c1 -> p - b -> p );
	c2 -> p = elc.alp2 * ( elc.ap * ( 2 * elc.exp_ap + b1 -> p ) - c -> p );
}

// ptr   - &sh_bc[0]
// ptr_1 - &sh_bc[-1]
// ptr_2 - &sh_bc[-2]
template<class T>
void recursion<T>::run_m3_leq_3( recursion_elem<T> const & elc )
{
	recursion_bc<T> * ptr_1 = &sh_bc[-1], *ptr = ptr_1 + 1, * ptr_2 = ptr_1 - 1;
	plus_minus<T> *b = &(ptr -> b), *c = &(ptr -> c), *b1 = &(ptr_1 -> b), *c1 = &(ptr_1 -> c), *b2 = &(ptr_2 -> b), *c2 = &(ptr_2 -> c);
	//T const & ka_p_kb = elc.ka_p_kb, & abs__ka_m_kb = elc.abs__ka_m_kb, & alp2 = elc.alp2;
	T const & ka_p_kb = elc.ka_p_kb, & alp2 = elc.alp2;

	//b1 -> m = T(_2_sqrtPi) * elc.erfh_am;
	//c1 -> m = T(_2_sqrtPi) * elc.daw_am;
	b1 -> p = T(_2_sqrtPi) * elc.erfh_ap;
	c1 -> p = T(_2_sqrtPi) * elc.daw_ap;

	//b2 -> m = alp2 * ( elc.am * c1 -> m - b -> m );
	//c2 -> m = alp2 * ( elc.am * ( 2 * elc.exp_am + b1 -> m ) - c -> m );
	b2 -> p = alp2 * ( elc.ap * c1 -> p - b -> p );
	c2 -> p = alp2 * ( elc.ap * ( 2 * elc.exp_ap + b1 -> p ) - c -> p );
	// P.S.
	// ka_p_kb = elc.alp2 * elc.ap;
	// ka_m_kb = elc.alp2 * elc.am;
	T pow_ka_p_kb_x2 = 2 * ka_p_kb;
	//T pow_ka_m_kb_x2 = 2 * abs__ka_m_kb;
	int fac_im1 = 1, im1 = 2;
	for( int i = 3; i <= -_min_n; ++i )
	{
		fac_im1 *= im1;
		pow_ka_p_kb_x2 *= ka_p_kb;
		//pow_ka_m_kb_x2 *= abs__ka_m_kb;
		ptr = &sh_bc[-i];// &sh_bc[-3]
		ptr_1 = ptr + 1;//  &sh_bc[-2]
		ptr_2 = ptr + 2;//  &sh_bc[-1]

		b = &(ptr -> b);
		c = &(ptr -> c);
		b1 = &(ptr_1 -> b);
		c1 = &(ptr_1 -> c);
		b2 = &(ptr_2 -> b);
		c2 = &(ptr_2 -> c);

		if( im1%2 == 0 )
		{
			//b -> m = ( pow_ka_m_kb_x2/ fac_im1 * elc.exp_am + abs__ka_m_kb * c1->m - alp2 * b2->m ) / im1 ;
			//c -> m = ( abs__ka_m_kb * b1->m - alp2 * c2->m ) / im1 ;
			b -> p = ( pow_ka_p_kb_x2/ fac_im1 * elc.exp_ap + ka_p_kb * c1->p - alp2 * b2->p ) / im1 ;
			c -> p = ( ka_p_kb * b1->p - alp2 * c2->p ) / im1 ;
		}
		else
		{
			//b -> m = alp2 * ( elc.am * c1->m - b2->m ) / im1 ;
			//c -> m = ( pow_ka_m_kb_x2/ fac_im1 * elc.exp_am + abs__ka_m_kb * b1 -> m - alp2 * c2 -> m ) / im1 ;
			b -> p = alp2 * ( elc.ap * c1->p - b2->p ) / im1 ;
			c -> p = ( pow_ka_p_kb_x2/ fac_im1 * elc.exp_ap + ka_p_kb * b1->p - alp2 * c2->p ) / im1 ;
		}
		im1 = i;
	}
}

template<class T>
void recursion<T>::run_3( recursion_elem<T> const & elc )
{
#ifdef  __recursion_log__
	std::cout << "-------- recursion<T>::run(recursion_elem<T> const &) [" << this << "] --------" << std::endl;
#endif
	run_0_eq_3( elc );// sh_bc[0] initialization
	if( _max_n > 0 )
		run_p1_geq_3( elc );// _max_n >= 1
	if( _min_n == -1 )
		return run_m1_eq_3( elc );
	if( _min_n == -2 )
		return run_m2_eq_3( elc );
	if( _min_n < -2 )
		return run_m3_leq_3( elc );
}

// sh_bc[0]
template<class T>
void recursion<T>::run_0_eq( recursion_elem<T> const & elc )
{
	recursion_bc<T> * ptr = &sh_bc[0];
	plus_minus<T> * b = &(ptr -> b), *c = &(ptr -> c);
	_g.m = T(1) / (elc.exp_am * T(4));
	_g.p = T(1) / (elc.exp_ap * T(4));

	b = &(ptr -> b);
	c = &(ptr -> c);
	b -> m = elc.sqrt_pialp;
	b -> p = elc.sqrt_pialp;

	c -> m = elc.sqrt_pialp * elc.erf_am;
	c -> p = elc.sqrt_pialp * elc.erf_ap;
}

// _max_n >= 1
template<class T>
void recursion<T>::run_p1_geq( recursion_elem<T> const & elc )
{
	recursion_bc<T> * ptr = &sh_bc[1], *ptr_1 = ptr - 1, * ptr_2;
	plus_minus<T> * b = &(ptr -> b), *c = &(ptr -> c), *c1 = &( ptr_1 -> c ), *b1, *b2, *c2;

	b = &(ptr -> b);
	c = &(ptr -> c);

	b -> m = elc.exp_am / elc.alp + elc.am * c1 -> m;
	c -> m = elc.am * elc.sqrt_pialp;

	b -> p = elc.exp_ap / elc.alp + elc.ap * c1 -> p;
	c -> p = elc.ap * elc.sqrt_pialp;

	// ptr ---- &sh_bc[1]
	int im1 = 1;
	for( int i = 2; i <= _max_n; ++i )
	{
		ptr_1 = ptr;// &sh_bc[1]
		ptr_2 = ptr_1 - 1;// &sh_bc[0]

		++ptr;// &sh_bc[2]
		b = &(ptr -> b);
		c = &(ptr -> c);
		b1 = &(ptr_1 -> b);
		c1 = &(ptr_1 -> c);
		b2 = &(ptr_2 -> b);
		c2 = &(ptr_2 -> c);

		b -> m = im1 / elc.alp2 * b2 -> m + elc.am * c1 -> m;
		c -> m = im1 / elc.alp2 * c2 -> m + elc.am * b1 -> m;

		b -> p = im1 / elc.alp2 * b2 -> p + elc.ap * c1 -> p;
		c -> p = im1 / elc.alp2 * c2 -> p + elc.ap * b1 -> p;
		im1 = i;
	}
}

// ptr_1 - &sh_bc[-1]
template<class T>
void recursion<T>::run_m1_eq( recursion_elem<T> const & elc )
{
	recursion_bc<T> * ptr_1 = &sh_bc[-1];
	plus_minus<T> * b1 = &(ptr_1 -> b), *c1 = &(ptr_1 -> c);
	//
	b1 -> m = T(_2_sqrtPi) * elc.erfh_am;
	c1 -> m = T(_2_sqrtPi) * elc.daw_am;
	b1 -> p = T(_2_sqrtPi) * elc.erfh_ap;
	c1 -> p = T(_2_sqrtPi) * elc.daw_ap;
}

// ptr   - &sh_bc[0]
// ptr_1 - &sh_bc[-1]
// ptr_2 - &sh_bc[-2]
template<class T>
void recursion<T>::run_m2_eq( recursion_elem<T> const & elc )
{
	recursion_bc<T> * ptr_1 = &sh_bc[-1], * ptr = ptr_1 + 1, * ptr_2 = ptr_1 - 1;
	plus_minus<T> *b = &(ptr -> b), *c = &(ptr -> c), *b1 = &(ptr_1 -> b), *c1 = &(ptr_1 -> c), *b2 = &(ptr_2 -> b), *c2 = &(ptr_2 -> c);

	b1 -> m = T(_2_sqrtPi) * elc.erfh_am;
	c1 -> m = T(_2_sqrtPi) * elc.daw_am;
	b1 -> p = T(_2_sqrtPi) * elc.erfh_ap;
	c1 -> p = T(_2_sqrtPi) * elc.daw_ap;

	b2 -> m = elc.alp2 * ( elc.am * c1 -> m - b -> m );
	c2 -> m = elc.alp2 * ( elc.am * ( 2 * elc.exp_am + b1 -> m ) - c -> m );
	b2 -> p = elc.alp2 * ( elc.ap * c1 -> p - b -> p );
	c2 -> p = elc.alp2 * ( elc.ap * ( 2 * elc.exp_ap + b1 -> p ) - c -> p );
}

// ptr   - &sh_bc[0]
// ptr_1 - &sh_bc[-1]
// ptr_2 - &sh_bc[-2]
template<class T>
void recursion<T>::run_m3_leq( recursion_elem<T> const & elc )
{
	recursion_bc<T> * ptr_1 = &sh_bc[-1], *ptr = ptr_1 + 1, * ptr_2 = ptr_1 - 1;
	plus_minus<T> *b = &(ptr -> b), *c = &(ptr -> c), *b1 = &(ptr_1 -> b), *c1 = &(ptr_1 -> c), *b2 = &(ptr_2 -> b), *c2 = &(ptr_2 -> c);
	T const & ka_p_kb = elc.ka_p_kb, & abs__ka_m_kb = elc.abs__ka_m_kb, & alp2 = elc.alp2;

	b1 -> m = T(_2_sqrtPi) * elc.erfh_am;
	c1 -> m = T(_2_sqrtPi) * elc.daw_am;
	b1 -> p = T(_2_sqrtPi) * elc.erfh_ap;
	c1 -> p = T(_2_sqrtPi) * elc.daw_ap;

	b2 -> m = alp2 * ( elc.am * c1 -> m - b -> m );
	c2 -> m = alp2 * ( elc.am * ( 2 * elc.exp_am + b1 -> m ) - c -> m );
	b2 -> p = alp2 * ( elc.ap * c1 -> p - b -> p );
	c2 -> p = alp2 * ( elc.ap * ( 2 * elc.exp_ap + b1 -> p ) - c -> p );
	// P.S.
	// ka_p_kb = elc.alp2 * elc.ap;
	// ka_m_kb = elc.alp2 * elc.am;
	T pow_ka_p_kb_x2 = 2 * ka_p_kb;
	T pow_ka_m_kb_x2 = 2 * abs__ka_m_kb;
	int fac_im1 = 1, im1 = 2;
	for( int i = 3; i <= -_min_n; ++i )
	{
		fac_im1 *= im1;
		pow_ka_p_kb_x2 *= ka_p_kb;
		pow_ka_m_kb_x2 *= abs__ka_m_kb;
		ptr = &sh_bc[-i];// &sh_bc[-3]
		ptr_1 = ptr + 1;//  &sh_bc[-2]
		ptr_2 = ptr + 2;//  &sh_bc[-1]

		b = &(ptr -> b);
		c = &(ptr -> c);
		b1 = &(ptr_1 -> b);
		c1 = &(ptr_1 -> c);
		b2 = &(ptr_2 -> b);
		c2 = &(ptr_2 -> c);

		if( im1%2 == 0 )
		{
			b -> m = ( pow_ka_m_kb_x2/ fac_im1 * elc.exp_am + abs__ka_m_kb * c1->m - alp2 * b2->m ) / im1 ;
			c -> m = ( abs__ka_m_kb * b1->m - alp2 * c2->m ) / im1 ;
			b -> p = ( pow_ka_p_kb_x2/ fac_im1 * elc.exp_ap + ka_p_kb * c1->p - alp2 * b2->p ) / im1 ;
			c -> p = ( ka_p_kb * b1->p - alp2 * c2->p ) / im1 ;
		}
		else
		{
			b -> m = alp2 * ( elc.am * c1->m - b2->m ) / im1 ;
			c -> m = ( pow_ka_m_kb_x2/ fac_im1 * elc.exp_am + abs__ka_m_kb * b1 -> m - alp2 * c2 -> m ) / im1 ;
			b -> p = alp2 * ( elc.ap * c1->p - b2->p ) / im1 ;
			c -> p = ( pow_ka_p_kb_x2/ fac_im1 * elc.exp_ap + ka_p_kb * b1->p - alp2 * c2->p ) / im1 ;
		}
		im1 = i;
	}
}

template<class T>
void recursion<T>::run( recursion_elem<T> const & elc )
{
#ifdef  __recursion_log__
	std::cout << "-------- recursion<T>::run(recursion_elem<T> const &) [" << this << "] --------" << std::endl;
#endif
	run_0_eq( elc );// sh_bc[0] initialization
	if( _max_n > 0 )
		run_p1_geq( elc );// _max_n >= 1
	if( _min_n == -1 )
		return run_m1_eq( elc );
	if( _min_n == -2 )
		return run_m2_eq( elc );
	if( _min_n < -2 )
		return run_m3_leq( elc );
}

#ifdef  __recursion_log__
template<class T>
void recursion<T>::log(std::string const & s)const
{
	std::cout << "[" << this << "] recursion<T>::" << s << std::endl;
}
#endif

#ifdef  __recursion_print__
template<class T>
void recursion<T>::print(std::ostream & out)const
{
	recursion_bc<T> const * p_i;
	plus_minus<T> const * b, * c;
	int w = 26, p = 16;
	std::string s = " recursion ";
	int t_w = 4*w + 4 - s.size();
	int w_ = t_w/2;
	out.setf( std::ios::scientific );
	out.precision( p );
	for(int i = 0; i < w_; ++i) out << '-'; out << s; for(int i = 0; i < t_w -w_; ++i) out << '-'; out << std::endl;
	out <<  std::setw(4) << "i" <<
		std::setw(w) << "B[+]" << std::setw(w) << "B[-]" <<
		std::setw(w) << "C[+]" << std::setw(w) << "C[-]" << std::endl << std::endl;
	for(int i = _min_n; i <= _max_n; ++i)
	{
		p_i = &sh_bc[i];
		b = &(p_i -> b);
		c = &(p_i -> c);
		out <<  std::setw(4) << i << 
			std::setw( w ) << b -> p << std::setw( w ) << b -> m <<
			std::setw( w ) << c -> p << std::setw( w ) << c -> m << std::endl;
	}
	out << std::endl;
	out << std::setw(4) << "min" << " : " << std::setw(4) << this->_min_n << std::endl;
	out << std::setw(4) << "max" << " : " << std::setw(4) << this->_max_n << std::endl;
	out << std::endl;
	out.unsetf( std::ios::scientific );
}

template<class T>
void print(std::ostream & out, recursion<T> const & recurs)
{
	recurs.print(out);
}
#endif//__recursion_print__

#endif//RECURSION_H
