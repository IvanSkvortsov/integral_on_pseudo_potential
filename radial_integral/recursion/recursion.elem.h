#ifndef RECURSION_ELEMENTS_H
#define RECURSION_ELEMENTS_H
#include<cmath>// sqrt
#include"spec_func/erfh.h"// erfh : erfh_d erfh_ld erfh_mp
#include"spec_func/dawson.h"// dawson : dawson_d dawson_ld dawson_mp
#include"../../../lib_math/math.constants.h"// Pi, sqrt( Pi ), 2 * sqrt( Pi ) : Pi, sqrtPi, _2_sqrtPi
#if defined __MPREAL_H__
#include"spec_func/dawson_mp.cpp"
#include"spec_func/erfh_mp.cpp"
#endif//__MPREAL_H__

#define __recursion_elem_print__
#ifdef  __recursion_elem_print__
#include<iostream>
#include<iomanip>
#endif

// recursion_elem
//
// It contains elements that are used in Radial Integral reccursivly realizied calculation

using std::sqrt;
template<class T>
struct recursion_elem
{
	T alp;// alp := alp_i + alp_j + alp_k
	T alp2;// 2 * alp;
	T alp4;// 4 * alp;
	T sqrt_alp;// sqrt( alp );
	T ap_sqrtalp, am_sqrtalp;// a * sqrt( alp );
	T ap;// (ka + kb) / ( 2 * alp ); 
	T am;// |ka - kb| / ( 2 * alp ); 
	T sqrt_pialp;// sqrt( Pi / alp );
	T erf_ap, erf_am;//    erf( a * sqrt( alp ) );
	T erfh_ap, erfh_am;// erfh( a * sqrt( alp ) );
	T daw_ap, daw_am;// dawson( a * sqrt( alp ) );
	T exp_ap, exp_am;//    exp( -alp * a**2 );
	T ka_p_kb, abs__ka_m_kb;// (ka + kb), |ka - kb|
	int sgn__ka_m_kb;// sign( ka - kb )
	// constructor
	recursion_elem();
	// destructor
	~recursion_elem(){set_zero();}
	void init_k(T const & _kb )
	{
		ka_p_kb = _kb;
		abs__ka_m_kb = _kb;
		sgn__ka_m_kb = -1;
	}
	void init_k(T const & _ka_p_kb, T const & _abs__ka_m_kb, int const & _sgn__ka_m_kb)
	{
		//init_ka_x_kb(_ka_p_kb, _abs__ka_m_kb, _sgn__ka_m_kb);
		ka_p_kb =  _ka_p_kb;
		abs__ka_m_kb = _abs__ka_m_kb;
		sgn__ka_m_kb = _sgn__ka_m_kb;
	}
	void run_k(T const & _ka, T const & _kb)
	{
		//run_ka_x_kb( _ka, _kb );
		ka_p_kb = _ka + _kb;
		abs__ka_m_kb = _ka - _kb;
		sgn__ka_m_kb = ( (abs__ka_m_kb > 0) - (abs__ka_m_kb < 0) );
		if( sgn__ka_m_kb == -1 )
			abs__ka_m_kb = -abs__ka_m_kb;
	}
	void init_a(T const & _alp)
	{
		alp = _alp;
	}
	void run_2(){this->run_3();}
	void run_3()
	{
		alp2 = 2 * alp;
		alp4 = 4 * alp;
		sqrt_alp = sqrt( alp );
		ap = ka_p_kb / alp2;
		//am = abs__ka_m_kb / alp2;

		ap_sqrtalp = ap * sqrt_alp;
		//am_sqrtalp = am * sqrt_alp;
		sqrt_pialp = T(sqrtPi) / sqrt_alp;

		erf_ap = erf( ap_sqrtalp );
		//erf_am = erf( am_sqrtalp );

		erfh_ap = erfh( ap_sqrtalp );
		//erfh_am = erfh( am_sqrtalp );

		daw_ap = dawson( ap_sqrtalp );
		//daw_am = dawson( am_sqrtalp );

		exp_ap = exp(-ka_p_kb * ka_p_kb / alp4 );
		//exp_am = exp( -abs__ka_m_kb * abs__ka_m_kb / alp4 );
	}
	void run_1()
	{
		alp2 = 2 * alp;
		alp4 = 4 * alp;
		sqrt_alp = sqrt( alp );
		ap = ka_p_kb / alp2;
		am = abs__ka_m_kb / alp2;

		ap_sqrtalp = ap * sqrt_alp;
		am_sqrtalp = am * sqrt_alp;
		sqrt_pialp = T(sqrtPi) / sqrt_alp;

		erf_ap = erf( ap_sqrtalp );
		erf_am = erf( am_sqrtalp );

		erfh_ap = erfh( ap_sqrtalp );
		erfh_am = erfh( am_sqrtalp );

		daw_ap = dawson( ap_sqrtalp );
		daw_am = dawson( am_sqrtalp );

		exp_ap = exp(-ka_p_kb * ka_p_kb / alp4 );
		exp_am = exp( -abs__ka_m_kb * abs__ka_m_kb / alp4 );
	}
	void run()
	{
		this->run_1();
	}
	void reverse()
	{
		sgn__ka_m_kb = -sgn__ka_m_kb;
	}
	void init( T const & _alp );
	void init_ka_x_kb(T const & _ka_p_kb, T const & _abs__ka_m_kb, int const & _sgn__ka_m_kb);
	void run_ka_x_kb(T const & _ka, T const & _kb);
	void set_zero();
#ifdef  __recursion_elem_print__
	void print(std::ostream & out)const;
#endif
};

template<class T>
recursion_elem<T>::recursion_elem():alp(0), alp2(0), alp4(0), sqrt_alp(0), ap_sqrtalp(0), am_sqrtalp(0), ap(0), am(0), sqrt_pialp(0),
	erf_ap(0), erf_am(0), erfh_ap(0), erfh_am(0), daw_ap(0), daw_am(0), exp_ap(0), exp_am(0), ka_p_kb(0), abs__ka_m_kb(0), sgn__ka_m_kb(0){}

template<class T>
void recursion_elem<T>::init_ka_x_kb( T const & _ka_p_kb, T const & _abs__ka_m_kb, int const & _sgn__ka_m_kb)
{
	// ka_p_kb      := ka + kb;
	// abs__ka_m_kb := abs( ka - kb );
	// sgn__ka_m_kb := sgn( ka - kb );
	ka_p_kb =  _ka_p_kb;
	abs__ka_m_kb = _abs__ka_m_kb;
	sgn__ka_m_kb = _sgn__ka_m_kb;
}

template<class T>
void recursion_elem<T>::run_ka_x_kb( T const & _ka, T const & _kb)
{
	// ka_p_kb      := ka + kb;
	// abs__ka_m_kb := abs( ka - kb );
	// sgn__ka_m_kb := sgn( ka - kb );
	ka_p_kb =  _ka + _kb;
	abs__ka_m_kb = _ka - _kb;
	sgn__ka_m_kb = sgn( abs__ka_m_kb );
	if( sgn__ka_m_kb == -1 )
		abs__ka_m_kb = -abs__ka_m_kb;
}

template<class T>
void recursion_elem<T>::init( T const & _alp )
{
	alp = _alp;
	alp2 = 2 * alp;
	alp4 = 4 * alp;
	sqrt_alp = sqrt( alp );
	ap = ka_p_kb / alp2;
	am = abs__ka_m_kb / alp2;

	ap_sqrtalp = ap * sqrt_alp;
	am_sqrtalp = am * sqrt_alp;
	sqrt_pialp = T(sqrtPi) / sqrt_alp;

	erf_ap = erf( ap_sqrtalp );
	erf_am = erf( am_sqrtalp );

	erfh_ap = erfh( ap_sqrtalp );
	erfh_am = erfh( am_sqrtalp );

	daw_ap = dawson( ap_sqrtalp );
	daw_am = dawson( am_sqrtalp );

	exp_ap = exp(-ka_p_kb * ka_p_kb / alp4 );
	exp_am = exp( -abs__ka_m_kb * abs__ka_m_kb / alp4 );
}
template<class T>
void recursion_elem<T>::set_zero()
{
	alp = 0;
	alp2 = 0;
	alp4 = 0;
	sqrt_alp = 0;
	ap_sqrtalp = 0;
	am_sqrtalp = 0;
	ap = 0;
	am = 0;
	sqrt_pialp = 0;
	erf_ap = 0;
	erf_am = 0;
	erfh_ap = 0;
	erfh_am = 0;
	daw_ap = 0;
	daw_am = 0;
	exp_ap = 0;
	exp_am = 0;

	ka_p_kb = 0;
	abs__ka_m_kb = 0;
	sgn__ka_m_kb = 0;
}

#ifdef  __recursion_elem_print__
template<class T>
void recursion_elem<T>::print(std::ostream & out)const
{
	int p = 30, w = p + 10;
	out.setf(std::ios::scientific);
	out << std::setw(10) << "alp" << std::setw( w ) << std::setprecision( p ) << alp << std::endl;
	out << std::setw(10) << "alp2" << std::setw( w ) << std::setprecision( p ) << alp2 << std::endl;
	out << std::setw(10) << "alp4" << std::setw( w ) << std::setprecision( p ) << alp4 << std::endl;
	out << std::setw(10) << "sqrt_alp" << std::setw( w ) << std::setprecision( p ) << sqrt_alp << std::endl;
	out << std::setw(10) << "ap_sqrtalp" << std::setw( w ) << std::setprecision( p ) << ap_sqrtalp << std::endl;
	out << std::setw(10) << "am_sqrtalp" << std::setw( w ) << std::setprecision( p ) << am_sqrtalp << std::endl;
	out << std::setw(10) << "ap" << std::setw( w ) << std::setprecision( p ) << ap << std::endl;
	out << std::setw(10) << "am" << std::setw( w ) << std::setprecision( p ) << am << std::endl;
	out << std::setw(10) << "sqrt_pialp" << std::setw( w ) << std::setprecision( p ) << sqrt_pialp << std::endl;
	out << std::setw(10) << "erf_ap" << std::setw( w ) << std::setprecision( p ) << erf_ap << std::endl;
	out << std::setw(10) << "erf_am" << std::setw( w ) << std::setprecision( p ) << erf_am << std::endl;
	out << std::setw(10) << "erfh_ap" << std::setw( w ) << std::setprecision( p ) << erfh_ap << std::endl;
	out << std::setw(10) << "erfh_am" << std::setw( w ) << std::setprecision( p ) << erfh_am << std::endl;
	out << std::setw(10) << "daw_ap" << std::setw( w ) << std::setprecision( p ) << daw_ap << std::endl;
	out << std::setw(10) << "daw_am" << std::setw( w ) << std::setprecision( p ) << daw_am << std::endl;
	out << std::setw(10) << "exp_ap" << std::setw( w ) << std::setprecision( p ) << exp_ap << std::endl;
	out << std::setw(10) << "exp_am" << std::setw( w ) << std::setprecision( p ) << exp_am << std::endl;
	out << std::setw(10) << "ka+kb" << std::setw( w ) << std::setprecision( p ) << ka_p_kb << std::endl;
	out << std::setw(10) << "|ka-kb|" << std::setw( w ) << std::setprecision( p ) << abs__ka_m_kb << std::endl;
	out << std::setw(10) << "sgn(ka-kb)" << std::setw( w ) << std::setprecision( p ) << sgn__ka_m_kb << std::endl;
	out.unsetf( std::ios::scientific );
	//print<T>(out, *this);
	return;
}

template<class T>
void print(std::ostream & out, recursion_elem<T> const & elc)
{
	/*
	T alp;// alp := alp_i + alp_j + alp_k
	T alp2;// 2 * alp;
	T sqrt_alp;// std::sqrt( alp );
	T ap_sqrtalp, am_sqrtalp;// a * sqrt( alp );
	T ap;// (ka + kb) / ( 2 * alp ); 
	T am;// |ka - kb| / ( 2 * alp ); 
	T sqrt_pialp;// sqrt( Pi / alp );
	T erf_ap, erf_am;//    erf( a * sqrt( alp ) );
	T erfh_ap, erfh_am;// erfh( a * sqrt( alp ) );
	T daw_ap, daw_am;// dawson( a * sqrt( alp ) );
	T exp_ap, exp_am;//    exp( -alp * a**2 );
	*/

	int w = 16, p = 6;
	out.setf(std::ios::scientific);
	out << std::setw(10) << "alp" << std::setw( w ) << std::setprecision( p ) << elc.alp << std::endl;
	out << std::setw(10) << "alp2" << std::setw( w ) << std::setprecision( p ) << elc.alp2 << std::endl;
	out << std::setw(10) << "alp4" << std::setw( w ) << std::setprecision( p ) << elc.alp4 << std::endl;
	out << std::setw(10) << "sqrt_alp" << std::setw( w ) << std::setprecision( p ) << elc.sqrt_alp << std::endl;
	out << std::setw(10) << "ap_sqrtalp" << std::setw( w ) << std::setprecision( p ) << elc.ap_sqrtalp << std::endl;
	out << std::setw(10) << "am_sqrtalp" << std::setw( w ) << std::setprecision( p ) << elc.am_sqrtalp << std::endl;
	out << std::setw(10) << "ap" << std::setw( w ) << std::setprecision( p ) << elc.ap << std::endl;
	out << std::setw(10) << "am" << std::setw( w ) << std::setprecision( p ) << elc.am << std::endl;
	out << std::setw(10) << "sqrt_pialp" << std::setw( w ) << std::setprecision( p ) << elc.sqrt_pialp << std::endl;
	out << std::setw(10) << "erf_ap" << std::setw( w ) << std::setprecision( p ) << elc.erf_ap << std::endl;
	out << std::setw(10) << "erf_am" << std::setw( w ) << std::setprecision( p ) << elc.erf_am << std::endl;
	out << std::setw(10) << "erfh_ap" << std::setw( w ) << std::setprecision( p ) << elc.erfh_ap << std::endl;
	out << std::setw(10) << "erfh_am" << std::setw( w ) << std::setprecision( p ) << elc.erfh_am << std::endl;
	out << std::setw(10) << "daw_ap" << std::setw( w ) << std::setprecision( p ) << elc.daw_ap << std::endl;
	out << std::setw(10) << "daw_am" << std::setw( w ) << std::setprecision( p ) << elc.daw_am << std::endl;
	out << std::setw(10) << "exp_ap" << std::setw( w ) << std::setprecision( p ) << elc.exp_ap << std::endl;
	out << std::setw(10) << "exp_am" << std::setw( w ) << std::setprecision( p ) << elc.exp_am << std::endl;
	out.unsetf( std::ios::scientific );
}
#endif//__recursion_elem_print__

#endif//RECURSION_ELEMENTS_H
