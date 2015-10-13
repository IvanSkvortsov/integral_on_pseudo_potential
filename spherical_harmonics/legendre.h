#ifndef __LEGENDRE_H__
#define __LEGENDRE_H__
#include"polynomial.h"
#include"math.h"// pow
//#include"../../lib_math/math.functions.h"// NewtonC, par_fact
//#include"../../lib_math/math.constants.h"// Pi, sqrtPi, etc.
#include"../lib_math/math.functions.h"// NewtonC, par_fact
#include"../lib_math/math.constants.h"// Pi, sqrtPi, etc.

template<class T>
struct polynom_4
{
	typedef T value_type;
	T d;
	int ct, st, cp, sp;
	polynom_4():d(1),ct(0),st(0),cp(0),sp(0){}
	polynom_4(T const & __d, int __ct, int __st, int __cp, int __sp):d(__d), ct(__ct), st(__st), cp(__cp), sp(__sp){}
	polynom_4(T const & __d, int const * __x):d(__d), ct(*__x), st(*(__x+1)), cp(*(__x+2)), sp(*(__x+3)){}
	polynom_4<T> & operator*=(polynom_4<T> const & v)
	{
		d *= v.d;
		ct += v.ct;
		st += v.st;
		cp += v.cp;
		sp += v.sp;
		return *this;
	}
	polynom_4<T> & operator*=(T const & v)
	{
		d *= v;
		return *this;
	}
	int & operator[](int i){return (&(int &)ct)[i];}
	int const & operator[](int i)const{return (&(int const &)ct)[i];}
	void set_zero_x()
	{
		ct = 0; st = 0; cp = 0; sp = 0;
	}
	bool is_even_x()const
	{
		return ct%2==0 && st%2==0 && cp%2==0 && sp%2==0;
	}
	static bool is_equal_x(polynom_4<T> const & a, polynom_4<T> const & b)
	{
		return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
	}
	T calc_r(T const * X)const
	{
		return d * pow_int<T>(X[0], ct) * pow_int<T>(X[1], st) * pow_int<T>(X[2], cp) * pow_int<T>(X[3], sp);
	}
};

template<class T>
polynom_4<T> operator*(polynom_4<T> const & a, polynom_4<T> const & b)
{
	polynom_4<T> c( a );
	c *= b;
	return c;
}

template<class T>
bool operator==(polynom_4<T> const & a, polynom_4<T> const & b)
{
	return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
}

template<class T>
bool operator!=(polynom_4<T> const & a, polynom_4<T> const & b)
{
	return !(a==b);
}

#define __polynom_4_PRINT__
#ifdef  __polynom_4_PRINT__
#include<iostream>
#include<iomanip>
template<class T>
void print_polynom_4(std::ostream & out, polynom_4<T> const & pol, int const & p = 4, int const & w_x = 6)
{
	int w = p + 8;
	out.setf( std::ios::scientific );
	out.precision( p );
	out <<  std::setw( w ) << pol.d <<
		std::setw( w_x ) << pol.ct <<
		std::setw( w_x ) << pol.st <<
		std::setw( w_x ) << pol.cp <<
		std::setw( w_x ) << pol.sp << std::endl;
}
#endif//__polynom_4_PRINT__

#define __legendre_error_message_ON_
#ifdef  __legendre_error_message_ON_
#include<iostream>
#endif

template<class T, class pT = polynom_4<T> >
struct legendre: public polynomial<pT>
{
	typedef pT polynomial_type;
	int _n, _m;
	legendre():polynomial<pT>(), _n(0), _m(0){}
	legendre(legendre<T, pT> const & legen):polynomial<pT>( legen ), _n(legen._n), _m(legen._m){}
	legendre(int const & __size):polynomial<pT>( __size ), _n(0), _m(0){}
	legendre<T, pT> & operator=(legendre<T, pT> const & legen)
	{
		if( this == &legen ) return *this;
		this->polynomial<pT>::operator=( legen );
		_n = legen.n();
		_m = legen.m();
		return *this;
	}
	int const & n()const{return _n;}
	int const & m()const{return _m;}
	void run( int __n)
	{
		this->_n = __n;
		unsigned int _2pow_n = 1 << (this->_n<31?this->_n:30);
		T t_i2pow_n = T(_2pow_n);
		for(int i = 30; i < this->_n; ++i)
			t_i2pow_n *= 2;
		t_i2pow_n  = 1/t_i2pow_n;
		this->resize( this->_n/2 + 1 );
		pT * p_legen = this->_data;
		for(int k = 0; k < this->size(); ++k)
		{
			p_legen->d = NewtonC<T>(this->_n, k) * NewtonC<T>(2*this->_n - 2*k, this->_n) * (k%2?-1:1) * t_i2pow_n;
			p_legen->set_zero_x();
			p_legen->ct = this->_n - 2 * k;
			//p_legen->operator[](0) = _n - 2 * k;
			p_legen++;
		}
	}
	int run( int __n, int __m )
	{
		this->run( __n );
		return this->run_m( __m );
	}
	int run_m( int __m)
	{
		int abs_m = (__m < 0? -__m : __m);
		if( abs_m > this->_n )
		{
#ifdef  __legendre_error_message_ON_
			std::cerr << "Error: legendre<T, pT>::run_m(int __m)" << std::endl;
			std::cerr << "__m : " << __m << std::endl;
			std::cerr << "_n  : " << _n << std::endl;
#endif
			return 1;
		}
		_m = __m;
		if( !_m ) return 0;
		int ct_pow = 0, iter = 0;//, * zero_k = new int[this->_size];
		pT * p_legen = this->_data;
		for(int k = 0; k < this->_size; ++k)
		{
			ct_pow = p_legen->operator[](0);
			if( ct_pow < abs_m )
			{
				//zero_k[iter] = k;
				p_legen->d  = T(0);
				p_legen->set_zero_x();
				++iter;
			}
			else
			{
				p_legen->d *= par_fact<T>(ct_pow-abs_m+1, abs_m);//fact<T>(ct_pow)/fact<T>(ct_pow-m);
				//p_legen->operator[](0) = ct_pow - abs_m;
				//p_legen->operator[](1) = abs_m;
				p_legen->ct = ct_pow - abs_m;
				p_legen->st = abs_m;
			}
			++p_legen;
		}
		//if( !iter )
		//{
			//delete [] zero_k;
		//	return 0;
		//}
		// little optimization
		// fast variant
		//this->resize( this->_size - iter );
		this->_size -= iter;
		return 0;
		/*
		legendre<T, pT> tmp( this->_size - iter );
		p_legen = &tmp[0];
		int * p_zero = zero_k;
		for(int k = 0; k < this->_size; ++k)
		{
			if( *p_zero == k )
			{
				p_zero++;
				continue;
			}
			*p_legen++ = this->_data[k];
		}
		tmp._n = this->_n;
		tmp._m = this->_m;
		this->operator=(tmp);
		*/
		//delete [] zero_k;
		if( abs_m%2==0 )
			return 0;
		pT * p = this->data();
		for(int i = 0; i < this->size(); ++i)
		{
			p->d = -p->d;
			++p;
		}
		return 0;
	}
	void run_cos_sin( int __m)
	{
		_m = __m;
		if( !_m )
		{
			this->resize( 1 );
			this->_data->d = T(1);//T(isqrt_2);
			this->_data->set_zero_x();
			return;
		}
		int _is_m_neg = (_m < 0 ? 1 : 0), abs_m = ( _is_m_neg ? -_m : _m ), __new_size = 0;
		for(int i = _is_m_neg; i <= abs_m; i += 2) __new_size++;
		this->resize( (_is_m_neg ? (abs_m+1)/2 : (abs_m/2 + 1) ) );
		if( this->_size != __new_size )
		{
#ifdef  __legendre_error_message_ON_
			std::cerr << "Error: legendre<T, pT>::run_cos_sin(int __m)" << std::endl;
			std::cerr << "__m : " << __m << std::endl;
			std::cerr << "this->_size : " << this->_size << std::endl;
			std::cerr << " __new_size : " << __new_size << std::endl;
#endif
			exit(1);
		}
		polynom_4<T> * p_cs = this->_data;
		int __k = 0;
		for(int i = _is_m_neg; i <= abs_m; i += 2)
		{
			//p_cs->d = ( __k%2==0 ? NewtonC<T>(abs_m, i) : -NewtonC<T>(abs_m, i) );
			//p_cs->d = ( (i-_is_m_neg)%4 ? -NewtonC<T>(abs_m, i) : NewtonC<T>(abs_m, i) );
			p_cs->d = ( ((i-_is_m_neg)/2)%2 ? -NewtonC<T>(abs_m, i) : NewtonC<T>(abs_m, i) );
			p_cs->set_zero_x();
			//p_cs->operator[](2) = abs_m - i;
			//p_cs->operator[](3) = i;
			p_cs->cp = abs_m - i;
			p_cs->sp = i;
			++p_cs;
			++__k;
		}
		if( __k != this->_size )
		{
#ifdef  __legendre_error_message_ON_
			std::cerr << "Error: legendre<T, pT>::run_cos_sin(int __m)" << std::endl;
			std::cerr << "__m : " << __m << std::endl;
			std::cerr << "this->_size : " << this->_size << std::endl;
			std::cerr << "__k : " << __k << std::endl;
#endif
			exit(1);
		}
	}
};

#ifdef  __polynom_4_PRINT__

#define __legendre_PRINT__
#ifdef  __legendre_PRINT__
#include<iostream>
#include<iomanip>

#ifndef __smart_title_PRINT__
#define __smart_title_PRINT__
void title_print(std::ostream & out, int const & line_size, std::string const & title_str)
{
	for(int i = 0; i < line_size; ++i) out << '-';
	out << std::endl;
	int rem = line_size - title_str.size(), side = 0;
	if( rem < 0 ) exit(1);
	side = rem / 2 - 1;
	if( side < 0 ) exit(1);
	for(int i = 0; i < side; ++i) out << '-';
	out << ' ' << title_str << ' ';
	for(int i = 0; i < rem - side - 2; ++i) out << '-';
	out << std::endl;
}
#endif//__smart_title_PRINT__

template<class T>
void print_legendre(std::ostream & out, legendre<T> const & legen, int const & n, int const & p = 4)
{
	int w = p + 8, w_x = 6;
	out.setf( std::ios::scientific );
	out.precision( p );
	title_print(out, 8 + w + 4 * w_x, "legendre polynomial" );
	//out << "----------------------------" << std::endl;
	//out << "--- lengendre polynomial ---" << std::endl;
	out <<  std::setw( 4 ) << "n" <<
	        std::setw( 4 ) << "m" <<
		std::setw( w ) << "d" <<
		std::setw( w_x ) << "cos_t" <<
		std::setw( w_x ) << "sin_t" <<
		std::setw( w_x ) << "cos_p" <<
		std::setw( w_x ) << "sin_p" << std::endl << std::endl;
	out << std::setw( 4 ) << n << std::setw( 4 ) << legen.m();
	for(int i = 0; i < legen.size(); ++i)
	{
		if( i ) out << std::setw( 8 ) << "";
		print_polynom_4<T>(out, legen[i], p, w_x);
	}
	out << std::endl;
	out.unsetf( std::ios::scientific );
}
#endif//__legendre_PRINT__
#endif//__polynom_4_PRINT__


#endif//__LEGENDRE_H__
