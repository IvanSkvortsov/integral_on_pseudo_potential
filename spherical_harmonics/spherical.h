#ifndef __SPHERICAL_H__
#define __SPHERICAL_H__
#include"polynomial.h"
#include"legendre.h"
#include"../angular_integral/angular.omega.xyz.h"
#include"math.h"// sqrt

template<class T>
struct polynom_xyz;

template<class T, class pT = polynom_xyz<T> >
struct spherical;

template<class T>
T Sigma_sph(angular_omega_xyz<T> const & sigma_ijk_, spherical<T> const & spher)
{
	if( !spher.is_even() )
		return T(0);
	T value = T(0);
	typename spherical<T>::polynomial_type const * p_pol = &spher[0];
	for(int i = 0; i < spher.size(); ++i)
	{
		value += sigma_ijk_(p_pol->x, p_pol->y, p_pol->z) * p_pol->d;
		++p_pol;
	}
	return value;
}

template<class T>
struct polynom_xyz
{
	typedef T value_type;
	T d;
	int x, y, z;
	polynom_xyz():d(1), x(0), y(0), z(0){}
	polynom_xyz(T const & __d, int __x, int __y, int __z):d(__d), x(__x), y(__y), z(__z){}
	polynom_xyz(T const & __d, int const * __x):d(__d), x(*__x), y(*(__x+1)), z(*(__x+2)){}
	int & operator[](int i){return ((int *)&x)[i];}
	int const & operator[](int i)const{return (&(int const &)x)[i];}
	int sum_x()const{return (x + y + z);}
	void set_x(int const * __x)
	{
		x = __x[0]; y = __x[1]; z = __x[2];
	}
	void set_x(int const & __x, int const & __y, int const & __z)
	{
		x = __x; y = __y; z = __z;
	}
	void set_zero_x()
	{
		x = 0; y = 0; z = 0;
	}
	polynom_xyz<T> & operator*=(polynom_xyz<T> const & v)
	{
		d *= v.d;
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	polynom_xyz<T> & operator*=(T const & v)
	{
		d *= v;
		return *this;
	}
	bool is_even_x()const
	{
		return x%2==0 && y%2==0 && z%2==0;
	}
	static bool is_equal_x(polynom_xyz<T> const & a, polynom_xyz<T> const & b)
	{
		return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
	}
	T calc_r(T const * X)const
	{
		return this->d * pow( X[0], this->x ) * pow( X[1], this->y ) * pow( X[2], this->z );
	}
};

template<class T>
polynom_xyz<T> operator*(polynom_xyz<T> const & a, polynom_xyz<T> const & b)
{
	polynom_xyz<T> c( a );
	c *= b;
	return c;
}

#define __polynom_xyz_PRINT__
#ifdef  __polynom_xyz_PRINT__
#include<iostream>
#include<iomanip>
#include<fstream>
template<class T>
void print_polynom_xyz(std::ostream & out, polynom_xyz<T> const & pol, int const & p = 4, int const & w_x = 6)
{
	int w = p + 8;
	//angular_omega_xyz<T> sigma_;
	//sigma_.run(15, 15, 15);
	out.setf( std::ios::scientific );
	out.precision( p );
	out <<  std::setw( w ) << pol.d <<
		std::setw( w_x ) << pol.x <<
		std::setw( w_x ) << pol.y <<
		std::setw( w_x ) << pol.z <<
		//std::setw( w ) << pol.d * sigma_(pol.x, pol.y, pol.z) <<
		//std::setw( w ) << pol.d * sigma_(pol.x+1, pol.y, pol.z) <<
		//std::setw( w ) << pol.d * sigma_(pol.x, pol.y+1, pol.z) <<
		//std::setw( w ) << pol.d * sigma_(pol.x, pol.y, pol.z+1) <<
		std::setw(w_x) << ( pol.d == T(0) ? '0' : '+') <<
		std::endl;
}
#endif//__polynom_xyz_PRINT__

#define __spherical_error_message_ON__
#ifdef  __spherical_error_message_ON__
#include<iostream>
#endif


int __3pow_n(int const & n)
{
	if( n > 19 )
	{
#ifdef  __spherical_error_message_ON__
		std::cerr << "Error: [__3pow_n(int const & n)]" << std::endl;
		std::cerr << "n : " << n << std::endl;
#endif
		exit(1);
	}
	int res = 1;
	for(int i = 0; i < n; ++i)
		res *= 3;
	return res;
}
//#define __spherical_log__
#ifdef  __spherical_log__
#include<iostream>
#endif
//template<class T, class pT = polynom_xyz<T> >
// TODO: in method 'optimize_less' decrease memory allocation/deallocation
// TODO: int _n, _m -- make some class that's parent for 'spherical' and 'legendre' (reason is that both classes have _n, _m member-elements)
template<class T, class pT >
struct spherical : public polynomial<pT>
{
	typedef pT polynomial_type;
	int _n, _m;
	spherical():polynomial<pT>(), _n(0), _m(0){}
	spherical(spherical<T, pT> const & spher):polynomial<pT>( spher ), _n(spher._n), _m(spher._m){}
	spherical(int const & __size):polynomial<pT>( __size ), _n(0), _m(0){}
	spherical<T, pT> & operator=(spherical<T, pT> const & spher)
	{
		if( this == &spher ) return *this;
		this->polynomial<pT>::operator=( spher );
		_n = spher.n();
		_m = spher.m();
		return *this;
	}
	int n(int const & __n){return _n = __n;}
	int m(int const & __m){return _m = __m;}
	int const & n()const{return _n;}
	int const & m()const{return _m;}
	int run_r2()
	{
		if( this->_size != 3 )
			this->resize( 3 );
		pT * p = this->_data;
		for(int i = 0; i < 3; ++i)
		{
			p->set_zero_x();
			p->d = T(1);
			p->operator[](i) = 2;
			++p;
		}
		return 0;
	}
	int run(int const & __n, int const & __m, int to_norma = 1)
	{
#ifdef  __spherical_log__
		this->log(std::cout, "run(int const & , int const & , int )");
#endif
		legendre<T> P_nm, cos_sin;
		if( P_nm.run( __n, __m ) )
		{
#ifdef  __spherical_error_message_ON__
			std::cerr << "Error: from spherical<T, pT>::run(int const & __n, int const & __m)" << std::endl;
			std::cerr << "__n : " << __n << std::endl;
			std::cerr << "__m : " << __m << std::endl;
#endif
			std::exit(1);
		}
		cos_sin.run_cos_sin( __m );
		P_nm *= cos_sin;
		this->run(P_nm);
		if( to_norma == 1 ) this->normalize();
		else if( to_norma >= 2 ) this->normalize_4pi();
		return 0;
	}
	int run(legendre<T> const & cs__P_nm)
	{
		this->resize( cs__P_nm.size() );
		_n = cs__P_nm.n();
		_m = cs__P_nm.m();
		typename spherical<T>::polynomial_type * p_xyz = this->_data;
		typename legendre<T>::polynomial_type const * p_4 = &cs__P_nm[0];
		for(int i = 0; i < this->_size; ++i)
		{
			p_xyz->d = p_4->d;
			p_xyz->operator[](0) = p_4->cp;// [2]
			p_xyz->operator[](1) = p_4->sp;// [3]
			p_xyz->operator[](2) = p_4->ct;// [0]
			//p_xyz->x = p_4->cp;
			//p_xyz->y = p_4->sp;
			//p_xyz->z = p_4->ct;
			++p_xyz;
			++p_4;
		}
		return 0;
	}
	T norm()const
	{
		if( this->_n == 0 ) return T(1);
		T norma_ = T(this->_n * 2 + 1);
		if( this->_m == 0 )
		{
			norma_ /= (4 * T(Pi));
			norma_ = sqrt( norma_ );
			return norma_;
		}
		int abs_m = (this->_m < 0 ? -this->_m : this->_m);
		norma_ *= fact<T>(this->_n - abs_m) / (fact<T>( this->_n + abs_m) * 2 * T(Pi));
		//norma_ = sqrt(norma_);
		return norma_;
	}
	T norm_4pi()const// norm multiplied by sqrt(4 * Pi)
	{
		if( this->_n == 0 ) return T(1);
		T norma_ = T(this->_n * 2 + 1);
		if( this->_m == 0 )
		{
			norma_ = sqrt( norma_ );
			return norma_;
		}
		int abs_m = (this->_m < 0 ? -this->_m : this->_m);
		norma_ *= fact<T>(this->_n - abs_m) * T(2) / fact<T>( this->_n + abs_m);
		//norma_ = sqrt(norma_);
		return norma_;
	}
	T normalize()
	{
		if( this->_n == 0 ) return T(1)/sqrt(4 * T(Pi));
		T norma_ = T(this->_n * 2 + 1);
		if( this->_m == 0 )
		{
			norma_ /= (4 * T(Pi));
			norma_ = sqrt( norma_ );
			*this *= norma_;
			//this->polynomial<pT>::operator*=(norma_);
			return norma_;
		}
		int abs_m = (this->_m < 0 ? -this->_m : this->_m);
		norma_ *= fact<T>(this->_n - abs_m) / (fact<T>( this->_n + abs_m) * 2 * T(Pi));
		norma_ = sqrt(norma_);
		*this *= norma_;
		//this->polynomial<pT>::operator*=(norma);
		return norma_;
	}
	T normalize_4pi()// norm multiplied by sqrt(4 * Pi)
	{
		if( this->_n == 0 ) return T(1);
		T norma_ = T(this->_n * 2 + 1);
		if( this->_m == 0 )
		{
			norma_ = sqrt( norma_ );
			*this *= norma_;
			//this->polynomial<pT>::operator*=(norma_);
			return norma_;
		}
		int abs_m = (this->_m < 0 ? -this->_m : this->_m);
		norma_ *= fact<T>(this->_n - abs_m) * T(2) / fact<T>( this->_n + abs_m);
		norma_ = sqrt(norma_);
		*this *= norma_;
		//this->polynomial<pT>::operator*=(norma);
		return norma_;
	}
	void optimize_less()
	{
		int sz = this->_size, * lack_n = new int[sz], * lack_pow = new int[sz], * lack_sz = new int[sz], new_sz = 0, tmp_si = 0;
		polynom_xyz<T> * p_xyz = this->_data, * p_b = 0, * tmp_p = 0;
		for(int i = 0; i < sz; ++i)
		{
			if( p_xyz->sum_x() > this->_n )
			{
#ifdef  __spherical_error_message_ON__
				std::cerr << "Error: from spherical<T, pT>::optimize_less()" << std::endl;
				std::cerr << "n : " << _n << std::endl;
				std::cerr << "m : " << _m << std::endl;
				std::cerr << "sum_x : " << p_xyz->sum_x() << std::endl;
				std::cerr << "x : " << p_xyz->x << std::endl;
				std::cerr << "y : " << p_xyz->y << std::endl;
				std::cerr << "z : " << p_xyz->z << std::endl;
#endif
				exit(1);
			}
			lack_n[i] = this->_n - p_xyz++->sum_x();
		}
		for(int i = 0; i < sz; ++i)
		{
			tmp_si = lack_n[i];
			if( tmp_si%2 )
			{
#ifdef  __spherical_error_message_ON__
				std::cerr << "Error: from spherical<T, pT>::optimize_less()" << std::endl;
				std::cerr << "n : " << _n << std::endl;
				std::cerr << "m : " << _m << std::endl;
				std::cerr << "tmp_si : " << tmp_si << std::endl;
#endif
				exit(1);
			}
			tmp_si /= 2;
			lack_pow[i] = tmp_si;
			new_sz += ( lack_sz[i] = __3pow_n( tmp_si ) );
		}
		if( new_sz == sz )
		{
			delete [] lack_n;
			delete [] lack_pow;
			delete [] lack_sz;
			return;
		}
		spherical<T, pT> spher_new( new_sz ), spher_buf, spher_r2(3);
		spher_r2.run_r2();
		p_b = &spher_new[0];
		int _lack_sz = 0;
		p_xyz = this->_data;
		for(int i = 0; i < sz; ++i)
		{
			*p_b = *p_xyz;
			_lack_sz = lack_sz[i];
			spher_buf.set( 1, _lack_sz, p_b);// dangerous magic: explicit set of member elements (do not repeat!)
			for(int j = 0; j < lack_pow[i]; ++j)
				spher_buf *= spher_r2;
			p_b += _lack_sz;
			spher_buf.zero();                // dangerous magic: zeroize all member elements (do not repeat!)
			++p_xyz;
		}
		this->polynomial<pT>::operator=(spher_new);
		delete [] lack_n;
		delete [] lack_pow;
		delete [] lack_sz;
	}
	void optimize_ez()
	{
		this->optimize_equal();
		this->optimize_zero();
	}
	void optimize()
	{
		this->optimize_equal();
		this->optimize_zero();
		this->optimize_less();
		this->optimize_equal();
		this->optimize_zero();
	}
	static void generate(scp_array<spherical<T, pT> > & v_spher, int __n )
	{
		legendre<T> P_n, cs__P_nm, cos_sin;
		P_n.run( __n );
		v_spher.resize(2 * __n + 1);
		spherical<T, pT> * p_spher = &v_spher[0];
		for(int __m = -__n; __m <= __n; ++__m)
		{
			cs__P_nm = P_n;
			cs__P_nm.run_m( __m );
			cos_sin.run_cos_sin( __m );
			cs__P_nm *= cos_sin;
			p_spher->n( __n );
			p_spher->m( __m );
			p_spher->run( cs__P_nm );
			p_spher++->optimize_ez();
		}
	}
private:
#ifdef  __spherical_log__
	void log(std::ostream & out, std::string s)const
	{
		out << "[" << this << "] spherical<T>::" << s << std::endl;
	}
#endif
};


#ifdef  __polynom_xyz_PRINT__
#define __spherical_PRINT__
#ifdef  __spherical_PRINT__
#include<iostream>
#include<iomanip>
#include<fstream>

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
void print_spherical(std::ostream & out, spherical<T> const & spher, int const & n, int const & m)
{
	if( n != spher.n() || m != spher.m() )
	{
#ifdef  __spherical_error_message_ON__
		std::cerr << "Error : [print_spherical]" << std::endl;
		std::cerr << "n : " << n << std::endl;
		std::cerr << "m : " << m << std::endl;
		std::cerr << "s.n() : " << spher.n() << std::endl;
		std::cerr << "s.m() : " << spher.m() << std::endl;
#endif
		//return;
	}
	int p = 4, w = p + 8, w_x = 6, sz = 0;
	out.setf( std::ios::scientific );
	out.precision( p );
	sz = 8 + w + 3 * w_x;
	title_print(out, sz, "spherical harmonic");
	//out << "--------------------------" << std::endl;
	//out << "--- spherical harmonic ---" << std::endl;
	out <<  std::setw( 4 ) << "n" << std::setw( 4 ) << "m" <<
		std::setw( w ) << "d" <<
		std::setw( w_x ) << "x" <<
		std::setw( w_x ) << "y" <<
		std::setw( w_x ) << "z" << std::endl << std::endl;
	out << std::setw( 4 ) << n << std::setw( 4 ) << m;
	for(int i = 0; i < spher.size(); ++i)
	{
		if( i ) out << std::setw(8) << "";
		print_polynom_xyz<T>(out, spher[i], p, w_x);
	}
	out << std::endl;
	out.unsetf( std::ios::scientific );
}
#endif//__spherical_PRINT__
#endif//__polynom_xyz_PRINT__

#endif//__SPHERICAL_H__
