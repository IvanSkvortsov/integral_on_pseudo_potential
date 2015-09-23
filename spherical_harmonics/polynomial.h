#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__
#include"../scp.array.h"

#ifndef __zero_float
#define __zero_float 1.0e-10
#endif

#define __ERROR_polynomial__  "Error: polynomial<pT>::"
#ifdef  __ERROR_polynomial__
#include<iostream>
#endif

template<class T>
bool is_zero_f(T const & x)
{
	if( x == T(0) )
		return true;
	if( x < T(0) )
		if( -x < T(__zero_float) ) return true;
		else return false;
	else if( x < T(__zero_float) )
		return true;
	return false;
}

bool is_equal_x(int const * a, int const * b)
{
	return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

template<class pT>
struct polynomial: public scp_array<pT>
{
	typedef typename pT::value_type value_type;
	typedef pT polynomial_type;
	polynomial():scp_array<pT>(){}
	polynomial(polynomial<pT> const & v):scp_array<pT>( v ){}
	polynomial(int const & __size):scp_array<pT>( __size ){}
	virtual ~polynomial(){}
	polynomial<pT> & operator=(polynomial<pT> const & v)
	{
		if( this == &v ) return *this;
		this->scp_array<pT>::operator=( v );
		return *this;
	}
	polynomial<pT> & operator+=(polynomial<pT> const & v)
	{
		int v_size = v.size(), old_size = this->_size;
		if( !v_size )
			return *this;
		polynomial<pT> tmp( *this );
		this->resize( old_size + v_size );
		this->operator=( tmp );
		this->resize( old_size + v_size );
		for(int i = 0; i < v_size; ++i)
			this->_data[i+old_size] = v[i];
		return *this;
	}
	polynomial<pT> & operator*=(polynomial<pT> const & v)
	{
		if( !this->_size ) return *this;
		else if( !v.size() )
		{
			this->resize( 0 );
			return *this;
		}
		int v_size = v.size(), old_size = this->_size, new_size = this->_size * v_size;
		polynomial<pT> tmp( *this );
		pT * pol_a, * pol_c;
		this->resize( new_size );
		pol_c = this->_data;
		for(int i = 0; i < old_size; ++i)
		{
			pol_a = &tmp[i];
			for(int j = 0; j < v_size; ++j)
				*pol_c++ = *pol_a * v[j];
		}
		return *this;
	}
	polynomial<pT> & operator+=(pT const & pol_b)
	{
		this->append( pol_b );
		return *this;
	}
	polynomial<pT> & operator*=(pT const & pol_b)
	{
		pT * p_pol = this->_data;
		for(int i = 0; i < this->_size; ++i)
		{
			*p_pol *= pol_b;
			++p_pol;
		}
		return *this;
	}
	polynomial<pT> & operator*=(value_type const & val)
	{
		pT * p_pol = this->_data;
		for(int i = 0; i < this->_size; ++i)
		{
			*p_pol++ *= val;
		}
		return *this;
	}
	value_type calc_r(value_type const * r)const
	{
		value_type value = value_type(0);
		pT const * p_pol = this->_data;
		for(int i = 0; i < this->size(); ++i)
			value += p_pol++->calc_r( r );
		return value;
	}
	int is_even()const
	{
		pT const * p_pol = this->_data;
		for(int i = 0; i < this->_size; ++i)
			if( p_pol++->is_even_x() )
				return 1;
		return 0;
	}
	void optimize_equal()
	{
		pT * p_a = this->_data, * p_b = 0;
		for(int i = 0; i < this->_size; ++i)
		{
			p_b = p_a + 1;
			for(int j = i+1; j < this->_size; ++j)
			{
				//if( is_equal_x( &p_a->x, &p_b->x) )
				if( pT::is_equal_x( *p_a, *p_b) )
				{
					p_a->d += p_b->d;
					p_b->d = value_type(0);
				}
				++p_b;
			}
			++p_a;
		}
	}
	void optimize_zero()
	{
		int iter = 0;
		pT * p_a = this->_data, * p_b = p_a + this->_size - 1, tmp_v;
		for(int i = 0; i < this->_size; ++i)
		{
			if( is_zero_f<value_type>( p_a->d ) )
			{
				if( p_a > p_b )
					break;
				for(int j = i+1; j < this->_size; ++j)
				{
					if( !is_zero_f<value_type>( p_b->d ) )
					{
						tmp_v = *p_a;
						*p_a = *p_b;
						*p_b = tmp_v;
						--p_b;
						break;
					}
					else
						++iter;
					--p_b;
				}
				++iter;
			}
			++p_a;
		}
		if( iter >= this->_size )
		{
#ifdef  __ERROR_polynomial__
			std::cerr << "Error: [optimize_zero]" << std::endl;
			std::cerr << "size : " << this->_size << std::endl;
			std::cerr << "iter : " << iter << std::endl;
#endif
			exit(1);
		}
		this->_size -= iter;
	}
	void set(int const & __size, int const & __capacity, pT * p)
	{
		this->_size = __size;
		this->_capacity = __capacity;
		this->_data = p;
	}
	void zero()
	{
		this->_size = 0;
		this->_capacity = 0;
		this->_data = 0;
	}
};

template<class pT>
polynomial<pT> operator+(polynomial<pT> const & a, polynomial<pT> const & b)
{
	polynomial<pT> c( a );
	c += b;
	return c;
}

template<class pT>
polynomial<pT> operator+(polynomial<pT> const & a, pT const & b)
{
	polynomial<pT> c( a );
	c += b;
	return c;
}

template<class pT>
polynomial<pT> operator+(pT const & a, polynomial<pT> const & b)
{
	polynomial<pT> c( b );
	c += a;
	return c;
}

template<class pT>
polynomial<pT> operator*(polynomial<pT> const & a, polynomial<pT> const & b)
{
	polynomial<pT> c( a );
	c *= b;
	return c;
}

template<class pT>
polynomial<pT> operator*(polynomial<pT> const & a, pT const & b)
{
	polynomial<pT> c( a );
	c *= b;
	return c;
}

template<class pT>
polynomial<pT> operator*(pT const & a, polynomial<pT> const & b)
{
	polynomial<pT> c( b );
	c *= a;
	return c;
}

#endif//__POLYNOMIAL_H__
