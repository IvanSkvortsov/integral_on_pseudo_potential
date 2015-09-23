#ifndef __SCP_ARRAY_H__
#define __SCP_ARRAY_H__
#include<stdlib.h>

#define __scp_error_message_ON__

#ifdef  __scp_error_message_ON__
#include<iostream>
#endif

template<class T>
struct scp_array
{
	int _size, _capacity;
	T * _data;
	scp_array():_size(0), _capacity(0), _data(0){}
	scp_array(scp_array<T> const & v):_size(v._size), _capacity(_size), _data(new T[_size])
	{
		for(int i = 0; i < _size; ++i)
			_data[i] = v._data[i];
	}
	scp_array(int const & __size):_size(__size), _capacity(_size), _data(new T[_size]){}
	virtual ~scp_array()
	{
		if( _data )
			delete [] _data;
		_size = 0;
		_capacity = 0;
		_data = 0;
	}
	scp_array<T> & operator=(scp_array<T> const & v)
	{
		if( this == &v ) return *this;
		this->resize( v.size() );
		for(int i = 0; i < _size; ++i)
			_data[i] = v[i];
		return *this;
	}
	T & operator[](int i){return _data[i];}
	T const & operator[](int i)const{return _data[i];}
	int const & size()const{return _size;}
	int const & capacity()const{return _capacity;}
	void append(T const & v)
	{
		if( _capacity > _size )
		{
			_data[_size] = v;
			++_size;
		}
		else
		{
			scp_array<T> tmp( *this );
			this->resize( _size + 1 );
			for(int i = 0; i < tmp._size; ++i)
				_data[i] = tmp._data[i];
			_data[_size-1] = v;
		}
	}
	void resize(int __size)
	{
		if( __size < 0 )
		{
#ifdef  __scp_error_message_ON__
			std::cerr << "Error: scp_array<T>::resize(int __size)" << std::endl;
			std::cerr << "__size : " << __size << std::endl;
#endif
			exit(1);
		}
		if( _capacity < __size )
		{
			if( _data ) delete [] _data;
			_size = __size;
			if( _capacity == 0 )
				_capacity = __size;
			else
				while( _capacity < __size ) _capacity *= 2;
			if( _capacity ) _data = new T[_capacity];
		}
		else _size = __size;
	}
	void reserve(int __capacity)
	{
		if( __capacity < 0 )
		{
#ifdef  __scp_error_message_ON__
			std::cerr << "Error: scp_array<T>::reserve(int __capacity)" << std::endl;
			std::cerr << "__capacity : " << __capacity << std::endl;
#endif
			exit(1);
		}
		if( _capacity < __capacity )
		{
			if( _data ) delete [] _data;
			if( _capacity == 0 )
				_capacity = __capacity;
			else
				while( _capacity < __capacity ) _capacity *= 2;
			if( _capacity ) _data = new T[_capacity];
		}
	}
};

#endif//__SCP_ARRAY_H__
