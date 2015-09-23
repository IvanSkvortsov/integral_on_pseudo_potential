#ifndef __TAUS_H__
#define __TAUS_H__
#include"recursion.h"

template<class T>
void taus_run_1(std::vector<T> & v, int & begin, int & end, T * & tau, T * & rho, T * & sigma, T * & _sigma,
		recursion<T> const & recurs, int sign__ka_m_kb)
{
	int size = recurs.size();
	recursion_bc<T> const * ptr;
	plus_minus<T> const * b, * c;
	T const & gm = recurs.gm(), & gp = recurs.gp();
	// set begin, end
	begin = recurs.begin();
	end = recurs.end();
	// set size
	v.resize( 4 * size );
	// set pointers
	tau = v.data() - begin;
	rho = tau + size;
	sigma = rho + size;
	_sigma = sigma + size;
	// calculate by using recursion functions
	for(int i = begin; i <= end; ++i)
	{
		ptr = recurs.p_bc(i);
		b = &(ptr->b);
		c = &(ptr->c);
		tau[i] = gp * b->p + gm * b->m;
		rho[i] = gp * b->p - gm * b->m;
		 sigma[i] = gp * c->p + sign__ka_m_kb * gm * c->m;
		_sigma[i] = gp * c->p - sign__ka_m_kb * gm * c->m;
	}
}

template<class T>
void taus_run_2(std::vector<T> & v, int & begin, int & end, T * & kappa, T * & etta, recursion<T> const & recurs)
{
	int size = recurs.size();
	recursion_bc<T> const * ptr;
	plus_minus<T> const * b, * c;
	T const & gp = recurs.gp();
	// set begin, end
	begin = recurs.begin();
	end = recurs.end();
	// set size
	v.resize( 2 * size );
	kappa = v.data() - begin;
	etta = kappa + size;
	// calculate by using recursion functions
	for(int i = begin; i <= end; ++i)
	{
		ptr = recurs.p_bc(i);
		b = &(ptr->b);
		c = &(ptr->c);
		kappa[i] = 2 * gp * b->p;
		etta[i]  = 2 * gp * c->p;
	}
}

#endif//__TAUS_H__
