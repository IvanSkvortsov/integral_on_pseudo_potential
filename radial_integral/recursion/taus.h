#ifndef __TAUS_H__
#define __TAUS_H__
#include"recursion.h"

#define __taus_print__
#ifdef  __taus_print__
#include<iostream>
#include<iomanip>
#endif
template<class T>
struct taus
{
	std::vector<T> v;
	T * tau, * sigma, * _sigma, * rho;
	T * kappa, * eta;
	int min_n, max_n;
	taus();
	void run_1b(recursion<T> const & recurs);
	void run_2b(recursion<T> const & recurs, int sgn__ka_m_kb);
#ifdef  __taus_print__
	void print(      std::ostream & out = std::cout)const;
	//
	void tau_print(  std::ostream & out, int prec)const;
	void kappa_print(std::ostream & out, int prec)const;
#endif
};

template<class T>
taus<T>::taus():v(), tau(0), sigma(0), _sigma(0), rho(0), kappa(0), eta(0), min_n(0), max_n(0){}

template<class T>
void taus<T>::run_1b(recursion<T> const & recurs)
{
	min_n = recurs.min_n();
	max_n = recurs.max_n();
	int sz = max_n - min_n + 1;
	v.resize( 2 * sz );
	for(int i = 0; i < v.size(); ++i) v[i] = T(0);
	this->kappa = &v[0] - min_n;
	this->eta = &v[0] + (sz - min_n);
	//
	recursion_bc<T> const * ptr;
	plus_minus<T> const * b, * c;
	T const & gp = recurs.gp();
	// calculate by using recursion functions
	for(int i = min_n; i <= max_n; ++i)
	{
		ptr = recurs.p_bc(i);
		b = &(ptr->b);
		c = &(ptr->c);
		kappa[i] = 2 * gp * b->p;
		eta[i]  = 2 * gp * c->p;
	}
	this->tau = 0; this->sigma = 0; this->_sigma = 0; this->rho = 0;
}

template<class T>
void taus<T>::run_2b(recursion<T> const & recurs, int sgn__ka_m_kb)
{
	min_n = recurs.min_n();
	max_n = recurs.max_n();
	int sz = max_n - min_n + 1;
	v.resize( 4 * sz );
	for(int i = 0; i < v.size(); ++i) v[i] = T(0);
	this->tau = &v[0] - min_n;
	this->sigma = &v[0] + (sz - min_n);
	this->_sigma = &v[0] + (2*sz - min_n);
	this->rho = &v[0] + (3*sz - min_n);
	//
	recursion_bc<T> const * ptr;
	plus_minus<T> const * b, * c;
	T const & gm = recurs.gm(), & gp = recurs.gp();
	// calculate by using recursion functions
	for(int i = min_n; i <= max_n; ++i)
	{
		ptr = recurs.p_bc(i);
		b = &(ptr->b);
		c = &(ptr->c);
		//tau[i] = b->p + b->m;
		//rho[i] = b->p - b->m;
		// sigma[i] = c->p + sign__ka_m_kb * c->m;
		//_sigma[i] = c->p - sign__ka_m_kb * c->m;
		tau[i] = gp * b->p + gm * b->m;
		rho[i] = gp * b->p - gm * b->m;
		 sigma[i] = gp * c->p + sgn__ka_m_kb * gm * c->m;
		_sigma[i] = gp * c->p - sgn__ka_m_kb * gm * c->m;
	}
	this->kappa = 0; this->eta = 0;
}


#ifdef  __taus_print__
template<class T>
void taus<T>::tau_print(std::ostream & out, int prec)const
{
	int w = prec + 8;
	std::string s = " taus "; 
	int t_width = 4 * w + 4 - s.size(), w_ = t_width/2;
	// title
	for(int i = 0; i < w_; ++i) out << '-'; out << s; for(int i = 0; i < t_width-w_; ++i) out << '-'; out << std::endl;
	out <<  std::setw(4) << "i" << std::setw(w) << "tau" << std::setw(w) << "sigma" <<
		std::setw(w) << "_sigma" << std::setw(w) << "rho" <<
		std::endl << std::endl;
	for(int i = this->min_n; i <= this->max_n; ++i)
		out <<  std::setw(4) << i <<
			std::setw(w) << tau[i] << std::setw(w) << sigma[i] <<
			std::setw(w) << _sigma[i] << std::setw(w) << rho[i] <<
			std::endl;
	out << std::endl;
	out << std::setw(4) << "max" << " : " << std::setw(4) << this->max_n << std::endl;
	out << std::setw(4) << "min" << " : " << std::setw(4) << this->min_n << std::endl;
	out << std::endl;
}

template<class T>
void taus<T>::kappa_print(std::ostream & out, int prec)const
{
	int w = prec + 8;
	std::string s = " taus "; 
	int t_width = 2 * w + 4 - s.size(), w_ = t_width/2;
	// title
	for(int i = 0; i < w_; ++i) out << '-'; out << s; for(int i = 0; i < t_width-w_; ++i) out << '-'; out << std::endl;
	out << std::setw(4) << "i" << std::setw(w) << "eta" << std::setw(w) << "kappa" << std::endl << std::endl;
	for(int i = this->min_n; i <= this->max_n; ++i)
		out << std::setw(4) << i << std::setw(w) << eta[i] << std::setw(w) << kappa[i] << std::endl;
	out << std::endl;
	out << std::setw(4) << "max" << " : " << std::setw(4) << this->max_n << std::endl;
	out << std::setw(4) << "min" << " : " << std::setw(4) << this->min_n << std::endl;
	out << std::endl;
}
template<class T>
void taus<T>::print(std::ostream & out)const
{
	int prec = 16;
	out.setf( std::ios::scientific );
	out.precision( prec );
	if( tau ) this->tau_print(out, prec);
	if( kappa ) this->kappa_print(out, prec);
	out.unsetf( std::ios::scientific );
}
#endif//__taus_print__

template<class T>
void taus_run_2b(std::vector<T> & v, int & min_n, int & max_n, T * & tau, T * & rho, T * & sigma, T * & _sigma,
		recursion<T> const & recurs, int sign__ka_m_kb)
{
	int size = recurs.size();
	recursion_bc<T> const * ptr;
	plus_minus<T> const * b, * c;
	T const & gm = recurs.gm(), & gp = recurs.gp();
	// set min_n, max_n
	min_n = recurs.min_n();
	max_n = recurs.max_n();
	// set size
	v.resize( 4 * size );
	// set pointers
	tau = &v[0] - min_n;
	rho = &v[0] + (size - min_n);
	sigma = &v[0] + (2*size - min_n);
	_sigma = &v[0] + (3*size - min_n);
	// calculate by using recursion functions
	for(int i = min_n; i <= max_n; ++i)
	{
		ptr = recurs.p_bc(i);
		b = &(ptr->b);
		c = &(ptr->c);
		//tau[i] = b->p + b->m;
		//rho[i] = b->p - b->m;
		// sigma[i] = c->p + sign__ka_m_kb * c->m;
		//_sigma[i] = c->p - sign__ka_m_kb * c->m;
		tau[i] = gp * b->p + gm * b->m;
		rho[i] = gp * b->p - gm * b->m;
		 sigma[i] = gp * c->p + sign__ka_m_kb * gm * c->m;
		_sigma[i] = gp * c->p - sign__ka_m_kb * gm * c->m;
	}
}

template<class T>
void taus_run_1b(std::vector<T> & v, int & min_n, int & max_n, T * & kappa, T * & eta, recursion<T> const & recurs)
{
	int size = recurs.size();
	recursion_bc<T> const * ptr;
	plus_minus<T> const * b, * c;
	T const & gp = recurs.gp();
	// set min_n, max_n
	min_n = recurs.min_n();
	max_n = recurs.max_n();
	// set size
	v.resize( 2 * size );
	kappa = &v[0] - min_n;
	eta = &v[0] + (size - min_n);
	// calculate by using recursion functions
	for(int i = min_n; i <= max_n; ++i)
	{
		ptr = recurs.p_bc(i);
		b = &(ptr->b);
		c = &(ptr->c);
		kappa[i] = 2 * gp * b->p;
		eta[i]  = 2 * gp * c->p;
	}
}

// syntactic sugar, and mainly for compatibility with older code
template<class T>
void taus_run_1(std::vector<T> & v, int & min_n, int & max_n, T * & tau, T * & rho, T * & sigma, T * & _sigma,
		recursion<T> const & recurs, int sign__ka_m_kb)
{
	return taus_run_2b<T>(v, min_n, max_n, tau, rho, sigma, _sigma, recurs, sign__ka_m_kb);
}
template<class T>
void taus_run_2(std::vector<T> & v, int & min_n, int & max_n, T * & kappa, T * & eta, recursion<T> const & recurs)
{
	return taus_run_1b<T>(v, min_n, max_n, kappa, eta, recurs);
}
#endif//__TAUS_H__
