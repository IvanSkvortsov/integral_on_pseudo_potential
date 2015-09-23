#ifndef __RECURSION_TYPE_H__
#define __RECURSION_TYPE_H__
#include<iostream>
#include<fstream>
#include<iomanip>
#include"../../math/l_mymath.h"

// struct recursion_t
//
// alp          = alp_A_i + alp_B_j + alp_C_lk
// alp_A_i      - power of gaussian exponent of i'th function of A basis set
// alp_B_j      - power of gaussian exponent of j'th function of B basis set
// alp_C_lk     - power of gaussian exponent of pseudo-potential k'th function with angular momentum l
//
// d            - coefficient of gaussian exponent of pseudo-potential (consistent with alp_C_lk)
//
// alp_i        - see alp_A_i
// alp_j        - see alp_B_j
// 
// ca_len       - sqrt( CA[0] * CA[0] + CA[1] * CA[1] + CA[2] * CA[2] )
// cb_len       - sqrt( CB[0] * CB[0] + CB[1] * CB[1] + CB[2] * CB[2] )
// CA           - (C - A)
// CB           - (C - B)
// C            - center of pseudo-potential
// A            - center of basis set A
// B            - center of basis set B
//
// ka_len       - ca_len * 2 * alp_A_i
// kb_len       - cb_len * 2 * alp_B_j
//
// *_p_*        - here 'p' means 'plus'
// *_m_*        - here 'm' means 'minus'
// ka_p_kb      - ka_len + kb_len
// abs__ka_m_kb - abs( ka_len - kb_len )
// sgn__ka_m_kb - sign( ka_len - kb_len )
//
// nk           - power of r in pseudo-potential
//
// lmax         - maximum angular momentum of orbitals in core space of atom (pseudo-potential)
// l            - angular momentum of pseudo-potential function
// la           - angular momentum of A basis set's orbital
// lb           - angular momentum of B basis set's orbital
//
// l_p_la       = l + la
// l_p_lb       = l + lb
// lmax_p_la    = lmax + la
// lmax_p_lb    = lmax + lb
// la_p_lb      = la + lb

template<class T>
struct recursion_t
{
	T alp, d, alp_i, alp_j, ka_p_kb, abs__ka_m_kb, ka_len, kb_len, ca_len, cb_len;
	int nk, la, lb, l, lmax, sgn__ka_m_kb, l_p_la, l_p_lb, lmax_p_la, lmax_p_lb;
	int la_p_lb, begin, end, N_size, lmb_a_size, lmb_b_size;
	//---- constructor -----//
	recursion_t();
	//----- destructor -----//
	virtual ~recursion_t(){set_zero();}
	void set_zero();
	T init_geometry_a(T const & __ca_len){return ca_len = __ca_len;}
	T init_geometry_b(T const & __cb_len){return cb_len = __cb_len;}
	T init_alp(T const & __alp){return alp = __alp;}
	T run_alp_i(T const & __alp_i){alp_i = __alp_i; return ka_len = 2 * ca_len * alp_i;}
	T run_alp_j(T const & __alp_j){alp_j = __alp_j; return kb_len = 2 * cb_len * alp_j;}
	void run_kx()
	{
		ka_p_kb = ka_len + kb_len;
		abs__ka_m_kb = ka_len - kb_len;
		sgn__ka_m_kb = sgn( abs__ka_m_kb );
		if( sgn__ka_m_kb == -1 )
			abs__ka_m_kb = -abs__ka_m_kb;
	}
	int init_lmax(int __lmax){return lmax = __lmax;}
	int init_momentum_a(int __la){lmax_p_la = lmax + __la;return la = __la;}
	int init_momentum_b(int __lb){lmax_p_lb = lmax + __lb;return lb = __lb;}
	int run_la_p_lb()
	{
		la_p_lb = la + lb;
		N_size = la_p_lb + 1;
		return la_p_lb;
	}
	int run_la_p_lb(int __la, int __lb)
	{
		if( la != __la )
			init_momentum_a(__la);
		if( lb != __lb )
			init_momentum_b(__lb);
		la_p_lb = __la + __lb;
		N_size = la_p_lb + 1;
		return la_p_lb;
	}
	int run_momentum_pp(int __l)
	{
		l_p_la = __l + la;
		l_p_lb = __l + lb;
		lmb_a_size = l_p_la + 1;
		lmb_b_size = l_p_lb + 1;
		return l = __l;
	}
	void run_nk(int __nk)
	{
		// begin = nk - 2 - (l + la) - (l + lb)
		// end   = nk - 2 + la + lb
		nk = __nk;
		end = nk - 2;
		begin = end - (2 * l + la_p_lb);
		end += la_p_lb;
	}
	void reverse();
	recursion_t<T> & operator=(recursion_t<T> const & rec_t);
};

// reverse
template<class T>
void recursion_t<T>::reverse()
// la, lb, l_p_la, l_p_lb, lmb_a_size, lmb_b_size
// ka_len, kb_len, ca_len, cb_len
// alp_i, alp_j
// sgn__ka_m_kb
{
	int tmp = 0;
	// la     <-> lb
	tmp = la;
	la = lb;
	lb = tmp;
	// l_p_la <-> l_p_lb
	tmp = l_p_la;
	l_p_la = l_p_lb;
	l_p_lb = tmp;
	// lmb_a_size <-> lmb_b_size
	tmp = lmb_a_size;
	lmb_a_size = lmb_b_size;
	lmb_b_size = lmb_a_size;
	T tmp_f = 0;
	// ka_len <-> kb_len
	tmp_f = ka_len;
	ka_len = kb_len;
	kb_len = tmp_f;
	// ca_len <-> cb_len
	tmp_f = ca_len;
	ca_len = cb_len;
	cb_len = tmp_f;
	// alp_i  <-> alp_j
	tmp_f = alp_i;
	alp_i  = alp_j;
	alp_j  = tmp_f;
	// sgn__ka_m_kb
	sgn__ka_m_kb = -sgn__ka_m_kb;
}

// set_zero
template<class T>
void recursion_t<T>::set_zero()
{
	alp = 0;
	alp_i = 0;
	alp_j = 0;
	ka_p_kb = 0;
	abs__ka_m_kb = 0;
	sgn__ka_m_kb = 0;
	ka_len = 0;
	kb_len = 0;
	ca_len = 0;
	cb_len = 0;
	nk = 0;
	d = 0;
	la = 0;
	lb = 0;
	l = 0;
	lmax = 0;
	lmax_p_la = 0;
	lmax_p_lb = 0;
	//
	l_p_la = 0;
	l_p_lb = 0;
	la_p_lb = 0;
	begin = 0;
	end = 0;
	N_size = 0;
	lmb_a_size = 0;
	lmb_b_size = 0;
	return;
}
// constructor
template<class T>
recursion_t<T>::recursion_t(): lmax(0), lmax_p_la(0), lmax_p_lb(0), l(0), la(0), lb(0),
	l_p_la(0), l_p_lb(0), la_p_lb(0),
	d(0), alp(0), alp_i(0), alp_j(0), nk(0),
	ka_p_kb(0), abs__ka_m_kb(0), sgn__ka_m_kb(0), ka_len(0), kb_len(0), ca_len(0), cb_len(0),
	begin(0), end(0),
	N_size(0), lmb_a_size(0), lmb_b_size(0){}

// operator=
template<class T>
recursion_t<T> & recursion_t<T>::operator=(recursion_t<T> const & rec_t)
{
	alp = rec_t.alp;
	alp_i = rec_t.alp_i;
	alp_j = rec_t.alp_j;
	ka_p_kb = rec_t.ka_p_kb;
	abs__ka_m_kb = rec_t.abs__ka_m_kb;
	sgn__ka_m_kb = rec_t.sgn__ka_m_kb;
	ka_len = rec_t.ka_len;
	kb_len = rec_t.kb_len;
	ca_len = rec_t.ca_len;
	cb_len = rec_t.cb_len;
	nk = rec_t.nk;
	d = rec_t.d;
	la = rec_t.la;
	lb = rec_t.lb;
	l = rec_t.l;
	l_p_la = rec_t.l_p_la;
	l_p_lb = rec_t.l_p_lb;
	lmax = rec_t.lmax;
	lmax_p_la = rec_t.lmax_p_la;
	lmax_p_lb = rec_t.lmax_p_lb;
	//
	l_p_la = rec_t.l_p_la;
	l_p_lb = rec_t.l_p_lb;
	la_p_lb = rec_t.la_p_lb;
	begin = rec_t.begin;
	end = rec_t.end;
	N_size = rec_t.N_size;
	lmb_a_size = rec_t.lmb_a_size;
	lmb_b_size = rec_t.lmb_b_size;
	return *this;
}
// print
// T alp, alp_i, alp_j, ka_p_kb, abs__ka_m_kb, ka_len, kb_len, ca_len, cb_len;
// int nk, d, la, lb, l, lmax, sgn__ka_m_kb, l_p_la, l_p_lb;
// qu
// int la_p_lb, begin, end, N_size, lmb_a_size, lmb_b_size;
template<class T>
void print( ostream & out, recursion_t<T> const & rec_t, int w = 14, int p = 6 )
{
	out << "---- RECURSION_T_H ----" << endl;
	out << setw( w ) << "l" << " : " << rec_t.l << endl;
	out << setw( w ) << "la" << " : " << rec_t.la << endl;
	out << setw( w ) << "lb" << " : " << rec_t.lb << endl;
	out << setw( w ) << "lmax" << " : " << rec_t.lmax << endl;
	out << setw( w ) << "lmax_p_la" << " : " << rec_t.lmax_p_la << endl;
	out << setw( w ) << "lmax_p_lb" << " : " << rec_t.lmax_p_lb << endl;
	out << setw( w ) << "nk" << " : " << rec_t.nk << endl;
	//
	out << setw( w ) << "la_p_lb" << " : " << rec_t.la_p_lb << endl;
	out << setw( w ) << "begin" << " : " << rec_t.begin << endl;
	out << setw( w ) << "end" << " : " << rec_t.end << endl;
	out << setw( w ) << "N_size" << " : " << rec_t.N_size << endl;
	out << setw( w ) << "lmb_a_size" << " : " << rec_t.lmb_a_size << endl;
	out << setw( w ) << "lmb_b_size" << " : " << rec_t.lmb_b_size << endl;
	//
	out.setf( ios::fixed );
	out << setw( w ) << "d" << " : " << setprecision( p ) << rec_t.d << endl;
	out << setw( w ) << "alp" << " : " << setprecision( p ) << rec_t.alp << endl;
	out << setw( w ) << "alp_i" << " : " << setprecision( p ) << rec_t.alp_i << endl;
	out << setw( w ) << "alp_j" << " : " << setprecision( p ) << rec_t.alp_j << endl;
	out << setw( w ) << "ka_len" << " : " << setprecision( p ) << rec_t.ka_len << endl;
	out << setw( w ) << "kb_len" << " : " << setprecision( p ) << rec_t.kb_len << endl;
	out << setw( w ) << "ca_len" << " : " << setprecision( p ) << rec_t.ca_len << endl;
	out << setw( w ) << "cb_len" << " : " << setprecision( p ) << rec_t.cb_len << endl;
	out << setw( w ) << "ka_p_kb" << " : " << setprecision( p ) << rec_t.ka_p_kb << endl;
	out << setw( w ) << "abs__ka_m_kb" << " : " << setprecision( p ) << rec_t.abs__ka_m_kb << endl;
	out << setw( w ) << "sgn__ka_m_kb" << " : " << rec_t.sgn__ka_m_kb << endl;
	out << "--------------------" << endl;
	out.unsetf( ios::fixed );
}
#endif//__RECURSION_TYPE_H__
