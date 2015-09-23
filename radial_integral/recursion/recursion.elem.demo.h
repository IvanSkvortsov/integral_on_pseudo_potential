#include"spec_func/mpreal.h"
#include"recursion.elem.h"
#include<fstream>
#include<string>

using namespace std;
using mpfr::mpreal;

template<class T>
struct param_set
{
	param_set()
	{
		set_zero();
	}
	virtual ~param_set()
	{
		set_zero();
	}
	T alp, ka_p_kb, abs__ka_m_kb;
	void set_zero()
	{
		//alp = 0; ka_p_kb = 0; abs__ka_m_kb = 0;
	}
};

template<class lT, class rT>
void assign(param_set<lT> & lset, param_set<rT> const & rset)
{
	lset.alp = rset.alp;
	lset.ka_p_kb = rset.ka_p_kb;
	lset.abs__ka_m_kb = rset.abs__ka_m_kb;
}

template<class T>
void demo_dread(T & alp, T & ka_p_kb, T & abs__ka_m_kb)
{
	alp = 10.5;
	ka_p_kb = 3.5;
	abs__ka_m_kb = 1.5;
}

void demo_sread(mpreal & alp, mpreal & ka_p_kb, mpreal & abs__ka_m_kb)
{
	alp = "10.5";
	ka_p_kb = "3.5";
	abs__ka_m_kb = "1.5";
}

template<class T>
void run_recursion_elem(T const & alp, T const & ka_p_kb, T const & abs__ka_m_kb)
{
	recursion_elem<T> elc;
	// init
	elc.init_a(alp);
	elc.init_ka_x_kb(ka_p_kb, abs__ka_m_kb, 1);
	// run
	elc.run();
	// print
	elc.print(cout);
}

void recursion_elem_demo()
{
	mpfr::mpreal::set_default_prec( mpfr::digits2bits(35) );
	mpreal alp_mp, ka_p_kb_mp, abs__ka_m_kb_mp;
	demo_sread(alp_mp, ka_p_kb_mp, abs__ka_m_kb_mp);

	cout << "Double test : " << endl;
	double alp_d, ka_p_kb_d, abs__ka_m_kb_d;
	demo_dread<double>(alp_d, ka_p_kb_d, abs__ka_m_kb_d);
	run_recursion_elem<double>(alp_d, ka_p_kb_d, abs__ka_m_kb_d);
	cout << endl;

	cout << "LDouble test : " << endl;
	long double alp_ld, ka_p_kb_ld, abs__ka_m_kb_ld;
	alp_ld = alp_mp.toLDouble();
	ka_p_kb_ld = ka_p_kb_mp.toLDouble();
	abs__ka_m_kb_ld = abs__ka_m_kb_mp.toLDouble();
	run_recursion_elem<long double>(alp_ld, ka_p_kb_ld, abs__ka_m_kb_ld);
	cout << endl;
	
	cout << "MPReal test : " << endl;
	run_recursion_elem<mpreal>(alp_mp, ka_p_kb_mp, abs__ka_m_kb_mp);
	cout << endl;
}
