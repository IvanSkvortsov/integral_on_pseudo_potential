#include"ecp.integral.h"

using namespace std;

template<class T>
void print_basis_prim(std::ostream & out, Primitive<T> const & a, int prec = 10)
{
	int w = prec + 8;
	for(int i = 0; i < 3; ++i) out << setw(4) << a.i[i];
	for(int i = 0; i < 3; ++i) out << setw(w) << a.r[i];
	out << setw( w ) << a.alp;
	out << endl;
}

template<class T>
void print_ecp_prim(std::ostream & out, ecp_primitive<T> const & ecp, int prec = 10)
{
	int w = prec + 8;
	out << setw(4) << ecp.l << setw(4) << ecp.n;
	out << setw(4) << "";
	for(int i = 0; i < 3; ++i) out << setw(w) << ecp.r[i];
	out << setw( w ) << ecp.alp;
	out << endl;
}

template<class T>
void ecp_demo(char const * file)
{
	Primitive<T> a, b;
	vector<ecp_primitive<T> > v_ecp;
	read_task<T>( a, b, v_ecp, file );
	cout << "task read from file \"" << file << "\"" << endl;
	T value = T(0);
	cout.setf( ios::scientific );
	int prec = 16, w = prec + 8;
	cout.precision( prec );
	cout << setw(4) << "i" << setw(4) << "j" << setw(4) << "k" << 
		setw(w) << "x" << setw(w) << "y" << setw(w) << "z" << setw(w) << "alpha" << endl;
	print_basis_prim<T>( cout, a, prec );
	print_basis_prim<T>( cout, b, prec );
	cout << endl;
	cout << setw(4) << "l" << setw(4) << "n" << setw(4) << "" << 
		setw(w) << "x" << setw(w) << "y" << setw(w) << "z" << setw(w) << "alpha" << endl;
	print_ecp_prim<T>( cout, v_ecp[0], prec );
	//
	for(int i = 0; i < v_ecp.size(); ++i)
	{
		value = ecp_semi_local_integral<T>( a, b, v_ecp[i] );
		cout << setw(14) << "ecp_integral" << " : " << setw( w ) << value << setw(4) << v_ecp[i].l << endl;// << endl;
		//cout << setw(14) << "ecp_integral" << " : " << setw( w ) << value / (4 * T(Pi)) << endl << endl;
		//cout << setw(14) << "ecp_integral" << " : " << setw( w ) << value / (16 * T(Pi) * T(Pi)) << endl << endl;
		//cout << setw(14) << "ecp_integral" << " : " << setw( w ) << value / (64 * T(Pi) * T(Pi) * T(Pi)) << endl << endl;
		//cout << setw(14) << "ecp_integral" << " : " << setw( w ) << value / (256 * T(Pi) * T(Pi) * T(Pi) * T(Pi)) << endl << endl;
		//print_ecp_prim<T>( cout, v_ecp[i], prec );
	}
}

//using mpfr::mpreal;

void ecp_demo_d(char const * file )
{
	//mpfr::mpreal::set_default_prec( mpfr::digits2bits(30) );
	//ecp_demo<mpreal>( file );
	ecp_demo<double>( file );
}
