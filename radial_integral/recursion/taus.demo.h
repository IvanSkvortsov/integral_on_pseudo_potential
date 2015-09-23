#include"taus.h"

using namespace std;

template<class T>
void taus_print(std::ostream & out, T const * t, int begin, int end, int prec = 10)
{
	out.setf( std::ios::scientific );
	out.precision( prec );
	int w = prec + 8;
	for(int i = begin; i <= end; ++i)
		out << setw(4) << i << setw(w) << t[i] << endl;
	out << endl;
}

template<class T>
void taus_demo()
{
	recursion_elem<T> elc;
	elc.run_k( 2.5, 1.0 );
	elc.init_a( 10.5 );
	elc.run();
	//
	recursion<T> recurs;
	recurs.resize( 3, 8 );
	recurs.run( elc );
	//
	vector<T> v;
	T * tau, * rho, * sigma, * _sigma;
	int begin, end;
	taus_run_1<T>( v, begin, end, tau, rho, sigma, _sigma, recurs, 1);
	taus_print<T>( cout, tau, begin, end );
}

void taus_demo_d()
{
	taus_demo<double>();
}
