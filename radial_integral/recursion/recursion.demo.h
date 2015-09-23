#include"recursion.h"
using namespace std;

template<class T>
void print( ostream & out, recursion<T> const & recurs, int begin, int end )
{
	cout << "begin : " << begin << endl;
	cout << "end   : " << end << endl;
	cout << endl;
	cout << "recursion::begin " << recurs.begin() << endl;
	cout << "recursion::end   " << recurs.end() << endl;
	cout << endl;
}

template<class T>
void test_recursion(int begin, int end)
{
	cout << "-----------" << endl;
	cout << "begin : " << begin << endl;
	cout << "end   : " << end << endl;
	cout << endl;
	recursion<T> recurs;
	recurs.begin( begin );
	recurs.end( end );
	if( recurs.resize() < 0 )
	{
		cerr << "Error: from test_recursion(int, int)" << endl;
		return;
	}
	cout << "-----------" << endl;
	cout << "recursion::begin " << recurs.begin() << endl;
	cout << "recursion::end   " << recurs.end() << endl;
	cout << "recursion::run(recursion_elem<T> const &)" << endl;
	cout << "-----------" << endl;
	recursion_elem<T> elc;
	elc.run_k( 2.5, 1.0 );
	elc.init_a( 10.5 );
	elc.run();
	//elc.print(cout);
	recurs.run( elc );
	cout << endl;
	cout << "-----------" << endl;
}

template<class T>
void test_recurs()
{
	int begin, end;

	begin = -7;
	end   = 7;
	test_recursion<T>( begin, end );

	begin = 3;
	end   = 7;
	test_recursion<T>( begin, end );

	begin = -7;
	end   = -4;
	test_recursion<T>( begin, end );

	begin = 3;
	end   = 37;
	test_recursion<T>( begin, end );

	begin = 3;
	end   = 1;
	test_recursion<T>( begin, end );

	begin = 0;
	end   = 0;
	test_recursion<T>( begin, end );

	begin = 0;
	end   = 5;
	test_recursion<T>( begin, end );

	begin = -7;
	end   = -7;
	test_recursion<T>( begin, end );

	begin = 7;
	end   = 7;
	test_recursion<T>( begin, end );

}
