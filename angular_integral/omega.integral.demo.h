#include"omega.integral.h"
#include<iostream>
#include<iomanip>
#include<fstream>

using namespace std;

template<class T>
int read_index(omega_index & idx, T * r1, T * r2, char const * file)
{
	ifstream inp( file );
	if( !inp.is_open() )
	{
		cerr << "Error: [read_index(omega_idex & idx, char const * file)] can't open file \"" << file << "\"" << endl;
		return 1;
	}
	int a, b, c;
	inp >> a >> b >> c;
	idx.set_indexA( a, b, c );
	inp >> a >> b >> c;
	idx.set_indexB( a, b, c );
	inp >> idx.l >> idx.lmb1 >> idx.lmb2;
	for(int i = 0; i < 3; ++i) inp >> r1[i];
	for(int i = 0; i < 3; ++i) inp >> r2[i];
	return 0;
}


template<class T>
void print_index(std::ostream & out, omega_index const & idx, T const * r1, T const * r2)
{
	out << "-- index --" << endl;
	out << setw(4) << "abc" << " : " << setw(4) << idx.a << setw(3) << idx.b << setw(3) << idx.c << endl;
	out << setw(4) << "def" << " : " << setw(4) << idx.d << setw(3) << idx.e << setw(3) << idx.f << endl;
	//out << setw(4) << "a" << " : " << setw(4) << idx.a << endl;
	//out << setw(4) << "b" << " : " << setw(4) << idx.b << endl;
	//out << setw(4) << "c" << " : " << setw(4) << idx.c << endl;
	//out << setw(4) << "d" << " : " << setw(4) << idx.d << endl;
	//out << setw(4) << "e" << " : " << setw(4) << idx.e << endl;
	//out << setw(4) << "f" << " : " << setw(4) << idx.f << endl;
	out << setw(4) << "l" << " : " << setw(4) << idx.l << endl;
	out << setw(4) << "lmb1" << " : " << setw(4) << idx.lmb1 << endl;
	out << setw(4) << "lmb2" << " : " << setw(4) << idx.lmb2 << endl;
	int prec = 12, w = prec + 8;
	out.setf( ios::scientific );
	out.precision( prec );
	out << setw(4) << "r1" << " : " << 
		setw(w) << r1[0] <<
		setw(w) << r1[1] <<
		setw(w) << r1[2] << endl;
	out << setw(4) << "r2" << " : " << 
		setw(w) << r2[0] <<
		setw(w) << r2[1] <<
		setw(w) << r2[2] << endl;
	out << "-----------" << endl;
	out.unsetf( ios::scientific );
}

template<class T>
int demo_omega_int(char const * file)
{
	omega_index idx;
	T r1[3], r2[3];
	if( read_index( idx, r1, r2, file ) )
		exit(1);
	print_index(cout, idx, r1, r2);
	int maxl = idx.max_l();
	if( maxl == 0 )
	{
		cerr << "Warning: maximum angular momentum is zero [nothing to do here]" << endl;
		return 1;
	}
	vector<vector<spherical<T> > > v2sph( maxl+1 );
	for(int i = 0; i < v2sph.size(); ++i)
	{
		v2sph[i].resize( 2 * i + 1 );
		for(int j = 0; j < v2sph[i].size(); ++j)
		{
			v2sph[i][j].run(i, j-i);
		}
	}
	for(int m = -idx.lmb1; m <= idx.lmb1; ++m)
		print_spherical<T>(std::cout, v2sph[idx.lmb1][idx.lmb1 + m], idx.lmb1, m);
	for(int m = -idx.lmb2; m <= idx.lmb2; ++m)
		print_spherical<T>(std::cout, v2sph[idx.lmb2][idx.lmb2 + m], idx.lmb2, m);
	angular_omega_xyz<T> omega_xyz_;
	omega_xyz_.run( 50, 50, 50 );
	spherical<T> sph_buf;
	T value = omega_sph3(idx, r1, r2, v2sph, omega_xyz_, sph_buf);
	int prec = 16, w = prec + 8;
	cout << endl;
	cout.precision( prec );
	cout.setf( ios::scientific );
	cout << setw( w ) << value << endl;
}

template<class T>
int demo_omega_integral(char const * file)
{
	omega_index idx;
	T r1[3], r2[3];
	if( read_index( idx, r1, r2, file ) )
		exit(1);
	print_index(cout, idx, r1, r2);
	int maxl = idx.max_l();
	if( maxl == 0 )
	{
		cerr << "Warning: maximum angular momentum is zero [nothing to do here]" << endl;
		return 1;
	}
	omega_integral<T> o_i;
	_abc_ la( idx.a, idx.b, idx.c ), lb( idx.d, idx.e, idx.f );
	angular_omega_xyz<T> omega_xyz_;
	omega_xyz_.run( 10, 10, 10 );
	o_i.run( idx.l, la, lb, r1, r2, omega_xyz_ );
	char file_out[] = "omega.integral.log";
	ofstream out( file_out );
	print_omega_integral(out, o_i);
	cout << "results've been written into \"" << file_out << "\"" << endl;
}
int demo_omega_integral_d(char const * file)
{
	//return demo_omega_int<double>( file );
	return demo_omega_integral<double>( file );
}
