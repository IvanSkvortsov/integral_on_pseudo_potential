#include"omega.integral.demo.h"

int main(int argc, char *argv[])
{
	if( argc != 2 )
	{
		cerr << "Error: [main] usage './omega.integral.exe file.inp'" << endl;
		return 1;
	}
	demo_omega_integral_d(argv[1]);
	return 0;
}
