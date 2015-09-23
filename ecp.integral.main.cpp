#include"ecp.integral.demo.h"

int main(int argc, char ** argv)
{
	if( argc != 2 )
	{
		cerr << "Error : [main] usage './main.exe file'" << endl;
		return 1;
	}
	ecp_demo_d(argv[1]);
	return 0;
}
