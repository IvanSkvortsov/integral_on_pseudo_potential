all:
	make recursion.elem.exe
	make recursion.exe

recursion.elem.exe: recursion.elem.main.cpp recursion.elem.h ./spec_func/erfh.o ./spec_func/dawson.o
	g++ recursion.elem.main.cpp spec_func/erfh.o spec_func/dawson.o -lgmp -lmpfr -o $@
	./$@

recursion.exe: recursion.main.cpp recursion.h recursion.elem.h spec_func/erfh.o spec_func/dawson.o
	g++ recursion.main.cpp spec_func/erfh.o spec_func/dawson.o -o $@
	./$@

taus.exe: taus.h taus.demo.h taus.main.cpp recursion.h recursion.elem.h recursion.bc.h spec_func/erfh.o spec_func/dawson.o
	g++ taus.main.cpp spec_func/erfh.o spec_func/dawson.o -o $@
	./$@

spec_func/daswon.o: spec_func/daswon.cpp spec_func/daswon_d.cpp spec_func/daswon_ld.cpp 
	g++ -c spec_func/daswon.cpp -o $@

spec_func/erfh.o: spec_func/erfh.cpp spec_func/erfh_d.cpp spec_func/erfh_ld.cpp
	g++ -c spec_func/erfh.cpp -o $@

clean:
	rm -f *.o *.exe
	rm -f spec_func/*.o spec_func/*.exe
