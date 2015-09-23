
all:
	make dawson_d.o
	make dawson_ld.o
	make dawson.o
	make erfh_d.o
	make erfh_ld.o
	make erfh.o

dawson_d.o: dawson_d.cpp dawson_d.h
	g++ -c dawson_d.cpp -o $@

dawson_ld.o: dawson_ld.cpp dawson_ld.h
	g++ -c dawson_ld.cpp -o $@

dawson.o: dawson.cpp dawson_d.cpp dawson_ld.cpp
	g++ -c dawson.cpp -o $@

erfh_d.o: erfh_d.cpp erfh_d.h
	g++ -c erfh_d.cpp -o $@

erfh_ld.o: erfh_ld.cpp erfh_ld.h
	g++ -c erfh_ld.cpp -o $@

erfh.o: erfh.cpp erfh_d.cpp erfh_ld.cpp
	g++ -c erfh.cpp -o $@

clean:
	rm -f *.exe *.o
