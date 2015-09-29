
ecp.integral.exe: ecp.integral.h ecp.integral.demo.h ecp.integral.main.cpp \
	radial_integral/recursion/spec_func/dawson_d.o \
	radial_integral/recursion/spec_func/dawson_ld.o \
	radial_integral/recursion/spec_func/erfh_d.o\
	radial_integral/recursion/spec_func/erfh_ld.o \
	radial_integral/recursion/recursion.h radial_integral/recursion/recursion.elem.h radial_integral/recursion/recursion.bc.h \
	radial_integral/recursion/taus.h \
	radial_integral/lib_qu/qu.h radial_integral/lib_qu/hankel.h \
	radial_integral/radial.h \
	angular_integral/angular.omega.xyz.h angular_integral/omega.integral.h \
	spherical_harmonics/spherical.h spherical_harmonics/legendre.h spherical_harmonics/polynomial.h spherical_harmonics/scp.array.h
	g++ ecp.integral.main.cpp \
	radial_integral/recursion/spec_func/dawson_d.o \
	radial_integral/recursion/spec_func/dawson_ld.o \
	radial_integral/recursion/spec_func/erfh_d.o \
	radial_integral/recursion/spec_func/erfh_ld.o \
	-o $@
	./$@ ecp.integral.inp

clean:
	rm ecp.integral.exe
