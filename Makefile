c_warnings_off = incompatible-pointer-types parentheses-equality shorten-64-to-32 \#warnings
f_warnings_off = maybe-uninitialized
blddir = build
CFLAGS ?= $(addprefix -Wno-,${c_warnings_off})
F90FLAGS ?= $(addprefix -Wno-,${f_warnings_off})
LDFLAGS ?= -lscalapack
FC ?= gfortran
MPIFC ?= mpifort

all: mbd.so
	mpiexec -n 2 python pymbd.py

mbd.so: mbd.f90 ${blddir}/mbd_interface.o
	mkdir -p ${blddir}
	CFLAGS="${CFLAGS}" f2py -c --build-dir ${blddir} --f90exec=$(shell which ${FC}) --f90flags="${F90FLAGS}" -m $(basename $@) ${LDFLAGS} $^
	rsync -a ${blddir}/*.mod .

${blddir}/%.o: %.f90
	mkdir -p ${blddir}
	${MPIFC} -c -J${blddir} -o $@ $<

build:
	mkdir build
	
clean:
	-rm -r build

distclean: clean
	-rm mbd.so
	-rm *.mod
