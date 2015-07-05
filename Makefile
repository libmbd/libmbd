include system.mk

blddir = build

all: mbd.so
	mpiexec -n 2 python pymbd.py

mbd.so: mbd_interface.f90 mbd.f90
	mkdir -p ${blddir}
	CFLAGS="${CFLAGS}" f2py -c --build-dir ${blddir} --fcompiler=${FVENDOR} \
		   --f90exec=${FC} --f90flags="${FFLAGS}" --compiler=${CVENDOR} \
		   -m $(basename $@) ${LDFLAGS} $^
	rsync -a ${blddir}/*.mod .
	rm -r ${blddir}

distclean:
	-rm mbd.so
	-rm *.mod
