include system.mk
blddir = build
sources = mbd_interface.f90 mbd.f90

all: mbd

mbd:
	mkdir -p ${blddir}
	CFLAGS="${CFLAGS}" f2py -c --build-dir ${blddir} --fcompiler=${FVENDOR} \
		   --f90exec=${FC} --f90flags="${FFLAGS}" --compiler=${CVENDOR} \
		   -m $@ ${LDFLAGS} ${sources}
	rsync -a ${blddir}/*.mod .

test:
	@${MAKE} -C tests

clean:
	-rm *.mod
	-rm -r build

distclean: clean
	-rm mbd.so
