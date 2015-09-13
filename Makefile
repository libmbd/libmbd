include system.mk
blddir = build
sources = mbd_interface.f90 mbd.f90

all: mbd.so

mbd.so: ${sources}
	mkdir -p ${blddir}
	CFLAGS="${CFLAGS}" f2py -c --build-dir ${blddir} --fcompiler=${FVENDOR} \
		   --f90exec=${FC} --f90flags="${FFLAGS}" --compiler=${CVENDOR} \
		   -m $(basename $@) ${LDFLAGS} $^
	rsync -a ${blddir}/*.mod .

test:
	@${MAKE} -C tests

clean:
	-rm *.mod
	-rm -r build

distclean: clean
	-rm mbd.so
