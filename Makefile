include system.mk
sources = mbd_interface.f90 mbd.f90

blddir = build

all: mbd.so

mbd.so: ${sources}
	mkdir -p ${blddir}
	CFLAGS="${CFLAGS}" f2py -c --build-dir ${blddir} --fcompiler=${FVENDOR} \
		   --f90exec=${FC} --f90flags="${FFLAGS}" --compiler=${CVENDOR} \
		   -m $(basename $@) ${LDFLAGS} $^
	rsync -a ${blddir}/*.mod .
	rm -r ${blddir}

test:
	python test.py

test2:
	python -m doctest -v README.md 

clean:
	-rm *.mod

distclean: clean
	-rm mbd.so
