include system.mk
blddir = build
extern = mbd_interface.f90 mbd_helper.f90
FFLAGS = -Og -fcheck=all

all: mbd

mbd: $(addprefix ${blddir}/,${extern:.f90=.o})
	@mkdir -p ${blddir}
	CFLAGS="${CFLAGS}" f2py -c --build-dir ${blddir} --fcompiler=${FVENDOR} \
		   --f90exec=${FC} --f90flags="${FFLAGS}" --compiler=${CVENDOR} \
		   $(addprefix ${blddir}/,mbd_interface.o mbd_helper.o) -m $@ ${LDFLAGS} mbd.f90
	rsync -a ${blddir}/*.mod .

${blddir}/%.o: %.f90
	@mkdir -p ${blddir}
	${FC} -c -fPIC -J ${blddir} ${FFLAGS} -o $@ $^

test:
	@${MAKE} -C tests

clean:
	-rm *.mod
	-rm -r build

distclean: clean
	-rm mbd.so
