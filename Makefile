SRCDIR ?= $(PWD)
BUILDS ?= Sbuild-serial Sbuild-serial-4.9 Mbuild-mpi

all:

test: $(BUILDS)

S%:
	make -C $* check

M%: S%
	mpirun -n 2 $*/src/mbd_tests
	mpirun -n 2 $*/src/mbd_api_tests

doc:
	cd docs && doxygen

setup: $(addprefix setup_,$(BUILDS))

setup_Sbuild-serial: DIR_serial
	cd -P build-serial && cmake $(SRCDIR)

setup_Sbuild-serial-4.9: DIR_serial-4.9
	cd -P build-serial-4.9 && cmake $(SRCDIR) -DCMAKE_Fortran_COMPILER=gfortran-4.9

setup_Mbuild-mpi: DIR_mpi
	cd -P build-mpi && cmake $(SRCDIR) -DENABLE_SCALAPACK_MPI=ON

DIR_%:
	mkdir -p build-$*

clean: $(addprefix clean_,$(BUILDS))

clean_S% clean_M%:
	$(MAKE) -C $* clean

distclean: $(addprefix distclean_,$(BUILDS))

distclean_S% distclean_M%:
	rm -rf $*/*
