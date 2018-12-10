SRCDIR ?= $(PWD)
MAKEENV ?= serial-build

all:

test: $(MAKEENV)

serial-%:
	make -C $* check

mpi-%: serial-%
	mpirun -n 2 $*/tests/mbd_unit_tests
	mpirun -n 2 $*/tests/mbd_api_tests

doc:
	cd docs && doxygen

setup: $(addprefix setup_,$(MAKEENV))

setup_serial-build: | build
	cd -P $| && cmake $(SRCDIR)

setup_serial-build-gfortran49: | build-gfortran49
	cd -P $| && cmake $(SRCDIR) -DCMAKE_Fortran_COMPILER=gfortran-4.9

setup_serial-build-gfortran5: | build-gfortran5
	cd -P $| && cmake $(SRCDIR) -DCMAKE_Fortran_COMPILER=gfortran-5

setup_mpi-build-mpi: | build-mpi
	cd -P $| && cmake $(SRCDIR) -DENABLE_SCALAPACK_MPI=ON

build build-gfortran49 build-mpi:
	mkdir -p $@

clean: $(addprefix clean_,$(MAKEENV))

clean_serial-% clean_mpi-%:
	$(MAKE) -C $* clean

distclean: $(addprefix distclean_,$(MAKEENV))

distclean_serial-% distclean_mpi-%:
	rm -rf $*/*
