SRCDIR ?= $(PWD)
MAKEENV ?= serial
BLDDIR ?= build
MPIFORT ?= gfortran
MPI_NODES ?= 2

all:

test: $(MAKEENV)

serial: serial-.
mpi: mpi-.

serial-%:
	make -C $(BLDDIR)/$* check

mpi-%: serial-%
	env OMP_NUM_THREADS=1 mpirun -n $(MPI_NODES) $(BLDDIR)/$*/tests/mbd_unit_tests
	env OMP_NUM_THREADS=1 mpirun -n $(MPI_NODES) $(BLDDIR)/$*/tests/mbd_api_tests

setup: $(addprefix setup_,$(MAKEENV))

setup_serial-default: | $(BLDDIR)/default
	cd -P $| && cmake $(SRCDIR)

setup_serial-gfortran49: | $(BLDDIR)/gfortran49
	cd -P $| && cmake $(SRCDIR) -DCMAKE_Fortran_COMPILER=gfortran-4.9

setup_serial-gfortran5: | $(BLDDIR)/gfortran5
	cd -P $| && cmake $(SRCDIR) -DCMAKE_Fortran_COMPILER=gfortran-5

setup_mpi-mpi: | $(BLDDIR)/mpi
	cd -P $| && cmake $(SRCDIR) -DENABLE_SCALAPACK_MPI=ON -DCMAKE_Fortran_COMPILER=$(MPIFORT)

setup_mpi-elsi: | $(BLDDIR)/elsi
	cd -P $| && cmake $(SRCDIR) -DENABLE_SCALAPACK_MPI=ON -DENABLE_ELSI=ON -DCMAKE_Fortran_COMPILER=$(MPIFORT)

$(BLDDIR)/%:
	mkdir -p $(BLDDIR)/$*

clean: $(addprefix clean_,$(MAKEENV))

clean_serial-% clean_mpi-%:
	$(MAKE) -C $(BLDDIR)/$* clean

distclean:
	rm -rf $(wildcard $(BLDDIR)/*)
