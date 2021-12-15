ifndef VIRTUAL_ENV
$(error Must be run inside a Python virtual environment)
endif

BLDDIR ?= $(CURDIR)/build
SRCDIR = $(CURDIR)
export LIBMBD_PREFIX = $(VIRTUAL_ENV)
ifdef MPI_NODES
RUN_CMD = env OMP_NUM_THREADS=1 mpiexec -n $(MPI_NODES)
endif
PYMBD_EXTRAS = test
ifneq (,$(findstring ENABLE_SCALAPACK_MPI=ON,$(CMAKE_ARGS)))
PYMBD_EXTRAS += mpi
endif
EMPTY =
SPACE = $(EMPTY) $(EMPTY)
COMMA = ,

all: install_editable test

install_libmbd:
	cmake -S $(SRCDIR) -B $(BLDDIR) -DCMAKE_INSTALL_PREFIX=$(LIBMBD_PREFIX) $(CMAKE_ARGS)
	make -C $(BLDDIR) all install

install_editable: install_libmbd
	poetry install $(foreach ext,$(PYMBD_EXTRAS),-E $(ext))

install: install_libmbd
	pip install cffi
	poetry build
	pip install pymbd[$(subst $(SPACE),$(COMMA),$(PYMBD_EXTRAS))] --pre -f ./dist

test_libmbd:
	ctest --test-dir $(BLDDIR) --output-on-failure

test: test_libmbd
	$(RUN_CMD) pytest -v --durations=3 $(PYTEST_FLAGS)

doc:
	pip install sphinx toml git+https://github.com/libmbd/ford@7b44574da7ec20f4ab4b1842ec7561de2a601930
	ford -I. doc/libmbd.md -o build
	sphinx-build -W -d $(BLDDIR)/doctrees doc doc/build/pymbd
	touch doc/build/.nojekyll

distclean:
	-rm -r $(BLDDIR)/*
	-rm src/pymbd/_libmbd.*.so
