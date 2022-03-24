ifndef VIRTUAL_ENV
$(error Must be run inside a Python virtual environment)
endif

BLDDIR ?= $(CURDIR)/build
export LIBMBD_PREFIX = $(VIRTUAL_ENV)
ifdef MPI_NODES
override RUN_CMD := env OMP_NUM_THREADS=1 mpiexec $(MPIEXEC_EXTRA_FLAGS) -n $(MPI_NODES) $(RUN_CMD)
endif
PYMBD_EXTRAS = test
ifneq (,$(findstring ENABLE_SCALAPACK_MPI=ON,$(CMAKE_ARGS)))
PYMBD_EXTRAS += mpi
endif
EMPTY =
SPACE = $(EMPTY) $(EMPTY)
COMMA = ,

all: install_editable test

run_cmake:
	cmake -B $(BLDDIR) -DCMAKE_INSTALL_PREFIX=$(LIBMBD_PREFIX) $(CMAKE_ARGS)

build_libmbd: run_cmake
	make -C $(BLDDIR) all

install_libmbd: build_libmbd
	make -C $(BLDDIR) install

install_editable: install_libmbd
	pip install -e .[$(subst $(SPACE),$(COMMA),$(PYMBD_EXTRAS))]

install: install_libmbd
	pip install .[$(subst $(SPACE),$(COMMA),$(PYMBD_EXTRAS))]

test_libmbd:
	ctest --test-dir $(BLDDIR) --output-on-failure

test: test_libmbd
	$(RUN_CMD) pytest -v --durations=3

build_doc:
	pip install "sphinx<3" "Jinja2<3.1" toml git+https://github.com/libmbd/ford@7b44574da7ec20f4ab4b1842ec7561de2a601930
	ford -I. doc/libmbd.md -o build
	sphinx-build -W -d $(BLDDIR)/doctrees doc doc/build/pymbd
	touch doc/build/.nojekyll

distclean:
	-rm -r $(BLDDIR)/*
	-rm src/pymbd/_libmbd.*.so
