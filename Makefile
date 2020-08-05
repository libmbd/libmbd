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
ifneq (,$(findstring ENABLE_SCALAPACK_MPI=ON,$(CMAKE_FLAGS)))
PYMBD_EXTRAS += mpi
endif
EMPTY =
SPACE = $(EMPTY) $(EMPTY)
COMMA = ,

all: install_editable test

install_libmbd:
	mkdir -p $(BLDDIR)
	cd $(BLDDIR) && cmake $(SRCDIR) -DCMAKE_INSTALL_PREFIX=$(LIBMBD_PREFIX) $(CMAKE_FLAGS)
	make -C $(BLDDIR) all install

install_editable: install_libmbd
	poetry install $(foreach ext,$(PYMBD_EXTRAS),-E $(ext))

install: install_libmbd
	pip install cffi
	poetry build
	python -m pip install pymbd[$(subst $(SPACE),$(COMMA),$(PYMBD_EXTRAS))] -f ./dist

test_libmbd:
	cd $(BLDDIR) && ctest --output-on-failure

test: test_libmbd
	$(RUN_CMD) pytest -v --durations=3 $(PYTEST_FLAGS)

doc:
	python -m pip install sphinx toml git+https://github.com/jhrmnn/ford@7b44574da7ec20f4ab4b1842ec7561de2a601930
	ford -I. docs/libmbd.md -o build
	sphinx-build -W -d $(BLDDIR)/doctrees docs docs/build/pymbd
	touch docs/build/.nojekyll
