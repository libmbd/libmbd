SOURCES = mbd.f90 mbd_api.f90 mbd_common.f90 mbd_defaults.f90 mbd_linalg.f90 mbd_parallel.f90 mbd_types.f90 mbd_vdw_param.f90

VPATH = ../../src

all: serial_dist scalapack_dist

serial_dist: $(addprefix serial/,$(SOURCES))

scalapack_dist: $(addprefix scalapack/,$(SOURCES))

serial/%.f90: %.f90 | serial
	cp $^ $@

serial/%.f90: %.F90 | serial
	gfortran -E -P $^ -o $@

scalapack/%.f90: %.f90 | scalapack
	cp $^ $@

scalapack/%.f90: %.F90 | scalapack
	gfortran -E -P $^ -o $@ -DWITH_SCALAPACK

serial scalapack:
	mkdir $@

clean:
	rm -rf serial scalapack
