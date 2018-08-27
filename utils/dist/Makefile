SOURCES = mbd.f90 mbd_api.f90 mbd_common.f90 mbd_constants.f90 mbd_lapack.f90 mbd_linalg.f90 mbd_matrix_type.f90 mbd_system_type.f90 mbd_ts.f90 mbd_vdw_param.f90

VPATH = ../../src

all: serial_dist scalapack_mpi_dist

serial_dist: $(addprefix serial/,$(SOURCES))

scalapack_mpi_dist: $(addprefix scalapack_mpi/,$(SOURCES))

serial/%.f90: %.f90 | serial
	cp $^ $@

serial/%.f90: %.F90 | serial
	gfortran -E -P $^ -o $@

scalapack_mpi/%.f90: %.f90 | scalapack_mpi
	cp $^ $@

scalapack_mpi/%.f90: %.F90 | scalapack_mpi
	gfortran -E -P $^ -o $@ -DWITH_MPI -DWITH_SCALAPACK

serial scalapack_mpi:
	mkdir $@

clean:
	rm -rf serial scalapack_mpi
