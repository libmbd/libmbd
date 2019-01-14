LIB = libmbd.a

$(LIB): mbd.o mbd_blacs.o mbd_c_api.o mbd_constants.o mbd_coulomb.o mbd_damping.o mbd_dipole.o mbd_formulas.o mbd_geom.o mbd_gradients.o mbd_hamiltonian.o mbd_lapack.o mbd_linalg.o mbd_matrix.o mbd_methods.o mbd_mpi.o mbd_rpa.o mbd_scalapack.o mbd_scs.o mbd_ts.o mbd_utils.o mbd_vdw_param.o
	ar -r $@ $^

%.o: %.f90
	$(FXX) $(FXXOPT) -c $<

%.o: %.F90
	$(FXX) $(FXXOPT) -DWITH_SCALAPACK -DWITH_MPI -c $<

mbd.o: mbd_constants.o mbd_damping.o mbd_formulas.o mbd_geom.o mbd_gradients.o mbd_methods.o mbd_ts.o mbd_utils.o mbd_vdw_param.o
mbd_blacs.o: mbd_constants.o
mbd_c_api.o: mbd_constants.o mbd_coulomb.o mbd_damping.o mbd_dipole.o mbd_geom.o mbd_gradients.o mbd_matrix.o mbd_methods.o mbd_ts.o mbd_utils.o
mbd_constants.o: 
mbd_coulomb.o: mbd_constants.o mbd_damping.o mbd_dipole.o mbd_geom.o mbd_lapack.o mbd_linalg.o mbd_matrix.o
mbd_damping.o: mbd_constants.o mbd_gradients.o mbd_utils.o
mbd_dipole.o: mbd_constants.o mbd_damping.o mbd_geom.o mbd_gradients.o mbd_lapack.o mbd_linalg.o mbd_matrix.o mbd_utils.o
mbd_formulas.o: mbd_constants.o mbd_gradients.o mbd_utils.o
mbd_geom.o: mbd_blacs.o mbd_constants.o mbd_lapack.o mbd_mpi.o mbd_utils.o
mbd_gradients.o: mbd_constants.o
mbd_hamiltonian.o: mbd_constants.o mbd_damping.o mbd_dipole.o mbd_geom.o mbd_gradients.o mbd_matrix.o mbd_utils.o
mbd_lapack.o: mbd_constants.o mbd_utils.o
mbd_linalg.o: mbd_constants.o
mbd_matrix.o: mbd_blacs.o mbd_constants.o mbd_lapack.o mbd_scalapack.o mbd_utils.o
mbd_methods.o: mbd_blacs.o mbd_constants.o mbd_damping.o mbd_formulas.o mbd_geom.o mbd_gradients.o mbd_hamiltonian.o mbd_lapack.o mbd_rpa.o mbd_scs.o mbd_utils.o
mbd_mpi.o: 
mbd_rpa.o: mbd_constants.o mbd_damping.o mbd_dipole.o mbd_formulas.o mbd_geom.o mbd_matrix.o mbd_utils.o
mbd_scalapack.o: mbd_blacs.o mbd_constants.o mbd_lapack.o mbd_utils.o
mbd_scs.o: mbd_constants.o mbd_damping.o mbd_dipole.o mbd_formulas.o mbd_geom.o mbd_gradients.o mbd_matrix.o mbd_utils.o
mbd_ts.o: mbd_constants.o mbd_damping.o mbd_geom.o mbd_utils.o
mbd_utils.o: mbd_constants.o mbd_mpi.o
mbd_vdw_param.o: mbd_constants.o mbd_utils.o

.PHONY: clean distclean
clean:
	rm -f *.o

distclean: clean
	rm -f *.mod
	rm -f $(LIB)
