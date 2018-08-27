LIB = libmbd.a

$(LIB): mbd.o mbd_api.o mbd_c_api.o mbd_common.o mbd_constants.o mbd_coulomb.o mbd_damping_type.o mbd_dipole.o mbd_gradients_type.o mbd_lapack.o mbd_linalg.o mbd_matrix_type.o mbd_system_type.o mbd_ts.o mbd_vdw_param.o
	ar -r $@ $^

%.o: %.f90
	$(FXX) $(FXXOPT) -c $<

%.o: %.F90
	$(FXX) $(FXXOPT) $(MACROS) -c $<

mbd.o: mbd_common.o mbd_constants.o mbd_damping_type.o mbd_dipole.o mbd_gradients_type.o mbd_lapack.o mbd_matrix_type.o mbd_system_type.o
mbd_api.o: mbd.o mbd_common.o mbd_constants.o mbd_damping_type.o mbd_gradients_type.o mbd_system_type.o mbd_ts.o mbd_vdw_param.o
mbd_c_api.o: mbd.o mbd_constants.o mbd_coulomb.o mbd_damping_type.o mbd_dipole.o mbd_gradients_type.o mbd_matrix_type.o mbd_system_type.o mbd_ts.o
mbd_common.o: mbd_constants.o
mbd_constants.o: 
mbd_coulomb.o: mbd_constants.o mbd_damping_type.o mbd_dipole.o mbd_lapack.o mbd_linalg.o mbd_matrix_type.o mbd_system_type.o
mbd_damping_type.o: mbd_common.o mbd_constants.o mbd_gradients_type.o
mbd_dipole.o: mbd_common.o mbd_constants.o mbd_damping_type.o mbd_gradients_type.o mbd_lapack.o mbd_linalg.o mbd_matrix_type.o mbd_system_type.o
mbd_gradients_type.o: mbd_constants.o
mbd_lapack.o: mbd_common.o mbd_constants.o
mbd_linalg.o: mbd_constants.o
mbd_matrix_type.o: mbd_common.o mbd_constants.o mbd_lapack.o
mbd_system_type.o: mbd_common.o mbd_constants.o mbd_lapack.o mbd_matrix_type.o
mbd_ts.o: mbd_common.o mbd_constants.o mbd_damping_type.o mbd_system_type.o
mbd_vdw_param.o: mbd_common.o mbd_constants.o

.PHONY: clean distclean
clean:
	rm -f *.o

distclean: clean
	rm -f *.mod
	rm -f $(LIB)
