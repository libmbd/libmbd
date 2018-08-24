LIB = libmbd.a

$(LIB): mbd_api.o mbd.o mbd_types.o mbd_vdw_param.o mbd_common.o mbd_linalg.o mbd_parallel.o mbd_defaults.o
	ar -r $@ $^

%.o: %.f90
	$(FXX) $(FXXOPT) -c $<

%.o: %.F90
	$(FXX) $(FXXOPT) $(MACROS) -c $<

mbd_api.o: mbd_common.o mbd_vdw_param.o mbd_defaults.o mbd.o
mbd.o: mbd_common.o mbd_linalg.o mbd_types.o mbd_parallel.o mbd_defaults.o
mbd_linalg.o: mbd_common.o mbd_types.o mbd_parallel.o
mbd_types.o: mbd_common.o mbd_parallel.o
mbd_defaults.o: mbd_common.o
mbd_parallel.o: mbd_common.o
mbd_vdw_param.o: mbd_common.o

.PHONY: clean distclean
clean:
	rm -f *.o

distclean: clean
	rm -f *.mod
	rm -f $(LIB)
