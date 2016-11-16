FVENDOR = gnu95
CVENDOR = unix
FC = mpifort
LDFLAGS = --link-lapack_opt
FFLAGS = -Og -fcheck=all
F2PY = f2py
BLDDIR = build
-include system.mk

all: mbd

mbd: mbd.f90 $(addprefix $(BLDDIR)/,mbd_interface.o mbd_helper.o)
	@mkdir -p $(BLDDIR)
	$(F2PY) -c --build-dir $(BLDDIR) --fcompiler=$(FVENDOR) \
		   --f90exec=$(FC) --f90flags="$(FFLAGS)" --compiler=$(CVENDOR) \
		   $(wordlist 2,3,$^) -m $@ $(LDFLAGS) $<
	rsync -a $(BLDDIR)/*.mod .

mbd_math: mbd_math.f90 $(addprefix $(BLDDIR)/,mbd.o mbd_interface.o mbd_helper.o)
	@mkdir -p $(BLDDIR)
	$(F2PY) -c --build-dir $(BLDDIR) --fcompiler=$(FVENDOR) \
		   --f90exec=$(FC) --f90flags="$(FFLAGS)" --compiler=$(CVENDOR) \
		   $(filter-out $<,$^) -m $@ $(LDFLAGS) $<
	rsync -a $(BLDDIR)/*.mod .

$(BLDDIR)/%.o: %.f90
	@mkdir -p $(BLDDIR)
	$(FC) -c -fPIC -J $(BLDDIR) $(FFLAGS) -o $@ $^

clean:
	rm -f *.mod
	rm -rf build

distclean: clean
	rm -f mbd.*so
	rm -rf mbd.*dSYM
