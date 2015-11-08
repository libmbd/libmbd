tests = argon-dimer argon-double-chain

.PHONY: ${tests}

all: ${tests}

${tests}:
	@echo "Testing $@ in serial..."
	@python3 $@/test.py >/dev/null
	@echo "Testing $@ in parallel..."
	@mpiexec -n 2 python3 $@/test.py >/dev/null
