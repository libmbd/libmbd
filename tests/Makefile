tests = $(wildcard */)

.PHONY: ${tests}

all: ${tests}

${tests}:
	@echo "Testing $@ in serial..."
	@../pymbd.py <$@/input.json | python3 compare.py $@/ref.json
	@echo "Testing $@ in parallel..."
	@mpiexec -n 2 ../pymbd.py $@/input.json | python3 compare.py $@/ref.json

