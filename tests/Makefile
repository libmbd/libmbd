tests = $(wildcard */)

.PHONY: ${tests}

all: ${tests}

${tests}:
	@echo "Testing $@ in serial..."
	@../pymbd.py <$@/input.json 2>/dev/null | python compare.py $@/ref.json
	@echo "Testing $@ in parallel..."
	@mpiexec -n 2 ../pymbd.py $@/input.json 2>/dev/null | python compare.py $@/ref.json

