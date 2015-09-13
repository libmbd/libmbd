tests = $(wildcard */)

.PHONY: ${tests}

all: ${tests}

${tests}:
	@echo "Testing $@..."
	@../pymbd.py $@/input.json 2>/dev/null | python compare.py $@/ref.json

