#!/usr/bin/env python3
import re

tests = []
with open('tests/mbd_grad_test_cases.f90') as f:
    for l in f:
        m = re.match(r' *subroutine *test_(\w+) *\(', l)
        if m:
            tests.append(f'grad/{m.group(1)}')
with open('tests/mbd_api_tests.F90') as f:
    for l in f:
        m = re.match(r' *case *\( *["\'](\w+)["\'] *\)', l)
        if m:
            tests.append(f'api/{m.group(1)}')
print(';'.join(tests))
