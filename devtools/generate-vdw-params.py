#!/usr/bin/env python3
import csv
from importlib.resources import open_text

with open_text('pymbd', 'vdw-params.csv') as f:
    reader = csv.DictReader(f, quoting=csv.QUOTE_NONNUMERIC)
    rows = list(reader)
rows = [row for row in rows if row['symbol'][-1] not in '+-']
n = len(rows)
print(f'real(dp), parameter :: ts_vdw_params(3, {n}) = reshape([ &')
for i, row in enumerate(rows):
    alpha_0 = row['alpha_0(TS)'] or 0
    C6 = row['C6(TS)'] or 0
    R_vdw = row['R_vdw(TS)'] or 0
    print(
        f'    {alpha_0}d0, {C6}d0, {R_vdw}d0{"," if i < n-1 else ""} &'
        f'  ! {row["symbol"]}'
    )
print(f'], [3, {n}])')
print()
print(f'real(dp), parameter :: tssurf_vdw_params(3, {n}) = reshape([ &')
for i, row in enumerate(rows):
    alpha_0 = row['alpha_0(TSsurf)'] or row['alpha_0(TS)'] or 0
    C6 = row['C6(TSsurf)'] or row['C6(TS)'] or 0
    R_vdw = row['R_vdw(TSsurf)'] or row['R_vdw(TS)'] or 0
    print(
        f'    {alpha_0}d0, {C6}d0, {R_vdw}d0{"," if i < n-1 else ""} &'
        f'  ! {row["symbol"]}'
    )
print(f'], [3, {n}])')
