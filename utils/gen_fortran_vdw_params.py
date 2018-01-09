import csv


with open('pymbd/vdw-params.csv') as f:
    reader = csv.DictReader(f, quoting=csv.QUOTE_NONNUMERIC)
    rows = list(reader)
n = len(rows)
print(f'real(dp), parameter :: default_vdw_params(3, {n}) = reshape([ &')
for i, row in enumerate(rows):
    print('    {0[alpha_0]}d0, {0[C6]}d0, {0[R_vdw]}d0{1} &  ! {0[species]}'.format(
        row, ',' if i < n-1 else ''
    ))
print(f'], [3, {n}])')
