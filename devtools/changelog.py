#!/usr/bin/env python3
import re
import subprocess
from string import Template

END = 'END'
s = subprocess.run(
    'git tag -l "*.*.*"'
    ' --format "%(refname:strip=2),%(taggerdate:short),'
    f'%(contents:subject)\n\n%(contents:body){END}"'
    ' --sort=-creatordate',
    shell=True,
    capture_output=True,
).stdout.decode()
data = [
    (v, d, m.strip() if m[0] == '#' else None)
    for v, d, m in re.findall(rf'(\d+\.\d+\.\d+),([\d-]+),(.*?){END}', s, re.DOTALL)
]
with open('devtools/CHANGELOG.md.in') as f:
    t = Template(f.read())
changes = [f'## [{v}] - {d}\n\n{m}' for v, d, m in data if m]
link = '[{0}]: https://github.com/libmbd/libmbd/{1}'
vs = [x[0] for x in data]
links = [
    *(
        link.format(
            vs[i - 1] if i else 'unreleased',
            f'compare/{vs[i]}...{vs[i-1] if i else "HEAD"}',
        )
        for i in range(len(vs))
    ),
    link.format(vs[-1], f'releases/tag/{vs[-1]}'),
]
with open('CHANGELOG.md', 'w') as f:
    f.write(t.substitute({'changes': '\n\n'.join(changes), 'links': '\n'.join(links)}))
