from pathlib import Path
import json
import cffi

from configure import conf_file

if Path(conf_file).exists():
    with open(conf_file) as f:
        blddir = json.load(f)['blddir']
root = Path(__file__).parent

header = 'mbd.h'
ffibuilder = cffi.FFI()
ffibuilder.set_source(
    'pymbd._libmbd',
    f'#include "{header}"',
    include_dirs=[str(root/'src')],
    libraries=['mbd'],
    library_dirs=[blddir],
)
ffibuilder.cdef((Path('src')/header).read_text())

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
