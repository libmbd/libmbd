set -v

mkdir build
(cd build && cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_INSTALL_LIBDIR=lib)
make -C build
make -C build install
env LIBMBD_PREFIX=$PREFIX poetry build
$PYTHON -m pip install pymbd -f ./dist --no-deps --ignore-installed --no-cache-dir -vvv
