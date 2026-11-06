#!/bin/sh -e

build_base=${1:-build}
srcdir=$build_base/src
builddir=$build_base/temp.local
installdir=$build_base/local
mkdir -p $srcdir $builddir $installdir


version=1.26
libdeflate_src="https://github.com/ebiggers/libdeflate/releases/download/v$version/libdeflate-$version.tar.gz"
test -d $srcdir/libdeflate-$version || curl -RL $libdeflate_src | (cd $srcdir; tar xz)

set -x
cmake -S $srcdir/libdeflate-$version -B $builddir/libdeflate \
	-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$installdir \
	-DLIBDEFLATE_BUILD_STATIC_LIB=0 -DLIBDEFLATE_BUILD_GZIP=0

cmake --build $builddir/libdeflate
cmake --install $builddir/libdeflate
set +x
