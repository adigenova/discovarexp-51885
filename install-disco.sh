ROOT=$PWD
CORES=1

if [ "$(uname)" == "Darwin" ]; then
    # Do something under Mac OS X platform        
    CORES=$(sysctl -n hw.ncpu)
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # Do something under GNU/Linux platform
    CORES=$(nproc)
fi


cd ${ROOT}/deps/jemalloc
./autogen.sh
#lib jmallloc
./configure --prefix=$PWD/ljemalloc
make -j ${CORES}
make install
#to link static the library
cp $PWD/ljemalloc/lib/libjemalloc.a ${ROOT}/src/
#we comeback to root directory
cd ${ROOT}
#dynamic link to jemalloc
#./configure --prefix=$PWD/dnovo --with-jemalloc=${ROOT}/deps/jemalloc/ljemalloc/lib
#static link to jemalloc
./configure --prefix=$PWD/dnovo
make -j ${CORES}
make install
#strip the code
strip $PWD/dnovo/bin/*
