WORKDIR=$HOME
export PATH=$PATH:$WORKDIR/bin
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$WORKDIR/include
export CPP_INCLUDE_PATH=$CPP_INCLUDE_PATH:$WORKDIR/include
export LIBRARY_PATH=$LIBRARY_PATH:$WORKDIR/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WORKDIR/lib
# Build Lowdin
./configure -p $WORKDIR/bin -s /tmp -l "-lblas -llapack"
make
make install
