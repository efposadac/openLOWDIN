# openLOWDIN Quantum chemistry package #

This program has been developed under direction of:

* [Prof. A REYES' Lab](http://www.qcc.unal.edu.co). Universidad Nacional de Colombia.

[![Makefile CI](https://github.com/efposadac/openLOWDIN/actions/workflows/makefile.yml/badge.svg)](https://github.com/efposadac/openLOWDIN/actions/workflows/makefile.yml)

© All rights reserved, 2023.

Welcome to openLOWDIN Quantum Chemistry Package.

Installation notes.
=============

### Prerequisites: ###

* A standard FORTRAN compiler. gfortran and intel FORTRAN compiler have been tested.
* Lapack or MKL libraries.
* Arpack library.
* LIBINT library version 2.0, and version 1.1.5
* LIBXC library

NOTE: If you have the libraries in your own path please be sure to export the LIBRARY_PATH environment variable. ie:
`export LIBRARY_PATH=$LIBRARY_PATH:[your library path]`

### Compile: (see below for a step-by-step example) ###

* run `./configure` in LOWDIN root directory. Be sure that you have permissions to write in the installation directory and have properly exported the `$PATH` environment.

* run `make`

### Install: ###

* run `make install`

### Uninstall ###

* run `make uninstall`

### Documentation ###

* run `make doc`

The `make doc` command produces both latex and html documentation using doxygen program. Be sure you have installed doxygen, for instance in a debian-based distribution run:

`# apt-get install doxygen graphviz`

To use latex documentation in doc/latex folder, run command:

`pdflatex refman.tex`

To visualize the html documentation use:

`<web browser> doc/html/index.html`

### Clean the project ###

* run `make clean` and then `make distclean`

---

### Step-by-step  installation example: (replace apt-get with your preferred package manager) ###

        sudo apt-get update
        sudo apt-get -y install wget git build-essential liblapack-dev libblas-dev libgsl0-dev autotools-dev automake libtool gfortran python3 gawk libeigen3-dev libgmp-dev libboost-all-dev libarpack2-dev
        # Define ENV Variables
        export WORKDIR=$PWD/dependencies
        export PATH=$PATH:$WORKDIR/bin
        export C_INCLUDE_PATH=$C_INCLUDE_PATH:$WORKDIR/include:$WORKDIR/include/libint2:/usr/include/eigen3
        export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$WORKDIR/include:$WORKDIR/include/libint2:/usr/include/eigen3
        export LIBRARY_PATH=$LIBRARY_PATH:$WORKDIR/lib
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WORKDIR/lib
        # Create work directories
        mkdir $WORKDIR
        mkdir $WORKDIR/bin
        mkdir $WORKDIR/lib
        cd $WORKDIR

	# Libint2
        # If you have Ubuntu, you can get this precompiled Libint2 library
        wget https://www.dropbox.com/s/d3d44j238lkfwcr/libint-master-SEP052019.tgz
        tar xzvf libint-master-SEP052019.tgz
	# Otherwise, download and compile with minimal (default am), G12, fPIC options (libint2 commit 668b10c4bdca5876984058742d4212675eb93f3f)
        cd -
	
        # Libint1
        git clone https://github.com/evaleev/libint.git
        cd libint
        git checkout v1
        aclocal -I lib/autoconf
        autoconf
        ./configure --prefix=$WORKDIR
        make -j 4
        make install
        make clean
        make distclean
        cd -

	# Libxc
        cd $WORKDIR        
        # If you have Ubuntu, you can get this precompiled Libxc library
        wget https://www.dropbox.com/s/6cja3zzhl1cq46i/libxc-master-MAY242023.tgz
        tar xzvf libxc-master-MAY242023.tgz
	# Otherwise, download and compile with default options (libxc commit 4bd0e1e36347c6d0a4e378a2c8d891ae43f8c951)
        cd ..
	
        # Configure Lowdin
        ./configure -p $WORKDIR/bin -s /tmp -l "-lblas -llapack -larpack"
        # Build Lowdin
        make -j 4
        # Install Lowdin
        make install
        # Run Tests
        make test

### Further info: ###
fernando.posada@temple.edu
felix.moncada@fysik.su.se
