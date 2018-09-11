Bootstrap: docker
From: ubuntu:16.04


%files

    md_meta_EBMetaD/ 	/opt


%environment

    plumedir=/builds/md_meta_EBMetaD
    PATH=/usr/local/gromacs/bin:${PATH}

    export plumedir PATH


%labels

   AUTHOR jmh5sf@virginia.edu


%post

    apt-get update && apt-get -y install wget mpich gcc-5 g++-5 libfftw3-dev make patch

    mkdir /builds
    cd /builds
    mv /opt/md_meta_EBMetaD /builds

    # Get gromacs 4.5.1
    wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-4.5.1.tar.gz
    tar xf gromacs-4.5.1.tar.gz
    cd gromacs-4.5.1
    ./configure --enable-mpi

    # Patch with plumed
    export plumedir=/builds/md_meta_EBMetaD
    cp ${plumedir}/patches/plumedpatch_gromacs_4.5.1.sh .
    bash plumedpatch_gromacs_4.5.1.sh -patch

    # Now make
    make -j8; make install
