#!/bin/bash

if [ "$#" -eq 0 ];
then
 echo "USAGE :"
 echo "$NAME  (-patch) (-revert)   "
 echo " -compile  : create links and compile " 
 echo " -clean    : make clean and delete links "
fi

case "$1" in
(-clean) 
   echo "doing make clean"
   make clean
   echo "deleing sym links"
   rm -f PLUMED_PATCH biasexchange.c  hills.c recon* metadyn.* restraint* ptmetad.c  read_restraint.c  testderivatives.c
;;
(-compile)
if [ -e PLUMED_PATCH ] ; then
   echo "sym links in place"
   # Just make the program as everything is linked
   make
else
  # Link recon files
  ln -s $plumedir/recon_src/recon_basins.cpp . 
  ln -s $plumedir/recon_src/recon_basins.h .
  ln -s $plumedir/recon_src/recon_utils.h .
  ln -s $plumedir/recon_src/recon_utils.cpp . 
  ln -s $plumedir/recon_src/recon_types.h .
  ln -s $plumedir/recon_src/recon_cbind.h .
  ln -s $plumedir/recon_src/recon_bespoke.h .
  ln -s $plumedir/recon_src/recon_bespoke.cpp .
  ln -s $plumedir/recon_src/recon_cbind.cpp

  ln -s $plumedir/common_files/*.c .
  ln -s $plumedir/common_files/*.h .

  echo "created sym links"   

  # Write some crap
  echo "LINKED ALL THE FILES" > PLUMED_PATCH

  make
fi

esac
