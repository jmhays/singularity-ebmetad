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
   rm -f PLUMED_PATCH recon* restraint* metadyn* biasexchange.c hills.c ptmetad.c read_restraint.c testderivatives.c plumed.inc
;;
(-compile)
if [ -e PLUMED_PATCH ] ; then
   echo "sym links in place"
   # Just make the program as everything is linked
   make
else
  # Link recon files
  ln -s $plumedir/recon_src/*.cpp .
  ln -s $plumedir/recon_src/*.h .

  # Link plumed files
  ln -s $plumedir/common_files/*.c .
  ln -s $plumedir/common_files/*.h . 

  echo "created sym links"   

  # Write some crap
  echo "LINKED ALL THE FILES" > PLUMED_PATCH

  # This create the include file for the makefile
  {
   echo -n "PLUMED_OBJECTS=" 
   for file in $plumedir/common_files/*.c
     do f=${file##*/}
     echo " \\"
     echo -n " ${f%.c}.o"
   done
   for file in $plumedir/recon_src/*.cpp
     do f=${file##*/}
     echo " \\"
     echo -n " ${f%.cpp}.o"
   done
  } > plumed.inc

  make
fi

esac
