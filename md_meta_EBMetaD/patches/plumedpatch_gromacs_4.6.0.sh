#!/bin/bash
# PATCH SCRIPT FOR GROMACS
#

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# definitions specific to this code
CODE="gromacs-4.6.0"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
CAMSHIFT_LINKS="$plumedir/camshift_src/*.cpp $plumedir/camshift_src/*.h"
RECON_LINKS="$plumedir/recon_src/*.cpp $plumedir/recon_src/*.h"
WHERE_LINKS="./src/kernel/"
# Leave this line here blank for normal patching of plumed
#RECON_LIBS= 
# CVS version only : for reconnaissance metadynamics enter here information on the
# locations of the lstdc++ libraries.
RECON_LIBS="-lstdc++"
# Leave this line here blank for normal patching of plumed
#CAMSHIFT_LIBS= 
CAMSHIFT_LIBS="-L$almost2dir/src/lib -lAlm -L$almost2dir/lib/sqlite-3.6.23.1/build/lib -lsqlite3 -lz -lbz2 -L$almost2dir/src/forcefield -lnbimpl"
CAMSHIFT_INCL="-I$almost2dir/include/almost -I$almost2dir/lib/sqlite-3.6.23.1"

function to_do_before_patch () {
  echo > /dev/null
}

function to_do_after_patch () {
  {
    echo "set(PLUMED_SRC"
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo "	${f%.c}.c"
    done
    if [ "$RECON_LIBS" != "" ] ; then
      for file in $plumedir/recon_src/*.cpp
        do f=${file##*/}
        echo " ${f%.cpp}.cpp"
      done
    fi
    if [ "$CAMSHIFT_LIBS" != "" ] ; then
      for file in $plumedir/camshift_src/*.cpp
        do f=${file##*/}
        echo " ${f%.cpp}.cpp"
      done
    fi
    echo ")"
  } > ./src/kernel/plumed.inc
  if [ "$RECON_LIBS" = "" ] ; then
    { 
      echo "set(RECON_FLAGS )" 
    } > ./src/kernel/recon_patch.inc
  else
    {
      echo "set(RECON_FLAGS -DRECONMETAD)" 
      echo "add_definitions(\${RECON_FLAGS})"
    } > ./src/kernel/recon_patch.inc
  fi
  echo "set(CAMSHIFT_LIBS \"$CAMSHIFT_LIBS\")" > ./src/kernel/camshift_patch.inc
  if [ "$CAMSHIFT_LIBS" = "" ] ; then
    { 
      echo "" 
    } >> ./src/kernel/camshift_patch.inc
  else
    {
      echo "set(CAMSHIFT_FLAGS -DHAVE_ALMOST)" 
      echo "set(CAMSHIFT_INCL \"$CAMSHIFT_INCL\")"
      echo "add_definitions(\${CAMSHIFT_FLAGS})"
      echo "add_definitions(\${CAMSHIFT_INCL})"
    } >> ./src/kernel/camshift_patch.inc
  fi
}

function to_do_before_revert () {
  rm ./src/kernel/camshift_patch.inc
  rm ./src/kernel/recon_patch.inc
  rm ./src/kernel/plumed.inc
}

function to_do_after_revert () {
  echo > /dev/null
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

