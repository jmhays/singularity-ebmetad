#!/bin/bash
# PATCH SCRIPT FOR OPEP 
#

# this script needs to be launched from the root directory of the host code
destination="$PWD"

# definitions specific to this code
CODE="OPEP"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
WHERE_LINKS="./"


function to_do_before_patch () {
  echo > /dev/null
}

function to_do_after_patch () {
  echo > /dev/null 
}

function to_do_before_revert () {
  echo > /dev/null 
}

function to_do_after_revert () {
  echo > /dev/null
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

