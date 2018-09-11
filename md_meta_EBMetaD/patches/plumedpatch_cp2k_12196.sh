#!/bin/bash
# PATCH SCRIPT FOR CP2K 
#

# this script needs to be launched from the root directory of the host code
# to be run in cp2k directory
destination="$PWD"

# definitions specific to this code
CODE="cp2k_12196"
LINKED_FILES="$plumedir/common_files/*.h $plumedir/common_files/*.c"
RECON_LINKS="$plumedir/recon_src/*.cpp $plumedir/recon_src/*.h"
WHERE_LINKS="./src/plumed"
# Leave this line here blank for normal patching of plumed
RECON_LIBS= 
# CVS version only : for reconnaissance metadynamics enter here information on the
# locations of the lstdc++ libraries.
#RECON_LIBS="-lstdc++"
# We must define what c++ compiler we are using to compile the reconnaissance routines
RECON_CPP="mpicxx"

function to_do_before_patch () {
#  echo > /dev/null
# here put the interface files
echo "CREATING ${destination}/src/plumed_methods.F"
cat >${destination}/src/plumed_methods.F<<EOF  
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief The plumed wrapper 
!> \param the plumed wrapper 
!> \author DB 6/11/2011 
! *****************************************************************************
MODULE plumed_methods
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type,&
                                              use_mixed_force
   USE cell_types,                      ONLY: cell_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_p_type,&
                                              cp_subsys_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
   USE simpar_types,                    ONLY: simpar_type
   USE particle_list_types,             ONLY: particle_list_p_type,&
                                             particle_list_type
   USE particle_types,                  ONLY: particle_type
   ! beware: boltzman is in kjoule : needs a fact of 4.184_dp
   USE physcon,                         ONLY: boltzmann,&
                                              femtoseconds,&
                                              kcalmol,&
                                              angstrom,&
                                              bohr,&
                                              kelvin,& 
                                              massunit
   USE plumed_types,                    ONLY: plumed_type

   USE kinds,                           ONLY: dp

#include "cp_common_uses.h"
  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'plumed_methods'

 CONTAINS  
   SUBROUTINE plumed_wrapper(force_env,myplumed,dt,nsteps,fe_section,error)
    CHARACTER(len=*), PARAMETER :: routineN = 'plumed_wrapper', &
     routineP = moduleN//':'//routineN
    TYPE(plumed_type), POINTER               :: myplumed
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simpar_type), POINTER               :: simpar
    CHARACTER(LEN=50)                        :: plumed_file
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER                                  :: natom,iw,i,stat,nsteps,j 
    TYPE(section_vals_type), POINTER         :: fe_section
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
    POINTER                                :: my_particles
    LOGICAL                                  :: failure
    REAL (KIND=dp) :: dt ,ene 
    failure=.FALSE.
 
    NULLIFY(subsys,cell,logger,particles_i,my_particles)
    logger => cp_error_get_logger(error)

 
    CALL force_env_get(force_env, subsys=subsys, cell=cell, error=error)
    ! a subsys contains what I need
    CALL cp_subsys_get(subsys,natom=natom,particles=particles_i,error=error)
    my_particles => particles_i%els
    ! cell                 has cell : see cell_types.F   
    ! timestep
    ! my_particles(i)%f(:) has force
    ! my_particles(i)%r(:) has position   in bohrs
    ! my_particles(i)%v(:) has velocity    
    ! my_particles(i)%atomic_kind%mass has mass
 
    myplumed%istep=myplumed%istep+1
    ! now try a dump
    if (myplumed%istep.eq.1) then
       NULLiFY(myplumed%x,myplumed%v,myplumed%f,myplumed%mass) 
       ALLOCATE(myplumed%x(3*natom),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ALLOCATE(myplumed%v(3*natom),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ALLOCATE(myplumed%f(3*natom),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ALLOCATE(myplumed%mass(natom),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ALLOCATE(myplumed%charge(natom),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    endif
    if(myplumed%istep.eq.1)then
      do i=1,natom
        myplumed%mass(i)=my_particles(i)%atomic_kind%mass/massunit
      enddo
    endif
    do i=1,natom
      myplumed%x(3*(i-1)+1)=my_particles(i)%r(1)*angstrom
      myplumed%x(3*(i-1)+2)=my_particles(i)%r(2)*angstrom
      myplumed%x(3*(i-1)+3)=my_particles(i)%r(3)*angstrom
    enddo
    myplumed%f=0.
    myplumed%cell(1)=cell%hmat(1,1) 
    myplumed%cell(2)=cell%hmat(2,2) 
    myplumed%cell(3)=cell%hmat(3,3) 
    ! the timestep ?
    iw = cp_print_key_unit_nr(logger,fe_section,"FREE_ENERGY_INFO",&
           extension=".free_energy",error=error)
!    IF (iw>0) THEN
!          WRITE(iw,'("PLUMED: FILE IS  ",A)'),myplumed%plumed_file  
!          WRITE(iw,'("PLUMED: NATOMS ",I5)'),natom  
!          WRITE(iw,'("PLUMED: ISTEP ",I5)'),myplumed%istep
!          WRITE(iw,'("PLUMED: NSTEP ",I10)'),nsteps
!          WRITE(iw,'("PLUMED: DT ",F12.6)'),dt*femtoseconds
!          if (cell%orthorhombic) then
!          WRITE(iw,'("PLUMED: ORTHO ")')
!          endif
!          WRITE(iw,'("PLUMED: HMAT IS  ",5F12.6)'),cell%hmat(1,1),cell%hmat(1,2),cell%hmat(1,3)
!          WRITE(iw,'("PLUMED: HMAT IS  ",5F12.6)'),cell%hmat(2,1),cell%hmat(2,2),cell%hmat(2,3)
!          WRITE(iw,'("PLUMED: HMAT IS  ",5F12.6)'),cell%hmat(3,1),cell%hmat(3,2),cell%hmat(3,3)
!          do i=1,natom
!             WRITE(iw,'("PLUMED: ",I5," POS ",3F8.3," VEL ",3F8.3," FORCE ",3F8.3," MASS ",F8.3)'),i,&   
!                      myplumed%x(3*(i-1)+1), myplumed%x(3*(i-1)+2), myplumed%x(3*(i-1)+3),&
!                      my_particles(i)%v(1),my_particles(i)%v(2),my_particles(i)%v(3),&
!                      my_particles(i)%f(1),my_particles(i)%f(2),my_particles(i)%f(3),&
!                      myplumed%mass(i)
!             !WRITE(iw,'("ATOM  ",I5,"  C   ALA     2    ",3F8.3,"  1.00  1.00")')i,&   
!             !         myplumed%x(3*(i-1)+1), myplumed%x(3*(i-1)+2), myplumed%x(3*(i-1)+3)
!          enddo
!    ENDIF
    if (myplumed%istep.eq.1)then
       call init_metadyn(natom, dt*femtoseconds ,myplumed%mass, & 
        myplumed%charge,1,1,trim(adjustl(myplumed%plumed_file))//char(0)) 
    endif
 
 
    call meta_force_calculation(myplumed%cell, myplumed%istep, myplumed%x,  myplumed%x,   myplumed%x, &
                                                               myplumed%f,  myplumed%f,   myplumed%f, ene)

!    do i=1,natom
!        WRITE(iw,'("PLUMED: ",I5," FORCE ",3F8.3)'),i,&   
!                 myplumed%f(3*(i-1)+1), myplumed%f(3*(i-1)+2), myplumed%f(3*(i-1)+3)
!    enddo


    ! convert forces into internal units .....still to work on it 
    !  f== e(kcal/mol)*(au/kcalmol)/x(ang)(bohr/ang)=  
    !  efact=au/kcalmol  = 1/kcalmol
    !  rfact=bohr/ang =  bohr
    !myplumed%f=myplumed%f/(kcalmol*bohr)
    myplumed%f=myplumed%f/(kcalmol)

    ! f=m*a=
    !myplumed%f=myplumed%f*bohr/kcalmol

    ! put derivatives 
    do i=1,natom
      particles_i%els(i)%f(1)=particles_i%els(i)%f(1)+myplumed%f(3*(i-1)+1)
      particles_i%els(i)%f(2)=particles_i%els(i)%f(2)+myplumed%f(3*(i-1)+2)
      particles_i%els(i)%f(3)=particles_i%els(i)%f(3)+myplumed%f(3*(i-1)+3)
    enddo

   END SUBROUTINE plumed_wrapper
END MODULE plumed_methods
EOF
echo "CREATING ${destination}/src/plumed_types.F"
cat >${destination}/src/plumed_types.F<<EOF  
MODULE plumed_types

  USE kinds,                           ONLY: dp
 
  ! a container type for passing the position and other stuff  
 PUBLIC :: plumed_type

 TYPE plumed_type
    INTEGER                                      :: mysize,istep
    CHARACTER(LEN=50)                            :: plumed_file
    REAL(KIND=dp), DIMENSION(:), POINTER         :: x
    REAL(KIND=dp), DIMENSION(:), POINTER         :: v
    REAL(KIND=dp), DIMENSION(:), POINTER         :: f
    REAL(KIND=dp), DIMENSION(:), POINTER         :: mass
    REAL(KIND=dp), DIMENSION(:), POINTER         :: charge 
    REAL(KIND=dp), DIMENSION(3)                  :: cell 
 END TYPE plumed_type 
 
END MODULE plumed_types
EOF
cat >${destination}/src/PLUMEDDEFS<<EOF
OBJECTS_PLUMED=\
 plumed_types.o\
 plumed_methods.o
EOF
if [ ! -e ${destination}/src/plumed ]; then
	mkdir ${destination}/src/plumed
fi

}

function to_do_after_patch () {
  {
    echo -n "PLUMED_OBJECTS="
    for file in $plumedir/common_files/*.c
      do f=${file##*/}
      echo " \\"
      echo -n "	${f%.c}.o"
    done
    if [ "$RECON_LIBS" != "" ] ; then
      for file in $plumedir/recon_src/*.cpp
        do f=${file##*/}
        echo " \\"
        echo -n " ${f%.cpp}.\$(OBJEXT)"
      done
    fi
    echo
  } > ./src/plumed/plumed.inc
  echo "RECON_LIBS=" $RECON_LIBS > ./src/recon_patch.inc
  if [ "$RECON_LIBS" = "" ] ; then
    { 
      echo "RECON_FLAGS=" 
      echo "RECON_CPP=\$(CC)"
    } >> ./src/plumed/recon_patch.inc
  else
    {
      echo "RECON_FLAGS=-DRECONMETAD" 
      echo "RECON_CPP="$RECON_CPP
    } >> ./src/plumed/recon_patch.inc
  fi
}

function to_do_before_revert () {
  rm ./src/plumed/recon_patch.inc
  rm ./src/plumed/plumed.inc
}

function to_do_after_revert () {
  if [ -e "${destination}/src/plumed" ]; then 
	rm -rf ${destination}/src/plumed
  fi 
}

#########

NAME="$0"
source $plumedir/patches/patch_tool.sh

