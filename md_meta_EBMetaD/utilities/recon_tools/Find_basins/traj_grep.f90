program traj_grep
implicit none
character(len=80) :: plumedfile, trajinp, string
real*8 :: lenunits
integer :: ftype,natoms,imcon
integer :: i,parse_pdb

plumedfile="plumed.dat"
ftype=0
lenunits=1.0

if (iargc().eq.0) then
   write(6,*)"USAGE:"
   write(6,*)" -plumed    plumed input filename. default plumed.dat"
   write(6,*)" -ihistory  read dlpoly history file"
   write(6,*)" -ipdb      read pdb file"
   write(6,*)" -nm        in plumed.dat lengths are in nm because of gromacs"
   stop
end if

do i=1,iargc()
   call getarg(i, string)
   if(INDEX(string,'-ihistory').ne.0) then
      ftype=1
      call getarg(i+1,trajinp)
   end if
   if(INDEX(string,'-ipdb').ne.0) then
      ftype=2
      call getarg(i+1,trajinp)
   end if
   if(INDEX(string,'-plumed').ne.0) call getarg(i+1,plumedfile)
   if(INDEX(string,'-nm').ne.0) lenunits=10
end do

if (ftype.eq.0) then
   write(0,*)"You must enter a trajectory file"
   stop
else if(ftype.eq.1) then
   open(10,file=trajinp,status='old')
   call parse_history(10, natoms, imcon)
else if(ftype.eq.2) then
   open(10,file=trajinp,status='old')
   natoms=parse_pdb(10)
   write(6,*)"Number of atoms in file",natoms
   imcon=0
end if

call search_traj( ftype, natoms, imcon, lenunits, 10, trim(plumedfile)//char(0) ) 

end program traj_grep
! ============================================================================
integer function parse_pdb(ifile)
implicit none
integer,intent(in) :: ifile
character(len=80) :: line
character(len=8) :: str
integer :: ios

ios=0
parse_pdb=0
do while(ios==0)
   read(ifile,*,iostat=ios)line
   if (ios.ne.0) then
      write(0,*)"!ERROR - problem in reading pdb file"
      stop
   end if
   read(line,*)str 
   if(index(str,"ATOM").ne.0) parse_pdb = parse_pdb + 1
   if(index(str,"HETATM").ne.0) parse_pdb = parse_pdb + 1
   if(index(str,"END").ne.0) exit
   if(index(str,"ENDMDL").ne.0) exit
end do
rewind(ifile)

end function parse_pdb
! ============================================================================
subroutine parse_history(ifile, natoms, imcon)
implicit none
integer,intent(in) :: ifile
integer, intent(out) :: imcon,natoms
integer :: htype

read(ifile,*)
read(ifile,*)htype,imcon,natoms

end subroutine parse_history
! ==============================================================================
subroutine readhistoryframe(ifile,natoms,element,box,pos,ios)
implicit none
integer,intent(in) :: ifile,natoms
character(len=5),dimension(1:natoms),intent(out) :: element
real*8,dimension(1:3*natoms),intent(out) :: pos
real*8,dimension(1:3),intent(out) :: box
real*8,dimension(1:9) :: cell
integer,intent(out) :: ios
character(len=8) :: timestep
integer :: imcon,htype,step,nat
integer :: j


ios=0
read(ifile,*,iostat=ios) timestep,step,nat,htype,imcon
if (ios.ne.0) return

cell(:)=0
if(imcon.gt.0) then
   read(ifile,*)cell(1),cell(2),cell(3)
   read(ifile,*)cell(4),cell(5),cell(6)
   read(ifile,*)cell(7),cell(8),cell(9)
   box(1)=cell(1)
   box(2)=cell(5)
   box(3)=cell(9)
else
   box(1)=0.0
   box(2)=0.0
   box(3)=0.0
end if

do j=1,natoms
   read(ifile,*)element(j)
   read(ifile,*)pos(j),pos(j+natoms),pos(j+2*natoms)
   if(htype.gt.0) read(ifile,*)
   if(htype.gt.1) read(ifile,*)
end do

end subroutine readhistoryframe
! ==============================================================================
subroutine readpdbframe(ifile,natoms,element,resname,resno,box,pos,lenunits,ios)
implicit none
integer,intent(in) :: ifile,natoms
real*8,intent(in) :: lenunits
character(len=5),dimension(1:natoms),intent(out) :: element
character(len=4),dimension(1:natoms),intent(out) :: resname
integer,dimension(1:natoms),intent(out) :: resno
real*8,dimension(1:3*natoms),intent(out) :: pos
real*8,dimension(1:3),intent(out) :: box
integer,intent(out) :: ios
real*8 :: occu,charge
character(len=1) :: ics, altloc, icode
character(len=6) :: lname
character(len=80) :: line
integer :: num, kk, j

11   format(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,f8.3,f8.3,f8.3,f6.2,f6.2)

! This is a cheat - maybe I will sort this out one day
box(1)=0.0
box(2)=0.0
box(3)=0.0

kk=0
ios=0
do while(ios==0)
   read (ifile,'(a80)',iostat=ios) line
   if(ios/=0) exit

   if(index(line,"END").ne.0) exit
   if(index(line,"ENDMDL").ne.0) exit

   if(index(line,"ATOM").ne.0) then
      kk=kk+1
      read (line,11)lname,num,element(kk),altloc,resname(kk),ics,resno(kk),icode,pos(kk),pos(kk+natoms),pos(kk+2*natoms),occu,charge
   endif
enddo

if (ios.eq.0) then
   do while(ios==0)
      read (ifile,'(a80)',iostat=ios) line
      if(ios/=0) exit
      if(index(line,"ATOM").ne.0) exit  
   end do
end if

if (ios.ne.0) return

! Go back one in ifile so that we are on the first atom again
backspace(ifile)

if (kk.ne.natoms) then
   write(0,*)"ERROR - problems reading pdb file"
   stop
end if

! Convert pos to appropriate length units
if (lenunits.ne.1.0) then
   do j=1,3*natoms
      pos(j)=pos(j)/lenunits
   end do
end if 

end subroutine readpdbframe
