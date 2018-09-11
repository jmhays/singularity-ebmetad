program network 
implicit none
character(len=80) :: plumedfile, cvfile, string
integer :: i,ineighbours,outdim,eps

plumedfile="plumed.dat"
cvfile="no_cvs"
ineighbours=0
outdim=2
eps=0

if (iargc().eq.0) then
   write(6,*)"USAGE:"
   write(6,*)" -neigh     number of neighbours to use in making connections in basins"
   write(6,*)" -expo      exponent for calculation of embedded coordinates"
   write(6,*)" -colvars   file containing colvars."
   write(6,*)" -plumed    OPT. plumed input filename. default plumed.dat"
   write(6,*)" -dim       OPT. dimensionality we are reducing to. default 2"
   stop
end if


do i=1,iargc()
   call getarg(i, string)
   if(INDEX(string,'-neigh').ne.0) then
     call getarg(i+1,string)
     read(string,*)ineighbours
   end if
   if(INDEX(string,'-expo').ne.0) then
     call getarg(i+1,string)
     read(string,*)eps
   end if
   if(INDEX(string,'-dim').ne.0) then
     call getarg(i+1,string)
     read(string,*)outdim
   end if
   if(INDEX(string,'-plumed').ne.0) call getarg(i+1,plumedfile)
   if(INDEX(string,'-colvars').ne.0) call getarg(i+1,cvfile)
end do

if(cvfile.eq."no_cvs") then
  write(0,*)"ERROR - you have not specified the name of the file containing the colvars"
  stop
else if(ineighbours.eq.0) then
   write(0,*)"ERROR - you have not specified a number of neighbours - use -neigh"
   stop
else if(eps.eq.0) then
   write(0,*)"ERROR - you have not specified the power for embedding of cvs - use -expo"
   stop
end if

call isomap( ineighbours, eps, outdim, trim(plumedfile)//char(0), trim(cvfile)//char(0) ) 

end program network 
