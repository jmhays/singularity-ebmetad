program generate_torsions
implicit none
character( len=80 ) :: ifilename, string
integer :: i, ftype, lev

lev=-1
ftype=0

if( iargc().eq.0 ) then
  write(6,*)"USAGE:"
  write(6,*)"  -all       generate all backbone torsions + sidechain torsions"
  write(6,*)"  -backbone  generate only the backbone torsions"
  write(6,*)"  -chi1      generate the backbone torsions + all chi1 angles"
  write(6,*)"  -igro      read in the gromacs input named"
  stop
end if

do i=1,iargc()
   call getarg(i, string)
   if(INDEX(string,"-backbone").ne.0) lev=0
   if(INDEX(string,"-chi1").ne.0) lev=1
   if(INDEX(string,"-all").ne.0) lev=2
   if(INDEX(string,"-igro").ne.0) then
      ftype=1
      call getarg(i+1,ifilename)
   endif
end do

if (lev.lt.0) then
  write(0,*)"ERROR - no information on what torsions you would like use -backbone/-chi1/-all"
  stop
end if

if(ftype.eq.0) then
  write(0,*)"ERROR - no input file specified use -gro"
  stop
end if

call gen_torsions( lev, ftype, trim(ifilename)//char(0) )

end program generate_torsions
