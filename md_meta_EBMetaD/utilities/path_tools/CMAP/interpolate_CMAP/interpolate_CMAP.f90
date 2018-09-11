 program interpolate_CMAP 
 implicit none

 integer :: i,j,k,kk
 integer :: argcount,IARGC
 integer :: iostat,npts,ncon,nframes
 integer, allocatable :: map_par(:,:)
 real*8  :: cat_init(4,3)
 real*8  :: cat_end(4,3)
 real*8  :: cat_rom(4,4)
 real*8  :: tau,integral
 real*8  :: val,u,du,tmp2
 real*8, allocatable  :: cmatrix(:,:),cmap(:,:)
 character*80 :: line,w_char,cmapvalue,outfile
 character*80 :: complex_name,numb_to_char
 character*15 :: str

 logical :: write_map

13   FORMAT(i4,1x,i4,1x,i4,4x,f12.8)
14   FORMAT(a3)

 integral  = 0.d0
 tau       = 0.5d0
 npts      = 1000
 outfile   = "CMAP.interpol"
 cmapvalue = ""
 write_map = .false.

 write(*,*)
 write(*,*) " --> CMAP_TOOLS :: interpolate_CMAP <--"
 write(*,*)

 argcount=IARGC()
 if(argcount.eq.0)then
   write(*,*)"USAGE :"
   write(*,*)"./interpolate_CMAP -cmap CMAPVALUES (-npts 1000 -tau 0.5 -out interpol)"
   write(*,*) 
   write(*,*)"-cmap     : name of CMAP file with ALL the frames that define the path concatenated"
   write(*,*)"-npts     : number of interpolated points           (default=1000)" 
   write(*,*)"-tau      : parameter tau for catmull interpolation (default=0.5)" 
   write(*,*)"-out      : interpolated CMAP output file           (default=CMAP.interpol)"
   write(*,*)"-writemap : write single-frame CMAP files for the original CMAP path" 
   write(*,*)
  stop
 endif
 do i=1,argcount
   call getarg(i,w_char) 
   if (INDEX(w_char,'-cmap').NE.0)then
     call getarg(i+1,w_char)
     read(w_char,*)cmapvalue
   endif
   if (INDEX(w_char,'-npts').NE.0)then
     call getarg(i+1,w_char)
     read(w_char,*)npts
   endif
   if (INDEX(w_char,'-tau').NE.0)then
     call getarg(i+1,w_char)
     read(w_char,*)tau
   endif
   if (INDEX(w_char,'-out').NE.0)then
     call getarg(i+1,w_char)
     read(w_char,*)outfile
   endif
   if (INDEX(w_char,'-writemap').NE.0)then
     call getarg(i+1,w_char)
     write_map=.true. 
   endif
 enddo

 if(cmapvalue=="") then
  write(*,*) "ERROR! You need to specify the name of the CMAP file with ALL the frames that define the path concatenated"
  write(*,*)
  stop
 endif

 write(*,*) "CMAP input file               :: ", cmapvalue
 write(*,*) "Output traj file              :: ",outfile
 write(*,*) "Number of interpolated frames :: ",npts
 write(*,*) "Tau for Catmull-Rom           :: ",tau
 if(write_map) write(*,*) "** Writing single-frame CMAP files for the original CMAP path (useful for choose_best_CMAP!)"
!! allocating stuff
 open(200,file=cmapvalue,form="formatted")
 open(300,file=outfile,form="formatted") 
 ncon=0
 nframes=0

 do
  read(200,'(a80)',iostat=iostat) line
  if(iostat/=0) exit
  read(line,*)str
  if(index(str,"END").ne.0) then 
    nframes=nframes+1
    cycle
  endif
  if(nframes==1) ncon=ncon+1
 enddo

 write(*,*) "Number of contacts            :: ",ncon
 write(*,*) "Number of frames              :: ",nframes

 allocate(cmap(nframes,ncon))
 allocate(map_par(2,ncon))
 allocate(cmatrix(4,ncon))

!! rewind file, reading again and printing out individual MAP file!!

 rewind(200)
 
 do i=1,nframes
 
  if(write_map) then
   write(numb_to_char,'(i15)')i
   numb_to_char=(adjustl(numb_to_char))
   complex_name="MAP."//trim(numb_to_char)
   open(100,file=complex_name,form="formatted")
  endif 

   do j=1,ncon
    read(200,13) k,map_par(1,j),map_par(2,j),cmap(i,j)
    if(write_map) write(100,13) k,map_par(1,j),map_par(2,j),cmap(i,j)
   enddo
   read(200,14)
   if(write_map) then 
    write(100,14) "END" 
    close(100)
   endif
 
 enddo
 
 close(200)

!! catmull rom matrixes
 cat_init(1,1) = 1.d0
 cat_init(1,2) = 0.d0
 cat_init(1,3) = 0.d0
 
 cat_init(2,1) = -tau 
 cat_init(2,2) =  tau 
 cat_init(2,3) = 0.d0  

 cat_init(3,1) = 3.d0*tau-3d0
 cat_init(3,2) = 3.d0-2.d0*tau
 cat_init(3,3) = -tau

 cat_init(4,1) = 2.d0-2.d0*tau
 cat_init(4,2) = tau-2.d0
 cat_init(4,3) = tau

 cat_end(1,1) = 0.d0 
 cat_end(1,2) = 1.d0 
 cat_end(1,3) = 0.d0
 
 cat_end(2,1) = -tau
 cat_end(2,2) = 0.d0
 cat_end(2,3) =  tau

 cat_end(3,1) = 2.d0*tau
 cat_end(3,2) = tau-3.d0
 cat_end(3,3) = 3.d0-3.d0*tau

 cat_end(4,1) = -tau
 cat_end(4,2) = 2.d0-tau
 cat_end(4,3) = 2*tau-2

 cat_rom(1,1) = 0.d0 
 cat_rom(1,2) = 1.d0
 cat_rom(1,3) = 0.d0 
 cat_rom(1,4) = 0.d0 

 cat_rom(2,1) = -tau
 cat_rom(2,2) = 0.d0
 cat_rom(2,3) = tau 
 cat_rom(2,4) = 0.d0

 cat_rom(3,1) = 2.d0*tau
 cat_rom(3,2) = tau-3.d0
 cat_rom(3,3) = 3.d0-2.d0*tau
 cat_rom(3,4) = -tau

 cat_rom(4,1) = -tau
 cat_rom(4,2) = 2.d0-tau
 cat_rom(4,3) = tau-2.d0
 cat_rom(4,4) =  tau

!! interpolation
 cmatrix=0.d0
 do k=1,4 ! u components 
   do j=1,3   !nnodes 
     do i=1,ncon ! ncon
      cmatrix(k,i)=cmatrix(k,i)+cat_init(k,j)*cmap(j,i) 
     enddo
   enddo
 enddo

 do i=0,int((npts-1)/(nframes-1))-1
     u=real(i)/(real(npts-1.d0)/real(nframes-1.d0))
     du=1.d0/(real(npts-1.d0)/real(nframes-1.d0))
     do j=1,ncon
       val=cmatrix(1,j)+cmatrix(2,j)*u+cmatrix(3,j)*u*u+cmatrix(4,j)*u*u*u
       tmp2=(cmatrix(2,j)+2*cmatrix(3,j)*u+3.d0*cmatrix(4,j)*u*u )**2.d0
       integral=integral+sqrt(tmp2)*du
       write(300,13) j,map_par(1,j),map_par(2,j),val
     enddo 
     write(300,'("END")')
 enddo 

! intermediate interpolation:loop
 do kk=1,nframes-3
   cmatrix=0.d0
   do k=1,4 ! u components 
     do j=1,4   !nnodes 
       do i=1,ncon ! ncon
        cmatrix(k,i)=cmatrix(k,i)+cat_rom(k,j)*cmap(j+kk-1,i) 
       enddo
     enddo
   enddo

   do i=0,int((npts-1)/(nframes-1))-1
       u=real(i)/(real(npts-1.d0)/real(nframes-1.d0))
       du=1.d0/(real(npts-1.d0)/real(nframes-1.d0))
       do j=1,ncon
           val=cmatrix(1,j)+cmatrix(2,j)*u+cmatrix(3,j)*u*u+cmatrix(4,j)*u*u*u
            tmp2=(cmatrix(2,j)+2*cmatrix(3,j)*u+3.d0*cmatrix(4,j)*u*u )**2.d0 
            integral=integral+sqrt(tmp2)*du
       write(300,13) j,map_par(1,j),map_par(2,j),val
       enddo 
       write(300,'("END")')
   enddo 
 enddo 

! Final step
 cmatrix=0.d0
 do k=1,4 ! u components 
   do j=1,3   !nnodes 
     do i=1,ncon ! ncon
      cmatrix(k,i)=cmatrix(k,i)+cat_end(k,j)*cmap(nframes-3+j,i) 
     enddo
   enddo
 enddo

 do i=0,int((npts-1)/(nframes-1))-1
     u=real(i)/(real(npts-1.d0)/real(nframes-1.d0))
     du=1.d0/(real(npts-1.d0)/real(nframes-1.d0))
     do j=1,ncon
       tmp2=0.d0  
       val=cmatrix(1,j)+cmatrix(2,j)*u+cmatrix(3,j)*u*u+cmatrix(4,j)*u*u*u
       tmp2=(cmatrix(2,j)+2*cmatrix(3,j)*u+3.d0*cmatrix(4,j)*u*u )**2.d0
       integral=integral+sqrt(tmp2)*du
       write(300,13) j,map_par(1,j),map_par(2,j),val
     enddo
     write(300,'("END")') 
 enddo
 write(*,*) 
 write(*,'(1x,a14,1x,f10.4,1x,a8)') "PATH LENGTH ::",integral,"CONTACTS"
 write(*,*)

 end 
