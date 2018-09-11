 PROGRAM matrix_metrix 

 IMPLICIT none

 integer      :: i,ii,j,jj,l,ll,m
 integer      :: iostat,nframes
 integer      :: argcount,IARGC
 integer      :: tot_contact 
 character*8  :: str
 character*45 :: matrixfile,inputfile 
 character*80 :: wq_char,line
 logical      :: go
 real*8, allocatable :: cmap(:,:)
 real*8       :: err, mean, mean2, stddev


11   FORMAT(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,f8.3,f8.3,f8.3,f6.2,f6.2)
13   FORMAT(i4,1x,i4,1x,i4,4x,f12.8)
14   FORMAT(a3)

 matrixfile="METRIX_MATRIX" 
 argcount = IARGC()
 go=.false.

 write(*,*)
 write(*,*) " --> CMAP_TOOLS :: metrix_matrix <--"
 write(*,*)

 do i=1,argcount
   call getarg(i,wq_char)
   if (INDEX(wq_char,'-cmap').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)inputfile
     go=.true.
   endif
   if (INDEX(wq_char,'-out').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)matrixfile
   endif
 enddo

 if(.not.(go)) then
   write(*,*)"USAGE : "
   write(*,*)"./metrix_matrix -cmap CMAPVALUES (-out METRIX_MATRIX.gplt)" 
   write(*,*) 
   write(*,*)"-cmap        : input CMAP file" 
   write(*,*)"-out         : output metric matrix filename (optional)"
   write(*,*)
   stop
 endif

 open (UNIT=20,file=matrixfile,status="unknown")
 open (UNIT=10,FILE=inputfile,FORM="FORMATTED")


!!  Count number of FRAMES and CONTACTS 

 tot_contact  = 0
 nframes      = 0
 do
  read(10,'(a80)',iostat=iostat) line
  if(iostat/=0) exit
  read(line,*)str
  if(index(str,"END").ne.0) then
   nframes=nframes+1
  else
   if(nframes==0) tot_contact=tot_contact+1
  endif 
 enddo

 write(*,*) "CMAP   file        : ",inputfile
 write(*,*) "METRIX file        : ",matrixfile
 write(*,*) "Number of Contacts : ",tot_contact
 write(*,*) "Number of Frames   : ",nframes
 rewind(10)

 allocate(cmap(nframes,tot_contact))

 do i=1,nframes
  do ii=1,tot_contact
   read(10,13) ll,l,m,cmap(i,ii)
  enddo       
  read(10,14) str
 enddo 

 mean=0.d0
 mean2=0.d0

!! evaluate matrix and print out
 do i=1,nframes
  do j=1,nframes
   err=sum((cmap(i,:)-cmap(j,:))**2)
   write(20,*) i,j,err
   if((j-i).eq.1) then
    mean=mean+err
    mean2=mean2+err**2
   endif
  enddo
  write(20,*) 
 enddo

 mean2  = mean2/dble(nframes-1)
 mean   = mean/dble(nframes-1)
 stddev = sqrt((mean2-mean**2)/dble(nframes-1))

 write(*,*)
 write(*,'(a29,f11.6)') " Average Distance i, i+1  :: ", mean
 write(*,'(a29,f11.6)') " Standard Deviation       :: ", stddev 
 write(*,'(a29,f11.6)') " Percentage               :: ", stddev/mean*dble(100)
 write(*,'(a29,f11.6)') "** SUGGESTED LAMBDA       :: ", 2.3d0/mean
 write(*,*) 
 end
