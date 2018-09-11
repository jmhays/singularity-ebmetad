 program choose_best_CMAP 

 implicit none

 integer :: tot_atom,iostat,iter
 integer :: found,minima,flag
 integer :: tot_con,tot_acon,tot_gcon
 integer :: iframe, cell_logic, tot_frames
 integer :: i,ii,iii,j,jj,k,kk,l,ll,m,mm
 integer :: argcount,IARGC,nframes
 integer :: indice
 integer, allocatable :: index_frame(:)
 integer :: idum = -123 
 real*8, external    :: ran1
 real*8  :: dist,minerr,beta,fakex,fakey,fakez
 real*8  :: err1,errn,pippo,pippo1,pippo2,pippo3
 
 character*80 :: cmapgroup,cmapvalue,line,cmapindex
 character*7  :: contact
 character*80 :: complex_name,numb_to_char,pdb,dcd,traj,wq_char
 character*15 :: str

!! PDB stuff
 character*6, allocatable ::   recordname(:)
 character*4, allocatable ::   atom(:),atmname(:)
 character*1, allocatable ::   altloc(:),chainid(:),icode(:)
 character*3, allocatable ::   resname(:)
 integer, allocatable     ::   atmnum(:),resnum(:)
 real*8, allocatable      ::   occ(:),temp(:)

!! CMAP variables
 real*8, allocatable ::  cmap_coord1(:),cmap_tmp(:),cmap_frame(:,:)
 real*8, allocatable ::  cmap_traj(:),cmap_coord2(:),cmap_run(:)
 integer, pointer    ::  map_par1(:,:)
 real*8, pointer     ::  map_par2(:,:)
 real*4, pointer     ::  x4(:)
 real*4, pointer     ::  y4(:)
 real*4, pointer     ::  z4(:)

 integer             :: ngroup_max, max_atom_group
 parameter(ngroup_max=10,max_atom_group=50)
!! group variables
 integer             :: ngroup,g_number(ngroup_max),g_index(ngroup_max,max_atom_group)
 real*8              :: rcm(3,ngroup_max)

!! DCD stuff
 character*4  :: car4
 integer      :: ntap,nstart,nsanc,nset,ntitle,namin,i5(5),i9(9)
 character*50 :: filein
 character*80 :: car(10)
 real*8       :: delta,cell(6)

!! MC stuff 
 real*8,pointer :: err_table(:,:)
 real*8         :: err
 real*8         :: upot_new,upot,upot_old
 real*8         :: offset
 real*8         :: errmin,xp,redux

!! logical
 logical              :: dcd_on,traj_on,firstneigh 
 logical, allocatable :: done(:)
 logical              :: exist,go1,go2

12   format(a7,1x,i4,1x,i4,1x,i4,1x,f5.1,1x,i2,1x,i2,1x,f5.1,1x,f8.5) 
13   format(i4,1x,i4,1x,i4,4x,f12.8)
14   format(a3)


 cell_logic  = 0
 xp          = 2.d0
 offset      = 0.d0
 redux       = 1.0d0
 firstneigh  = .false.
 beta        = 10000
 go1         = .false.
 go2         = .false.
 cmapindex   = ""
 cmapgroup   = ""
 dcd_on      = .false. 
 traj_on     = .false.
 iter        = 1
 tot_con     = 0
 tot_acon    = 0
 tot_gcon    = 0
 dcd         = ""
 pdb         = ""
 traj        = ""

 write(*,*)
 write(*,*) " --> CMAP_TOOLS :: choose_best_CMAP <--"
 write(*,*) 


 argcount = IARGC()
 do i=1,argcount
   call getarg(i,wq_char)
   if (INDEX(wq_char,'-cmap').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) cmapvalue 
     go1=.true.
   endif
   if (INDEX(wq_char,'-pdb').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) pdb 
     dcd_on=.true.
   endif
   if (INDEX(wq_char,'-dcd').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) dcd 
     dcd_on=.true.
   endif
   if (INDEX(wq_char,'-traj').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) traj 
     traj_on=.true.
   endif
   if (INDEX(wq_char,'-xp').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) xp 
   endif
   if (INDEX(wq_char,'-offset').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)offset 
   endif
   if (INDEX(wq_char,'-beta').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)beta  
   endif
   if (INDEX(wq_char,'-firstneigh').NE.0)then
     firstneigh=.true.
   endif
   if (INDEX(wq_char,'-nframes').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)nframes
     go2=.true. 
   endif 
   if (INDEX(wq_char,'-iter').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)iter
   endif
   if (INDEX(wq_char,'-redux').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)redux 
   endif
   if (INDEX(wq_char,'-index').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)cmapindex
     dcd_on=.true.
   endif
   if (INDEX(wq_char,'-group').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)cmapgroup
   endif
   if (INDEX(wq_char,'-cell').NE.0)then
     cell_logic=1 
   endif
 enddo

 if(.not.(go1.and.go2))then
    write(*,*)"USAGE : "
    write(*,*)"./choose_best_CMAP -cmap MAP. -nframes 20 -index CMAPINDEX -cmapgroup CMAPGROUP"
    write(*,*)"                   -pdb protein.pdb -dcd protein.dcd -traj interpol.CMAP -iter 10 -cell -beta 100"
    write(*,*)"                   -redux 1.0 -firstneigh -xp 2.0 -offset 0.0"
    write(*,*)  
    write(*,*)"-cmap       : basename for CMAP files"
    write(*,*)"-nframes    : number of frames"
    write(*,*)"-index      : CMAP index file                            "
    write(*,*)"-group      : CMAP group file                            (optional)"
    write(*,*)"-pdb        : pdb connected to dcd                       (optional)"
    write(*,*)"-dcd        : dcd to scan                                (optional)"
    write(*,*)"-traj       : trajectory file in CMAP format             (optional)"
    write(*,*)"-iter       : number of mc iterations                    (default=10)"
    write(*,*)"-cell       : if dcd with cell information               (optional)"
    write(*,*)"-beta       : mc temperature                             (default=10000)"    
    write(*,*)"-redux      : scaling factor in the potential definition (default=1.0)"
    write(*,*)"-firstneigh : only first neighbour frames                (optional)"
    write(*,*)"-xp         : exponential in the potential definition    (default=2.0)"
    write(*,*)"-offset     : offset in the potential definition         (default=0.0)"
    write(*,*)
    stop
  endif 

!! various controls
 if(dcd_on.and.(cmapindex=="".or.pdb=="".or.dcd=="")) then
  write(*,*) "ERROR! You have chosen to extract maps from a DCD file but you are missing something:"
  write(*,*) "       the index file (-index), the pdb (-pdb) or the dcd (-dcd) ??"
  write(*,*)
  stop
 endif

 if(.not.(dcd_on.or.traj_on)) then
  write(*,*) "ERROR! You have to specify either a DCD file or a TRAJ file in CMAP format."
  write(*,*)
  stop
 endif

!! info printout
 write(*,*) "CMAP basename            :: ",cmapvalue
 write(*,*) "Number of frames         :: ",nframes

 if(dcd_on) then
  write(*,*) "** Extracting intermediate frames from a DCD trajectory."
  write(*,*) "PDB file name            :: ",pdb
  write(*,*) "DCD file name            :: ",dcd
  write(*,*) "CMAP index file          :: ",cmapindex
 else
  write(*,*) "** Extracting intermediate frames from a CMAP trajectory."
  write(*,*) "TRAJ file name           :: ",traj
 endif

 write(*,*) "Number of MC iterations :: ",iter
 write(*,*) "Beta MC temperature     :: ",beta
 write(*,*) "Redux                   :: ",redux
 write(*,*) "Offset                  :: ",offset
 write(*,*) "Exponent                :: ",xp

 if(firstneigh) write(*,*) "Using only first neighbour frames"
 if(dcd_on) then
  if(cell_logic.eq.1) then
   write(*,*) "**** WARNING!! Reading dcd WITH cell information"
  else
    write(*,*) "**** WARNING!! Reading dcd WITHOUT cell information"
  endif 
 endif
 write(*,*)

!! allocation
 allocate(done(nframes), index_frame(nframes))
 allocate(err_table(nframes,nframes))

 index_frame = 0

!!! CASE of PDB+DCD+CMAP index file

 if(dcd_on) then 

!! reading PDB for number of atoms
  open(120,file=pdb,form="formatted")
  tot_atom=0 
  do
     read(120,'(a80)',iostat=iostat)line
     if(iostat/=0)exit
     read(line,*)str
     if(index(str,"END").ne.0)exit
     if(index(str,"ATOM").ne.0) tot_atom=tot_atom+1 
  enddo
  rewind(120)

!! allocation
  allocate(x4(tot_atom))
  allocate(y4(tot_atom))
  allocate(z4(tot_atom))
  allocate(recordname(tot_atom))
  allocate(atom(tot_atom))
  allocate(atmnum(tot_atom))
  allocate(atmname(tot_atom))
  allocate(altloc(tot_atom))
  allocate(chainid(tot_atom))
  allocate(icode(tot_atom))
  allocate(resname(tot_atom))
  allocate(resnum(tot_atom))
  allocate(occ(tot_atom))
  allocate(temp(tot_atom))

!! reading PDB
  kk=0
  do
   read(120,'(a80)',iostat=iostat)line
   if(iostat/=0) exit
   read(line,*)str
   if(index(str,"END").ne.0)exit
   if(index(str,"ATOM").ne.0)then
    kk=kk+1
     read(line,'(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,f8.3,f8.3,f8.3,f6.2,f6.2)')&
      &recordname(kk),atmnum(kk),atmname(kk),&
      &altloc(kk),resname(kk),chainid(kk),&
      &resnum(kk),icode(kk),fakex,fakey,fakez,&
      &occ(kk),temp(kk)
   endif
  enddo
  close(120)
 

!! reading CMAP index file for allocation...
  open(unit=55,file=cmapindex,form="formatted")
  kk=0
  ii=0
  do
      read(55,'(a80)',iostat=iostat)line
      if(iostat/=0)exit
      read(line,*)str
      if(index(str,"CONTACT").ne.0) kk=kk+1
      if(index(str,"GROUP").ne.0) ii=ii+1
  enddo
  rewind(55)

  write(*,'(1x,a5,1x,i4,1x,a19,1x,i3,1x,a14)') "Found",kk,"atomic contacts and",ii,"group contacts" 

  tot_acon=kk
  tot_gcon=ii
  tot_con=tot_acon+tot_gcon

  allocate(map_par1(tot_con,4))
  allocate(map_par2(tot_con,3))
  allocate(cmap_tmp(tot_con))
  allocate(cmap_coord1(tot_con))
  allocate(cmap_coord2(tot_con))
  allocate(cmap_run(tot_con))
  allocate(cmap_frame(tot_con,nframes))

!! ...real reading 
  do ii=1,tot_con
   read(55,12) contact,i,map_par1(ii,1),map_par1(ii,2),map_par2(ii,1),map_par1(ii,3),&
                 &map_par1(ii,4),map_par2(ii,2),map_par2(ii,3)
  enddo
  close(55)

  if((tot_gcon.ne.0).and.(cmapgroup.eq.""))then
   write(*,*) "No CMAPGROUP file specified. Exiting...."
   stop 
  endif

!! reading CMAP group file
  if(tot_gcon.ne.0)then 
   open(3000,file=cmapgroup,status="unknown")
   i=1
   do
    read(3000,*,iostat=iostat) contact,mm,g_number(i),(g_index(i,ii),ii=1,g_number(i))
    if(iostat/=0) exit
    i=i+1
   enddo
   ngroup=i-1 
   write(*,*) "Number of group found: ",ngroup
   close(3000)
  endif !! case of GROUP file

!! reading DCD trajectory

!!   First scanning of dcd for allocation       
  write(*,*) "First scanning of DCD for allocation..."
  open (unit=11,file=dcd,form="unformatted")
  read (11)car4,nset, nstart,nsanc,i5,namin,delta,i9
  read (11)ntitle,(car(i),i=1,ntitle)
  read (11)ntap
  write(*,*)"END read DCD HEADER ",ntap," ATOMS "

  if(tot_atom.ne.ntap) then
   write(*,*) "** ERROR!! Mismatch between pdb and dcd number of atoms" 
   stop
  endif
  
  tot_frames=0
  do
   if(cell_logic.eq.1) then
     read(11,iostat=iostat)cell
     if (iostat/=0) exit
   endif
   read(11,iostat=iostat)(x4(m),m=1,tot_atom)
   if (iostat/=0) exit
   read(11,iostat=iostat)(y4(m),m=1,tot_atom)
   if (iostat/=0) exit
   read(11,iostat=iostat)(z4(m),m=1,tot_atom)
   if (iostat/=0) exit
   tot_frames=tot_frames+1 
  enddo
 
  write(*,*) "** DCD has ",tot_frames," frames"
  allocate(cmap_traj(tot_frames*tot_con))
  rewind(11)
 
!!   Second scanning: building the running contact maps array
  read (11)car4,nset, nstart,nsanc,i5,namin,delta,i9
  read (11)ntitle,(car(i),i=1,ntitle)
  read (11)ntap

  do i=1,tot_frames
   if(cell_logic.eq.1) read(11)cell
   read(11)(x4(m),m=1,tot_atom)
   read(11)(y4(m),m=1,tot_atom)
   read(11)(z4(m),m=1,tot_atom)

   call cmap_running(tot_atom,atmnum,tot_acon,tot_gcon,x4,y4,z4,cmap_run,map_par1,map_par2,ngroup,g_number,g_index)
   do j=1,tot_con
    cmap_traj(j+tot_con*(i-1))=cmap_run(j)
   enddo 
  enddo
  close(11)

 endif !! end case DCD, PDB, INDEX case 

!! case of TRAJ input file
 if(traj_on) then

  open(240,file=traj,form="formatted")
 
!! first reading for allocation
  tot_frames = 0
  tot_con    = 0
  do
   read(240,'(a80)',iostat=iostat)line
   if(iostat/=0)exit
   read(line,*)str
   if(index(str,"END").ne.0) then 
    tot_frames=tot_frames+1 
    cycle
   endif
   if(tot_frames==1) tot_con=tot_con+1
  enddo 

  write(*,*) "Total number of contacts          :: ",tot_con
  write(*,*) "Total number of CMAP in TRAJ file :: ",tot_frames

  allocate(map_par1(tot_con,4))
  allocate(map_par2(tot_con,3))
  allocate(cmap_tmp(tot_con))
  allocate(cmap_coord1(tot_con))
  allocate(cmap_coord2(tot_con))
  allocate(cmap_run(tot_con))
  allocate(cmap_frame(tot_con,nframes))
  allocate(cmap_traj(tot_frames*tot_con))
  rewind(240)

!! reading for storing CMAPs
  do i=1,tot_frames
   do j=1,tot_con
     read(240,13) ll,map_par1(j,1),map_par1(j,2),cmap_traj(j+tot_con*(i-1))
   enddo
   read(240,*) 
  enddo
  close(240)

 endif !! end of input TRAJ case

!! reading input CMAP
 do i=1,nframes

  done(i)=.false. 
  write(numb_to_char,'(i15)')i
  numb_to_char=(adjustl(numb_to_char))
  complex_name=trim(cmapvalue)//trim(numb_to_char)
  inquire(file=complex_name,EXIST=exist)
  if((.not.exist).and.((i==1).or.(i==nframes))) then
     write(*,*)"ERROR!! ",complex_name, " is missing. You need these files for first and last frame!"
     stop
  endif

!! if CMAP is found...
  if(exist) then 
   write(*,*)"** FOUND CMAP file ",adjustl(trim(complex_name))," for frame ",i
!! reading CMAP file
   open(unit=55,file=complex_name,form="formatted")
   do kk=1,tot_con
    read(55,13) ll,l,m,cmap_frame(kk,i)
   enddo
   close(55) 
  
!! done initialization
   done(i)=.true.
  endif

 enddo  !! loop on CMAP input reading

 write(*,*) 
 write(*,*) "** Starting FRAMES initialization..."

!! For the frames already initialized, look for the closer frame in cmap_traj
!! This is useful when optimizing path with SM-like procedure
 do i=2,nframes-1
  
  if(done(i)) then
   errmin=1.d40 
   do kk=1,tot_frames
!! copy map on a tmp array
    do ii=1,tot_con
     cmap_run(ii)=cmap_traj(ii+tot_con*(kk-1))
    enddo
!!
    dist=sum((cmap_frame(:,i)-cmap_run(:))**2)
    if(dist.le.errmin) then
      errmin   = dist
      mm       = kk
      cmap_tmp = cmap_run
    endif
   enddo
   cmap_frame(:,i) = cmap_tmp
   index_frame(i)  = mm
  endif
 enddo

!! initializing the other frames 
 do iframe=1,tot_frames
 
  do j=1,tot_con
   cmap_run(j)=cmap_traj(j+tot_con*(iframe-1))
  enddo

  err1=sum((cmap_run(:)-cmap_frame(:,1))**2)
  errn=sum((cmap_run(:)-cmap_frame(:,nframes))**2)
  indice=int( (err1**2/(err1**2+errn**2))*real(nframes-1) +1.d0 ) 

  if(indice==1) indice = 2
  if(indice==nframes) indice=nframes-1
 
  if(.not.done(indice)) then
   index_frame(indice)=iframe
   done(indice)=.true.
   write(*,*)"DONE INIT FRAME ",indice
   do j=1,tot_con
     cmap_frame(j,indice)=cmap_run(j)
   enddo 
  endif

 enddo

!! check if all the frames have been initialized

 write(*,*) "** Initialization FINISHED **"
 do i=2,nframes-1

  write(*,*) " Frame ",i," assigned to traj frame # ",index_frame(i)
  if(.not.done(i)) then
   write(*,*) " WARNING!! Frame ",i," NOT ASSIGNED "
   do j=i+1,nframes
   if (done(j)) exit 
   enddo 
   index_frame(i)=int((index_frame(j)-index_frame(i-1))/(j-i+1))+index_frame(i-1)  
   done(i)=.true.
   do j=1,tot_con
     cmap_frame(j,i)=cmap_traj(j+tot_con*(index_frame(i)-1))
   enddo
   write(*,*) " GUESSING... Frame ",i," assigned to traj frame # ",index_frame(i) 
  endif

 enddo 
   
!! starting choose best frame 

 upot=1.d40

!! loop on MC iterations 
 do ll=1,iter 

       found=0

       write(*,*) "ITER ",ll

!! loop on CMAP traj
       do iframe=1,tot_frames
        
!! storing running CMAP
        do j=1,tot_con
           cmap_run(j)=cmap_traj(j+tot_con*(iframe-1))
        enddo

!! find the candidate for substitution
        minerr=10000.d0
        do k=2,nframes-1
         err=sum((cmap_frame(:,k)-cmap_run(:))**2) 
         if(err.lt.minerr) then
          indice=k
          minerr=err
         endif
        enddo 

!! build err table
        do i=1,nframes-1 
          err_table(i,i)=0.d0
          if(i.eq.indice) then
           cmap_coord1(:)=cmap_run(:) 
          else
           cmap_coord1(:)=cmap_frame(:,i)
          endif 
          
          do j=i+1,nframes
            err_table(j,j)=0.d0
            if(j.eq.indice) then
             cmap_coord2(:)=cmap_run(:) 
            else
             cmap_coord2(:)=cmap_frame(:,j)
            endif 
          
            err=sum((cmap_coord1(:)-cmap_coord2(:))**2)
 
            err=err/redux  !units avoid sublinear behaviourr
            err_table(i,j)=err
            err_table(j,i)=err
          enddo
        enddo
       
!! new potential energy after MC move
        upot_new=0.d0
        pippo1=0.d0
        pippo2=0.d0
        pippo3=0.d0

        do i=1,nframes-1
            upot_new=upot_new+2.d0*((err_table(i,i+1) -offset )**xp)
            pippo1=pippo1+2.d0*((err_table(i,i+1) -offset )**xp)
        enddo
        if(.not.firstneigh) then 
         do i=1,nframes-2
            do j=i+1,nframes-1
              upot_new=upot_new+ (abs(err_table(i,j)-err_table(i,j+1)) - offset )**xp
              pippo2=pippo2+ (abs(err_table(i,j)-err_table(i,j+1)) - offset )**xp
            enddo
         enddo
         do j=3,nframes
            do i=2,j-1
              upot_new=upot_new+ (abs(err_table(i,j)-err_table(i-1,j)) - offset )**xp
              pippo3=pippo3+(abs(err_table(i,j)-err_table(i-1,j)) - offset )**xp 
            enddo
         enddo
        endif

!! MC move
        if(upot_new.lt.upot)then
              found=found+1
              cmap_frame(:,indice)=cmap_run(:)
              index_frame(indice)=iframe
              upot=upot_new
        else
           if(ran1(idum).lt.exp(-beta*(upot_new-upot) ) )then
              found=found+1
              upot=upot_new
              cmap_frame(:,indice)=cmap_run(:) 
              index_frame(indice)=iframe
              write(*,'(a5,1x,f18.6,i3,3(1x,f18.6))')"UPOT ",upot,indice,pippo1,pippo2,pippo3
           endif
        endif   

 
       enddo   !! end loop on CMAP traj 


       if(found.eq.nframes-2) then
        write(*,*) "DONE!!! No more optimal frames found..."
        goto 700
       endif

 enddo  ! ITER loop


 700  continue

 write(*,*)
 write(*,*) "*** Finishing... Writing final.CMP"
 do i=2,nframes-1
  write(*,*) "Frame ",i, " comes from TRAJ frame # ",index_frame(i)
 enddo

!! Writing final stuff....
 open(file="final.CMAP",unit=49,form="formatted")

 do m=1,nframes
  do k=1,tot_con
     write(49,13) k,map_par1(k,1),map_par1(k,2),cmap_frame(k,m)
  enddo
  write(49,14) "END"
 enddo

 if(dcd_on) then 
  write(*,*) "Writing final.pdb with all the atoms of DCD !!"
  open(3000,file="final.pdb",status="unknown")

  do j=2,nframes-1

     open (unit=11,file=dcd,form="unformatted") 
     read (11)car4,nset, nstart,nsanc,i5,namin,delta,i9
     read (11)ntitle,(car(i),i=1,ntitle)
     read (11)ntap
     i=0

     do
       if(cell_logic.eq.1) then 
        read(11,iostat=iostat)cell
        if(iostat/=0) exit
       endif
       read(11,iostat=iostat)(x4(m),m=1,tot_atom)
       if(iostat/=0) exit
       read(11,iostat=iostat)(y4(m),m=1,tot_atom)
       if(iostat/=0) exit
       read(11,iostat=iostat)(z4(m),m=1,tot_atom)
       if(iostat/=0) exit
       i=i+1
       if(i.ne.index_frame(j)) cycle

       write(3000,'(a16,i3,a17,i8)') "REMARK. FRAME # ",j, " from TRAJ frame ",index_frame(j)
       do ii=1,tot_atom
           write(3000,'(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,f8.3,f8.3,f8.3,f6.2,f6.2)') &
             &recordname(ii),atmnum(ii),atmname(ii),&
             &altloc(ii),resname(ii),chainid(ii),&
             &resnum(ii),icode(ii),x4(ii),y4(ii),z4(ii),&
             &occ(ii),temp(ii)
       enddo
       write(3000,'(a3)') "END"
 
     enddo !! loop on DCD

     close(11)

  enddo !! loop on frames

 endif !! end of if dcd_on

 end  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE cmap_running(n,num,tot_acon,tot_gcon,x,y,z,cmap,map1,map2,ngroup,g_number,g_index)

  IMPLICIT NONE

  integer      n,tot_acon,tot_gcon
  real*4       x(n),y(n),z(n)
  integer      num(n) 
  real*8       dist,cmap(tot_acon+tot_gcon),map2(tot_acon+tot_gcon,3)
  integer      i,ii,j,jj,kk,kkk
  integer      i1,i2
  integer      map1(tot_acon+tot_gcon,4)
  integer      ngroup_max, max_atom_group
  parameter(ngroup_max=10,max_atom_group=50)
  integer      ngroup,g_number(ngroup_max),g_index(ngroup_max,max_atom_group)
  real*8       rcm(3,ngroup_max)


  do i=1,tot_acon

   jj=map1(i,1)
   kk=map1(i,2)

!! find right indexes
   do ii=1,n
    if(num(ii)==jj) i1=ii
    if(num(ii)==kk) i2=ii
   enddo

   dist=sqrt(real((x(i1)-x(i2))**2+(y(i1)-y(i2))**2+(z(i1)-z(i2))**2)) 

   if(dist.gt.map2(i,2)) then 
    cmap(i)=0.d0
   else
    if(dist.ne.map2(i,1)) then
     cmap(i)=map2(i,3)*(1.d0-(dist/map2(i,1))**map1(i,3))/(1.d0-(dist/map2(i,1))**map1(i,4))
    else
     cmap(i)=map2(i,3)*map1(i,3)/map1(i,4)
    endif
   endif

  enddo


  if(tot_gcon.ne.0) then
!! calculating center of mass

    do j=1,ngroup

     rcm(:,j)=0.d0

     do i=1,g_number(j)

      ii=g_index(j,i)

!! find right index
      do kkk=1,n
       if(num(kkk)==ii) i1=kkk
      enddo

      rcm(1,j)=rcm(1,j)+real(x(i1))
      rcm(2,j)=rcm(2,j)+real(y(i1))
      rcm(3,j)=rcm(3,j)+real(z(i1)) 

     enddo
     rcm(:,j)=rcm(:,j)/real(g_number(j))

    enddo

!! evaluating contacts
    do j=1,tot_gcon

     i=j+tot_acon
     jj=map1(i,1)
     kk=map1(i,2)

     dist=sqrt((rcm(1,jj)-rcm(1,kk))**2+(rcm(2,jj)-rcm(2,kk))**2+(rcm(3,jj)-rcm(3,kk))**2)

     if(dist.gt.map2(i,2)) then 
      cmap(i)=0.d0
     else
      if(dist.ne.map2(i,1)) then
       cmap(i)=map2(i,3)*(1.d0-(dist/map2(i,1))**map1(i,3))/(1.d0-(dist/map2(i,1))**map1(i,4))
      else
       cmap(i)=map2(i,3)*map1(i,3)/map1(i,4)
      endif
     endif

    enddo

  endif

 return

    end

