  PROGRAM create_CMAP

  implicit none

  integer      :: mm,ll,l,i,j,ii,jj,numatom,kk,aa,bb
  integer      :: argcount,iostat,tot_con,tot_gcon,new_con
  integer      :: g1, g2
  integer, allocatable :: index_1(:), index_2(:)
  integer      :: flag,ngroup,nn,nd,exclude

  integer      :: ngroup_max, max_atom_group
  parameter(ngroup_max=10,max_atom_group=50)
  integer      :: g_number(ngroup_max),g_index(ngroup_max,max_atom_group)
  real*8       :: rcm(3,ngroup_max), dist, cmap, dcut, d0 
  character*1  :: ics
  character*8  :: str
  character*45 :: cmapgroup,cmapvalue,cmapindex,pdbname
  character*80 :: line,wq_char
  logical      :: go,from_scratch
  logical      :: par1,par2,par3,par4
  logical      :: log_all, log_select, log_ofr, log_nm

!! PDB stuff
  character*7              :: dummy
  character*6, allocatable :: name(:)
  character*5, allocatable :: atom(:)
  character*3, allocatable :: residue(:)
  character*1, allocatable :: icode(:),altloc(:)
  integer, allocatable     :: num(:),map_par1(:,:)
  integer, allocatable     :: numres(:)
  real*8, allocatable      :: x(:),y(:),z(:)
  real*8, allocatable      :: occu(:),charge(:),map_par2(:,:)


 11   FORMAT(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,f8.3,f8.3,f8.3,f6.2,f6.2)
 12   FORMAT(a7,1x,i4,1x,i4,1x,i4,1x,f5.2,1x,i2,1x,i2,1x,f5.1,1x,f8.5) 
 13   FORMAT(i4,1x,i4,1x,i4,4x,f12.8)
 14   FORMAT(a3)
     
  argcount = IARGC()

!! initializing stuff
  go         = .false.
  par1       = .false.
  par2       = .false.
  par3       = .false.
  par4       = .false. 
  log_all    = .false.
  log_select = .false.
  log_ofr    = .false.
  log_nm     = .false.
  exclude    = 0
  cmapvalue  = ""
  cmapindex  = ""
  cmapgroup  = ""

  write(*,*)
  write(*,*) " --> CMAP_TOOLS :: create_CMAP <--"
  write(*,*)

  do i=1,argcount

    call getarg(i,wq_char)
    if (INDEX(wq_char,'-pdb').NE.0)then
      call getarg(i+1,wq_char)
      read(wq_char,*)pdbname
      go=.true.
    endif
    if (INDEX(wq_char,'-cmap').NE.0)then
      call getarg(i+1,wq_char)
      read(wq_char,*)cmapvalue
    endif
    if (INDEX(wq_char,'-index').NE.0)then
      call getarg(i+1,wq_char)
      read(wq_char,*)cmapindex
    endif
    if (INDEX(wq_char,'-group').NE.0)then
      call getarg(i+1,wq_char)
      read(wq_char,*)cmapgroup
    endif
    if (INDEX(wq_char,'-d0').NE.0)then
      call getarg(i+1,wq_char)
      read(wq_char,*) d0
      par1=.true.
    endif
    if (INDEX(wq_char,'-nn').NE.0)then
      call getarg(i+1,wq_char)
      read(wq_char,*) nn
      par2=.true.
    endif
    if (INDEX(wq_char,'-nd').NE.0)then
      call getarg(i+1,wq_char)
      read(wq_char,*) nd
      par3=.true.
    endif
    if (INDEX(wq_char,'-dcut').NE.0)then
      call getarg(i+1,wq_char)
      read(wq_char,*) dcut
      par4=.true.
    endif
    if (INDEX(wq_char,'-exclude').NE.0)then
      call getarg(i+1,wq_char)
      read(wq_char,*) exclude
    endif
    if (INDEX(wq_char,'-all').NE.0)then
      call getarg(i+1,wq_char)
      log_all=.true. 
    endif
    if (INDEX(wq_char,'-select').NE.0)then
      call getarg(i+1,wq_char)
      log_select=.true. 
    endif
    if (INDEX(wq_char,'-rmofr').NE.0)then
      log_ofr=.true. 
    endif
    if (INDEX(wq_char,'-nm').NE.0)then
      log_nm=.true. 
    endif
  enddo

  if(.not.go) then
     write(*,*)"USAGE : "
     write(*,*)"./create_CMAP -pdb frame.pdb (-index CMAPINDEX -group CMAPGROUP -cmap CMAPVALUE &
                &-d0 8.5 -dcut 20.0 -nn 6 -nd 10 -exclude 1 -all/-select -rmofr)"
     write(*,*)                                                
     write(*,*)"-pdb     : input pdb file  " 
     write(*,*)"-index   : index file       (optional)"
     write(*,*)"-group   : group file       (optional)"
     write(*,*)"-cmap    : CMAP output file (optional)"
     write(*,*)"-d0      : switching param. (optional)"
     write(*,*)"-dcut    : switching param. (optional)"
     write(*,*)"-nn      : switching param. (optional)"
     write(*,*)"-nd      : switching param. (optional)"
     write(*,*)"-exclude : switching param. (optional)"
     write(*,*)"-all     : all contacts between PDB atoms"
     write(*,*)"-select  : contact between two groups of selected atoms"
     write(*,*)"-rmofr   : remove those contact farer than dcut"
     write(*,*)"-nm      : use nanometers as output length unit"
     write(*,*)
     stop
  endif


 if(cmapindex.eq."") then
!! no CMAP index in input
  from_scratch=.true.
  write(*,*) "** No index file in input"
  if (.not.(par1.and.par2.and.par3.and.par4)) then
   write(*,*) "ERROR! If no index file is specified, you have to provide the parameters for the &
               &switching function on the command line, ex:"
   write(*,*)
   write(*,*) "./create_CMAP -pdb frame.pdb -d0 8.5 -dcut 20.0 -nn 6 -nd 10 (-exclude 2) -all (-rmofr)"
   write(*,*)
   stop
  endif
  if (par1.and.par2.and.par3.and.par4) then
   write(*,*) "...using parameters provided on the command line..."
   cmapindex="CMAPINDEX"
   if (cmapvalue.eq."") cmapvalue="CMAPVALUES"
  endif
  if(.not.(log_all.or.log_select)) then
   write(*,*) "ERROR! If no index file is specified, you have to choose between: "
   write(*,*) "       -all    :: all the possible contacts" 
   write(*,*) "       -select :: contacts between the two groups specified in the occupancy column"
   write(*,*)
   stop
  endif
 else
!! CMAP index provided in input
  from_scratch=.false.
  write(*,*) "** Index file in input: evaluating parameters from file"
   if (cmapvalue.eq."") cmapvalue="CMAPVALUES"
 endif

  write(*,*)
  write(*,*) "PDB   file    : ",pdbname
  write(*,*) "index file    : ",cmapindex
  write(*,*) "group file    : ",cmapgroup 
  write(*,*) "CMAP  file    : ",cmapvalue
  if(from_scratch) then
   write(*,*) "d0            : ",d0
   write(*,*) "dcut          : ",dcut
   write(*,*) "nn            : ",nn
   write(*,*) "nd            : ",nd
   write(*,*) "exclude       : ",exclude
   write(*,*)
  endif

 open(2000,file=cmapvalue,form="formatted")
 open(1000,file=cmapindex,form="formatted")

!! open PDB and counting number of atoms
  open (UNIT=10,FILE=pdbname,FORM="FORMATTED")

  numatom=0
  do
   read(10,'(a80)',iostat=iostat)line
   if(iostat/=0) exit
   read(line,*)str
   if(index(str,"ATOM").ne.0) numatom=numatom+1
  enddo

  rewind(10)
  write(*,*) "The number of atoms in pdb is :: ",numatom

  allocate(name(numatom))
  allocate(atom(numatom))
  allocate(residue(numatom))
  allocate(num(numatom))
  allocate(numres(numatom))
  allocate(x(numatom))
  allocate(y(numatom))
  allocate(z(numatom))
  allocate(occu(numatom))
  allocate(charge(numatom))
  allocate(altloc(numatom))
  allocate(icode(numatom))

!! read pdb
 kk = 0 
 g1 = 0
 g2 = 0
 do
  read (10,'(a80)',iostat=iostat) line
  if(iostat/=0) exit
  if(index(line,"ATOM").ne.0) then
   kk=kk+1
   read (line,11) name(kk),num(kk),atom(kk),altloc(kk),residue(kk),ics,numres(kk),icode(kk),x(kk),y(kk),z(kk),occu(kk),charge(kk)
   if(int(occu(kk))==1) g1 = g1 + 1
   if(int(occu(kk))==2) g2 = g2 + 1
  endif
  if((index(line,"END").ne.0) .or. (index(line,"TER").ne.0)) exit
 enddo 

 if(log_select) then
  write(*,*) "Group 1 number of atoms       :: ", g1
  write(*,*) "Group 2 number of atoms       :: ", g2
  if((g1==0).or.(g2==0)) then
   write(*,*) "ERROR! You have selected the option contacts between two groups, but you have not"
   write(*,*) "       specified the two groups in the occupancy column of the PDB"
   write(*,*) 
   stop
  endif
!! allocating and creating index vector
  allocate(index_1(g1), index_2(g2))
  jj=0
  kk=0
  do i=1,numatom
   if(occu(i)==1) then
    jj=jj+1
    index_1(jj)=i
   endif
   if(occu(i)==2) then
    kk=kk+1
    index_2(kk)=i
   endif 
  enddo
 endif !! end log_select case

 if(from_scratch) then

  if(log_all) then 
     tot_con=numatom*(numatom-1)/2
  else
     tot_con=g1*g2
  endif

  tot_gcon=0
  allocate(map_par1(tot_con,4))
  allocate(map_par2(tot_con,3))

!! creating CMAP index in case of log_all=.true.
  if(log_all) then
   jj=1
   do i=1,numatom
    do j=i+1,numatom
     map_par1(jj,1)=num(i)
     map_par1(jj,2)=num(j)
     map_par1(jj,3)=nn
     map_par1(jj,4)=nd
     if(abs(numres(i)-numres(j)).le.exclude) then 
      map_par2(jj,3)=0.d0
     else
      map_par2(jj,3)=1.d0
     endif
     map_par2(jj,1)=d0
     map_par2(jj,2)=dcut
     jj=jj+1   
    enddo
   enddo 
  else
!! creating CMAP index in case log_select=.true.
  jj=1
   do i=1,g1
    mm=index_1(i)
    do j=1,g2 
     ll=index_2(j)
     map_par1(jj,1)=num(mm)
     map_par1(jj,2)=num(ll)
     map_par1(jj,3)=nn
     map_par1(jj,4)=nd
     if(abs(numres(mm)-numres(ll)).le.exclude) then
      map_par2(jj,3)=0.d0
     else
      map_par2(jj,3)=1.d0
     endif
     map_par2(jj,1)=d0
     map_par2(jj,2)=dcut
     jj=jj+1
    enddo
   enddo
  endif
   
!! common printout
  new_con=0
  do jj=1,tot_con
    if(log_ofr) then
      aa=map_par1(jj,1)
      bb=map_par1(jj,2) 
      do l=1,numatom
        if(num(l).eq.aa) then 
          ll=l
          flag=flag+1
        endif
        if(num(l).eq.bb) then
          mm=l
          flag=flag+1
        endif
      enddo    
      dist=sqrt((x(ll)-x(mm))**2+(y(ll)-y(mm))**2+(z(ll)-z(mm))**2) 
      if(dist.le.map_par2(jj,2).and.map_par2(jj,3).gt.0.d0) then
        if(log_nm) then
          write(1000,12) "CONTACT",new_con+1,map_par1(jj,1),map_par1(jj,2),map_par2(jj,1)/10,map_par1(jj,3),map_par1(jj,4),&
                      map_par2(jj,2)/10,map_par2(jj,3)
        else 
          write(1000,12) "CONTACT",new_con+1,map_par1(jj,1),map_par1(jj,2),map_par2(jj,1),map_par1(jj,3),map_par1(jj,4),&
                      map_par2(jj,2),map_par2(jj,3)
        endif
        new_con=new_con+1
      endif
    else
        if(log_nm) then 
              write(1000,12) "CONTACT",jj,map_par1(jj,1),map_par1(jj,2),map_par2(jj,1)/10,map_par1(jj,3),map_par1(jj,4),&
                    map_par2(jj,2)/10,map_par2(jj,3)
        else 
              write(1000,12) "CONTACT",jj,map_par1(jj,1),map_par1(jj,2),map_par2(jj,1),map_par1(jj,3),map_par1(jj,4),&
                    map_par2(jj,2),map_par2(jj,3)
        endif
        new_con=new_con+1
    endif
  enddo
  write(*,*) "Total number of contacts      :: ",new_con

 else  !! case of CMAP index file

!! checking number of contacts in CMAP index file

  j=0
  i=0
  do
   read(1000,'(a80)',iostat=iostat)line
   if(iostat/=0) exit
   read(line,*)str
   if(index(str,"CONTACT").ne.0) j=j+1
   if(index(str,"GROUP").ne.0) i=i+1
  enddo

  tot_con=j
  tot_gcon=i
  rewind(1000)

  write(*,'(1x,a5,1x,i4,1x,a17,1x,i3,1x,a14)') "Found",tot_con,"atomic contacts /",tot_gcon,"group contacts" 

!! reading CMAP group file 
  if((tot_gcon.ne.0).and.(cmapgroup.eq.""))then
   write(*,*) "ERROR! Found GROUP contact in CMAP index file, but no CMAPGROUP file specified."
   stop 
  endif
  if(tot_gcon.ne.0)then
   open(3000,file=cmapgroup,status="unknown")
   i=1
   do
    read(3000,*,iostat=iostat) dummy,mm,g_number(i),(g_index(i,ii),ii=1,g_number(i))
    if(iostat/=0) exit
    i=i+1
   enddo
   g_index=g_index
   ngroup=i-1 
   write(*,*) "Number of group found: ",ngroup
  endif

!! allocating stuff 
  allocate(map_par1(tot_con+tot_gcon,4))
  allocate(map_par2(tot_con+tot_gcon,3))
!! reading cmapindex file

  do ii=1,tot_con+tot_gcon
   read(1000,12) dummy,i,map_par1(ii,1),map_par1(ii,2),map_par2(ii,1),map_par1(ii,3),map_par1(ii,4),map_par2(ii,2),map_par2(ii,3)
  enddo       

 endif !! end case CMAP index file

!! writing CMAP file
!! 1) atomic contacts
   new_con=0
 do i=1,tot_con
 
  jj=map_par1(i,1)
  kk=map_par1(i,2) 
!! find the indexes
  flag=0 
  do l=1,numatom
   if(num(l).eq.jj) then 
    ll=l
    flag=flag+1
   endif
   if(num(l).eq.kk) then
    mm=l
    flag=flag+1
   endif
  enddo    

  if(flag.ne.2) then
     write(*,*) "*** Warning!! Some atoms needed for atomic contacts are not in the pdb file! "
     stop
  endif

  dist=sqrt((x(ll)-x(mm))**2+(y(ll)-y(mm))**2+(z(ll)-z(mm))**2) 
  
  if(dist.gt.map_par2(i,2)) then 
   cmap=0.d0
  else
   cmap=map_par2(i,3)*(1.d0-(dist/map_par2(i,1))**map_par1(i,3))/(1.d0-(dist/map_par2(i,1))**map_par1(i,4))
  endif

   if(log_ofr) then
      if(cmap.gt.0.d0.and.map_par2(i,3).gt.0.d0) then
        write(2000,13) new_con+1,jj,kk,cmap
        new_con=new_con+1
      endif
   else
     write(2000,13) i,jj,kk,cmap
   endif

 enddo
       
!! 2) Group contacts

!! Calculating center of mass
 if(tot_gcon.ne.0) then

  do i=1,ngroup !! cycle on the number of groups
   rcm(:,i)=0.d0

   do j=1,g_number(i) !! cycle on the number of atoms of the group
    jj=g_index(i,j)
    flag=0
    do l=1,numatom
     if(num(l).eq.jj) then
      ll=l
      flag=1
     endif
    enddo
    if(flag.eq.0) then
     write(*,*) "*** Warning!! Some atoms needed for group contacts are not in the pdb file! "
     stop
    endif

    rcm(1,i)=rcm(1,i)+x(ll)
    rcm(2,i)=rcm(2,i)+y(ll)
    rcm(3,i)=rcm(3,i)+z(ll)
   enddo

   rcm(:,i)=rcm(:,i)/real(g_number(i))

  enddo

  do i=1,tot_gcon !! cycle on GROUP contacts

   j=i+tot_con
   jj=map_par1(j,1)
   kk=map_par1(j,2)

   dist=sqrt((rcm(1,jj)-rcm(1,kk))**2+(rcm(2,jj)-rcm(2,kk))**2+(rcm(3,jj)-rcm(3,kk))**2)

   if(dist.gt.map_par2(j,2)) then
    cmap=0.d0
   else
    cmap=map_par2(j,3)*(1.d0-(dist/map_par2(j,1))**map_par1(j,3))/(1.d0-(dist/map_par2(j,1))**map_par1(j,4))
   endif

   write(2000,13) j,jj,kk,cmap

  enddo

 endif !! if GROUP contacts

 write(2000,14) "END"
 write(*,*)    

 end
