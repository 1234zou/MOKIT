! written by jxzou at 20201208: copy from bas_gms2molcas.f90 and do modifications
! This is a subroutine used to transform the basis sets in GAMESS format to those in Molpro format

! Note: Currently isotopes are not tested.
program main
 use pg, only: iout
 implicit none
 integer :: i
 character(len=4) :: str
 character(len=240) :: fname
 ! fname: input file contains basis sets and Cartesian coordinates in GAMESS format
 logical :: alive, spherical

 i = iargc()
 if(i<1 .or. i>2) then
  write(iout,'(/,A,/)') 'Example1: bas_gms2molpro a.inp (generating a Molpro a.com file)'
  write(iout,'(A,/)') "Example2: bas_gms2molpro a.inp -sph (without 'Cartesian')"
  stop
 end if

 str = ' '
 fname = ' '
 spherical = .false.
 call getarg(1,fname)

 inquire(file=TRIM(fname),exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine bas_gms2molpro: file does not exist!'
  write(iout,'(A)') 'Filename='//TRIM(fname)
  stop
 end if

 if(i == 2) then
  call getarg(2,str)

  if(str == '-sph') then
   spherical = .true.
  else
   write(iout,'(A)') 'ERROR in subroutine bas_gms2molpro: wrong command line arguments!'
   write(iout,'(A)') 'str='//str
   stop
  end if

 end if

 call bas_gms2molpro(fname, spherical)
 stop
end program main

! Transform the basis sets in GAMESS format to those in Molpro format
subroutine bas_gms2molpro(fort7, spherical)
 use pg, only: iout, natom, ram, ntimes, coor, elem, prim_gau, all_ecp, ecp_exist
 implicit none
 integer :: i, j, k, m, n, nline, rc, rel
 integer :: fid1, fid2
 integer :: charge, mult
 real(kind=8) :: rtmp
 real(kind=8) :: exp1   ! the first exponent
 character(len=7) :: str
 character(len=240), intent(in) :: fort7
 character(len=240) :: buf, input ! input is the Molpro .com file
 character(len=240) :: orbfile, orbfile2 ! Molpro MOs, plain text file
 character(len=1) :: stype0, stype
 logical, intent(in) :: spherical
 logical :: uhf, X2C

 ! initialization
 buf = ' '
 input = ' '
 orbfile = fort7
 orbfile2 = fort7
 call convert2molpro_fname(orbfile, '.a')
 call convert2molpro_fname(orbfile2, '.b')
 i = index(fort7, '.', back=.true.)
 input = fort7(1:i)//'com'

 open(newunit=fid2,file=TRIM(input),status='replace')
 write(fid2,'(A)') '***,auto-generated file by bas_gms2molpro of MOKIT'
 write(fid2,'(A)') 'symmetry,nosym,noorient'
 write(fid2,'(A)') 'geometry={'

 call read_natom_from_gms_inp(fort7, natom)
 allocate(elem(natom), ram(natom), coor(3,natom), ntimes(natom))
 call read_elem_nuc_coor_from_gms_inp(fort7, natom, elem, ram, coor)
 call calc_ntimes(natom, elem, ntimes)

 do i = 1, natom, 1
  str = ' '
  write(str,'(A,I0)') TRIM(elem(i)), ntimes(i)
  write(fid2,'(A7,1X,3F18.8)') str, coor(1:3,i)
 end do ! for i
 deallocate(coor)

 write(fid2,'(A)') '}'
 if(.not. spherical) write(fid2,'(A)') 'Cartesian'
 write(fid2,'(A)') 'basis={'

 uhf = .false.
 charge = 0; mult = 1

 open(newunit=fid1,file=TRIM(fort7),status='old',position='rewind')
 read(fid1,'(A)') buf
 call upper(buf)
 if(index(buf,'UHF') > 0) uhf = .true.

 i = index(buf,'ICHAR')
 if(i > 0) read(buf(i+7:),*) charge
 i = index(buf,'MULT')
 if(i > 0) read(buf(i+5:),*) mult

 ! check if any ECP data exists
 ecp_exist = .false.
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf(3:5))
  if(buf(2:5) == '$ECP') then
   ecp_exist = .true.
   exit
  end if
 end do ! for while

 if(ecp_exist) then ! if exists, read
  allocate(all_ecp(natom))

  do i = 1, natom, 1
   read(fid1,'(A)') buf
   call upper(buf)
   if(index(buf,'NONE') /= 0) cycle

   all_ecp(i)%ecp = .true.
   k = index(buf,'GEN')
   read(buf(k+3:),*) all_ecp(i)%core_e, m
   all_ecp(i)%highest = m
   allocate(all_ecp(i)%potential(m+1))

   do j = 1, m+1, 1
    read(fid1,'(A)') buf
    if(j == 1) then
     k = index(buf,'ul')
     if(k == 0) then
      write(iout,'(A)') "ERROR in subroutine bas_gms2molpro: ECP/PP not starts with '-ul potential'."
      write(iout,'(A)') 'You should check the format of ECP/PP data in file '//TRIM(fort7)//'.'
      stop
     end if
    end if

    read(buf,*) n
    all_ecp(i)%potential(j)%n = n
    allocate(all_ecp(i)%potential(j)%col1(n), source=0.0d0)
    allocate(all_ecp(i)%potential(j)%col2(n), source=0)
    allocate(all_ecp(i)%potential(j)%col3(n), source=0.0d0)
    do k = 1, n, 1
     read(fid1,*) all_ecp(i)%potential(j)%col1(k), all_ecp(i)%potential(j)%col2(k), &
     all_ecp(i)%potential(j)%col3(k)
    end do ! for k
   end do ! for j
  end do ! for i
 end if
 ! read ECP/PP done

 ! rewind and find the $DATA section
 rewind(fid1)
 do while(.true.)
  read(fid1,'(A)',iostat=rc) buf
  if(rc /= 0) exit
  call upper(buf(3:6))
  if(buf(2:6) == '$DATA') exit
 end do
 if(rc /= 0) then
  write(iout,'(A)') 'ERROR in subroutine bas_gms2molpro: No $DATA section found&
                   & in file '//TRIM(fort7)
  close(fid1)
  stop
 end if

 ! skip 2 lines: the Title line and the Point Group line
 read(fid1,'(A)') buf
 read(fid1,'(A)') buf

 ! initialization: clear all primitive gaussians
 call clear_prim_gau()

 ! read and print basis set data
 do i = 1, natom, 1
  read(fid1,'(A)',iostat=rc) buf
  ! 'buf' contains the element, ram and coordinates
  if(rc /= 0) exit
  if(buf(2:2) == '$') call upper(buf(3:5))
  if(buf(2:5) == '$END') exit
  write(fid2,'(A)') '! '//TRIM(elem(i))

  ! deal with primitive gaussians
  stype0 = ' '
  do while(.true.)
   read(fid1,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit

   read(buf,*) stype, nline
   if(stype0 == ' ') then
    !------------------------------------------------------
    ! the following 5 lines are added to determine whether
    ! an angular momentum occurs more than once
    call stype2itype(stype, k)
    if( allocated(prim_gau(k)%coeff) ) then
     stype0 = stype
     BACKSPACE(fid1)
    else
    !------------------------------------------------------
     call read_prim_gau1(stype, nline, fid1)
     stype0 = stype
    end if
   else ! stype0 /= ' '
    if(stype == stype0) then
     exp1 = 0.0d0
     read(fid1,*) j, exp1
     BACKSPACE(fid1)
     call read_prim_gau2(stype, nline, fid1, exp1)
    else ! stype /= stype0
     stype0 = ' ' ! reset stype0
     BACKSPACE(fid1)
    end if
   end if
  end do ! for while 2

  ! print basis sets and ECP/PP (if any) of this atom in Molcas format
  call prt_prim_gau_molpro(i, fid2)

  ! clear all primitive gaussians for next cycle
  call clear_prim_gau()
 end do ! for i

 deallocate(ram, ntimes, elem)
 if(allocated(all_ecp)) deallocate(all_ecp)
 close(fid1)

 if(natom==0 .or. rc/=0) then
  if(natom == 0) write(iout,'(A)') 'ERROR in subroutine bas_gms2molpro: zero atom found!'
  if(rc /= 0) write(iout,'(A)') "ERROR in subroutine bas_gms2molpro: it seems the '$DATA'&
                               & has no corresponding '$END'. EOF."
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(A)') '}'

 call check_DKH_in_gms_inp(fort7, rel)
 if(rel>-1 .and. rel<2) then
  write(iout,'(A)') 'ERROR in subroutine bas_gms2molpro: DKH0 or DKH1 not supported.'
  write(iout,'(A)') 'Please use DKH2 at least.'
  stop
 else if(rel > 1) then ! at least DKH2
  call check_X2C_in_gms_inp(fort7, X2C)
  if(X2C) then
   write(fid2,'(A)') 'SET,DKHO=101'
  else
   write(fid2,'(A,I0)') 'SET,DKHO=', rel
  end if
 end if

 write(fid2,'(2(A,I0))') 'wf,charge=',charge,',spin=',mult-1
 write(fid2,'(A)') '{matrop;'
 write(fid2,'(A)') 'read,mo,ORB,file='//TRIM(orbfile)//';'

 if(uhf) then
  write(fid2,'(A)') 'save,mo,2200.2,ORBITALS,ALPHA;'
  write(fid2,'(A)') 'read,mo,ORB,file='//TRIM(orbfile2)//';'
  write(fid2,'(A)') 'save,mo,2200.2,ORBITALS,BETA}'
  write(fid2,'(A3)') 'UHF'
 else
  write(fid2,'(A)') 'save,mo,2100.2,ORBITALS}'
  write(fid2,'(A2)') 'HF'
 end if

 write(fid2,'(A)',advance='no') '{put,xml'
 if(spherical) write(fid2,'(A)',advance='no') ';keepspherical'
 write(fid2,'(A)') '}'

 close(fid2)
 return
end subroutine bas_gms2molpro

! print primitive gaussians
subroutine prt_prim_gau_molpro(iatom, fid)
 use pg, only: prim_gau, elem, ntimes, all_ecp, ecp_exist
 implicit none
 integer :: i, j, k, m, n, nline, ncol
 integer, intent(in) :: iatom, fid
 integer, allocatable :: list(:)
 character(len=1), parameter :: am(0:5) = ['s','p','d','f','h','h']

 call get_highest_am()

 do i = 1, 7, 1
  if(.not. allocated(prim_gau(i)%coeff)) cycle

  write(fid,'(A,I0)',advance='no') &
   prim_gau(i)%stype//', '//TRIM(elem(iatom)), ntimes(iatom)

  nline = prim_gau(i)%nline
  ncol = prim_gau(i)%ncol
  do j = 1, nline-1, 1
   write(fid,'(A1,ES15.8)',advance='no') ',',prim_gau(i)%coeff(j,1)
  end do ! for j
  write(fid,'(A1,ES15.8)') ',',prim_gau(i)%coeff(nline,1)

  do j = 2, ncol, 1
   write(fid,'(A,I0)',advance='no') 'c, 1.', nline

   do k = 1, nline-1, 1
    write(fid,'(A1,ES15.8)',advance='no') ',',prim_gau(i)%coeff(k,j)
   end do ! for k
   write(fid,'(A1,ES15.8)') ',',prim_gau(i)%coeff(nline,j)
  end do ! for j
 end do ! for i

 if(.not. ecp_exist) return

 if(all_ecp(iatom)%ecp) then
  m = all_ecp(iatom)%highest
  write(fid,'(A,I0,1X,I0,1X,I0)') 'ECP '//TRIM(elem(iatom)),ntimes(iatom),all_ecp(iatom)%core_e,m
  allocate(list(m+1))
  list(1) = m
  forall(i=2:m+1) list(i) = i-2
  do i = 1, m+1, 1
   n = all_ecp(iatom)%potential(i)%n
   if(i == 1) then
    write(fid,'(I0,A)') n,'; !  ul potential'
   else ! i > 1
    write(fid,'(I0,A)') n,'; !  '//am(list(i))//'-ul potential'
   end if

   do j = 1, n, 1
    write(fid,'(I0,2(1X,F16.8))') all_ecp(iatom)%potential(i)%col2(j), &
     all_ecp(iatom)%potential(i)%col3(j), all_ecp(iatom)%potential(i)%col1(j)
   end do ! for j
  end do ! for i
 end if

 return
end subroutine prt_prim_gau_molpro

