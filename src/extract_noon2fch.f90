! written by jxzou at 20200503: extract Natural Orbital Occupation Numbers (NOONs)
!  from (1) .out file of PySCF or (2) .dat file of GAMESS, and write into a given .fch file

! updated by jxzou at 20200517: fix the bug(nmo=6*N) in subroutine read_noon_from_pyout
! updated by jxzou at 20201212: modify read_noon_from_gmsgms to read more digits of NOONs

! Note: [idx1,idx2] contains the singly occupied MOs and GVB pairs, npair=(idx2-idx1+1-nopen)/2

program main
 implicit none
 integer :: i, j, idx1, idx2, nopen
 integer, parameter :: iout = 6
 character(len=10) :: buf
 character(len=240) :: outname, fchname
 logical :: gau_order

 nopen = 0

 i = iargc()
 if(i<4 .or. i>6) then
  write(iout,'(/,A)') ' ERROR in subroutine extract_noon2fch: wrong command line arguments.'
  write(iout,'(A)')   ' Format: extract_noon_2fch outname fchname idx1 idx2 nopen [-gau]'
  write(iout,'(/,A)') ' Example 1(PySCF CASCI): extract_noon2fch a.out a.fch 19 24'
  write(iout,'(A)')   ' Example 2 (GAMESS GVB): extract_noon2fch a.dat a.fch 19 24 0'
  write(iout,'(A)')   ' Example 3 (GAMESS GVB): extract_noon2fch a.dat a.fch 19 24 0 -gau'
  write(iout,'(A)')   ' Example 4 (GAMESS CAS): extract_noon2fch a.gms a.fch 19 24'
  write(iout,'(A)')   ' Example 5 (ORCA CAS)  : extract_noon2fch a.out a.fch 19 24'
  write(iout,'(A,/)') ' Example 6 (PSI4 CAS)  : extract_noon2fch a.out a.fch 19 24'
  stop
 end if

 call getarg(1, outname)
 call require_file_exist(outname)

 call getarg(2, fchname)
 call require_file_exist(fchname)

 call getarg(3, buf)
 read(buf,*) idx1
 call getarg(4, buf)
 read(buf,*) idx2

 gau_order = .false.
 if(i > 4) then
  call getarg(5, buf)
  read(buf,*,iostat=j) nopen
  if(j /= 0) then
   write(iout,'(A)') 'ERROR in subroutine extract_noon2fch: wrong command line&
                    & arguments. Failed to read integer nopen.'
   stop
  end if
  if(nopen < 0) then
   write(iout,'(A)') 'ERROR in subroutine extract_noon2fch: wrong command line&
                    & arguments. nopen<0.'
   stop
  end if
 end if

 if(i == 6) then
  call getarg(6, buf)
  if(index(buf,'-gau') /= 0) then
   gau_order = .true.
  else
   write(iout,'(A)') 'ERROR in subroutine extract_noon2fch: wrong command line arguments.'
   write(iout,'(A)') 'The 5th argument='//TRIM(buf)
   stop
  end if
 end if

 i = idx2 - idx1
 if(i < 1) then
  write(iout,'(A)') 'ERROR in subroutine extract_noon2fch: wrong input indices.'
  write(iout,'(A,3I5)') 'idx1, idx2, nopen=', idx1, idx2, nopen
  stop
 end if

 if(index(outname,'.dat',back=.true.) /= 0) then
  if(MOD(i+1-nopen,2) /= 0) then
   write(iout,'(A)') 'ERROR in subroutine extract_noon2fch: wrong input indices.&
                    & In this case idx2-idx1+1-nopen must be an even integer.'
   write(iout,'(A,3I5)') 'idx1, idx2, nopen=', idx1, idx2, nopen
   stop
  end if
 end if

 call extract_noon2fch(outname, fchname, idx1, idx2, nopen, gau_order)
 stop
end program main

! extract NOONs from PySCF .out file, and print it into .fch(k) file
subroutine extract_noon2fch(outname, fchname, idx1, idx2, nopen, gau_order)
 implicit none
 integer :: i, fid1, fid2, nmo, nif, itype, RENAME
 integer, intent(in) :: idx1, idx2, nopen
 integer, parameter :: iout = 6
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: outname, fchname
 real(kind=8), allocatable :: noon(:), ev(:)
 logical, intent(in) :: gau_order

 nmo = idx2 - idx1 + 1   ! Note: here nmo is 2*npair + nopen for GVB
 allocate(noon(nmo))

 if(index(outname,'.dat',back=.true.) /= 0) then
  call read_noon_from_dat(nmo, noon, outname, nopen, gau_order) ! GVB NOONs
 else if(index(outname,'.gms',back=.true.) /= 0) then
  call read_noon_from_gmsgms(idx1, nmo, noon, outname) ! CASCI/CASSCF NOONs
 else
  call identify_itype_of_out(outname, itype)
  select case(itype)
  case(1)
   call read_noon_from_pyout(nmo, noon, outname)
  case(2)
   call read_noon_from_orca_out(nmo, noon, outname)
  case(3)
   call read_noon_from_psi4_out(nmo, noon, outname)
  case default
   write(iout,'(A,I0)') 'ERROR in subroutine extract_noon2fch: invalid itype=',itype
   write(iout,'(A)') TRIM(outname)
   stop
  end select
 end if

 ! read nif in file fchname, and copy into file fchname1 by the way
 open(newunit=fid1,file=TRIM(fchname),status='old',position='rewind')
 i = index(fchname, '.fch')
 fchname1 = fchname(1:i-1)//'_noon.fch'
 open(newunit=fid2,file=TRIM(fchname1),status='replace')
 do while(.true.)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  if(buf(1:13) == 'Number of ind') exit
 end do
 BACKSPACE(fid1)
 read(fid1,'(A49,2X,I10)') buf, nif

 ! process to 'Alpha Orbital Energies'
 do while(.true.)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  if(buf(1:8) == 'Alpha Or') exit
 end do

 ! print NOONs in to file fchname1
 allocate(ev(nif), source=0d0)
 if(idx1 > 1) ev(1:idx1-1) = 2d0
 ev(idx1:idx2) = noon
 deallocate(noon)
 write(fid2,'(5(1X,ES15.8))') (ev(i),i=1,nif)
 deallocate(ev)

 ! skip the 'Alpha Orbital Energies' in file fchname
 do while(.true.)
  read(fid1,'(A)') buf
  if(index(buf,'=') /= 0) exit
 end do

 ! copy remaining content
 BACKSPACE(fid1)
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
end subroutine extract_noon2fch

! read GVB NOONs from .dat file of GAMESS
subroutine read_noon_from_dat(nmo, noon, datname, nopen, gau_order)
 implicit none
 integer :: i, j, k, m, npair, fid
 integer, intent(in) :: nmo ! must be an even integer
 integer, intent(in) :: nopen ! number of singly occupied orbitals
 integer, parameter :: iout = 6
 real(kind=8), allocatable :: cicoeff(:)
 real(kind=8), intent(out) :: noon(nmo)
 character(len=240) :: buf
 character(len=240), intent(in) :: datname
 logical, intent(in) :: gau_order

 npair = (nmo - nopen)/2
 if(npair == 0) then
  write(iout,'(A)') 'Warning in subroutine read_noon_from_dat: npair=0. High&
                   & spin ROHF wfn assumed.'
  noon = 1d0
  return
 end if
 allocate(cicoeff(npair), source=0d0)

 open(newunit=fid,file=TRIM(datname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(index(buf,'$SCF') /= 0) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_noon_from_dat: no '$SCF'&
                   & found in file "//TRIM(datname)
  close(fid)
  stop
 end if
 k = 1
 if(index(buf,'CICOEF') /= 0) then
  i = index(buf,'='); j = index(buf,',')
  read(buf(i+1:j-1),*) cicoeff(1)
  k = k + 1
 end if

 do i = k, npair, 1
  read(fid,'(A)') buf
  j = index(buf,'='); m = index(buf,',')
  read(buf(j+1:m-1),*) cicoeff(i)
 end do ! for i
 close(fid)

 ! turn GVB CICOEFF into NOONs
 forall(i = 1:npair) cicoeff(i) = 2d0*cicoeff(i)*cicoeff(i)

 if(gau_order) then
  forall(i = 1:npair) noon(i) = cicoeff(i)
  if(nopen > 0) noon(npair+1:npair+nopen) = 1d0
  forall(i = 1:npair) noon(nmo-i+1) = 2d0 - cicoeff(i)
 else ! in GAMESS order
  if(nopen > 0) noon(1:nopen) = 1d0
  forall(i = 1:npair)
   noon(nopen+2*i-1) = cicoeff(i)
   noon(nopen+2*i) = 2d0 - cicoeff(i)
  end forall
 end if

 deallocate(cicoeff)
 return
end subroutine read_noon_from_dat

! read CASSCF NOONs from .gms file of GAMESS
! Note: CASCI NOONs are recorded in .dat file and no need to handle here
subroutine read_noon_from_gmsgms(idx1, nmo, noon, gmsname)
 implicit none
 integer :: i, j, k, fid, n, nmo1
 integer, intent(in) :: idx1, nmo
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: noon(nmo)
 real(kind=8), allocatable :: on(:)
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname

 open(newunit=fid,file=TRIM(gmsname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(6:13) == 'ATOMIC M') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_noon_from_gmsgms: no 'ATOMIC M'&
                   & found in file "//TRIM(gmsname)
  stop
 end if

 read(fid,'(A)') buf ! skip one line

 nmo1 = idx1 - 1 + nmo
 n = nmo1/5
 if(nmo1-5*n > 0) n = n + 1
 allocate(on(nmo1), source=0d0)

 do i = 1, n, 1
  do j = 1, 3
   read(fid,'(A)') buf
  end do ! for j
  k = min(5*i,nmo1)
  read(buf,*) on(5*i-4:k)

  read(fid,'(A)') buf
  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit
  end do ! for while
 end do ! for i

 close(fid)
 noon = on(idx1:)
 return
end subroutine read_noon_from_gmsgms

! read NOONs from .out file of PySCF
subroutine read_noon_from_pyout(nmo, noon, outname)
 implicit none
 integer :: i, j, k, nline, fid
 integer, intent(in) :: nmo
 real(kind=8), intent(out) :: noon(nmo)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: alive

 ! if the .out file is still opened, close it
 inquire(file=TRIM(outname),opened=alive)
 if(alive) then
  inquire(file=TRIM(outname),number=fid)
  close(fid)
 end if

 ! Note: always read the last 'Natural occ'
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(index(buf,'Natural occ') /= 0) exit
 end do ! for while

 j = index(buf,'['); k = index(buf,']')
 if(k == 0) then
  read(buf(j+1:),*) noon(1:6)
 else
  read(buf(j+1:k-1),*) (noon(i), i=1, min(6,nmo))
 end if

 if(nmo > 6) then
  nline = nmo/6 - 1
  do i = 1, nline, 1
   read(fid,'(A)') buf
   j = index(buf,']')  ! in case nmo is 6N
   if(j == 0) j = LEN_TRIM(buf) + 1
   read(buf(1:j-1),*) noon(6*i+1 : 6*i+6)
  end do ! for i

  k = nmo - (nline+1)*6
  if(k > 0) then
   read(fid,'(A)') buf
   j = index(buf,']')
   read(buf(1:j-1),*) (noon(i), i=nmo-k+1, nmo)
  end if
 end if

 close(fid)
 return
end subroutine read_noon_from_pyout

! read CASCI/CASSCF NOONs from ORCA .out file
subroutine read_noon_from_orca_out(nmo, noon, outname)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nmo
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: noon(nmo)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 noon = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(4:10) == 'N(occ)=') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_noon_from_orca_out: no ''&
                   & found in file "//TRIM(outname)//'.'
  close(fid)
  stop
 end if
 close(fid)

 i = index(buf,'=')
 read(buf(i+1:),*) (noon(i),i=1,nmo)
 return
end subroutine read_noon_from_orca_out

! read CASCI/CASSCF NOONs from a PSI4 .out file
subroutine read_noon_from_psi4_out(nmo, noon, outname)
 implicit none
 integer :: i, nline, fid
 integer, intent(in) :: nmo
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: noon(nmo)
 character(len=1) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 noon = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(4:19) == 'Active Space Nat') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_noon_from_psi4_out: no 'Active&
                   & Space' found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 nline = nmo/3

 do i = 1, nline, 1
  read(fid,*) str,noon(3*i-2),str,noon(3*i-1),str,noon(3*i)
 end do ! for i

 if(nmo - 3*nline > 0) then
  read(fid,'(A)') buf
  i = index(buf,'A'); buf(i:i) = ' '
  i = index(buf,'A')
  if(i > 0) buf(i:i) = ' '
  read(buf,*) noon(3*nline+1:nmo)
 end if

 close(fid)
end subroutine read_noon_from_psi4_out

! identify whether this is a ORCA or PySCF output file
subroutine identify_itype_of_out(outname, itype)
 implicit none
 integer :: i, fid
 integer, intent(out) :: itype ! 1/2/3 for PySCF/ORCA/PSI4
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 itype = 1
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do i = 1, 5, 1
  read(fid,'(A)') buf

  if(index(buf,'O   R   C   A') > 0) then
   itype = 2
   exit
  else if(buf(11:23) == 'Psi4: An Open') then
   itype = 3
   exit
  end if
 end do ! for i

 return
end subroutine identify_itype_of_out

