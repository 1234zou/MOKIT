! written by jxzou at 20200503: extract Natural Orbital Occupation Numbers (NOONs)
!  from (1) .out file of PySCF or (2) .dat file of GAMESS, and write into a given .fch file

! updated by jxzou at 20200517: fix the bug(nmo=6*N) in subroutine read_noon_from_pyout
! updated by jxzou at 20201212: modify read_noon_from_gmsgms to read more digits of NOONs

! Note: [idx1,idx2] contains the singly occupied MOs and GVB pairs, npair=(idx2-idx1+1-nopen)/2

program main
 implicit none
 integer :: i, j, idx1, idx2, nopen
 character(len=10) :: buf
 character(len=240) :: outname, fchname
 logical :: gau_order

 i = iargc()
 if(i<2 .or. i>6) then
  write(6,'(/,A)') ' ERROR in subroutine extract_noon2fch: wrong command line arguments.'
  write(6,'(A)')   ' Syntax: extract_noon_2fch outname fchname idx1 idx2 nopen [-gau]'
  write(6,'(/,A)') ' Example 1(PySCF CASCI): extract_noon2fch a.out a.fch 19 24'
  write(6,'(A)')   ' Example 2 (GAMESS GVB): extract_noon2fch a.dat a.fch 19 24 0'
  write(6,'(A)')   ' Example 3 (GAMESS GVB): extract_noon2fch a.dat a.fch 19 24 0 -gau'
  write(6,'(A)')   ' Example 4 (GAMESS CAS): extract_noon2fch a.gms a.fch 19 24'
  write(6,'(A)')   ' Example 5 (GAMESS MP2): extract_noon2fch a.gms a.fch -mp2'
  write(6,'(A)')   ' Example 6 (ORCA CAS)  : extract_noon2fch a.out a.fch 19 24'
  write(6,'(A)')   ' Example 7 (PSI4 CAS)  : extract_noon2fch a.out a.fch 19 24'
  write(6,'(A,/)') ' Example 8 (Q-Chem GVB): extract_noon2fch a.out a.fch'
  stop
 end if

 nopen = 0
 gau_order = .false. ! initialization

 call getarg(1, outname)
 call require_file_exist(outname)

 call getarg(2, fchname)
 call require_file_exist(fchname)

 select case(i)
 case(2)
  ! Q-Chem case, fake two integers
  idx1 = 1; idx2 = 2
 case(3)
  ! GAMESS MP2 NOs case. Here we use the variable gau_order=.T. to represent the
  ! MP2 NOs case
  gau_order = .true.
  idx1 = 1
  call read_nif_from_fch(fchname, idx2)
 case default
  call getarg(3, buf)
  read(buf,*) idx1
  call getarg(4, buf)
  read(buf,*) idx2
 end select

 if(i > 4) then
  call getarg(5, buf)
  read(buf,*,iostat=j) nopen
  if(j /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine extract_noon2fch: wrong command line a&
                    &rguments. Failed'
   write(6,'(A)') 'to read integer nopen.'
   stop
  end if
  if(nopen < 0) then
   write(6,'(/,A)') 'ERROR in subroutine extract_noon2fch: wrong command line a&
                    &rguments. nopen<0.'
   stop
  end if
 end if

 if(i == 6) then
  call getarg(6, buf)
  if(INDEX(buf,'-gau') /= 0) then
   gau_order = .true.
  else
   write(6,'(/,A)') 'ERROR in subroutine extract_noon2fch: wrong command line a&
                    &rguments.'
   write(6,'(A)') 'The 5th argument='//TRIM(buf)
   stop
  end if
 end if

 j = idx2 - idx1
 if(j < 1) then
  write(6,'(/,A)') 'ERROR in subroutine extract_noon2fch: wrong input indices.'
  write(6,'(A,3I5)') 'idx1, idx2, nopen=', idx1, idx2, nopen
  stop
 end if

 if(INDEX(outname,'.dat',back=.true.)>0 .and. i/=3) then
  if(MOD(j+1-nopen,2) /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine extract_noon2fch: wrong input indices.&
                   & In this case'
   write(6,'(A)') 'idx2-idx1+1-nopen must be an even integer.'
   write(6,'(A,3I5)') 'idx1, idx2, nopen=', idx1, idx2, nopen
   stop
  end if
 end if

 call extract_noon2fch(outname, fchname, idx1, idx2, nopen, gau_order)
end program main

! extract NOONs from PySCF .out file, and print it into .fch(k) file
subroutine extract_noon2fch(outname, fchname, idx1, idx2, nopen, gau_order)
 implicit none
 integer :: i, fid1, fid2, nfocc, nmo, nif, itype, RENAME
 integer, intent(in) :: idx1, idx2, nopen
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: outname, fchname
 real(kind=8), allocatable :: noon(:), ev(:)
 logical, intent(in) :: gau_order

 nmo = idx2 - idx1 + 1 ! Note: here nmo is 2*npair + nopen for GVB
 allocate(noon(nmo))

 if(INDEX(outname,'.dat',back=.true.) > 0) then
  call read_noon_from_dat(outname, nmo, noon, nopen, gau_order) ! GVB NOONs
 else if(INDEX(outname,'.gms',back=.true.) > 0) then
  if(gau_order) then ! MP2 NOONs
   call read_mp2_on_from_gmsgms(outname, nmo, noon)
  else               ! CASCI/CASSCF NOONs
   call read_noon_from_gmsgms(outname, idx1, nmo, noon)
  end if
 else
  call identify_itype_of_out(outname, itype)
  select case(itype)
  case(1)
   call read_noon_from_pyout(outname, nmo, noon)
  case(2)
   call read_noon_from_orca_out(outname, nmo, noon)
  case(3)
   call read_noon_from_psi4_out(outname, nmo, noon)
  case(4)
   deallocate(noon)
   call read_nfocc_nif_from_qchem_gvb_out(outname, nfocc, nmo)
   allocate(noon(nmo), source=0d0)
   call read_noon_from_qchem_out(outname, nfocc, noon(1:nfocc))
  case default
   write(6,'(/,A,I0)') 'ERROR in subroutine extract_noon2fch: invalid itype=',&
                        itype
   write(6,'(A)') TRIM(outname)
   stop
  end select
 end if

 call read_nif_from_fch(fchname, nif)

 fchname1 = TRIM(fchname)//'.t'
 open(newunit=fid1,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(fchname1),status='replace')

 ! process to 'Alpha Orbital Energies'
 do while(.true.)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  if(buf(1:8) == 'Alpha Or') exit
 end do

 ! print NOONs in to file fchname1
 if(itype == 4) then ! Q-Chem case
  allocate(ev(nif), source=noon)
 else
  allocate(ev(nif), source=0d0)
  if(idx1 > 1) ev(1:idx1-1) = 2d0
  ev(idx1:idx2) = noon
 end if
 deallocate(noon)
 write(fid2,'(5(1X,ES15.8))') (ev(i),i=1,nif)
 deallocate(ev)

 ! skip the 'Alpha Orbital Energies' in file fchname
 do while(.true.)
  read(fid1,'(A)') buf
  if(INDEX(buf,'=') > 0) exit
 end do ! for while

 ! copy remaining content
 BACKSPACE(fid1)
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
end subroutine extract_noon2fch

! read GVB NOONs from .dat file of GAMESS
subroutine read_noon_from_dat(datname, nmo, noon, nopen, gau_order)
 implicit none
 integer :: i, j, k, m, npair, fid
 integer, intent(in) :: nmo   ! must be an even integer
 integer, intent(in) :: nopen ! number of singly occupied orbitals
 real(kind=8), allocatable :: cicoeff(:)
 real(kind=8), intent(out) :: noon(nmo)
 character(len=240) :: buf
 character(len=240), intent(in) :: datname
 logical, intent(in) :: gau_order

 npair = (nmo - nopen)/2
 if(npair == 0) then
  write(6,'(/,A)') 'Warning in subroutine read_noon_from_dat: npair=0. High spi&
                   &n ROHF wfn assumed.'
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
  write(6,'(/,A)') "ERROR in subroutine read_noon_from_dat: no '$SCF' found in &
                   &file "//TRIM(datname)
  close(fid)
  stop
 end if

 k = 1
 if(index(buf,'CICOEF') /= 0) then
  i = INDEX(buf,'='); j = INDEX(buf,',')
  read(buf(i+1:j-1),*) cicoeff(1)
  k = k + 1
 end if

 do i = k, npair, 1
  read(fid,'(A)') buf
  j = INDEX(buf,'='); m = INDEX(buf,',')
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
end subroutine read_noon_from_dat

! read CASSCF NOONs from .gms file of GAMESS
! Note: CASCI NOONs are recorded in .dat file and no need to be handled here
subroutine read_noon_from_gmsgms(gmsname, idx1, nmo, noon)
 implicit none
 integer :: i, j, k, fid, n, nmo1
 integer, intent(in) :: idx1, nmo
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
  write(6,'(/,A)') "ERROR in subroutine read_noon_from_gmsgms: no 'ATOMIC M' fo&
                   &und in file "//TRIM(gmsname)
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
end subroutine read_noon_from_gmsgms

! read MP2 NOONs from a GAMESS .gms file
subroutine read_mp2_on_from_gmsgms(gmsname, nif, noon)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nif
 real(kind=8), intent(out) :: noon(nif)
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname

 noon = 0d0
 open(newunit=fid,file=TRIM(gmsname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(2:8) == 'MP2 NAT') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mp2_on_from_gmsgms: no 'MP2 NAT' f&
                   &ound in file "//TRIM(gmsname)
  close(fid)
  stop
 end if

 read(fid,'(5(1X,ES15.8))') noon
 close(fid)
end subroutine read_mp2_on_from_gmsgms

! read NOONs from .out file of PySCF
subroutine read_noon_from_pyout(outname, nmo, noon)
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

 j = INDEX(buf,'['); k = INDEX(buf,']')
 if(k == 0) then
  read(buf(j+1:),*) noon(1:6)
 else
  read(buf(j+1:k-1),*) (noon(i), i=1, min(6,nmo))
 end if

 if(nmo > 6) then
  nline = nmo/6 - 1
  do i = 1, nline, 1
   read(fid,'(A)') buf
   j = INDEX(buf,']')  ! in case nmo is 6N
   if(j == 0) j = LEN_TRIM(buf) + 1
   read(buf(1:j-1),*) noon(6*i+1 : 6*i+6)
  end do ! for i

  k = nmo - (nline+1)*6
  if(k > 0) then
   read(fid,'(A)') buf
   j = INDEX(buf,']')
   read(buf(1:j-1),*) (noon(i), i=nmo-k+1, nmo)
  end if
 end if

 close(fid)
end subroutine read_noon_from_pyout

! read CASCI/CASSCF NOONs from ORCA .out file
subroutine read_noon_from_orca_out(outname, nmo, noon)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nmo
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
  write(6,'(A)') "ERROR in subroutine read_noon_from_orca_out: no ''&
                   & found in file "//TRIM(outname)//'.'
  close(fid)
  stop
 end if
 close(fid)

 i = INDEX(buf,'=')
 read(buf(i+1:),*) (noon(i),i=1,nmo)
end subroutine read_noon_from_orca_out

! read CASCI/CASSCF NOONs from a PSI4 .out file
subroutine read_noon_from_psi4_out(outname, nmo, noon)
 implicit none
 integer :: i, nline, fid
 integer, intent(in) :: nmo
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
  write(6,'(A)') "ERROR in subroutine read_noon_from_psi4_out: no 'Active&
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
  i = INDEX(buf,'A'); buf(i:i) = ' '
  i = INDEX(buf,'A')
  if(i > 0) buf(i:i) = ' '
  read(buf,*) noon(3*nline+1:nmo)
 end if

 close(fid)
end subroutine read_noon_from_psi4_out

! identify whether this is a ORCA or PySCF output file
subroutine identify_itype_of_out(outname, itype)
 implicit none
 integer :: i, fid
 integer, intent(out) :: itype ! 1/2/3/4 for PySCF/ORCA/PSI4/Q-Chem
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
  else if(buf(1:5) == 'qchem') then
   itype = 4
   exit
  end if
 end do ! for i
end subroutine identify_itype_of_out

! read nfocc and nif from a Q-Chem output file
! nfocc := doubly occupied + singly occupied + 2*npair
subroutine read_nfocc_nif_from_qchem_gvb_out(outname, nfocc, nif)
 implicit none
 integer :: i, fid, ndb, npair, nopen, nvir
 integer, intent(out) :: nfocc, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ndb = 0; npair = 0; nopen = 0; nvir = 0; nfocc = 0; nif = 0
 buf = ' '
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:7) == 'Alpha:') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_nfocc_nif_from_qchem_gvb_out: no 'Al&
                 &pha:' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(buf(8:),*) ndb
 i = INDEX(buf,',')
 read(buf(i+1:),*) npair
 read(fid,*) nvir
 close(fid)
 nfocc = ndb + nopen + 2*npair
 nif = nfocc + nvir
end subroutine read_nfocc_nif_from_qchem_gvb_out

! read NOONs from a Q-Chem output file
subroutine read_noon_from_qchem_out(outname, nfocc, noon)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nfocc
 real(kind=8) :: r1, r2
 real(kind=8), intent(out) :: noon(nfocc)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 noon = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(4:14) == 'Orbital occ') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_noon_from_qchem_gvb_out: no 'Orbit&
                   &al occ' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, nfocc, 1
  read(fid,*) j, r1, r2
  noon(i) = r1 + r2
 end do ! for i

 close(fid)
end subroutine read_noon_from_qchem_out

! read the number of MOs from a Gaussian .fch(k) file
subroutine read_nif_from_fch(fchname, nif)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nif
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 nif = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:13) == 'Number of ind') exit
 end do

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nif_from_fch: no 'Number of ind' f&
                   &ound in file "//TRIM(fchname)
  stop
 end if

 read(buf(45:),*) nif
end subroutine read_nif_from_fch

