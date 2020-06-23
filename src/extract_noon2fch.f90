! written by jxzou at 20200503: extract Natural Orbital Occupation Numbers (NOONs)
!  from (1) .out file of PySCF or (2) .dat file of GAMESS, and write into a given .fch file

! updated by jxzou at 20200517: fix the bug(nmo=6*N) in subroutine read_noon_from_pyout

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
  write(iout,'(/,A)') ' ERROR in subroutine extract_noon2fch: wrong&
                      & command line arguments.'
  write(iout,'(A)') ' Format: extract_noon_2fch outname fchname idx1 idx2 nopen [-gau]'
  write(iout,'(/,A)') ' Example1(PySCF CASCI): extract_noon2fch a.out a.fch 19 24'
  write(iout,'(/,A)') ' Example2 (GAMESS GVB): extract_noon2fch a.dat a.fch 19 24 0'
  write(iout,'(/,A)') ' Example3 (GAMESS GVB): extract_noon2fch a.dat a.fch 19 24 0 -gau'
  write(iout,'(/,A)') ' Example4 (GAMESS CAS): extract_noon2fch a.gms a.fch 19 24'
  stop
 end if

 call getarg(1, outname)
 call getarg(2, fchname)
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
 integer :: i, fid1, fid2, nmo, nif
 integer :: RENAME
 integer, intent(in) :: idx1, idx2, nopen
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: outname, fchname
 real(kind=8), allocatable :: noon(:), ev(:)
 logical, intent(in) :: gau_order

 nmo = idx2 - idx1 + 1
 ! Note: here nmo is 2*npair + nopen for GVB
 allocate(noon(nmo), source=0.0d0)

 if(index(outname,'.dat',back=.true.) /= 0) then
  call read_noon_from_dat(nmo, noon, outname, nopen, gau_order) ! GVB NOONs
 else if(index(outname,'.gms',back=.true.) /= 0) then
  call read_noon_from_gmsgms(idx1, nmo, noon, outname) ! CASCI/CASSCF NOONs
 else
  call read_noon_from_pyout(nmo, noon, outname) ! CASCI/CASSCF NOONs
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
 allocate(ev(nif), source=0.0d0)
 if(idx1 > 1) ev(1:idx1-1) = 2.0d0
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
 return
end subroutine extract_noon2fch

! read GVB NOONs from .dat file of GAMESS
subroutine read_noon_from_dat(nmo, noon, datname, nopen, gau_order)
 implicit none
 integer :: i, j, k, m, npair, fid
 integer, intent(in) :: nmo ! must be an even integer
 integer, intent(in) :: nopen ! number of singly occupied orbitals
 real(kind=8), allocatable :: cicoeff(:)
 real(kind=8), intent(out) :: noon(nmo)
 character(len=240) :: buf
 character(len=240), intent(in) :: datname
 logical, intent(in) :: gau_order

 npair = (nmo - nopen)/2
 allocate(cicoeff(npair), source=0.0d0)

 open(newunit=fid,file=TRIM(datname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(index(buf,'SCF') /= 0) exit
 end do ! for while

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
 forall(i = 1:npair) cicoeff(i) = 2.0d0*cicoeff(i)*cicoeff(i)

 if(gau_order) then
  forall(i = 1:npair) noon(i) = cicoeff(i)
  if(nopen > 0) noon(npair+1:npair+nopen) = 1.0d0
  forall(i = 1:npair) noon(nmo-i+1) = 2.0d0 - cicoeff(i)
 else ! in GAMESS order
  if(nopen > 0) noon(1:nopen) = 1.0d0
  forall(i = 1:npair)
   noon(nopen+2*i-1) = cicoeff(i)
   noon(nopen+2*i) = 2.0d0 - cicoeff(i)
  end forall
 end if

 deallocate(cicoeff)
 return
end subroutine read_noon_from_dat

! read NOONs from .gms file of GAMESS
subroutine read_noon_from_gmsgms(idx1, nmo, noon, gmsname) ! CASCI/CASSCF NOONs
 implicit none
 integer :: i, k, fid, n, nmo1
 integer, intent(in) :: idx1, nmo
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: noon(nmo)
 real(kind=8), allocatable :: on(:)
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname
 logical :: casci

 casci = .true.
 open(newunit=fid,file=TRIM(gmsname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(11:26) == 'NATURAL ORBITALS') exit       ! CASCI
  if(buf(11:32) == 'MCSCF NATURAL ORBITALS') then ! CASSCF
   casci = .false.
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_noon_from_gmsgms: no 'NATURAL&
   & ORBITALS' or 'MCSCF NATURAL ORBITALS' found in file "//TRIM(gmsname)
  stop
 end if

 read(fid,'(A)') buf ! skip one line

 nmo1 = idx1 - 1 + nmo
 n = nmo1/5
 if(nmo1-5*n > 0) n = n + 1
 allocate(on(nmo1), source=0.0d0)

 do i = 1, n, 1
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  k = min(5*i,nmo1)
  read(buf,*) on(5*i-4:k)

  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit
  end do
  BACKSPACE(fid)
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

