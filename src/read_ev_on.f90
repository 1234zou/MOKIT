! Create by jxzou at 20250522: moved subroutines about reading eigenvalues
!  and occupation numbers from rwwfn.f90 to this file.

! Read Alpha/Beta eigenvalues in a given .fch(k) file.
! The Alpha/Beta Orbital Energies in .fch(k) file can be either energy levels
!  or NOONs, depending on the job type.
subroutine read_eigenvalues_from_fch(fchname, nif, ab, noon)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nif
!f2py intent(in) :: nif
 real(kind=8), intent(out) :: noon(nif)
!f2py intent(out) :: noon
!f2py depend(nif) :: noon
 character(len=1), intent(in) :: ab
!f2py intent(in) :: ab
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 character(len=8) :: key
 character(len=8), parameter :: key1 = 'Alpha Or'
 character(len=7), parameter :: key2 = 'Beta Or'

 noon = 0d0
 key = key1
 if(ab/='a' .and. ab/='A') key = key2//'b'

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == key) exit
 end do

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_eigenvalues_from_fch: no '"//key//&
                   "' found in "
  write(6,'(A)') 'file '//TRIM(fchname)//'.'
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, i
 if(i /= nif) then
  write(6,'(/,A)') 'ERROR in subroutine read_eigenvalues_from_fch: i /= nif.'
  write(6,'(A)') 'Inconsistency found between input nif and that in file '//TRIM(fchname)
  stop
 end if

 read(fid,'(5(1X,ES15.8))') (noon(i),i=1,nif)
 close(fid)
end subroutine read_eigenvalues_from_fch

! read orbital energies from an ORCA .mkl file
subroutine read_ev_from_mkl(mklname, nmo, ab, ev)
 implicit none
 integer :: i, k, rc, nline, fid
 integer, intent(in) :: nmo
!f2py intent(in) :: nmo
 real(kind=8), intent(out) :: ev(nmo)
!f2py intent(out) :: ev
!f2py depend(nmo) :: ev
 character(len=1), intent(in) :: ab
!f2py intent(in) :: ab
 character(len=3) :: str = ' '
 character(len=8) :: key = ' '
 character(len=8), parameter :: key1 = '$COEFF_A'
 character(len=8), parameter :: key2 = '$COEFF_B'
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname
!f2py intent(in) :: mklname

 ev = 0d0
 if(ab=='b' .or. ab=='B') then
  key = key2
 else
  key = key1
 end if

 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == key) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_ev_from_mkl: no '"//key//"' found &
                   &in file "//TRIM(mklname)
  close(fid)
  stop
 end if

 nline = nmo/5
 if(nmo-5*nline > 0) nline = nline + 1
 read(fid,'(A)') buf

 do i = 1, nline, 1
  k = min(5*i,nmo)
  read(unit=fid,fmt=*,iostat=rc) ev(5*i-4:k)
  if(rc /= 0) exit
  if(i == nline) exit

  do while(.true.)
   read(fid,'(A)',iostat=rc) buf
   if(rc /= 0) exit
   read(buf,*,iostat=rc) str
   if(rc == 0) then
   if(str=='a1g' .or. str=='A1g' .or. str=='A1G') exit
   end if
  end do ! for while

  if(rc /= 0) exit
 end do ! for i

 close(fid)
 if(rc /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_ev_from_mkl: incomplete .mkl file.'
  write(6,'(A)') 'Filename='//TRIM(mklname)
  stop
 end if
end subroutine read_ev_from_mkl

! read Alpha/Beta eigenvalues (orbital energies) from a given orbital file of BDF
subroutine read_ev_from_bdf_orb(orbname, nif, ab, ev)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nif
!f2py intent(in) :: nif
 real(kind=8), intent(out) :: ev(nif)
!f2py intent(out) :: ev
!f2py depend(nif) :: ev
 character(len=1), intent(in) :: ab
 character(len=240) :: buf
 character(len=240), intent(in) :: orbname
!f2py intent(in) :: orbname

 ev = 0d0
 open(newunit=fid,file=TRIM(orbname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:9) == 'ORBITAL E') exit
 end do

 if(i /= 0) then
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') "Warning in subroutine read_ev_from_bdf_orb: 'ORBITAL E' not f&
                 &ound in file "//TRIM(orbname)//'.'
  write(6,'(A)') 'Failed to read orbital energies. Set all orbital energies to &
                 &be zero. This does not'
  write(6,'(A)') 'affect subsequent computations in AutoMR of MOKIT. So continu&
                 &e. If you know your'
  write(6,'(A)') 'subsequent computations will explicitly use orbital energies &
                 &in fch(k) file, please'
  write(6,'(A)') 'stop subsequent computations.'
  write(6,'(A)') REPEAT('-',79)
  close(fid)
  return
 end if

 read(fid,*) (ev(i),i=1,nif)
 if(ab=='b' .or. ab=='B') read(fid,*) (ev(i),i=1,nif)
 close(fid)
end subroutine read_ev_from_bdf_orb

! read Alpha/Beta orbital energies from a specified .amo file
subroutine read_ev_from_amo(amoname, nif, ab, ev)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nif
!f2py intent(in) :: nif
 real(kind=8), intent(out) :: ev(nif)
!f2py intent(out) :: ev
!f2py depend(nif) :: ev
 character(len=1), intent(in) :: ab
!f2py intent(in) :: ab
 character(len=10) :: buf ! long buf is not needed here
 character(len=240), intent(in) :: amoname
!f2py intent(in) :: amoname

 ev = 0d0
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')

 select case(ab)
 case('a') ! Alpha MOs
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:3)=='En:' .or. buf(1:4)=='EnA:') exit
  end do ! for while
 case('b') ! Beta MOs
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:4) == 'EnB:') exit
  end do ! for while
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_ev_from_amo: MO type cannot be rec&
                   &ognized. ab='//ab
  close(fid)
  stop
 end select

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_ev_from_amo: specified orbital eig&
                   &envalue section'
  write(6,'(A)') 'not found in file '//TRIM(amoname)
  close(fid)
  stop
 end if

 read(fid,*,iostat=i) ev
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_ev_from_amo: failed to read orbita&
                   &l eigenvalues.'
  write(6,'(A)') 'This file seems problematic: '//TRIM(amoname)
  stop
 end if
end subroutine read_ev_from_amo

! read Alpha/Beta occupation numbers from a given .orb file of (Open)Molcas
subroutine read_on_from_orb(orbname, nif, ab, on)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nif
!f2py intent(in) :: nif
 real(kind=8), intent(out) :: on(nif)
!f2py intent(out) :: on
!f2py depend(nif) :: on
 character(len=1), intent(in) :: ab
!f2py intent(in) :: ab
 character(len=240) :: buf
 character(len=240), intent(in) :: orbname
!f2py intent(in) :: orbname
 character(len=5) :: key
 character(len=4), parameter :: key1 = '#OCC'
 character(len=5), parameter :: key2 = '#UOCC'

 on = 0d0
 key = key1//' '
 if(ab/='a' .and. ab/='A') key = key2

 open(newunit=fid,file=TRIM(orbname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,TRIM(key)) /= 0) exit
 end do

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_on_from_orb: no '"//TRIM(key)//&
                  &"' found in"
  write(6,'(A)') 'file '//TRIM(orbname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(fid,'(5(1X,ES21.14))') (on(i),i=1,nif)
 close(fid)
end subroutine read_on_from_orb

! read $OCCNO from a given GAMESS .dat file
subroutine read_on_from_gms_dat(datname, nmo, on, alive)
 implicit none
 integer :: i, k, nline, fid
 integer, intent(in) :: nmo
!f2py intent(in) :: nmo
 real(kind=8), intent(out) :: on(nmo)
!f2py intent(out) :: on
!f2py depend(nmo) :: on
 character(len=240) :: buf
 character(len=240), intent(in) :: datname
!f2py intent(in) :: datname
 logical, intent(out) :: alive ! whether $OCCNO exists

 on = 0d0; alive = .false.
 open(newunit=fid,file=TRIM(datname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:7) == '$OCCNO') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  return
 end if

 alive = .true.
 nline = nmo/5
 if(nmo-5*nline > 0) nline = nline + 1
 do i = 1, nline, 1
  k = min(5*i,nmo)
  read(fid,*) on(5*i-4:k)
 end do ! for i

 close(fid)
end subroutine read_on_from_gms_dat

! read occupation numbers from Molpro .xml file
subroutine read_on_from_xml(xmlname, nmo, ab, on)
 implicit none
 integer :: i, j, k, nline, fid
 integer, intent(in) :: nmo
!f2py intent(in) :: nmo
 real(kind=8), intent(out) :: on(nmo)
!f2py intent(out) :: on
!f2py depend(nmo) :: on
 character(len=240) :: buf
 character(len=1), intent(in) :: ab
!f2py intent(in) :: ab
 character(len=240), intent(in) :: xmlname
!f2py intent(in) :: xmlname

 on = 0d0
 nline = nmo/10
 if(nmo-nline*10 > 0) nline = nline + 1
 open(newunit=fid,file=TRIM(xmlname),status='old',position='rewind')

 if(ab=='a' .or. ab=='A') then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(INDEX(buf,"type=""ALPH")/=0 .or. INDEX(buf,"type=""CANO")/=0 .or. &
      INDEX(buf,"type=""NATU")/=0) exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_on_from_xml: none of 'ALPH', 'CAN&
                    &O' or 'NATU' is"
   write(6,'(A)') 'found in file '//TRIM(xmlname)
   close(fid)
   stop
  end if

 else ! UHF

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(INDEX(buf,"type=""BETA") /= 0) exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_on_from_xml: no type=""BETA found&
                    & in file "//TRIM(xmlname)
   close(fid)
   stop
  end if
 end if

 read(fid,'(A)') buf

 do i = 1, nmo, 1
  !read(fid,'(A)',iostat=k) buf
  !if(k /= 0) exit
  !j = INDEX(buf, """")
  !k = INDEX(buf(j+1:), """")
  !read(buf(j+1:j+k-1),*) on(i)
  do while(.true.)
   read(fid,'(A)') buf
   j = INDEX(buf, "occupation=""")
   if(j > 0) then
    k = INDEX(buf(j+12:), """")
    read(buf(j+12:j+k+10),*) on(i)
   end if
   k = LEN_TRIM(buf)
   if(buf(k:k) == '>') exit
  end do ! for while

  do j = 1, nline+1, 1
   read(fid,'(A)',iostat=k) buf
   if(k /= 0) exit
  end do ! for j
 end do ! for i

 close(fid)
 if(k /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_on_from_xml: not all occupation nu&
                   &mbers are found.'
  write(6,'(A)') "Did you forget to add '{put,xml}' in Molpro input file?"
  write(6,'(A)') 'Filname = '//TRIM(xmlname)
  stop
 end if
end subroutine read_on_from_xml

! read Alpha/Beta occupation numbers from a given orbital file of BDF
subroutine read_on_from_bdf_orb(orbname, nif, ab, on)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nif
!f2py intent(in) :: nif
 real(kind=8), intent(out) :: on(nif)
!f2py intent(out) :: on
!f2py depend(nif) :: on
 character(len=1), intent(in) :: ab
!f2py intent(in) :: ab
 character(len=240) :: buf
 character(len=240), intent(in) :: orbname
!f2py intent(in) :: orbname

 on = 0d0
 open(newunit=fid,file=TRIM(orbname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:10) == 'OCCUPATION') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "Warning in subroutine read_on_from_bdf_orb: no 'OCCUPATION' f&
                 &ound in file "//TRIM(orbname)//'.'
  write(6,'(A)') 'Failed to read CAS NOONs. This does not affect the CASCI/CASS&
                 &CF energy. So continue.'
  write(6,'(A)') 'To avoid this error, please use newer version of BDF.'
  close(fid)
  return
 end if

 read(fid,*) (on(i),i=1,nif)

 if(ab=='b' .or. ab=='B') then
  read(fid,'(A)') buf
  read(fid,*) (on(i),i=1,nif)
 end if

 close(fid)
end subroutine read_on_from_bdf_orb

! read CAS NOONs from a DALTON.MOPUN format file
subroutine read_on_from_dalton_mopun(orbname, nif, on)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nif
!f2py intent(in) :: nif
 real(kind=8), intent(out) :: on(nif)
!f2py intent(out) :: on
!f2py depend(nif) :: on
 character(len=240) :: buf
 character(len=240), intent(in) :: orbname
!f2py intent(in) :: orbname

 on = 0d0
 open(newunit=fid,file=TRIM(orbname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == '**NATOCC') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_on_from_dalton_mopun: no '**NATOCC&
                  &' found in file "//TRIM(orbname)//'.'
  close(fid)
  stop
 end if

 read(fid,*) on
 close(fid)
end subroutine read_on_from_dalton_mopun

! read CAS NOONs from a Dalton output file
! Note: only occupation numbers of active orbitals are printed in the Dalton output
subroutine read_on_from_dalton_out(outname, nacto, on)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nacto
!f2py intent(in) :: nacto
 real(kind=8), intent(out) :: on(nacto)
!f2py intent(out) :: on
!f2py depend(nacto) :: on
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname

 on = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:19) == 'Occupancies of nat') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_on_from_dalton_out: 'Occupancies o&
                   &f nat' not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, 4
  read(fid,'(A)') buf
 end do

 read(fid,*) on
 close(fid)
end subroutine read_on_from_dalton_out

! read occupation numbers from ORCA .mkl file
subroutine read_on_from_mkl(mklname, nmo, ab, on)
 implicit none
 integer :: i, k, nline, fid
 integer, intent(in) :: nmo
!f2py intent(in) :: nmo
 real(kind=8), intent(out) :: on(nmo)
!f2py intent(out) :: on
!f2py depend(nmo) :: on
 character(len=1), intent(in) :: ab
!f2py intent(in) :: ab
 character(len=6) :: key = ' '
 character(len=6), parameter :: key1 = '$OCC_A'
 character(len=6), parameter :: key2 = '$OCC_B'
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname
!f2py intent(in) :: mklname

 on = 0d0
 if(ab=='b' .or. ab=='B') then
  key = key2
 else
  key = key1
 end if

 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == key) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_on_from_mkl: no '"//key//"' found in&
                & file "//TRIM(mklname)
  close(fid)
  stop
 end if

 nline = nmo/5
 if(nmo-5*nline > 0) nline = nline + 1

 do i = 1, nline, 1
  k = min(5*i,nmo)
  read(fid,*) on(5*i-4:k)
 end do ! for i

 close(fid)
end subroutine read_on_from_mkl

! read occupation numbers from a specified MOLDEN format file
subroutine read_on_from_molden(molden, nmo, occ)
 implicit none
 integer :: i, j, k, ntag, fid
 integer, intent(in) :: nmo
!f2py intent(in) :: nmo
 real(kind=8), intent(out) :: occ(nmo)
!f2py intent(out) :: occ
!f2py depend(nmo) :: occ
 character(len=240) :: buf
 character(len=240), intent(in) :: molden
!f2py intent(in) :: molden

 occ = 0d0
 call find_ntag_before_mo_in_molden(molden, ntag)
 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4)=='[MO]' .or. buf(2:5)=='[MO]') exit
 end do ! for while

 do i = 1, nmo-1, 1
  do j = 1, ntag, 1
   read(fid,'(A)') buf
   if(buf(1:6)=='Occup=' .or. buf(2:7)=='Occup=') then
    k = INDEX(buf(1:7), '=')
    read(buf(k+1:),*) occ(i)
   end if
  end do ! for j

  do while(.true.)
   read(fid,'(A)') buf
   if(buf(1:3)=='Ene' .or. buf(2:4)=='Ene') then
    BACKSPACE(fid)
    exit
   end if
  end do ! for while
 end do ! for i

 do j = 1, ntag, 1
  read(fid,'(A)') buf
  if(buf(1:6)=='Occup=' .or. buf(2:7)=='Occup=') then
   k = INDEX(buf(1:7), '=')
   read(buf(k+1:),*) occ(nmo)
   exit
  end if
 end do ! for j

 close(fid)
end subroutine read_on_from_molden

! Read the number of Alpha MOs and Beta MOs from a specified .molden file.
! Note:
! 1) it is possible that only occupied MOs (or occupied plus some virtual
!  MOs) are printed in the .molden file.
! 2) it is possible that the number of alpha spin orbitals is larger than
!  that of beta spin orbitals.
subroutine read_nmo_from_molden(molden, nmo_a, nmo_b)
 implicit none
 integer :: i, nmo, fid
 integer, intent(out) :: nmo_a, nmo_b
!f2py intent(out) :: nmo_a, nmo_b
 character(len=6) :: key
 character(len=240) :: buf
 character(len=240), intent(in) :: molden
!f2py intent(in) :: molden

 nmo_a = 0; nmo_b = 0
 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4)=='[MO]' .or. buf(2:5)=='[MO]') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_nmo_from_molden: no [MO] found in &
                   &file '//TRIM(molden)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 buf = ADJUSTL(buf)
 i = INDEX(buf, '=')
 key = ADJUSTL(buf(1:i))
 BACKSPACE(fid)

 nmo = 0
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  if(buf(1:1)=='[' .or. buf(2:2)=='[') exit
  if(INDEX(buf(1:7),TRIM(key)) > 0) nmo = nmo + 1
  if(buf(1:5)=='Spin=' .or. buf(2:6)=='Spin=') then
   i = LEN_TRIM(buf)
   if(INDEX(buf(6:i),'Alpha') > 0) nmo_a = nmo_a + 1
  end if
 end do ! for while

 close(fid)
 nmo_b = nmo - nmo_a

 if(nmo_a < nmo_b) then
  write(6,'(/,A)') 'ERROR in subroutine read_nmo_from_molden: the number of Alp&
                   &ha spin orbitals'
  write(6,'(A)') 'is less than that of Beta spin orbitals. Not supported.'
  stop
 end if
end subroutine read_nmo_from_molden

! Read a square 2D array from an ORCA .json file. For example,
! key = 'S-Matrix' for AO overlap integral matrix
! key = 'autocipre' for AO density matrix
subroutine read_square_array_from_orca_json(json, key, n, a)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: n
!f2py intent(in) :: n
 real(kind=8), intent(out) :: a(n,n)
!f2py depend(n) :: a
!f2py intent(out) :: a
 character(len=240) :: buf
 character(len=240), intent(in) :: json
!f2py intent(in) :: json
 character(len=*), intent(in) :: key
!f2py intent(in) :: key

 a = 0d0
 open(newunit=fid,file=TRIM(json),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  k = LEN_TRIM(buf)
  if(INDEX(buf(1:k), ':') > 0) then
   if(INDEX(buf(1:k),TRIM(key)) > 0) exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_ao_dm_from_orca_json: "'//&
                    TRIM(key)//'" not found in'
  write(6,'(A)') 'file '//TRIM(json)
  close(fid)
  stop
 end if

 do i = 1, n, 1
  read(fid,'(A)') buf
  read(fid,*) a(:,i)
  read(fid,'(A)') buf
 end do ! for i

 close(fid)
end subroutine read_square_array_from_orca_json

! read MO coefficients and occupation numbers from an ORCA .json file
subroutine read_mo_and_on_from_orca_json(json, nbf, nif, mo, occ)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(out) :: mo(nbf,nif), occ(nif)
!f2py depend(nbf,nif) :: mo
!f2py depend(nif) :: occ
!f2py intent(out) :: mo, occ
 character(len=240) :: buf
 character(len=240), intent(in) :: json
!f2py intent(in) :: json

 mo = 0d0; occ = 0d0
 open(newunit=fid,file=TRIM(json),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(8:10) == 'MOs') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mo_from_orca_json: 'MOs' not found&
                   & in file "//TRIM(json)
  close(fid)
  stop
 end if

 do i = 1, nif, 1
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  read(fid,*) mo(:,i)
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  k = INDEX(buf, ':')
  read(buf(k+1:),*) occ(i)
  do k = 1, 4
   read(fid,'(A)') buf
  end do ! for k
 end do ! for i

 close(fid)
end subroutine read_mo_and_on_from_orca_json

! Write MO coefficients and occupation numbers into an ORCA .json file. This
! subroutine requires that the "MOs" section already exists in .json, and MO
! coefficients as well as occupation numbers will be updated/replaced.
subroutine write_mo_and_on_into_json(json, nbf, nif, mo, occ)
 implicit none
 integer :: i, j, fid, fid1, RENAME
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: mo(nbf,nif), occ(nif)
!f2py depend(nbf,nif) :: mo
!f2py depend(nif) :: occ
!f2py intent(in) :: mo, occ
 character(len=240) :: buf, new_json
 character(len=240), intent(in) :: json
!f2py intent(in) :: json

 call find_specified_suffix(json, '.json', i)
 new_json = json(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(json),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(new_json),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  if(buf(8:10) == 'MOs') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine write_mo_and_on_into_json: 'MOs' not fo&
                   &und in file "//TRIM(json)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 do i = 1, nif, 1
  write(fid1,'(8X,A)') '{'
  write(fid1,'(10X,A)') '"MOCoefficients": ['
  do j = 1, nbf-1, 1
   write(fid1,'(12X,ES19.12,A)') mo(j,i), ','
  end do ! for j
  write(fid1,'(12X,ES19.12)') mo(nbf,i)
  write(fid1,'(10X,A)') '],'
  write(fid1,'(10X,A,F12.8,A)') '"Occupancy": ',occ(i),','
  write(fid1,'(10X,A)') '"OrbitalEnergy": 0.0,'
  write(fid1,'(10X,A)') '"OrbitalSymLabel": "A",'
  write(fid1,'(10X,A)') '"OrbitalSymmetry": 0'
  if(i < nif) then
   write(fid1,'(8X,A)') '},'
  else
   write(fid1,'(8X,A)') '}'
  end if
 end do ! for i

 j = (8+nbf)*nif
 do i = 1, j, 1
  read(fid,'(A)') buf
 end do ! for i

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid)
 close(fid1)
 i = RENAME(TRIM(new_json), TRIM(json))
end subroutine write_mo_and_on_into_json

