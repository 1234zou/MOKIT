! written by jxzou at 20201209: move some subroutines from bas_gms2molcas.f90 to this file
! updated by jxzou at 20210222: merge subroutine read_prim_gau1 and read_prim_gau2

module pg
 implicit none
 integer :: natom     ! the number of atoms
 integer :: highest   ! highest angular momentum
 integer, allocatable :: ram(:), ntimes(:) ! ram: relative atomic mass
 ! I made a mistake, ram should be interpreted as atomic order
 real(kind=8), allocatable :: coor(:,:)    ! Cartesian coordinates
 character(len=2), allocatable :: elem(:)  ! elements

 ! 'L' will be divided into two parts: 'S' and 'P'
 type primitive_gaussian
  character(len=1) :: stype = ' ' ! 'S','P','D','F','G','H','I'
  integer :: nline = 0
  integer :: ncol  = 0
  real(kind=8), allocatable :: coeff(:,:)
 end type primitive_gaussian

 ! 7 for 'S', 'P', 'D', 'F', 'G', 'H', 'I'
 !        1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7
 type(primitive_gaussian) :: prim_gau(7)

 ! --- below are arrays for ECP/PP ---
 type ecp_potential
  integer :: n = 0   ! size for col1, 2, and 3
  real(kind=8), allocatable :: col1(:)
  integer     , allocatable :: col2(:)
  real(kind=8), allocatable :: col3(:)
 end type ecp_potential

 type ecp4atom
  logical :: ecp = .false. ! whether this atom has ECP/PP
  integer :: core_e  = 0   ! number of core electrons
  integer :: highest = 0   ! highest angular momentum
  type(ecp_potential), allocatable :: potential(:) ! size highest+1
 end type ecp4atom

 type(ecp4atom), allocatable :: all_ecp(:) ! size natom
 logical :: ecp_exist = .false.
end module pg

! check whether UHF in a given GAMESS .inp file
subroutine check_uhf_in_gms_inp(inpname, uhf)
 implicit none
 integer :: fid
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 logical, intent(out) :: uhf

 uhf = .false.
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  call upper(buf)
  if(index(buf,'SCFTYP=UHF') > 0) then
   uhf = .true.
   exit
  end if
  if(index(buf,'$END') > 0) exit
 end do ! for while

 close(fid)
end subroutine check_uhf_in_gms_inp

! find the number of atoms in GAMESS .inp file
subroutine read_natom_from_gms_inp(inpname, natom)
 implicit none
 integer :: i, fid, nline
 integer, intent(out) :: natom
 character(len=1) :: str
 character(len=240):: buf
 character(len=240), intent(in) :: inpname

 natom = 0
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf(2:6))
  if(buf(2:6) == '$DATA') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine read_natom_from_gms_inp: wrong format&
                   & in file '//TRIM(inpname)
  stop
 end if
 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do while(.true.)
  read(fid,'(A)') buf
  call upper(buf(3:5))
  if(buf(2:5) == '$END') exit

  do while(.true.)
   read(fid,'(A)') buf
   read(buf,*,iostat=i) str, nline
   if(i /= 0) exit

   do i = 1, nline, 1
    read(fid,'(A)') buf
   end do ! for i
  end do ! for while

  natom = natom + 1
 end do ! for while

 close(fid)
 if(natom == 0) then
  write(6,'(A)') 'ERROR in subroutine read_natom_from_gms_inp: zero atom&
                 & found in file '//TRIM(inpname)
  stop
 end if
end subroutine read_natom_from_gms_inp

! read charge, spin multiplicity, uhf(bool) bohrs(bool) from a given GAMESS
! .inp/.dat file
subroutine read_charge_and_mult_from_gms_inp(inpname, charge, mult, uhf, ghf, ecp)
 implicit none
 integer :: i, fid
 integer, intent(out) :: charge, mult
 character(len=240) :: buf
 character(len=480) :: buf1
 character(len=240), intent(in) :: inpname
 logical, intent(out) :: uhf, ghf, ecp

 charge = 0
 mult = 1
 uhf = .false.; ghf = .false.; ecp = .false.

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 read(fid,'(A)') buf
 read(fid,'(A)') buf1
 buf1 = TRIM(buf)//' '//TRIM(buf1)
 call upper(buf1)
 if(index(buf1,'UHF') > 0) uhf = .true.
 if(index(buf1,'GHF') > 0) ghf = .true.

 i = index(buf1,'ICHAR')
 if(i > 0) read(buf1(i+7:),*) charge
 i = index(buf1,'MULT')
 if(i > 0) read(buf1(i+5:),*) mult

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit

  if(buf(2:2) == '$') then
   call upper(buf(3:5))
   if(buf(2:5) == '$ECP') then
    ecp = .true.
    exit
   end if
  end if
 end do ! for while

 close(fid)
end subroutine read_charge_and_mult_from_gms_inp

! read elements, nuclear charges and Cartesian coordinates from a GAMESS .inp file
subroutine read_elem_nuc_coor_from_gms_inp(inpname, natom, elem, nuc, coor, ghost)
 use phys_cons, only: Bohr_const
 implicit none
 integer :: i, k, fid, nline
 integer, intent(in) :: natom
 integer, intent(out) :: nuc(natom)
 real(kind=8), allocatable :: nuc1(:)
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=1) :: str
 character(len=2), intent(out) :: elem(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 logical :: bohrs
 logical, intent(out) :: ghost(natom) ! size natom

 allocate(nuc1(natom), source=0d0)
 nuc = 0; coor = 0d0; elem = ' '; ghost = .false.

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 ! find in the first 6 lines whether the coordinates are in Angstrom or Bohr
 bohrs = .false.
 do i = 1, 6
  read(fid,'(A)') buf
  call upper(buf)
  if(index(buf,'UNITS=BOHR') /= 0) then
   bohrs = .true.
   exit
  end if
  if(index(buf,'$END') /= 0) exit
 end do ! for i
 ! Angstrom/Bohr determined

 rewind(fid)
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:6) == '$DATA') exit
 end do ! for while

 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do k = 1, natom, 1
  read(fid,*) elem(k), nuc1(k), coor(1:3,k)
  if(nuc1(k) < 0d0) then ! ghost atom (0 charge, has basis function)
   nuc1(k) = -nuc1(k)
   ghost(k) = .true.
  end if

  do while(.true.)
   read(fid,'(A)') buf
   read(buf,*,iostat=i) str, nline
   if(i /= 0) exit

   do i = 1, nline, 1
    read(fid,'(A)') buf
   end do ! for do
  end do ! for while

 end do ! for k

 close(fid)

 call standardize_elem(natom, elem)
 forall(i = 1:natom) nuc(i) = DNINT(nuc1(i))
 deallocate(nuc1)
 if(bohrs) coor = coor*Bohr_const
end subroutine read_elem_nuc_coor_from_gms_inp

! read nbf and nif from a given GAMESS .inp/.dat file
subroutine read_nbf_and_nif_from_gms_inp(inpname, nbf, nif)
 implicit none
 integer:: i, j, k, fid
 integer, intent(out) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 nbf = 0; nif = 0

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf)
  j = index(buf,'$GUESS'); k = index(buf,'NORB=')
  if(j/=0 .and. k/=0) then
   read(buf(k+5:),*) nif
   exit
  end if
 end do ! for while

 if(i /= 0) then
  close(fid)
  return
 end if

 rewind(fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf(3:6))
  if(buf(2:6)=='$DATA') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine bas_gms2py: No $DATA section found&
                 & in file '//TRIM(inpname)//'.'
  close(fid)
  stop
 end if

 ! read the Title Card line, find nbf in it
 read(fid,'(A)') buf
 i = index(buf,'nbf=')
 if(i == 0) then ! if not found, set to nif
  nbf = nif
 else
  read(buf(i+4:),*) nbf
 end if

 close(fid)
end subroutine read_nbf_and_nif_from_gms_inp

! read Cartesian-type nbf and nif from GAMESS .inp/.dat file
subroutine read_cart_nbf_nif_from_dat(datname, nbf, nif)
 implicit none
 integer :: i, j, k, fid
 integer, intent(out) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: datname

 nbf = 0; nif = 0
 open(newunit=fid,file=TRIM(datname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:2) == '$') call upper(buf(3:5))
  if(buf(2:5) == '$VEC') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_cart_nbf_nif_from_dat: no '$VEC'&
                & found in file "//TRIM(datname)
  close(fid)
  stop
 end if

 j = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:2) /= ' 1') exit
  j = j + 1
 end do ! for while

 k = j ! backup

 BACKSPACE(fid)
 BACKSPACE(fid)
 read(fid,'(A)') buf
 nbf = (j-1)*5 + (LEN_TRIM(buf)-5)/15

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:2) == '$') call upper(buf(3:5))
  if(buf(2:5) == '$END') exit
  j = j + 1
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_cart_nbf_nif_from_dat: no '$END'&
                & corresponds to '$VEC'"
  write(6,'(A)') 'in file '//TRIM(datname)
  stop
 end if

 nif = j/k
end subroutine read_cart_nbf_nif_from_dat

! read na, nb, nif and nbf from a given GAMESS .inp file
! Note: when spherical harmonic functions are used, the nbf here will <=
!  the number of basis functions in $VEC (where MOs are always expanded
!  on Cartesian functions)
subroutine read_na_nb_nif_nbf_from_gms_inp(inpname, na, nb, nif, nbf)
 implicit none
 integer :: i, fid
 integer, intent(out) :: na, nb, nif, nbf
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 na = 0; nb = 0; nif = 0; nbf = 0
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf(3:6))
  if(buf(2:6)=='$DATA') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine read_na_nb_nif_nbf_from_gms_inp: No&
                 & $DATA section found in'
  write(6,'(A)') 'file '//TRIM(inpname)//'.'
  close(fid)
  stop
 end if

 ! read the Title Card line
 read(fid,'(A)') buf
 close(fid)

 i = index(buf,'nbf=')
 read(buf(i+4:),*) nbf
 buf(i-1:) = ' '

 i = index(buf,'nif=')
 read(buf(i+4:),*) nif
 buf(i-1:) = ' '

 i = index(buf,'nb=')
 read(buf(i+3:),*) nb
 buf(i-1:) = ' '

 i = index(buf,'na=')
 read(buf(i+3:),*) na
end subroutine read_na_nb_nif_nbf_from_gms_inp

! read type all_ecp from a given GAMESS .inp/.dat file
subroutine read_all_ecp_from_gms_inp(inpname)
 use pg, only: natom, all_ecp, ecp_exist
 implicit none
 integer :: i, j, k, m, n, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 if(natom == 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_all_ecp_from_gms_inp: natom = 0.'
  write(6,'(A)') 'The variable natom should be initialized before calling this &
                 &subroutine.'
  stop
 end if

 allocate(all_ecp(natom))
 all_ecp(:)%ecp = .false.
 if(.not. ecp_exist) return

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:2) == '$') then
   call upper(buf(3:5))
   if(buf(3:5) == 'ECP') exit
  end if
 end do ! for while

 do i = 1, natom, 1
  read(fid,'(A)') buf
  call upper(buf)
  if(index(buf,'NONE') /= 0) cycle

  all_ecp(i)%ecp = .true.
  k = index(buf,'GEN')
  read(buf(k+3:),*) all_ecp(i)%core_e, m
  all_ecp(i)%highest = m
  allocate(all_ecp(i)%potential(m+1))

  do j = 1, m+1, 1
   read(fid,'(A)') buf
   if(j == 1) then
    k = index(buf,'ul')
    if(k == 0) then
     write(6,'(A)') "ERROR in subroutine read_all_ecp_from_gms_inp: ECP/PP does&
                    &not start with '-ul potential'."
     write(6,'(A)') 'You should check the format of ECP/PP data in file '//TRIM(inpname)
     stop
    end if
   end if

   read(buf,*) n
   all_ecp(i)%potential(j)%n = n
   allocate(all_ecp(i)%potential(j)%col1(n), source=0d0)
   allocate(all_ecp(i)%potential(j)%col2(n), source=0)
   allocate(all_ecp(i)%potential(j)%col3(n), source=0d0)
   do k = 1, n, 1
    read(fid,*) all_ecp(i)%potential(j)%col1(k), all_ecp(i)%potential(j)%col2(k), &
    all_ecp(i)%potential(j)%col3(k)
   end do ! for k
  end do ! for j
 end do ! for i

 close(fid)
end subroutine read_all_ecp_from_gms_inp

! deallocate the allocatable arrays in array prim_gau
subroutine clear_prim_gau()
 use pg, only: prim_gau
 implicit none
 integer :: i

 do i = 1, 7, 1
  prim_gau(i)%stype = ' '
  prim_gau(i)%nline = 0
  prim_gau(i)%ncol = 0
  if(allocated(prim_gau(i)%coeff)) deallocate(prim_gau(i)%coeff)
 end do ! for i
end subroutine clear_prim_gau

! read this type of primitive gaussians, i.e., 'S', 'L', etc.
subroutine read_prim_gau(stype, nline, fid)
 use pg, only: prim_gau
 implicit none
 integer :: i, j, k, itmp, ncol, ncol1, nline0, nline1
 integer, intent(in) :: nline, fid
 real(kind=8), parameter :: zero = 1d-6
 real(kind=8) :: exp0, rtmp, rtmp1
 real(kind=8), allocatable :: coeff(:,:)
 character(len=1), intent(in) :: stype
 logical :: share

 call stype2itype(stype, k)
 ! two cases: 'L', or not 'L'

 if(k /= 0) then ! 'S','P','D','F','G','H','I', not 'L'

  ! two subcases: whether or not this angular momentum has occurred
  if(allocated(prim_gau(k)%coeff)) then
   read(fid,*) itmp, exp0, rtmp
   BACKSPACE(fid)

   nline0 = prim_gau(k)%nline
   ncol = prim_gau(k)%ncol
   allocate(coeff(nline0,ncol), source=prim_gau(k)%coeff)
   deallocate(prim_gau(k)%coeff)

   ! further two sub-subcases: whether sharing the same exponents
   share = .false. ! initialization
   if(DABS(coeff(1,1)-exp0) < zero) then
    nline1 = max(nline0,nline)
    share = .true.
   end if
   if((.not.share) .and. nline==1 .and. DABS(coeff(nline0,1)-exp0)<zero) then
    ! this rare case occurs in nobasistransform cc-pVDZ for Zn
    nline1 = nline0
    share = .true.
   end if
   if(.not. share) nline1 = nline0 + nline

   allocate(prim_gau(k)%coeff(nline1,ncol+1),source=0d0)
   prim_gau(k)%coeff(1:nline0,1:ncol) = coeff
   deallocate(coeff)
   ncol = ncol + 1
   prim_gau(k)%ncol = ncol
   prim_gau(k)%nline = nline1

   do i = 1, nline, 1
    read(fid,*) itmp, exp0, rtmp
    if(nline == 1) rtmp = 1d0
    if(share) then
     if(i > nline0) then
      prim_gau(k)%coeff(i,1) = exp0
      prim_gau(k)%coeff(i,ncol) = rtmp
     else ! i <= nline0
      do j = 1, nline0, 1
       if(DABS(prim_gau(k)%coeff(j,1)-exp0) < zero) then
        prim_gau(k)%coeff(j,ncol) = rtmp
        exit
       end if
      end do ! for j
      if(j == nline0+1) then
       write(6,'(/,A)') "ERROR in subroutine read_prim_gau: I've never seen such basis."
       write(6,'(A)') "Did you forget to add 'int(nobasistransform)' in .gjf file?"
       write(6,'(3(A,I0))') 'stype='//stype//', i=',i,', nline=',nline,', nline0=',nline0
       write(6,'(2(A,E16.8))') 'exp0=', exp0, ', rtmp=', rtmp
       stop
      end if
     end if
    else ! not share
     prim_gau(k)%coeff(nline0+i,1) = exp0
     prim_gau(k)%coeff(nline0+i,ncol) = rtmp
    end if
   end do ! for i

  else ! never occurs before, read first time
   prim_gau(k)%stype = stype
   prim_gau(k)%nline = nline
   prim_gau(k)%ncol = 2
   allocate(prim_gau(k)%coeff(nline,2), source=0d0)
   do i = 1, nline, 1
    read(fid,*) itmp, prim_gau(k)%coeff(i,1), rtmp
    if(nline == 1) rtmp = 1d0
    prim_gau(k)%coeff(i,2) = rtmp
   end do ! for i
  end if

 else ! k == 0, this is 'L' for Pople-type basis sets
  ! sharing exponents for 'L' is unlikely in Pople basis sets, don't consider

  if(allocated(prim_gau(1)%coeff)) then ! 'S'
   nline0 = prim_gau(1)%nline
   ncol = prim_gau(1)%ncol
   prim_gau(1)%nline = nline0 + nline
   prim_gau(1)%ncol = ncol + 1
   allocate(coeff(nline0,ncol), source=prim_gau(1)%coeff)
   deallocate(prim_gau(1)%coeff)
   allocate(prim_gau(1)%coeff(nline0+nline,ncol+1),source=0d0)
   prim_gau(1)%coeff(1:nline0,1:ncol) = coeff
   deallocate(coeff)
   ncol = ncol + 1
  else
   prim_gau(1)%stype = 'S'
   prim_gau(1)%nline = nline
   prim_gau(1)%ncol = 2
   allocate(prim_gau(1)%coeff(nline,2),source=0d0)
   nline0 = 0; ncol = 2
  end if

  if(allocated(prim_gau(2)%coeff)) then ! 'P'
   nline1 = prim_gau(2)%nline
   ncol1 = prim_gau(2)%ncol
   prim_gau(2)%nline = nline1 + nline
   prim_gau(2)%ncol = ncol1 + 1
   allocate(coeff(nline1,ncol1), source=prim_gau(2)%coeff)
   deallocate(prim_gau(2)%coeff)
   allocate(prim_gau(2)%coeff(nline1+nline,ncol1+1),source=0d0)
   prim_gau(2)%coeff(1:nline1,1:ncol1) = coeff
   deallocate(coeff)
   ncol1 = ncol1 + 1
  else
   prim_gau(2)%stype = 'P'
   prim_gau(2)%nline = nline
   prim_gau(2)%ncol = 2
   allocate(prim_gau(2)%coeff(nline,2),source=0d0)
   nline1 = 0; ncol1 = 2
  end if

  do i = 1, nline, 1
   read(fid,*) itmp, exp0, rtmp, rtmp1
   if(nline == 1) then
    rtmp = 1d0; rtmp1 = 1d0
   end if
   prim_gau(1)%coeff(nline0+i,1) = exp0
   prim_gau(1)%coeff(nline0+i,ncol) = rtmp
   prim_gau(2)%coeff(nline1+i,1) = exp0
   prim_gau(2)%coeff(nline1+i,ncol1) = rtmp1
  end do ! for i
 end if

end subroutine read_prim_gau

! determine the highest angular momentum quantum number
subroutine get_highest_am()
 use pg, only: highest, prim_gau
 implicit none

 if(allocated(prim_gau(7)%coeff)) then
  highest = 6
 else if(allocated(prim_gau(6)%coeff)) then
  highest = 5
 else if(allocated(prim_gau(5)%coeff)) then
  highest = 4
 else if(allocated(prim_gau(4)%coeff)) then
  highest = 3
 else if(allocated(prim_gau(3)%coeff)) then
  highest = 2
 else if(allocated(prim_gau(2)%coeff)) then
  highest = 1
 else
  highest = 0
 end if
end subroutine get_highest_am

! print primitive gaussians
subroutine prt_prim_gau(iatom, fid)
 use pg, only: prim_gau, ram, highest, all_ecp, elem, ecp_exist
 implicit none
 integer :: i, j, k, m, n, nline, ncol
 integer, intent(in) :: iatom, fid
 integer, allocatable :: list(:)
 character(len=1), parameter :: am(0:6) = ['S','P','D','F','G','H','I']

 call get_highest_am()
 if(ecp_exist) then
  write(fid,'(5X,I0,A1,3X,I1)') ram(iatom)-all_ecp(iatom)%core_e,'.',highest
 else
  write(fid,'(5X,I0,A1,3X,I1)') ram(iatom), '.', highest
 end if

 do i = 1, 7, 1
  if(.not. allocated(prim_gau(i)%coeff)) cycle
  write(fid,'(A)') '* '//prim_gau(i)%stype//'-type functions'
  nline = prim_gau(i)%nline
  ncol = prim_gau(i)%ncol
  write(fid,'(2(1X,I4))') nline, ncol-1
  do j = 1, nline, 1
   write(fid,'(3X,ES16.9)') prim_gau(i)%coeff(j,1)
  end do ! for j
  do j = 1, nline, 1
   write(fid,'(10(ES16.9,2X))') (prim_gau(i)%coeff(j,k), k=2,ncol)
  end do ! for j
 end do ! for i

 if(.not. ecp_exist) return

 if(all_ecp(iatom)%ecp) then
  m = all_ecp(iatom)%highest
  write(fid,'(A2,1X,A2,1X,I0,1X,I0)') 'PP', elem(iatom), all_ecp(iatom)%core_e, m
  allocate(list(m+1))
  list(1) = m
  forall(i=2:m+1) list(i) = i-2

  do i = 1, m+1, 1
   n = all_ecp(iatom)%potential(i)%n
   if(i == 1) then
    write(fid,'(I0,A)') n,';!'//am(list(1))//' POTENTIAL'
   else ! i > 1
    write(fid,'(I0,A)') n,';!'//am(list(i))//'-'//am(list(1))//' POTENTIAL'
   end if

   do j = 1, n, 1
    write(fid,'(I0,2(1X,F16.8))') all_ecp(iatom)%potential(i)%col2(j),&
     all_ecp(iatom)%potential(i)%col3(j), all_ecp(iatom)%potential(i)%col1(j)
   end do ! for j
  end do ! for i
  write(fid,'(A1,/,A,/,A)') '*', 'Spectral', 'End of Spectral'
 end if
end subroutine prt_prim_gau

! update the number of times each atom occurred
subroutine update_ntimes(iatom)
 use pg, only: elem, ntimes
 implicit none
 integer :: i
 integer, intent(in) :: iatom
 character(len=2) :: ctmp

 ctmp = elem(iatom)

 do i = iatom-1, 1, -1
  if(ctmp == elem(i)) then
   ntimes(iatom) = ntimes(i) + 1
   exit
  end if
 end do ! for i
end subroutine update_ntimes

! update the number of times each atom occurred
subroutine calc_ntimes(natom, elem, ntimes)
 implicit none
 integer :: i, j
 integer, intent(in) :: natom
 integer, intent(out) :: ntimes(natom)
 character(len=2) :: tmp
 character(len=2), intent(in) :: elem(natom)

 ntimes = 1

 do i = 2, natom, 1
  tmp = elem(i)

  do j = i-1, 1, -1
   if(tmp == elem(j)) then
    ntimes(i) = ntimes(j) + 1
    exit
   end if
  end do ! for j

 end do ! for i
end subroutine calc_ntimes

! generate contracted string, e.g. 5s3p1d -> 3s2p1d
subroutine gen_contracted_string(nline, ncol, str1, str2)
 implicit none
 integer :: i
 integer, intent(in) :: nline(7), ncol(7)
 character(len=1), parameter :: am(7) = ['s','p','d','f','g','h','i']
 character(len=3) :: str
 character(len=21), intent(out) :: str1, str2
 !10s10p10d10f10g10h10i

 str1 = ' '; str2 = ' '

 do i = 1, 7, 1
  if(nline(i) > 0) then
   write(str,'(I0,A1)') nline(i), am(i)
   str1 = TRIM(str1)//TRIM(str)
   write(str,'(I0,A1)') ncol(i)-1,am(i)
   str2 = TRIM(str2)//TRIM(str)
  end if
 end do ! for i
end subroutine gen_contracted_string

! read Alpha or (both Alpha and Beta) MOs from a GAMESS .dat or .inp file
! Note: if you want to read both Alpha and Beta MOs, just double the variable
! nif
subroutine read_mo_from_dat(datname, nbf, nif, coeff)
 implicit none
 integer i, j, k, nline, nleft, fid
 integer, intent(in) :: nbf, nif
 character(len=5) :: str1
 character(len=30) :: str2
 character(len=240) :: buf
 character(len=240), intent(in) :: datname
 real(kind=8), intent(out) :: coeff(nbf,nif)

 coeff = 0d0
 open(newunit=fid,file=TRIM(datname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf(3:5))
  if(buf(2:5) == '$VEC') exit
 end do

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_mo_from_dat: No '$VEC' section in&
                 & file "//TRIM(datname)
  close(fid)
  stop
 end if

 nline = nbf/5
 nleft = nbf - nline*5

 do i = 1, nif, 1
  k = 1
  do j = 1, nline, 1
   read(fid,'(A)') buf
   buf = buf(6:)
   read(buf,'(5ES15.8)') coeff(k:k+4,i)
   k = k + 5
  end do ! for j

  if(nleft > 0) then
   read(fid,'(A)') buf
   buf = buf(6:)
   str1 = ' '
   write(str1,'(I5)') nleft
   str1 = ADJUSTL(str1)
   str2 = '('//TRIM(str1)//'ES15.8)'
   read(buf,TRIM(str2)) coeff(k:nbf,i)
  end if
 end do ! for i

 close(fid)
end subroutine read_mo_from_dat

! read the number of GVB pairs from a GAMESS .dat file
subroutine read_npair_from_dat(datname, npair)
 implicit none
 integer :: i, fid
 integer, intent(out) :: npair
 character(len=240) :: buf
 character(len=240), intent(in) :: datname

 npair = 0
 open(newunit=fid,file=TRIM(datname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:5) == '$SCF') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_npair_from_dat: no '$SCF' found&
                & in file "//TRIM(datname)
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 do while(.true.)
  read(fid,'(A)') buf
  i = index(buf, 'CICOEF')
  if(i > 0) npair = npair + 1
  if(index(buf,'END')>0 .or. index(buf,'VEC')>0) exit
 end do ! for while

 close(fid)
end subroutine read_npair_from_dat

! read CI coefficients from a GAMESS .dat or .inp file
subroutine read_ci_coeff_from_dat(fname, npair, coeff)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: npair
 character(len=240) :: buf
 character(len=240), intent(in) :: fname
 real(kind=8), intent(out) :: coeff(2,npair)

 buf = ' '; coeff = 0d0
 open(newunit=fid,file=TRIM(fname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  j = index(buf,'CICOEF(')
  if(j == 0) j = index(buf,'cicoef(')
  if(j /= 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine read_ci_coeff_from_dat: no GVB CI&
                & coefficients found in file '//TRIM(fname)
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 do i = 1, npair, 1
  read(fid,'(A)') buf
  j = index(buf,'=')
  k = index(buf,',')
  read(buf(j+1:k-1),*) coeff(1,i)
  read(buf(k+1:),*) coeff(2,i)
 end do ! for i

 close(fid)
end subroutine read_ci_coeff_from_dat

! print MOs into .dat file
! if replace is .true., the MOs in the original file will be replaced
! if replace is .false., a new *_new.dat file will be generated
subroutine write_mo_into_dat(datname, nbf, nif, coeff, replace)
 implicit none
 integer :: i, j, k, nline, nleft, fid1, fid2, RENAME
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: coeff(nbf,nif)
 character(len=240), intent(in) :: datname
 character(len=240) :: newdat, buf
 logical, intent(in) :: replace

 buf = ' '; newdat = ' '
 i = index(datname,'.dat',.true.)
 if(i == 0) i = index(datname,'.inp',.true.)
 newdat = datname(1:i-1)//'_new.dat'

 open(newunit=fid1,file=TRIM(datname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(newdat),status='replace')
 do while(.true.)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  call upper(buf(3:5))
  if(buf(2:5)=='$VEC') exit
 end do ! for while

 ! print MOs
 nline = nbf/5
 nleft = nbf - 5*nline

 do i = 1, nif, 1
  k = MOD(i,100)
  do j = 1, nline, 1
   write(fid2,'(I2,I3,5ES15.8)') k, MOD(j,1000), coeff(5*j-4:5*j,i)
  end do ! for j
  if(nleft > 0) then
   write(fid2,'(I2,I3,5ES15.8)') k, MOD(j,1000), coeff(5*j-4:nbf,i)
  end if
 end do ! for i
 write(fid2,'(A)') ' $END'
 ! print MOs done

 ! skip the MOs in datname
 do while(.true.)
  read(fid1,'(A)') buf
  call upper(buf(3:5))
  if(buf(2:5) == '$END') exit
 end do

 ! copy remaining contens
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while
 close(fid2)
 ! copy done

 if(replace) then
  close(fid1,status='delete')
  i = RENAME(TRIM(newdat), TRIM(datname))
 else
  close(fid1)
 end if
end subroutine write_mo_into_dat

! delete the first $VEC section in a specified .dat file
subroutine del_vec_in_dat(datname)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, datname1
 character(len=240), intent(in) :: datname

 datname1 = TRIM(datname)//'.t'
 open(newunit=fid,file=TRIM(datname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(datname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:5) == '$VEC') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:5) == '$END') exit
 end do ! for while

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(datname1), TRIM(datname))
end subroutine del_vec_in_dat

