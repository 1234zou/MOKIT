! written by jxzou at 20210115: move .mkl related subroutines/functions here

module mkl_content
 implicit none
 integer :: charge, mult, natom, ncontr, nbf, nif
 integer, allocatable :: nuc(:), shell_type(:), shl2atm(:)
 real(kind=8), allocatable :: alpha_coeff(:,:), beta_coeff(:,:)
 real(kind=8), allocatable :: coor(:,:), ev_a(:), ev_b(:)
 character(len=2), allocatable :: elem(:)

 type primitive_gaussian
  character(len=2) :: stype = ' ' ! S,P,SP,D,F,G,H,I
  integer :: nline = 0
  integer :: ncol  = 0
  real(kind=8), allocatable :: coeff(:,:)
 end type primitive_gaussian

 type pg4atom
  integer :: nc   ! size of array prim_gau
  type(primitive_gaussian), allocatable :: prim_gau(:)
 end type pg4atom

 type(pg4atom), allocatable :: all_pg(:) ! size is natom

contains

! check whether beta MOs exist in a given ORCA .mkl file
subroutine check_uhf_in_mkl(mklname, uhf)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname
 logical, intent(out) :: uhf

 uhf = .false.
 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit

  if(buf(1:8) == '$COEFF_B') then
   uhf = .true.
   close(fid)
   return
  end if
 end do ! for while

 close(fid)
end subroutine check_uhf_in_mkl

! read various information from a given .mkl file
subroutine read_mkl(mklname, uhf, read_mo)
 implicit none
 character(len=240), intent(in) :: mklname
 logical, intent(in) :: uhf, read_mo

 call read_natom_from_mkl(mklname, natom)
 allocate(elem(natom), nuc(natom), coor(3,natom))
 call read_elem_and_coor_from_mkl(mklname, natom, elem, nuc, coor, charge, mult)

 call read_ncontr_from_mkl(mklname, ncontr)
 allocate(shell_type(ncontr), shl2atm(ncontr))
 call read_shltyp_and_shl2atm_from_mkl(mklname, ncontr, shell_type, shl2atm)
 call read_all_pg(mklname)
 call un_normalized_all_pg()
 call merge_s_and_p_into_sp()

 if(.not. read_mo) return
 call read_nbf_and_nif_from_mkl(mklname, nbf, nif)
 allocate(alpha_coeff(nbf,nif), source=0d0)
 call read_mo_from_mkl(mklname, nbf, nif, 'a', alpha_coeff)
 if(uhf) then
  allocate(beta_coeff(nbf,nif), source=0d0)
  call read_mo_from_mkl(mklname, nbf, nif, 'b', beta_coeff)
 end if

 allocate(ev_a(nif))
 call read_ev_from_mkl(mklname, nif, 'a', ev_a)
 if(uhf) then
  allocate(ev_b(nif))
  call read_ev_from_mkl(mklname, nif, 'b', ev_b)
 end if
end subroutine read_mkl

! read the basis set data of all atoms
subroutine read_all_pg(mklname)
 implicit none
 integer :: i, j, k, nc, nline, fid
 integer, external :: detect_ncol_in_buf
 character(len=1) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname

 if(natom == 0) then
  write(6,'(A)') 'ERROR in subroutine read_all_pg: natom = 0.'
  write(6,'(A)') 'You must call subroutine read_natom_from_mkl to get the&
                   & value of natom, before calling this subroutine.'
  stop
 end if

 allocate(all_pg(natom))

 if(.not. allocated(shl2atm)) then
  write(6,'(A)') 'ERROR in subroutine read_all_pg: array shl2atm is not allocated.'
  write(6,'(A)') 'You must call subroutines read_ncontr_from_mkl and &
                   & read_shltyp_and_shl2atm_from_mkl'
  write(6,'(A)') 'to generate the array shl2atm, before calling this subroutine.'
  stop
 end if

 do i = 1, natom, 1
  nc = COUNT(shl2atm==i)
  all_pg(i)%nc = nc
  allocate(all_pg(i)%prim_gau(nc))
  all_pg(i)%prim_gau(:)%ncol = 2
  all_pg(i)%prim_gau(:)%nline = 0
 end do ! for i

 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == '$BAS') exit
 end do ! for while

 do i = 1, natom, 1
  nc = all_pg(i)%nc

  do j = 1, nc, 1
   read(fid,*) k, all_pg(i)%prim_gau(j)%stype

   nline = 0
   do while(.true.)
    read(fid,'(A)') buf
    if(LEN_TRIM(buf)==0 .or. index(buf(1:5),'$END')/=0) exit
    if(buf(1:2) == '$$') exit
    if(detect_ncol_in_buf(buf) == 3) then
     BACKSPACE(fid)
     exit
    end if
    nline = nline + 1
   end do ! for while

   all_pg(i)%prim_gau(j)%nline = nline
   allocate(all_pg(i)%prim_gau(j)%coeff(nline,2), source=0d0)
   if(buf(1:2) == '$$') exit
  end do ! for j
 end do ! for i

 rewind(fid)
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == '$BAS') exit
 end do ! for while

 do i = 1, natom, 1
  nc = all_pg(i)%nc

  do j = 1, nc, 1
   read(fid,*) k, str
   nline = all_pg(i)%prim_gau(j)%nline

   do k = 1, nline, 1
    read(fid,*) all_pg(i)%prim_gau(j)%coeff(k,1:2)
   end do ! for k
  end do ! for j
  read(fid,'(A)') buf
 end do ! for i

 close(fid)
end subroutine read_all_pg

subroutine un_normalized_all_pg()
 implicit none
 integer :: i, j, k, nc, nline, itype
 real(kind=8) :: norm_fac
 real(kind=8), allocatable :: coeff(:,:)

 do i = 1, natom, 1
  nc = all_pg(i)%nc

  do j = 1, nc, 1
   nline = all_pg(i)%prim_gau(j)%nline
   if(nline == 1) then
    all_pg(i)%prim_gau(j)%coeff(1,2) = 1d0
    cycle
   end if

   call stype2itype(all_pg(i)%prim_gau(j)%stype, itype)
   if(itype == 0) then
    write(6,'(A)') "ERROR in subroutine un_normalized_all_pg: 'L' type should&
                     & not be in the .mkl file."
    write(6,'(A)') 'This .mkl file violates the definition of mkl.'
    stop
   end if
   itype = itype - 1 ! 'S' begins with 0

   coeff = all_pg(i)%prim_gau(j)%coeff
   norm_fac = norm_fac_of_contract_gau(itype, nline, coeff)

   if(DABS(1D0-norm_fac) > 1D-4) then
    do k = 1, nline, 1
     coeff(k,2) = coeff(k,2)/norm_fac_of_prim_gau(itype, coeff(k,1))
    end do ! for k
    norm_fac = norm_fac_of_contract_gau(itype, nline, coeff)

    if(DABS(1D0-norm_fac) > 1D-3) then
     write(6,'(A)') 'ERROR in subroutine un_normalized_all_pg: the contraction&
                      & coefficients in .mkl file are wrong.'
     write(6,'(A)') 'They must be either normalized or un-normalized coefficients.'
    end if
   end if

   all_pg(i)%prim_gau(j)%coeff(:,2) = coeff(:,2)
   deallocate(coeff)
  end do ! for j
 end do ! for i

end subroutine un_normalized_all_pg

! calculate the normalization factor of a serial of primitive gaussians
function norm_fac_of_contract_gau(itype, nline, coeff) result(norm_fac)
 implicit none
 integer :: k, m, n
 integer, intent(in) :: itype, nline
 real(kind=8) :: norm_fac, rtmp
 real(kind=8), allocatable :: S(:)
 real(kind=8), intent(in) :: coeff(nline,2)

 norm_fac = 1d0
 n = (nline*nline-nline)/2
 allocate(S(n))

 n = 0
 do k = 1, nline-1, 1
  do m = k+1, nline, 1
   rtmp = coeff(k,1) + coeff(m,1)
   n = n + 1
   S(n) = (4d0*coeff(k,1)*coeff(m,1)/(rtmp*rtmp))**(0.25d0*DBLE(2*itype+3))
  end do ! for m
 end do ! for k

 norm_fac = 0d0
 do k = 1, nline, 1
  norm_fac = norm_fac + coeff(k,2)*coeff(k,2)
 end do ! for k

 n = 0
 do k = 1, nline-1, 1
  do m = k+1, nline, 1
   n = n + 1
   norm_fac = norm_fac + 2d0*coeff(k,2)*coeff(m,2)*S(n)
  end do ! for m
 end do ! for k

 deallocate(S)
end function norm_fac_of_contract_gau

! calculate the normalization factor of a primitive gaussian
function norm_fac_of_prim_gau(itype, alpha) result(norm_fac)
 implicit none
 integer, intent(in) :: itype
 real(kind=8) :: norm_fac
 real(kind=8), intent(in) :: alpha
 real(kind=8), parameter :: PI = 4d0*DATAN(1d0)

 norm_fac = 1d0

 if(itype == 0) then
  norm_fac = (2d0*alpha/PI)**(0.75d0)
 else
  norm_fac = (2d0*PI)**(0.5d0*DBLE(itype))
  norm_fac = norm_fac*(2d0*alpha/PI)**(0.25d0*DBLE(2*itype+3))
  ! the DSQRT(1*3*5...) should not be considered here
 end if
end function norm_fac_of_prim_gau

subroutine merge_s_and_p_into_sp()
 implicit none
 integer :: i, j, k, nc, nline
 real(kind=8), allocatable :: tmp_coeff(:,:), rtmp(:)

 do i = 1, natom, 1
  nc = all_pg(i)%nc
  if(nc == 1) cycle

  do j = 1, nc-1, 1
   nline = all_pg(i)%prim_gau(j)%nline
   if(nline /= all_pg(i)%prim_gau(j+1)%nline) cycle
   if(all_pg(i)%prim_gau(j)%stype /= 'S') cycle
   if(all_pg(i)%prim_gau(j+1)%stype /= 'P') cycle

   allocate(rtmp(nline))
   rtmp = all_pg(i)%prim_gau(j)%coeff(:,1) - all_pg(i)%prim_gau(j+1)%coeff(:,1)
   rtmp = DABS(rtmp)

   if( ALL( rtmp<1D-5 ) ) then
    all_pg(i)%prim_gau(j)%stype = 'SP'
    all_pg(i)%prim_gau(j)%ncol = 3
    allocate(tmp_coeff(nline,3))
    tmp_coeff(:,1:2) = all_pg(i)%prim_gau(j)%coeff
    deallocate(all_pg(i)%prim_gau(j)%coeff)
    allocate(all_pg(i)%prim_gau(j)%coeff(nline,3), source=tmp_coeff)
    deallocate(tmp_coeff)
    all_pg(i)%prim_gau(j)%coeff(:,3) = all_pg(i)%prim_gau(j+1)%coeff(:,2)
    deallocate(all_pg(i)%prim_gau(j+1)%coeff)
    all_pg(i)%prim_gau(j+1)%stype = ' '
    all_pg(i)%prim_gau(j+1)%ncol = 0
    all_pg(i)%prim_gau(j+1)%nline = 0

    if(j+1 < nc) then
     do k = j+1, nc-1, 1
      all_pg(i)%prim_gau(k)%stype = all_pg(i)%prim_gau(k+1)%stype
      all_pg(i)%prim_gau(k)%ncol  = all_pg(i)%prim_gau(k+1)%ncol
      all_pg(i)%prim_gau(k)%nline = all_pg(i)%prim_gau(k+1)%nline
      all_pg(i)%prim_gau(k)%coeff = all_pg(i)%prim_gau(k+1)%coeff
     end do ! for k
     all_pg(i)%prim_gau(nc)%stype = ' '
     all_pg(i)%prim_gau(nc)%ncol = 0
     all_pg(i)%prim_gau(nc)%nline = 0
     deallocate(all_pg(i)%prim_gau(nc)%coeff)
    end if

    nc = nc - 1
    all_pg(i)%nc = nc
   end if

   deallocate(rtmp)
  end do ! for j
 end do ! for i
end subroutine merge_s_and_p_into_sp

! find nbf and nif from a given ORCA .mkl file
subroutine read_nbf_and_nif_from_mkl(mklname, nbf, nif)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nbf, nif
 integer, external :: detect_ncol_in_buf
 real(kind=8) :: rtmp
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname

 nbf = 0; nif = 0
 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == '$COEFF_A') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine read_nbf_and_nif_from_mkl: incomplete file.'
  write(6,'(A)') 'Filename='//TRIM(mklname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 nif = detect_ncol_in_buf(buf)
 read(fid,'(A)') buf

 do while(.true.)
  read(fid,*,iostat=i) rtmp
  if(i /= 0) exit
  nbf = nbf + 1
 end do ! for while

 if(nbf == 0) then
  write(6,'(A)') 'ERROR in subroutine read_nbf_and_nif_from_mkl: zero basis function.'
  write(6,'(A)') 'There must be something wrong in file '//TRIM(mklname)//'.'
  close(fid)
  stop
 end if

 do while(.true.)
  nif = nif + detect_ncol_in_buf(buf)
  do i = 1, nbf+1, 1
   read(fid,'(A)') buf
  end do ! for i
  read(fid,'(A)') buf
  if(index(buf,'$END') /= 0) exit
 end do ! for while

 close(fid)
end subroutine read_nbf_and_nif_from_mkl

! find the array size of shell_type from a given .mkl file
subroutine read_ncontr_from_mkl(mklname, ncontr)
 implicit none
 integer :: i, fid
 integer, intent(out) :: ncontr
 integer, external :: detect_ncol_in_buf
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname

 ncontr = 0
 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == '$BAS') exit
 end do ! for while

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf)==0 .or. buf(1:4)=='$END') exit
  i = detect_ncol_in_buf(buf)
  if(i == 3) ncontr = ncontr + 1
 end do ! for while

 close(fid)
end subroutine read_ncontr_from_mkl

! read array shell_type from a given .mkl file
subroutine read_shltyp_and_shl2atm_from_mkl(mklname, ncontr, shltyp, shl2atm)
 implicit none
 integer :: i, j, k, fid, iatom
 integer, intent(in) :: ncontr
 integer, intent(out) :: shltyp(ncontr), shl2atm(ncontr)
 integer, external :: detect_ncol_in_buf
 character(len=1) :: ang
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname

 ang = ' '
 shltyp = 0; shl2atm = 0
 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == '$BAS') exit
 end do ! for while

 iatom = 1
 read(fid,'(A)') buf

 do i = 1, ncontr, 1
  shl2atm(i) = iatom
  read(buf,*,iostat=j) k, ang
  if(j /= 0) then
   write(6,'(A)') 'ERROR in subroutine read_shltyp_from_mkl: wrong format&
                    & in file '//TRIM(mklname)
   write(6,'(A)') 'buf='//TRIM(buf)
   stop
  end if

  select case(ang)
  case('S')
   shltyp(i) = 0
  case('P')
   shltyp(i) = 1
  case('D')
   shltyp(i) = -2
  case('F')
   shltyp(i) = -3
  case('G')
   shltyp(i) = -4
  case('H')
   shltyp(i) = -5
  case default
   write(6,'(A)') 'ERROR in subroutine read_shltyp_from_mkl: unsupported&
                    & angular momentum='//ang
   stop
  end select

  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf)==0 .or. index(buf,'$END')/=0) exit
   if(buf(1:2) == '$$') iatom = iatom + 1
   k = detect_ncol_in_buf(buf)
   if(k == 3) exit
  end do ! for while
 end do ! for i

 close(fid)
end subroutine read_shltyp_and_shl2atm_from_mkl

! read the number of atoms in .mkl file
subroutine read_natom_from_mkl(mklname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname

 natom = 0
 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5) == '$COOR') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_natom_from_mkl: no '$COOR' found&
                  & in file "//TRIM(mklname)
  close(fid)
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == '$END') exit
  natom = natom + 1
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_natom_from_mkl: no '$END' found&
                   & in file "//TRIM(mklname)
  close(fid)
  stop
 end if

 close(fid)
end subroutine read_natom_from_mkl

! read 3 arrays elem, nuc, coor, and the total charge as well as multiplicity
! from a given .mkl file
subroutine read_elem_and_coor_from_mkl(mklname, natom, elem, nuc, coor, charge, mult)
 use fch_content, only: nuc2elem
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 integer, intent(out) :: charge, mult, nuc(natom)
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=2), intent(out) :: elem(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname

 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')

 ! find charge and mult
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == '$CHAR_M') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_elem_and_coor_from_mkl: no&
                  & '$CHAR_M' found in file "//TRIM(mklname)
  close(fid)
  stop
 end if
 read(fid,*) charge, mult

 ! find and read coordinates
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5) == '$COOR') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_elem_and_coor_from_mkl: no&
                  & '$COOR' found in file "//TRIM(mklname)
  close(fid)
  stop
 end if

 do i = 1, natom, 1
  read(fid,*) nuc(i), coor(1:3,i)
 end do ! for i
 close(fid)

 forall(i=1:natom) elem(i) = nuc2elem(nuc(i))
end subroutine read_elem_and_coor_from_mkl

! read Alpha/Beta MO coefficients from ORCA .mkl file
subroutine read_mo_from_mkl(mklname, nbf, nif, ab, mo)
 implicit none
 integer :: i, j, k, fid, nbatch
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(out) :: mo(nbf,nif)
 character(len=240) :: buf
 character(len=11) :: key
 character(len=11), parameter :: key1 = '$COEFF_ALPH'
 character(len=11), parameter :: key2 = '$COEFF_BETA'
 character(len=1), intent(in) :: ab
 character(len=240), intent(in) :: mklname

 if(ab=='a' .or. ab=='A') then
  key = key1
 else if (ab=='b' .or. ab=='B') then
  key = key2
 else
  write(6,'(A)') 'ERROR in subroutine read_mo_from_mkl: invalid ab.'
  write(6,'(A)') 'ab = '//ab
  stop
 end if

 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == key) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_mo_from_mkl: no '"//TRIM(key)&
                    //"' found in file "//TRIM(mklname)//'.'
  stop
 end if

 nbatch = nif/5
 do i = 1, nbatch, 1
  read(fid,'(A)') buf
  read(fid,'(A)') buf

  do j = 1, nbf, 1
   read(fid,*) mo(j,5*i-4:5*i)
  end do ! for j
 end do ! for i

 k = nif - 5*nbatch
 if(k > 0) then
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  do j = 1, nbf, 1
   read(fid,*) mo(j,5*i-4:nif)
  end do ! for j
 end if

 close(fid)
end subroutine read_mo_from_mkl

! read orbital energies from an ORCA .mkl file
subroutine read_ev_from_mkl(mklname, nmo, ab, ev)
 implicit none
 integer :: i, k, rc, nline, fid
 integer, intent(in) :: nmo
 real(kind=8), intent(out) :: ev(nmo)
 character(len=1), intent(in) :: ab
 character(len=3) :: str = ' '
 character(len=8) :: key = ' '
 character(len=8), parameter :: key1 = '$COEFF_A'
 character(len=8), parameter :: key2 = '$COEFF_B'
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname

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
  write(6,'(A)') "ERROR in subroutine read_ev_from_mkl: no '"//key//"' found&
                   & in file "//TRIM(mklname)
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
  write(6,'(A)') 'ERROR in subroutine read_ev_from_mkl: incomplete .mkl file.'
  write(6,'(A)') 'Filename='//TRIM(mklname)
  stop
 end if
end subroutine read_ev_from_mkl

! read occupation numbers from ORCA .mkl file
subroutine read_on_from_mkl(mklname, nmo, ab, on)
 implicit none
 integer :: i, k, nline, fid
 integer, intent(in) :: nmo
 real(kind=8), intent(out) :: on(nmo)
 character(len=1), intent(in) :: ab
 character(len=6) :: key = ' '
 character(len=6), parameter :: key1 = '$OCC_A'
 character(len=6), parameter :: key2 = '$OCC_B'
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname

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
  write(6,'(A)') "ERROR in subroutine read_on_from_mkl: no '"//key//"' found&
                   & in file "//TRIM(mklname)
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

end module mkl_content

