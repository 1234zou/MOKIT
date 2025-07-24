! written by jxzou at 20210115: move .mkl related subroutines/functions here

module mkl_content
 implicit none
 integer :: charge, mult, natom, ncontr, nbf, nif
 integer, allocatable :: nuc(:), shell_type(:), shl2atm(:)
 real(kind=8), allocatable :: alpha_coeff(:,:), beta_coeff(:,:)
 real(kind=8), allocatable :: coor(:,:), ev_a(:), ev_b(:)
 character(len=2), allocatable :: elem(:)

 type :: primitive_gaussian
  character(len=2) :: stype = '  ' ! S,P,SP,D,F,G,H,I
  integer :: nline = 0
  integer :: ncol  = 0
  real(kind=8), allocatable :: coeff(:,:)
 end type primitive_gaussian

 type :: pg4atom
  integer :: nc   ! size of array prim_gau
  type(primitive_gaussian), allocatable :: prim_gau(:)
 end type pg4atom

 type(pg4atom), allocatable :: all_pg(:) ! size is natom

contains

subroutine copy_type_pg4atom(p_in, p_out)
 implicit none
 integer :: i, nc, nline, ncol
 type(pg4atom), intent(in) :: p_in
 type(pg4atom), intent(out) :: p_out

 nc = p_in%nc
 p_out%nc = nc

 if(nc > 0) then
  allocate(p_out%prim_gau(nc))
  do i = 1, nc, 1
   p_out%prim_gau(i)%stype = p_in%prim_gau(i)%stype
   nline = p_in%prim_gau(i)%nline
   ncol = p_in%prim_gau(i)%ncol
   p_out%prim_gau(i)%nline = nline
   p_out%prim_gau(i)%ncol = ncol
   if(nline>0 .and. ncol>0) then
    allocate(p_out%prim_gau(i)%coeff(nline,ncol))
    p_out%prim_gau(i)%coeff = p_in%prim_gau(i)%coeff
   end if
  end do ! for i
 end if
end subroutine copy_type_pg4atom

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
 integer :: i, k, nlin
 real(kind=8), parameter :: thres = 1d-6
 real(kind=8), allocatable :: r1(:), r2(:,:)
 character(len=240), intent(in) :: mklname
 logical, intent(in) :: uhf, read_mo

 call read_natom_from_mkl(mklname, natom)
 allocate(elem(natom), nuc(natom), coor(3,natom))
 call read_elem_and_coor_from_mkl(mklname, natom, elem, nuc, coor, charge, mult)

 call read_ncontr_from_mkl(mklname, ncontr)
 allocate(shell_type(ncontr), shl2atm(ncontr))
 call read_shltyp_and_shl2atm_from_mkl(mklname, ncontr, shell_type, shl2atm)
 call read_all_pg_from_mkl(mklname)
 call un_normalized_all_pg()
 call merge_s_and_p_into_sp()

 if(.not. read_mo) return
 call read_nbf_and_nif_from_mkl(mklname, nbf, nif)
 allocate(alpha_coeff(nbf,nif))
 call read_mo_from_mkl(mklname, nbf, nif, 'a', alpha_coeff)
 if(uhf) then
  allocate(beta_coeff(nbf,nif))
  call read_mo_from_mkl(mklname, nbf, nif, 'b', beta_coeff)
 end if

 allocate(ev_a(nif))
 call read_ev_from_mkl(mklname, nif, 'a', ev_a)
 if(uhf) then
  allocate(ev_b(nif))
  call read_ev_from_mkl(mklname, nif, 'b', ev_b)
 end if

 ! Check basis set linear dependency. ORCA set linear dependent MOs as zero in
 ! .mkl and .molden files.
 nlin = 0
 do i = nif, 1, -1
  if(DABS(ev_a(i)) < thres) then
   if(SUM(DABS(alpha_coeff(:,i))) < thres) nlin = nlin + 1
  else
   exit
  end if
 end do ! for i

 ! If there is linear dependency, remove MOs with zero coefficients and update
 ! nif.
 if(nlin > 0) then
  k = nif - nlin
  allocate(r1(k), source=ev_a(1:k))
  deallocate(ev_a)
  allocate(ev_a(k), source=r1)
  allocate(r2(nbf,k), source=alpha_coeff(:,1:k))
  deallocate(alpha_coeff)
  allocate(alpha_coeff(nbf,k), source=r2)
  if(uhf) then
   r1 = ev_b(1:k)
   deallocate(ev_b)
   allocate(ev_b(k), source=r1)
   r2 = beta_coeff(:,1:k)
   deallocate(beta_coeff)
   allocate(beta_coeff(nbf,k), source=r2)
  end if
  deallocate(r1, r2)
  nif = k
 end if
end subroutine read_mkl

! read the basis set data of all atoms from .mkl file
subroutine read_all_pg_from_mkl(mklname)
 implicit none
 integer :: i, j, k, nc, nline, fid
 integer, external :: detect_ncol_in_buf
 character(len=1) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname

 if(natom == 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_all_pg_from_mkl: natom = 0.'
  write(6,'(A)') 'You must call subroutine read_natom_from_mkl to get the value&
                 & of natom before'
  write(6,'(A)') 'calling this subroutine.'
  stop
 end if

 if(.not. allocated(shl2atm)) then
  write(6,'(/,A)') 'ERROR in subroutine read_all_pg_from_mkl: array shl2atm is &
                   &not allocated.'
  write(6,'(A)') 'You must call subroutines read_ncontr_from_mkl and read_shlty&
                 &p_and_shl2atm_from_mkl'
  write(6,'(A)') 'to generate the array shl2atm, before calling this subroutine.'
  stop
 end if
 allocate(all_pg(natom))

 do i = 1, natom, 1
  nc = COUNT(shl2atm == i)
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
 read(fid,'(A)') buf

 do i = 1, natom, 1
  nc = all_pg(i)%nc

  do j = 1, nc, 1
   read(buf,*) k, str
   all_pg(i)%prim_gau(j)%stype(1:1) = str

   nline = 0
   do while(.true.)
    read(fid,'(A)') buf
    if(LEN_TRIM(buf)==0 .or. INDEX(buf(1:5),'$END')>0) exit
    if(buf(1:2) == '$$') exit
    if(detect_ncol_in_buf(buf) == 3) exit
    nline = nline + 1
   end do ! for while

   all_pg(i)%prim_gau(j)%nline = nline
   allocate(all_pg(i)%prim_gau(j)%coeff(nline,2), source=0d0)
   if(buf(1:2) == '$$') then
    read(fid,'(A)') buf
    exit
   end if
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
    read(fid,*) all_pg(i)%prim_gau(j)%coeff(k,:)
   end do ! for k
  end do ! for j
  read(fid,'(A)') buf
 end do ! for i

 close(fid)
end subroutine read_all_pg_from_mkl

subroutine un_normalized_all_pg()
 implicit none
 integer :: i, j, k, nc, nline, itype
 real(kind=8) :: norm_fac0, norm_fac
 real(kind=8), allocatable :: coeff(:,:)

 do i = 1, natom, 1
  nc = all_pg(i)%nc

  do j = 1, nc, 1
   nline = all_pg(i)%prim_gau(j)%nline
   if(nline == 1) then
    all_pg(i)%prim_gau(j)%coeff(1,2) = 1d0
    cycle
   end if

   ! subroutine stype2itype accepts only one character length
   call stype2itype(all_pg(i)%prim_gau(j)%stype(1:1), itype)
   if(itype == 0) then
    write(6,'(/,A)') "ERROR in subroutine un_normalized_all_pg: 'L' type should&
                     & not be in the"
    write(6,'(A)') '.mkl file. It violates the definition of mkl.'
    stop
   end if
   itype = itype - 1 ! 'S' begins with 0, minus 1 is needed

   allocate(coeff(nline,2), source=all_pg(i)%prim_gau(j)%coeff)
   norm_fac = norm_fac_of_contract_gau(itype, nline, coeff)
   norm_fac0 = norm_fac

   if(DABS(1D0-norm_fac) > 1D-4) then
    do k = 1, nline, 1
     coeff(k,2) = coeff(k,2)/norm_fac_of_prim_gau(itype, coeff(k,1))
    end do ! for k
    norm_fac = norm_fac_of_contract_gau(itype, nline, coeff)
   end if

   if(DABS(1D0-norm_fac) > 1D-4) then
    coeff(:,2) = all_pg(i)%prim_gau(j)%coeff(:,2)/DSQRT(norm_fac0)
    norm_fac = norm_fac_of_contract_gau(itype, nline, coeff)
   end if

   if(DABS(1D0-norm_fac) > 1D-3) then
    write(6,'(/,A)') 'ERROR in subroutine un_normalized_all_pg: the contraction&
                     & coefficients'
    write(6,'(A)') 'are wrong. They must be either normalized or un-normalized &
                   &coefficients.'
    write(6,'(A,I0,A,F12.8)') 'itype=', itype, ', norm_fac=', norm_fac
    stop
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
 n = (nline*nline - nline)/2
 allocate(S(n), source=0d0)

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

 select case(itype)
 case(0)
  norm_fac = (2d0*alpha/PI)**(0.75d0)
 case(1,2,3,4,5,6)
  norm_fac = (2d0*PI)**(0.5d0*DBLE(itype))
  norm_fac = norm_fac*(2d0*alpha/PI)**(0.25d0*DBLE(2*itype+3))
  if(itype == 6) norm_fac = norm_fac/DSQRT(10395d0)
  ! 10395 = 3*5*7*9*11
 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine norm_fac_of_prim_gau: invalid itype&
                      &=', itype
  stop
 end select
end function norm_fac_of_prim_gau

subroutine merge_s_and_p_into_sp()
 implicit none
 integer :: i, j, k, m, nc, nline
 integer, allocatable :: itmp(:)
 character(len=2) :: str2(2)
 real(kind=8), allocatable :: tmp_coeff(:,:), rtmp(:)
 logical :: final_contr

 m = 0 ! an index in shell_type

 do i = 1, natom, 1
  nc = all_pg(i)%nc
  if(nc < 1) then
   write(6,'(/,A)') 'ERROR in subroutine merge_s_and_p_into_sp: nc<1.'
   write(6,'(A)') 'Internal inconsistency or problematic basis set data.'
   stop
  else if(nc == 1) then
   m = m + 1
   cycle
  end if
  final_contr = .false.

  ! the upper limit of the loop (i.e. nc-1) is determined when the loop starts,
  ! and this upper limit does not change even if nc is updated within any loop
  do j = 1, nc-1, 1
   m = m + 1
   if(j == nc) then ! this is possible for O 6-311G(d)
    final_contr = .true.
    exit
   end if
   nline = all_pg(i)%prim_gau(j)%nline
   if(nline /= all_pg(i)%prim_gau(j+1)%nline) cycle

   str2 = [all_pg(i)%prim_gau(j)%stype, all_pg(i)%prim_gau(j+1)%stype]
   if(.not. (str2(1)=='S ' .and. str2(2)=='P ')) cycle

   allocate(rtmp(nline))
   rtmp = DABS(all_pg(i)%prim_gau(j)%coeff(:,1) - all_pg(i)%prim_gau(j+1)%coeff(:,1))

   if( ALL( rtmp<1D-5 ) ) then
    all_pg(i)%prim_gau(j)%stype = 'SP'
    all_pg(i)%prim_gau(j)%ncol = 3
    allocate(tmp_coeff(nline,3), source=0d0)
    tmp_coeff(:,1:2) = all_pg(i)%prim_gau(j)%coeff
    deallocate(all_pg(i)%prim_gau(j)%coeff)
    allocate(all_pg(i)%prim_gau(j)%coeff(nline,3), source=tmp_coeff)
    deallocate(tmp_coeff)
    all_pg(i)%prim_gau(j)%coeff(:,3) = all_pg(i)%prim_gau(j+1)%coeff(:,2)
    deallocate(all_pg(i)%prim_gau(j+1)%coeff)
    all_pg(i)%prim_gau(j+1)%stype = '  '
    all_pg(i)%prim_gau(j+1)%ncol = 0
    all_pg(i)%prim_gau(j+1)%nline = 0

    if(j+1 < nc) then
     do k = j+1, nc-1, 1
      all_pg(i)%prim_gau(k)%stype = all_pg(i)%prim_gau(k+1)%stype
      all_pg(i)%prim_gau(k)%ncol  = all_pg(i)%prim_gau(k+1)%ncol
      all_pg(i)%prim_gau(k)%nline = all_pg(i)%prim_gau(k+1)%nline
      all_pg(i)%prim_gau(k)%coeff = all_pg(i)%prim_gau(k+1)%coeff
     end do ! for k
     all_pg(i)%prim_gau(nc)%stype = '  '
     all_pg(i)%prim_gau(nc)%ncol = 0
     all_pg(i)%prim_gau(nc)%nline = 0
     deallocate(all_pg(i)%prim_gau(nc)%coeff)
    end if

    nc = nc - 1
    all_pg(i)%nc = nc

    ! update integer arrays shell_type and shl2atm
    allocate(itmp(ncontr-1))
    itmp(1:m-1) = shell_type(1:m-1)
    itmp(m) = -1
    itmp(m+1:) = shell_type(m+2:)
    deallocate(shell_type)
    allocate(shell_type(ncontr-1), source=itmp)
    itmp(1:m) = shl2atm(1:m)
    itmp(m+1:) = shl2atm(m+2:)
    deallocate(shl2atm)
    allocate(shl2atm(ncontr-1), source=itmp)
    deallocate(itmp)
    ncontr = ncontr - 1
   end if

   deallocate(rtmp)
   if(j == nc) then ! this is possible for O 6-311G
    final_contr = .true.
    exit
   end if
  end do ! for j

  if(.not. final_contr) m = m + 1
 end do ! for i

 if(m /= ncontr) then
  write(6,'(/,A)') 'ERROR in subroutine merge_s_and_p_into_sp: m/=ncontr. Inter&
                   &nal inconsistency.'
  write(6,'(3(A,I0))') 'natom=', natom, ', m=', m, ', ncontr=', ncontr
  stop
 end if
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
  write(6,'(/,A)') 'ERROR in subroutine read_nbf_and_nif_from_mkl: incomplete f&
                   &ile.'
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
  write(6,'(/,A)') 'ERROR in subroutine read_nbf_and_nif_from_mkl: zero basis f&
                   &unction.'
  write(6,'(A)') 'There must be something wrong in file '//TRIM(mklname)//'.'
  close(fid)
  stop
 end if

 do while(.true.)
  if(INDEX(buf(1:5),'$END') > 0) exit
  nif = nif + detect_ncol_in_buf(buf)
  do i = 1, nbf+1, 1
   read(fid,'(A)') buf
  end do ! for i
  read(fid,'(A)') buf
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

 ang = ' '; shltyp = 0; shl2atm = 0
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
   write(6,'(/,A)') 'ERROR in subroutine read_shltyp_and_shl2atm_from_mkl: wron&
                    &g format in'
   write(6,'(A)') 'file '//TRIM(mklname)
   write(6,'(A)') 'buf='//TRIM(buf)
   close(fid)
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
  case('I')
   shltyp(i) = -6
  case default
   write(6,'(/,A)') 'ERROR in subroutine read_shltyp_and_shl2atm_from_mkl: unsu&
                    &pported angular'
   write(6,'(A)') 'momentum='//ang
   close(fid)
   stop
  end select

  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf)==0 .or. index(buf,'$END')>0) exit
   if(buf(1:2) == '$$') then
    iatom = iatom + 1
   else
    k = detect_ncol_in_buf(buf)
    if(k == 3) exit
   end if
  end do ! for while
 end do ! for i

 close(fid)
end subroutine read_shltyp_and_shl2atm_from_mkl

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

 mo = 0d0
 if(ab=='a' .or. ab=='A') then
  key = key1
 else if (ab=='b' .or. ab=='B') then
  key = key2
 else
  write(6,'(/,A)') 'ERROR in subroutine read_mo_from_mkl: invalid ab.'
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
  write(6,'(/,A)') "ERROR in subroutine read_mo_from_mkl: no '"//TRIM(key)//&
                   "' found in file "//TRIM(mklname)
  close(fid)
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

! delete a primitive exponent and its coefficient when the coefficient is 0
subroutine del_zero_coeff_in_prim_gau(pg)
 implicit none
 integer :: i, j, ncol, nline, ndel
 real(kind=8), parameter :: thres = 1d-9
 real(kind=8), allocatable :: rtmp(:,:)
 type(primitive_gaussian), intent(inout) :: pg

 ncol = pg%ncol
 if(ncol /= 2) return ! possible pg%ncol=3 in the future

 nline = pg%nline
 ndel = 0
 do i = 1, nline, 1
  if(DABS(pg%coeff(i,2)) < thres) ndel = ndel + 1
 end do ! for i

 if(ndel == nline) then
  write(6,'(/,A)') 'ERROR in subroutine del_zero_coeff_in_prim_gau: all exponen&
                   &ts in this sub-'
  write(6,'(A)') 'section will be deleted. Something must be wrong.'
  stop
 end if

 nline = nline - ndel
 allocate(rtmp(nline,ncol), source=0d0)
 j = 0
 do i = 1, pg%nline, 1
  if(DABS(pg%coeff(i,2)) < thres) cycle
  j = j + 1
  rtmp(j,:) = pg%coeff(i,:)
 end do ! for i

 deallocate(pg%coeff)
 allocate(pg%coeff(nline,ncol), source=rtmp)
 deallocate(rtmp)
 pg%nline = nline
end subroutine del_zero_coeff_in_prim_gau

! read the basis set data of all atoms from GAMESS .inp/.dat file
subroutine read_all_pg_from_gms_inp(inpname)
 implicit none
 integer :: i, j, k, m, nc, nline, fid
 integer, external :: detect_ncol_in_buf
 real(kind=8), allocatable :: coeff(:,:)
 character(len=1) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 if(natom == 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_all_pg_from_gms_inp: natom = 0.'
  stop
 end if

 if(.not. allocated(shl2atm)) then
  write(6,'(/,A)') 'ERROR in subroutine read_all_pg_from_gms_inp: array shl2atm&
                   & is not allocated.'
  stop
 end if

 if(allocated(all_pg)) deallocate(all_pg)
 allocate(all_pg(natom))

 call goto_data_section_in_gms_inp(inpname, fid)
 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,'(A)') buf
  nc = COUNT(shl2atm == i)
  all_pg(i)%nc = nc
  allocate(all_pg(i)%prim_gau(nc))
  j = 0

  do while(j <= nc)
   j = j + 1
   read(fid,*) str, nline

   if(str == 'L') then
    allocate(coeff(nline,3))
    do k = 1, nline, 1
     read(fid,*) m, coeff(k,1:3)
    end do ! for k
    ! split into S P, and check normalization in subroutine un_normalized_all_pg
    all_pg(i)%prim_gau(j:j+1)%stype = ['S ','P ']
    all_pg(i)%prim_gau(j:j+1)%nline = [nline,nline]
    all_pg(i)%prim_gau(j:j+1)%ncol = [2,2]
    allocate(all_pg(i)%prim_gau(j)%coeff(nline,2))
    allocate(all_pg(i)%prim_gau(j+1)%coeff(nline,2))
    all_pg(i)%prim_gau(j)%coeff = coeff(:,1:2)
    all_pg(i)%prim_gau(j+1)%coeff(:,1) = coeff(:,1)
    all_pg(i)%prim_gau(j+1)%coeff(:,2) = coeff(:,3)
    deallocate(coeff)
    j = j + 1
   else ! not 'L'
    all_pg(i)%prim_gau(j)%stype(1:1) = str
    all_pg(i)%prim_gau(j)%nline = nline
    all_pg(i)%prim_gau(j)%ncol = 2
    allocate(all_pg(i)%prim_gau(j)%coeff(nline,2), source=0d0)
    do k = 1, nline, 1
     read(fid,*) m, all_pg(i)%prim_gau(j)%coeff(k,1:2)
    end do ! for k
   end if
   if(j == nc) exit
  end do ! for while j

  read(fid,'(A)') buf ! skip the blank line
 end do ! for i

 close(fid)
end subroutine read_all_pg_from_gms_inp

! read the basis set data of all atoms from Amesp .amo file
subroutine read_all_pg_from_amo(amoname, ntype, prim_per_shell)
 implicit none
 integer :: i, j, k, nc, nline, fid
 integer, intent(in) :: ntype, prim_per_shell(ncontr)
 integer, allocatable :: angshl(:)
 character(len=2) :: str2 = '  '
 character(len=2), allocatable :: elem0(:)
 character(len=2), parameter :: i2s(0:6) = ['S ','P ','D ','F ','G ','H ','I ']
 character(len=240) :: buf
 character(len=240), intent(in) :: amoname
 type(pg4atom), allocatable :: all_pg0(:) ! size ntype

 if(natom==0 .or. ncontr==0) then
  write(6,'(/,A)') 'ERROR in subroutine read_all_pg_from_amo: natom or ncontr i&
                   &s 0.'
  stop
 end if

 if(.not. allocated(shl2atm)) then
  write(6,'(/,A)') 'ERROR in subroutine read_all_pg_from_amo: array shl2atm is &
                   &not allocated.'
  stop
 end if

 if(allocated(all_pg)) deallocate(all_pg)
 allocate(all_pg(natom))
 k = 0

 do i = 1, natom, 1
  nc = COUNT(shl2atm == i)
  all_pg(i)%nc = nc
  allocate(all_pg(i)%prim_gau(nc))
  do j = 1, nc, 1
   k = k + 1
   nline = prim_per_shell(k)
   all_pg(i)%prim_gau(j)%nline = nline
   all_pg(i)%prim_gau(j)%ncol = 2
   allocate(all_pg(i)%prim_gau(j)%coeff(nline,2), source=0d0)
  end do ! for j
 end do ! for i

 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:9) == 'MaxAtmBas') exit
 end do ! for while

 allocate(all_pg0(ntype), elem0(ntype))

 do i = 1, ntype, 1
  read(fid,*) elem0(i)
  elem0(i) = ADJUSTL(elem0(i))
  read(fid,'(A)') buf
  j = INDEX(buf,':')
  read(buf(j+1:),*) nc
  all_pg0(i)%nc = nc
  allocate(all_pg0(i)%prim_gau(nc))
  all_pg0(i)%prim_gau(:)%ncol = 2
  allocate(angshl(nc), source=0)
  read(fid,'(A)') buf
  read(fid,*) angshl
  read(fid,'(A)') buf
  read(fid,*) all_pg0(i)%prim_gau(:)%nline
  do j = 1, nc, 1
   nline = all_pg0(i)%prim_gau(j)%nline
   allocate(all_pg0(i)%prim_gau(j)%coeff(nline,2), source=0d0)
  end do ! for j
  forall(j=1:nc) all_pg0(i)%prim_gau(j)%stype = i2s(angshl(j))
  deallocate(angshl)

  read(fid,'(A)') buf
  do j = 1, nc, 1
   read(fid,*) all_pg0(i)%prim_gau(j)%coeff(:,1)
  end do ! for j

  read(fid,'(A)') buf
  do j = 1, nc, 1
   read(fid,*) all_pg0(i)%prim_gau(j)%coeff(:,2)
  end do ! for j

  read(fid,'(A)') buf
  read(buf(7:),*) k
  nline = k/5
  if(k - 5*nline > 0) nline = nline + 1
  do j = 1, nline, 1
   read(fid,'(A)') buf
  end do ! for j
 end do ! for i

 close(fid)
 do i = 1, natom, 1
  str2 = elem(i)
  do j = 1, ntype, 1
   if(elem0(j) == str2) exit
  end do ! for j
  call copy_type_pg4atom(all_pg0(j), all_pg(i))
 end do ! for i

 deallocate(elem0, all_pg0)
end subroutine read_all_pg_from_amo

end module mkl_content


module basis_data ! to store basis set data for an atom
 implicit none
 integer :: ncol(0:7), nline(0:7)

 type :: one_momentum
  real(kind=8), allocatable :: prim_exp(:) ! size nline(i)
  real(kind=8), allocatable :: coeff(:,:)  ! size (ncol(i),nline(i))
 end type one_momentum

 type(one_momentum) :: bas4atom(0:7) ! SPDFGHIJ

contains

! enlarge the arrays prim_exp and/or coeff in type bas4atom
subroutine enlarge_bas4atom(j, k, prim_exp, contr_coeff, enlarge_exp)
 implicit none
 integer :: p, q
 integer, intent(in) :: j, k
 real(kind=8), intent(in) :: prim_exp(k), contr_coeff(k)
 real(kind=8), allocatable :: r1(:), r2(:,:)
 logical, intent(in) :: enlarge_exp

 p = ncol(j); q = nline(j)

 if(enlarge_exp)then
  if(allocated(bas4atom(j)%prim_exp)) then
   allocate(r1(q-k), source=bas4atom(j)%prim_exp)
   deallocate(bas4atom(j)%prim_exp)
   allocate(bas4atom(j)%prim_exp(q))
   bas4atom(j)%prim_exp(1:q-k) = r1
   deallocate(r1)
   bas4atom(j)%prim_exp(q-k+1:q) = prim_exp
  else
   allocate(bas4atom(j)%prim_exp(k), source=prim_exp)
  end if

  if(allocated(bas4atom(j)%coeff)) then
   allocate(r2(p-1,q-k), source=bas4atom(j)%coeff)
   deallocate(bas4atom(j)%coeff)
   allocate(bas4atom(j)%coeff(p,q), source=0d0)
   bas4atom(j)%coeff(1:p-1,1:q-k) = r2
   deallocate(r2)
   bas4atom(j)%coeff(p,q-k+1:q) = contr_coeff
  else
   allocate(bas4atom(j)%coeff(1,k))
   bas4atom(j)%coeff(1,:) = contr_coeff
  end if

 else
  allocate(r2(p-1,q), source=bas4atom(j)%coeff)
  deallocate(bas4atom(j)%coeff)
  allocate(bas4atom(j)%coeff(p,q), source=0d0)
  bas4atom(j)%coeff(1:p-1,:) = r2
  deallocate(r2)
  bas4atom(j)%coeff(p,:) = contr_coeff
 end if
end subroutine enlarge_bas4atom

! enlarge the arrays prim_exp and/or coeff in type bas4atom, only for L or SP,
! i.e. Pople-type basis set
subroutine enlarge_bas4atom_sp(k, prim_exp, contr_coeff, contr_coeff_sp)
 implicit none
 integer :: i, p(0:1), q(0:1)
 integer, intent(in) :: k
 real(kind=8), intent(in) :: prim_exp(k), contr_coeff(k), contr_coeff_sp(k)
 real(kind=8), allocatable :: r1(:), r2(:,:), coeff(:,:)

 allocate(coeff(k,0:1))
 coeff(:,0) = contr_coeff; coeff(:,1) = contr_coeff_sp
 p = ncol(0:1); q = nline(0:1)

 do i = 0, 1 ! 0/1 for S/P, respectively
  if(allocated(bas4atom(i)%prim_exp)) then
   allocate(r1(q(i)-k), source=bas4atom(i)%prim_exp)
   deallocate(bas4atom(i)%prim_exp)
   allocate(bas4atom(i)%prim_exp(q(i)))
   bas4atom(i)%prim_exp(1:q(i)-k) = r1
   deallocate(r1)
   bas4atom(i)%prim_exp(q(i)-k+1:q(i)) = prim_exp
  else
   allocate(bas4atom(i)%prim_exp(k), source=prim_exp)
  end if

  if(allocated(bas4atom(i)%coeff)) then
   allocate(r2(p(i)-1,q(i)-k), source=bas4atom(i)%coeff)
   deallocate(bas4atom(i)%coeff)
   allocate(bas4atom(i)%coeff(p(i),q(i)), source=0d0)
   bas4atom(i)%coeff(1:p(i)-1,1:q(i)-k) = r2
   deallocate(r2)
   bas4atom(i)%coeff(p(i),q(i)-k+1:q(i)) = coeff(:,i)
  else
   allocate(bas4atom(i)%coeff(1,k))
   bas4atom(i)%coeff(1,:) = coeff(:,i)
  end if
 end do ! for i

 deallocate(coeff)
end subroutine enlarge_bas4atom_sp

! deallocate arrays in type bas4atom
subroutine clear_bas4atom
 implicit none
 integer :: i

 do i = 0, 7
  if(allocated(bas4atom(i)%prim_exp)) deallocate(bas4atom(i)%prim_exp)
  if(allocated(bas4atom(i)%coeff)) deallocate(bas4atom(i)%coeff)
 end do ! for i
end subroutine clear_bas4atom

end module basis_data


! print basis set and ECP(if any) data into GENBAS and ECPDATA
subroutine prt_cfour_genbas(ecp, mrcc)
 use fch_content
 use basis_data
 implicit none
 integer :: i, j, k, n, n1, n2, i1, i2, highest, iatom, fid
 character(len=1) :: str = ' '
 character(len=2) :: str2 = '  '
 character(len=1), parameter :: am_type1(0:6) = ['s','p','d','f','g','h','i']
 character(len=57), parameter :: c1 = 'generated by MOKIT(not necessarily PVTZ,&
                                      & borrowed string)'
 character(len=63), parameter :: c2 = 'generated by MOKIT(not necessarily ECP-1&
                                      &0-MDF, borrowed string)'
 logical :: cycle_atom
 logical, intent(in) :: ecp, mrcc

 open(newunit=fid,file='GENBAS',status='replace')
 iatom = 1; i1 = 1; i2 = 1 ! initialization

 do while(.true.)
  ncol = 0; nline = 0

  do i = i1, ncontr, 1
   if(shell2atom_map(i) == iatom+1) exit
   j = shell_type(i); k = prim_per_shell(i)

   if(j == -1) then ! L or SP
    ncol(0:1) = ncol(0:1) + 1
    nline(0:1) = nline(0:1) + k
    call enlarge_bas4atom_sp(k, prim_exp(i2:i2+k-1), contr_coeff(i2:i2+k-1), &
                                                  contr_coeff_sp(i2:i2+k-1))
   else ! j /= -1
    if(j < 0) j = -j ! -2,-3,... -> 2,3,...
    ncol(j) = ncol(j) + 1
    if(i == i1) then
     nline(j) = k
     allocate(bas4atom(j)%prim_exp(k), source=prim_exp(i2:i2+k-1))
     allocate(bas4atom(j)%coeff(1,k),source=0d0)
     bas4atom(j)%coeff(1,:) = contr_coeff(i2:i2+k-1)
    else ! i > i1
     if(k /= prim_per_shell(i-1)) then
      nline(j) = nline(j) + k
      call enlarge_bas4atom(j,k,prim_exp(i2:i2+k-1),contr_coeff(i2:i2+k-1),.true.)
     else ! k = prim_per_shell(i-1)
      if( ANY( DABS(prim_exp(i2:i2+k-1)-prim_exp(i2-k:i2-1)) >1d-5 ) ) then
       nline(j) = nline(j) + k
       call enlarge_bas4atom(j,k,prim_exp(i2:i2+k-1),contr_coeff(i2:i2+k-1),.true.)
      else ! identical primitive exponents
       call enlarge_bas4atom(j,k,prim_exp(i2:i2+k-1),contr_coeff(i2:i2+k-1),.false.)
      end if
     end if
    end if
   end if

   i2 = i2 + k
  end do ! for i

  ! One cannot use `highest = shell_type(i-1)` here, since it would lead to wrong
  ! result if diffuse functions are used, e.g. maug-cc-pVTZ. The shell_type(i-1)
  ! is not necessary the one which has the highest angular momentum.
  highest = MAXVAL(IABS(shell_type(i1:i-1)))
  i1 = i   ! remember to update i1
  if(highest > 7) then
   write(6,'(/,A)') 'ERROR in subroutine prt_cfour_genbas: angular momentum too&
                    & high. Not supported!'
   close(fid)
   stop
  end if

  cycle_atom = .false.
  if(iatom > 1) then
   if(ANY(elem(1:iatom-1) == elem(iatom))) cycle_atom = .true.
  end if

  if(.not. cycle_atom) then
   str2 = elem(iatom)
   if(mrcc) then ! MRCC
    write(fid,'(A)') TRIM(str2)//':PVTZ'
   else          ! CFOUR
    if(str2(2:2) /= ' ') call upper(str2(2:2))
    if(ecp) then
     if(LPSkip(iatom) == 0) then
      write(fid,'(A)') TRIM(str2)//':ECP-10-MDF'
     else
      write(fid,'(A)') TRIM(str2)//':PVTZ'
     end if
    else
     write(fid,'(A)') TRIM(str2)//':PVTZ'
    end if
   end if
   write(fid,'(A)') c1
   write(fid,'(/,I3)') highest+1
   write(fid,'(9I5)') (i, i=0,highest)
   write(fid,'(9I5)') ncol(0:highest)
   write(fid,'(9I5)') nline(0:highest)
   write(fid,'(/)',advance='no')

   do i = 0, highest, 1
    write(fid,'(5(1X,ES15.8))') bas4atom(i)%prim_exp
    write(fid,'(/)',advance='no')

    do j = 1, nline(i), 1
     write(fid,'(16(1X,ES15.8))') bas4atom(i)%coeff(1:ncol(i),j)
    end do ! for j
    write(fid,'(/)',advance='no')
   end do ! for i
  end if

  call clear_bas4atom()
  if(iatom == natom) exit
  iatom = iatom + 1
 end do ! for while

 close(fid)
 deallocate(ielem, prim_per_shell, prim_exp, contr_coeff)
 if(allocated(contr_coeff_sp)) deallocate(contr_coeff_sp)
 if(.not. ecp) return

 if(mrcc) then
  ! MRCC read ECP/PP data also from GENBAS
  open(newunit=fid,file='GENBAS',status='old',position='append')
 else
  open(newunit=fid,file='ECPDATA',status='replace')
 end if

 do i = 1, natom, 1
  if(LPSkip(i) /= 0) cycle
  str2 = elem(i)
  if(str2(2:2) /= ' ') call upper(str2(2:2))
  write(fid,'(A)') TRIM(str2)//':ECP-10-MDF'
  write(fid,'(A,/,A)') c2, '*'
  write(fid,'(4X,A,I3,4X,A,I2)') 'NCORE =', NINT(RNFroz(i)), 'LMAX =', LMax(i)
  str = am_type1(LMax(i))

  do j = 1, 10, 1
   n1 = KFirst(i,j); n2 = KLast(i,j)
   if(n1 == 0) exit
   if(j == 1) then
    write(fid,'(A)') str
   else
    write(fid,'(A)') am_type1(j-2)//'-'//str
   end if
   do n = n1, n2, 1
    write(fid,'(2X,ES15.8,4X,I0,2X,ES15.8)') CLP(n), NLP(n), ZLP(n)
   end do ! for n
  end do ! for j

  write(fid,'(A)') '*'
 end do ! for i

 close(fid)
 deallocate(KFirst, KLast, Lmax, LPSkip, NLP, RNFroz, CLP, CLP2, ZLP, elem)
end subroutine prt_cfour_genbas

! get/find basis function marks from the integer array shell_type
subroutine read_bas_mark_from_shltyp(ncontr,shell_type, nfmark, ngmark, nhmark,&
                                     nimark, f_mark, g_mark, h_mark, i_mark)
 implicit none
 integer :: i, k
 integer, intent(in) :: ncontr
 integer, intent(in) :: shell_type(ncontr)
 integer, intent(out) :: nfmark, ngmark, nhmark, nimark, f_mark(ncontr), &
  g_mark(ncontr), h_mark(ncontr), i_mark(ncontr)

 nfmark = 0; ngmark = 0; nhmark = 0; nimark = 0
 f_mark = 0; g_mark = 0; h_mark = 0; i_mark = 0
 k = 0

 do i = 1, ncontr, 1
  select case(shell_type(i))
  case(0) ! S
   k = k + 1
  case(1) ! P
   k = k + 3
  case(-1) ! L
   k = k + 4
  case(-2) ! D
   k = k + 5
  case(-3) ! F
   nfmark = nfmark + 1
   f_mark(nfmark) = k + 1
   k = k + 7
  case(-4) ! G
   ngmark = ngmark + 1
   g_mark(ngmark) = k + 1
   k = k + 9
  case(-5) ! H
   nhmark = nhmark + 1
   h_mark(nhmark) = k + 1
   k = k + 11
  case(-6) ! I
   nimark = nimark + 1
   i_mark(nimark) = k + 1
   k = k + 13
  end select
 end do ! for i
end subroutine read_bas_mark_from_shltyp

! Update MO coefficients using basis function marks since some MO coefficients
! in ORCA are negative to those in Gaussian
subroutine update_mo_using_bas_mark(nbf, nif, nfmark, ngmark, nhmark, nimark, &
                                ncontr, f_mark, g_mark, h_mark, i_mark, coeff)
 implicit none
 integer :: i, k
 integer, intent(in) :: nbf, nif, nfmark, ngmark, nhmark, nimark, ncontr
 integer, intent(in) :: f_mark(ncontr), g_mark(ncontr), h_mark(ncontr), &
  i_mark(ncontr)
 real(kind=8), intent(inout) :: coeff(nbf,nif)

 do i = 1, nfmark, 1
  k = f_mark(i)
  coeff(k+5,:) = -coeff(k+5,:)
  coeff(k+6,:) = -coeff(k+6,:)
 end do ! for i

 do i = 1, ngmark, 1
  k = g_mark(i)
  coeff(k+5,:) = -coeff(k+5,:)
  coeff(k+6,:) = -coeff(k+6,:)
  coeff(k+7,:) = -coeff(k+7,:)
  coeff(k+8,:) = -coeff(k+8,:)
 end do ! for i

 do i = 1, nhmark, 1
  k = h_mark(i)
  coeff(k+5,:) = -coeff(k+5,:)
  coeff(k+6,:) = -coeff(k+6,:)
  coeff(k+7,:) = -coeff(k+7,:)
  coeff(k+8,:) = -coeff(k+8,:)
 end do ! for i

 do i = 1, nimark, 1
  k = i_mark(i)
  coeff(k+5,:) = -coeff(k+5,:)
  coeff(k+6,:) = -coeff(k+6,:)
  coeff(k+7,:) = -coeff(k+7,:)
  coeff(k+8,:) = -coeff(k+8,:)
 end do ! for i
end subroutine update_mo_using_bas_mark

! find nprim from type all_pg
subroutine find_nprim_from_all_pg(ncontr, prim_per_shell, nprim, has_sp)
 use mkl_content, only: natom, all_pg
 implicit none
 integer :: i, j, k, nc, nline
 integer, intent(in) :: ncontr
 integer, intent(out) :: nprim, prim_per_shell(ncontr)
 logical, intent(out) :: has_sp

 has_sp = .false.
 nprim = 0; prim_per_shell = 0; k = 0

 do i = 1, natom, 1
  nc = all_pg(i)%nc
  do j = 1, nc, 1
   nline = all_pg(i)%prim_gau(j)%nline
   prim_per_shell(k+1) = nline
   k = k + 1
   nprim = nprim + nline
   if(.not. has_sp) then
    if(all_pg(i)%prim_gau(j)%ncol == 3) has_sp = .true.
   end if
  end do ! for j
 end do ! for i
end subroutine find_nprim_from_all_pg

subroutine all_pg2prim_exp_and_contr_coeff(has_sp)
 use fch_content, only: nprim, prim_exp, contr_coeff, contr_coeff_sp
 use mkl_content, only: natom, all_pg
 implicit none
 integer :: i, j, k, nc, nline, ncol
 logical, intent(in) :: has_sp

 allocate(prim_exp(nprim), source=0d0)
 allocate(contr_coeff(nprim), source=0d0)
 if(has_sp) allocate(contr_coeff_sp(nprim), source=0d0)

 ! convert type all_pg to arrays prim_exp, contr_coeff, and contr_coeff_sp
 k = 0
 do i = 1, natom, 1
  nc = all_pg(i)%nc
  do j = 1, nc, 1
   nline = all_pg(i)%prim_gau(j)%nline
   ncol = all_pg(i)%prim_gau(j)%ncol
   prim_exp(k+1:k+nline) = all_pg(i)%prim_gau(j)%coeff(:,1)
   contr_coeff(k+1:k+nline) = all_pg(i)%prim_gau(j)%coeff(:,2)
   if(ncol==3) contr_coeff_sp(k+1:k+nline) = all_pg(i)%prim_gau(j)%coeff(:,3)
   k = k + nline
  end do ! for j
 end do ! for i
end subroutine all_pg2prim_exp_and_contr_coeff

! print MO coefficients and orbital energies in a given file ID (.mkl file)
subroutine prt_mo_and_e_in_mkl(fid, nbf, nif, mo, ev)
 implicit none
 integer :: i, j, k, m
 integer, intent(in) :: fid, nbf, nif
 real(kind=8), intent(in) :: mo(nbf,nif), ev(nif)

 k = 0
 do while(.true.)
  if(k+1 > nif) exit
  if(k+5 > nif) then
   j = nif
  else
   j = k + 5
  end if
  write(fid,'(5(A4,1X))') (' a1g', i=k+1,j)
  write(fid,'(5(F14.8,1X))') (ev(i),i=k+1,j)
  do i = 1, nbf, 1
   write(fid,'(5(ES15.8,1X))') (mo(i,m),m=k+1,j)
  end do ! for i
  k = j
 end do ! for while
end subroutine prt_mo_and_e_in_mkl

! Convert an RHF .mkl file into a UHF one by copying alpha MOs to beta MOs.
! The input file will be updated. If brokensym is .True., the beta HOMO LUMO
! will be interchanged to make a broken symmetry guess.
subroutine mkl_r2u(mklname, brokensym)
 use mkl_content, only: read_nbf_and_nif_from_mkl, read_mo_from_mkl
 implicit none
 integer :: i, nbf, nif, ndb, fid
 real(kind=8), parameter :: thres = 1d-6
 real(kind=8), allocatable :: mo_a(:,:), on_a(:), e_a(:), rtmp(:)
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname
!f2py intent(in) :: mklname
 logical, intent(in) :: brokensym
!f2py intent(in) :: brokensym

 call read_nbf_and_nif_from_mkl(mklname, nbf, nif)
 allocate(mo_a(nbf,nif), on_a(nif), e_a(nif))
 call read_mo_from_mkl(mklname, nbf, nif, 'a', mo_a)
 call read_ev_from_mkl(mklname, nif, 'a', e_a)
 call read_on_from_mkl(mklname, nif, 'a', on_a)

 ndb = 0
 do i = 1, nif, 1
  if(DABS(on_a(i) - 2d0) < thres) then
   ndb = ndb + 1
  else if(DABS(on_a(i)) < thres) then
   exit
  end if
 end do ! for i

 if(brokensym) then
  allocate(rtmp(nbf), source=mo_a(:,ndb))
  mo_a(:,ndb) = mo_a(:,ndb+1)
  mo_a(:,ndb+1) = rtmp
  rtmp(1) = e_a(ndb)
  e_a(ndb) = e_a(ndb+1)
  e_a(ndb+1) = rtmp(1)
  deallocate(rtmp)
 end if
 on_a = on_a*0.5d0

 open(newunit=fid,file=TRIM(mklname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(1:6) == '$OCC_A') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') 'ERROR in subroutine mkl_r2u: problematic file '//&
                   TRIM(mklname)
  stop
 end if

 write(fid,'(5(F12.7,1X))') (on_a(i), i=1,nif)
 write(fid,'(A,/)') '$END'

 write(fid,'(A)') '$COEFF_BETA'
 call prt_mo_and_e_in_mkl(fid, nbf, nif, mo_a, e_a)
 write(fid,'(A)') '$END'
 write(fid,'(/,A)') '$OCC_BETA'
 write(fid,'(5(F12.7,1X))') (on_a(i), i=1,nif)
 write(fid,'(A,/)') '$END'
 close(fid)

 deallocate(mo_a, on_a, e_a)
end subroutine mkl_r2u

