! written by jxzou at 20230209: transfer MOs from Gaussian to CFOUR
! limitation: one element with different basis sets is not supported

! currently the generated internal coordinates and basis function order are
! correct, but need to use the rotated coordinates in CFOUR

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i
 character(len=240) :: fchname

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in subroutine fch2cfour: wrong command line argument!'
  write(6,'(A)')   ' Example: fch2cfour water.fch'
  write(6,'(A,/)') ' Note: fch2cfour can only be used for CFOUR >= v2.1'
  stop
 end if

 fchname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 call fch2cfour(fchname)
end program main

subroutine fch2cfour(fchname)
 use fch_content
 implicit none
 integer :: i, j, k, m, length, nbf0, nbf1, nif1, icart
 integer :: n5dmark, n7fmark, n9gmark, n11hmark
 integer :: n6dmark, n10fmark, n15gmark, n21hmark
 integer, allocatable :: idx(:), d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 character(len=240), intent(in) :: fchname
 real(kind=8), allocatable :: coeff(:,:), coeff2(:,:), norm(:)
 logical :: uhf, sph, ecp

 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)
 call find_irel_in_fch(fchname, irel)

 sph = .true.
 call find_icart_in_fch(fchname, .false., icart)
 if(icart == 2) sph = .false.

 uhf = .false.; ecp = .false.
 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 call read_fch(fchname, uhf)
 nbf0 = nbf ! make a copy of nbf
 if(LenNCZ > 0) ecp = .true.

 ! generate the file ZMAT
 call prt_cfour_zmat(natom, elem, coor, charge, mult, irel, uhf, sph, ecp)
 deallocate(coor)

 ! generate files GENBAS and ECPDATA(if needed)
 call prt_cfour_genbas(ecp, .false.)

 if(uhf) then ! UHF
  allocate(coeff(nbf,2*nif))
  coeff(:,1:nif) = alpha_coeff
  coeff(:,nif+1:2*nif) = beta_coeff
  deallocate(alpha_coeff, beta_coeff)
  nif1 = 2*nif
 else         ! R(O) HF
  allocate(coeff(nbf,nif), source=alpha_coeff)
  deallocate(alpha_coeff)
  nif1 = nif
 end if

 ! enlarge arrays shell_type and shell2atom_map, using d_mark as a tmp array
 k = ncontr
 allocate(d_mark(k), source=shell_type)
 deallocate(shell_type)
 allocate(shell_type(2*k), source=0)
 shell_type(1:k) = d_mark

 d_mark = shell2atom_map
 deallocate(shell2atom_map)
 allocate(shell2atom_map(2*k), source=0)
 shell2atom_map(1:k) = d_mark
 deallocate(d_mark)

! first we adjust the basis functions in each MO according to the Shell to atom map
! this is to ensure that D comes after L functions
 ! split the 'L' into 'S' and 'P'
 call split_L_func(k, shell_type, shell2atom_map, length)
 allocate(idx(nbf))
 forall(i = 1:nbf) idx(i) = i
 allocate(norm(nbf), source=1d0)

 ! sort the shell_type, shell_to_atom_map by ascending order
 ! MOs will be adjusted accordingly
 call sort_shell_and_mo_idx(length, shell_type, shell2atom_map, nbf, idx)
! adjust done

 ! record the indices of d, f, g and h functions
 k = length  ! update k
 allocate(d_mark(k), f_mark(k), g_mark(k), h_mark(k))

 if(sph) then
  call read_mark_from_shltyp_sph(k, shell_type, n5dmark, n7fmark, n9gmark, &
                                 n11hmark, d_mark, f_mark, g_mark, h_mark)
  ! adjust the order of 5d, 7f, etc. functions
  call fch2cfour_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, d_mark, &
                             f_mark, g_mark, h_mark, nbf, idx, norm)
 else
  call read_mark_from_shltyp_cart(k, shell_type, n6dmark, n10fmark, n15gmark,&
                                  n21hmark, d_mark, f_mark, g_mark, h_mark)
  ! adjust the order of 6d, 10f, etc. functions
  call fch2inporb_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, &
                           d_mark, f_mark, g_mark, h_mark, nbf, idx, norm)
 end if

 deallocate(d_mark, f_mark, g_mark, h_mark)

! move the 2nd, 3rd, ... Zeta basis functions forward
 i = 0
 nbf = 0
 do while(i < k)
  i = i + 1
  j = shell2atom_map(i)
  m = shell_type(i)
  nbf1 = nbf
  select case(m)
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   nbf = nbf + 3
  case(-2)   ! 5D
   nbf = nbf + 5
  case( 2)   ! 6D
   nbf = nbf + 6
  case(-3)   ! 7F
   nbf = nbf + 7
  case( 3)   ! 10F
   nbf = nbf + 10
  case(-4)   ! 9G
   nbf = nbf + 9
  case( 4)   ! 15G
   nbf = nbf + 15
  case(-5)   ! 11H
   nbf = nbf + 11
  case( 5)   ! 21H
   nbf = nbf + 21
  end select
  if(m == 0) cycle

  length = 1
  do while(i < k)
   i = i + 1
   if(shell_type(i) /= m) exit
   if(shell2atom_map(i) /= j) exit
   length = length + 1
  end do ! for while

  if(i < k) i = i - 1
  if(length > 1) then
   call zeta_mv_forwd_idx(nbf1, m, length, nbf0, idx, norm)
   nbf = nbf1 + length*(nbf-nbf1)
  end if
 end do ! for while

 deallocate(shell_type, shell2atom_map)
! move done

 nbf = nbf0
 allocate(coeff2(nbf,nif1), source=coeff)
 forall(i=1:nbf, j=1:nif1) coeff(i,j) = coeff2(idx(i),j)*norm(i)
 deallocate(norm, coeff2)

 allocate(alpha_coeff(nbf,nif))
 alpha_coeff = coeff(:,1:nif)
 if(uhf) then
  allocate(beta_coeff(nbf,nif))
  beta_coeff = coeff(:,nif+1:2*nif)
 end if
 deallocate(coeff)

 ! create/print CFOUR orbital file OLDMOS
 call prt_cfour_oldmos(nbf, nif, alpha_coeff, .false.)

 if(uhf) then
  call prt_cfour_oldmos(nbf, nif, beta_coeff, .true.)
 else
  if(mult /= 1) call prt_cfour_oldmos(nbf, nif, alpha_coeff, .true.)
 end if
end subroutine fch2cfour

! print/create/write the CFOUR input file ZMAT
subroutine prt_cfour_zmat(natom, elem, coor, charge, mult, irel, uhf, sph, ecp)
 use fch_content, only: LPSkip
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom, charge, mult, irel
 real(kind=8), intent(in) :: coor(3,natom)
 real(kind=8), external :: norm, ang, dih
 character(len=2) :: str
 character(len=2), intent(in) :: elem(natom)
 logical, intent(in) :: uhf, sph, ecp

 str = '  '
 open(newunit=fid,file='ZMAT',status='replace')
 write(fid,'(A)') 'Generated by fch2cfour of MOKIT'
 write(fid,'(A)') TRIM(elem(1))
 if(natom > 1) write(fid,'(A)') TRIM(elem(2))//' 1 B1'
 if(natom > 2) write(fid,'(A)') TRIM(elem(3))//' 1 B2 2 A1'

 do i = 4, natom, 1
  write(fid,'(3(A,I0))') TRIM(elem(i))//' 1 B',i-1,' 2 A',i-2,' 3 D',i-3
 end do ! for i

 write(fid,'(/)',advance='no')

 do i = 2, natom, 1
  write(fid,'(A,I0,A,F13.7)') 'B',i-1,'=',norm(coor(:,1)-coor(:,i))
 end do ! for i

 do i = 3, natom, 1
  write(fid,'(A,I0,A,F13.7)') 'A',i-2,'=',ang(coor(:,2), coor(:,1), coor(:,i))
 end do ! for i

 do i = 4, natom, 1
  write(fid,'(A,I0,A,F13.7)') 'D',i-3,'=',dih(coor(:,3),coor(:,1),coor(:,2),coor(:,i))
 end do ! for i

 write(fid,'(/,A)',advance='no') '*CFOUR(CALC=SCF,REF='
 if(uhf) then
  write(fid,'(A)',advance='no') 'U'
 else
  if(mult == 1) then
   write(fid,'(A)',advance='no') 'R'
  else
   write(fid,'(A)',advance='no') 'RO'
  end if
 end if
 write(fid,'(2(A,I0))') 'HF,SYM=OFF,CHARGE=', charge, ',MULTI=', mult

 if(.not. sph) write(fid,'(A)') 'SPHERICAL=OFF'
 select case(irel)
 case(-3)
  write(fid,'(A)') 'RELATIVISTIC=X2C1E'
 case(-1) ! do nothing
 case default
  write(6,'(/,A)') 'ERROR in subroutine prt_cfour_zmat: irel out of range!'
  write(6,'(A)') 'CFOUR does not support this type of Hamiltonian.'
  close(fid)
  stop
 end select

 if(ecp) then
  write(fid,'(A,/)') 'BASIS=SPECIAL,ECP=ON)'
 else
  write(fid,'(A,/)') 'BASIS=PVTZ)'
 end if

 if(ecp) then
  do i = 1, natom, 1
   if(LPSkip(i) == 0) then
    str = elem(i)
    if(str(2:2) /= ' ') call upper(str(2:2))
    write(fid,'(A)') TRIM(str)//':ECP-10-MDF'
   else
    write(fid,'(A)') TRIM(elem(i))//':PVTZ'
   end if
  end do ! for i

  write(fid,'(/)',advance='no')

  do i = 1, natom, 1
   if(LPSkip(i) == 0) then
    str = elem(i)
    if(str(2:2) /= ' ') call upper(str(2:2))
    write(fid,'(A)') TRIM(str)//':ECP-10-MDF'
   else
    write(fid,'(A)') TRIM(elem(i))//':NONE'
   end if
  end do ! for i

  write(fid,'(/)',advance='no')
 end if

 close(fid)
end subroutine prt_cfour_zmat

! print CFOUR orbital file OLDMOS
subroutine prt_cfour_oldmos(nbf, nif, coeff, append)
 implicit none
 integer :: i, j, k, nb, fid
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: coeff(nbf,nif)
 character(len=9) :: str
 logical, intent(in) :: append

 if(append) then ! usually for printing Beta MOs
  open(newunit=fid,file='OLDMOS',status='old',position='append')
 else            ! usually for printing Alpha MOs
  open(newunit=fid,file='OLDMOS',status='replace')
 end if

 nb = nif/4
 do i = 1, nb, 1
  write(fid,'(4E30.20)') ((coeff(j,k),k=4*i-3,4*i,1),j=1,nbf,1)
 end do ! for i

 i = nif - 4*nb
 if(i > 0) then
  str = ' '
  write(str,'(A1,I1,A7)') '(',i,'E30.20)'
  write(unit=fid,fmt=str) ((coeff(j,k),k=nif-i+1,nif,1),j=1,nbf,1)
 end if

 close(fid)
end subroutine prt_cfour_oldmos

subroutine fch2cfour_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, &
                           d_mark, f_mark, g_mark, h_mark, nbf, idx, norm)
 implicit none
 integer :: i, j
 integer, intent(in) :: n5dmark, n7fmark, n9gmark, n11hmark, k, nbf
 integer, intent(in) :: d_mark(k), f_mark(k), g_mark(k), h_mark(k)
 integer, intent(inout) :: idx(nbf)
 real(kind=8), intent(out) :: norm(nbf)

 do i = 1, n5dmark, 1
  j = d_mark(i)
  call fch2cfour_permute_5d(idx(j:j+4), norm(j:j+4))
 end do ! for i

 do i = 1, n7fmark, 1
  j = f_mark(i)
  call fch2cfour_permute_7f(idx(j:j+6), norm(j:j+6))
 end do ! for i

 do i = 1, n9gmark, 1
  j = g_mark(i)
  call fch2cfour_permute_9g(idx(j:j+8), norm(j:j+8))
 end do ! for i

 do i = 1, n11hmark, 1
  j = h_mark(i)
  call fch2cfour_permute_11h(idx(j:j+10), norm(j:j+10))
 end do ! for i
end subroutine fch2cfour_permute_sph

subroutine fch2cfour_permute_5d(idx, norm)
 use root_parameter, only: root3
 implicit none
 integer :: i, idx0(5)
 integer, parameter :: order(5) = [1,5,2,4,3]
 integer, intent(inout) :: idx(5)
 real(kind=8) :: norm0(5)
 real(kind=8), intent(inout) :: norm(5)
! From: the order of spherical d functions in Gaussian
! To: the order of spherical d functions in Molpro
! 1    2    3    4    5
! d0 , d+1, d-1, d+2, d-2
! d0 , d2-, d1+, d2+, d1-

 idx0 = idx
 norm(1) = 0.5d0*norm(1)/root3
 norm(4) = 0.5d0*norm(4)
 norm0 = norm

 forall(i = 1:5)
  idx(i) = idx0(order(i))
  norm(i) = norm0(order(i))
 end forall
end subroutine fch2cfour_permute_5d

subroutine fch2cfour_permute_7f(idx, norm)
 use root_parameter, only: root6, root10, root15
 implicit none
 integer :: i, idx0(7)
 integer, parameter :: order(7) = [2,3,1,6,5,7,4]
 integer, intent(inout) :: idx(7)
 real(kind=8) :: norm0(7)
 real(kind=8), intent(inout) :: norm(7)
 real(kind=8), parameter :: r1 = 2d0*root15
 real(kind=8), parameter :: r2 = 2d0*root10
 real(kind=8), parameter :: r3 = 2d0*root6
 real(kind=8), parameter :: ratio(7) = [r1, r2, r2, 2d0, 1d0, r3, r3]
! From: the order of spherical f functions in Gaussian
! To: the order of spherical f functions in Molpro
! 1    2    3    4    5    6    7
! f0 , f+1, f-1, f+2, f-2, f+3, f-3
! f1+, f1-, f0 , f3+, f2-, f3-, f2+

 idx0 = idx
 forall(i = 1:7) norm(i) = norm(i)/ratio(i)
 norm0 = norm

 forall(i = 1:7)
  idx(i) = idx0(order(i))
  norm(i) = norm0(order(i))
 end forall
end subroutine fch2cfour_permute_7f

subroutine fch2cfour_permute_9g(idx, norm)
 use root_parameter, only: root3, root6, root12, root21, root42, root105
 implicit none
 integer :: i, idx0(9)
 integer, parameter :: order(9) = [1,5,2,8,7,4,9,6,3]
 integer, intent(inout) :: idx(9)
 real(kind=8) :: norm0(9)
 real(kind=8), intent(inout) :: norm(9)
 real(kind=8), parameter :: r1 = 2d0*root42
 real(kind=8), parameter :: r2 = 2d0*root21
 real(kind=8), parameter :: r3 = 2d0*root6
 real(kind=8), parameter :: ratio(9) = [8d0*root105, r1, r1, 2d0*r2, r2, r3, &
                                        r3, 4d0*root12, 2d0*root3]
! From: the order of spherical g functions in Gaussian
! To: the order of spherical g functions in CFOUR
! 1    2    3    4    5    6    7    8    9
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4
! g0 , g2-, g1+, g4+, g3-, g2+, g4-, g3+, g1-

 idx0 = idx
 forall(i = 1:9) norm(i) = norm(i)/ratio(i)
 norm0 = norm

 forall(i = 1:9)
  idx(i) = idx0(order(i))
  norm(i) = norm0(order(i))
 end forall
end subroutine fch2cfour_permute_9g

subroutine fch2cfour_permute_11h(idx, norm)
 use root_parameter, only: root3, root6, root7, root30, root105
 implicit none
 integer :: i, idx0(11)
 integer, parameter :: order(11) = [2,3,4,6,9,11,8,7,1,10,5]
 integer, intent(inout) :: idx(11)
 real(kind=8) :: norm0(11)
 real(kind=8), intent(inout) :: norm(11)
 real(kind=8), parameter :: r1 = 24d0*root7
 real(kind=8), parameter :: r2 = 24d0*root6
 real(kind=8), parameter :: r3 = 2d0*root3
 real(kind=8), parameter :: r4 = 8d0*root30
 real(kind=8), parameter :: ratio(11) = [24d0*root105, r1, r1, 12d0, 6d0, r2, &
  r2, 4d0*r3, r3, r4, r4]
! From: the order of spherical h functions in Gaussian
! To: the order of spherical h functions in CFOUR
! 1    2    3    4    5    6    7    8    9    10   11
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5
! h1+, h1-, h2+, h3+, h4-, h5-, h4+, h3-, h0 , h5+, h2-

 idx0 = idx
 forall(i = 1:11) norm(i) = norm(i)/ratio(i)
 norm0 = norm

 forall(i = 1:11)
  idx(i) = idx0(order(i))
  norm(i) = norm0(order(i))
 end forall
end subroutine fch2cfour_permute_11h

! calculate the norm of a 3D vector
function norm(c)
 implicit none
 real(kind=8) :: norm
 real(kind=8), intent(in) :: c(3)

 norm = DSQRT(DOT_PRODUCT(c,c))
end function norm

! calculate bond angle
function ang(c1, c2, c3)
 use Sdiag_parameter, only: PI
 implicit none
 real(kind=8) :: v1(3), v2(3), ang
 real(kind=8), external :: norm
 real(kind=8), intent(in) :: c1(3), c2(3), c3(3)

 v1 = c1 - c2
 v2 = c3 - c2
 ang = DOT_PRODUCT(v1,v2)/(norm(v1)*norm(v2))
 ang = MAX(-1d0, MIN(ang, 1d0)) ! deal with numerical error
 ang = 180d0*DACOS(ang)/PI
end function ang

! calculate the 4,6,2,5,3ctor (normalized to 1) of two vectors
subroutine normal_vector(v1, v2, nv)
 implicit none
 real(kind=8), external :: norm
 real(kind=8), intent(in) :: v1(3), v2(3)
 real(kind=8), intent(out) :: nv(3)

 nv(1) = v1(2)*v2(3) - v1(3)*v2(2)
 nv(2) = v1(3)*v2(1) - v1(1)*v2(3)
 nv(3) = v1(1)*v2(2) - v1(2)*v2(1)
 nv = nv/norm(nv)
end subroutine normal_vector

function dih(c1, c2, c3, c4)
 use Sdiag_parameter, only: PI
 implicit none
 real(kind=8) :: nv1(3), nv2(3), dih
 real(kind=8), intent(in) :: c1(3), c2(3), c3(3), c4(3)

 call normal_vector(c3-c2, c1-c2, nv1)
 call normal_vector(c3-c2, c4-c2, nv2)
 dih = DOT_PRODUCT(nv1, nv2)
 dih = MAX(-1d0, MIN(dih,1d0)) ! deal with numerical error
 dih = 180d0*DACOS(dih)/PI
 if(DOT_PRODUCT(c1-c2, nv2) < 0d0) dih = -dih
end function dih

