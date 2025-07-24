! written by jxzou at 20220214: move math library wrappers into this file

! This is a function to replace zdotc(3,z,1,z,1). zdotc is unsupported when
! using f2py+gfortran+OpenBLAS.
function cmplx_v3_square(z) result(s)
 implicit none
 real(kind=8) :: s
 complex(kind=8), intent(in) :: z(3)

 s = REAL(z(1))*REAL(z(1)) + AIMAG(z(1))*AIMAG(z(1)) + &
     REAL(z(2))*REAL(z(2)) + AIMAG(z(2))*AIMAG(z(2)) + &
     REAL(z(3))*REAL(z(3)) + AIMAG(z(3))*AIMAG(z(3))
end function cmplx_v3_square

! This is a function to replace REAL(zdotc(3,y,1,z,1)). zdotc is unsupported
! when using f2py+gfortran+OpenBLAS.
function cmplx_v3_dot_real(z1, z2) result(s)
 implicit none
 real(kind=8) :: s
 complex(kind=8), intent(in) :: z1(3), z2(3)

 s = REAL(z1(1))*REAL(z2(1)) + AIMAG(z1(1))*AIMAG(z2(1)) + &
     REAL(z1(2))*REAL(z2(2)) + AIMAG(z1(2))*AIMAG(z2(2)) + &
     REAL(z1(3))*REAL(z2(3)) + AIMAG(z1(3))*AIMAG(z2(3))
end function cmplx_v3_dot_real

! Rotate two MOs using cos(alpha) and sin(alpha).
! The ifort directive "!dir$ ivdep" cannot be recognized by gfortran. So here we
!  use OpenMP SIMD, which is universal for various Fortran compilers.
subroutine rotate_mo_ij(cos_a, sin_a, nbf, mo_i, mo_j)
 implicit none
 integer :: k
 integer, intent(in) :: nbf
!f2py intent(in) :: nbf
 real(kind=8), intent(in) :: cos_a, sin_a
!f2py intent(in) :: cos_a, sin_a
 real(kind=8), intent(inout) :: mo_i(nbf), mo_j(nbf)
!f2py intent(in,out) :: mo_i, mo_j
!f2py depend(nbf) :: mo_i, mo_j
 real(kind=8), allocatable :: tmp_mo(:)

 allocate(tmp_mo(nbf), source=mo_i)

!$omp simd
 do k = 1, nbf, 1
  mo_i(k) = cos_a*tmp_mo(k) + sin_a*mo_j(k)
 end do ! for k
!$omp end simd

!$omp simd
 do k = 1, nbf, 1
  mo_j(k) = cos_a*mo_j(k) - sin_a*tmp_mo(k)
 end do ! for k
!$omp end simd

 deallocate(tmp_mo)
end subroutine rotate_mo_ij

! update gross(:,i,i), gross(:,j,j) and gross(:,j,i) involved in Jacobian 2*2
! rotation
subroutine update_gross_ii_jj_ji(cos_a, sin_a, natom, v, vdiff, gross_ii, &
                                 gross_jj, gross_ji)
 implicit none
 integer :: k
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 real(kind=8) :: cc, ss, cos_2a, sin_2a
 real(kind=8), intent(in) :: cos_a, sin_a, v(natom,3), vdiff(natom)
!f2py intent(in) :: cos_a, sin_a, v, vdiff
!f2py depend(natom) :: v, vdiff
 real(kind=8), intent(inout) :: gross_ii(natom),gross_jj(natom),gross_ji(natom)
!f2py intent(in,out) :: gross_ii, gross_jj, gross_ji
!f2py depend(natom) :: gross_ii, gross_jj, gross_ji

 cc = cos_a*cos_a; ss = sin_a*sin_a
 cos_2a = cc - ss; sin_2a = 2d0*sin_a*cos_a

!$omp simd
 do k = 1, natom, 1
  gross_ii(k) = cc*v(k,1) + ss*v(k,3) + sin_2a*v(k,2)
 end do ! for k
!$omp end simd

!$omp simd
 do k = 1, natom, 1
  gross_jj(k) = ss*v(k,1) + cc*v(k,3) - sin_2a*v(k,2)
 end do ! for k
!$omp end simd

!$omp simd
 do k = 1, natom, 1
  gross_ji(k) = cos_2a*v(k,2) - 0.5d0*sin_2a*vdiff(k)
 end do ! for k
!$omp end simd
end subroutine update_gross_ii_jj_ji

subroutine init_identity_mat(n, a)
 implicit none
 integer :: i
 integer, intent(in) :: n
 real(kind=8), intent(out) :: a(n,n)

 if(n == 1) then
  a(1,1) = 1d0
  return
 end if
 a = 0d0

!$omp simd
 do i = 1, n, 1
  a(i,i) = 1d0
 end do ! for n
!$omp end simd
end subroutine init_identity_mat

! initialize an RHF occupation matrix
subroutine init_rhf_occ_mat(ndb, nif, n)
 implicit none
 integer :: i
 integer, intent(in) :: ndb, nif
 real(kind=8), intent(out) :: n(nif,nif)

 n = 0d0
!$omp simd
 do i = 1, ndb, 1
  n(i,i) = 2d0
 end do ! for i
!$omp end simd
end subroutine init_rhf_occ_mat

! Reduce a fraction, where
! the product of nume_deno(1,:) is the numerator,
! the product of nume_deno(2,:) is the denominator.
! For example, (7*6*5)/(3*2) = 7*5/1
subroutine reduce_frac(n, nume_deno, deno_one)
 implicit none
 integer :: i, j, k1, k2
 integer, intent(in) :: n
 integer, intent(inout) :: nume_deno(2,n)
 logical :: changed
 logical, intent(in) :: deno_one

 do while(.true.)
  changed = .false.

  do i = 1, n, 1
   k1 = nume_deno(1,i)
   if(k1 == 1) cycle
   do j = 1, n-1, 1
    k2 = nume_deno(2,j)
    if(k2 == 1) cycle
    if(MOD(k1,k2) == 0) then
     k1 = k1/k2; nume_deno(1,i) = k1; nume_deno(2,j) = 1
     changed = .true.
    end if
   end do ! for j
  end do ! for i

  do i = 1, n-1, 1
   k1 = nume_deno(2,i)
   if(k1 == 1) cycle
   do j = 1, n, 1
    k2 = nume_deno(1,j)
    if(k2 == 1) cycle
    if(MOD(k1,k2) == 0) then
     k1 = k1/k2; nume_deno(2,i) = k1; nume_deno(1,j) = 1
     changed = .true.
    end if
   end do ! for j
  end do ! for i

  if(.not. changed) exit
 end do ! for while

 if(deno_one) then
  if(ANY( nume_deno(2,:)/=1 )) then
   write(6,'(/,A)') 'ERROR in subroutine reduce_frac: some integer in denominat&
                    &or is not reduced to 1.'
   write(6,'(A,I0)') 'n=', n
   write(6,'(20I3)') nume_deno(1,:)
   write(6,'(20I3)') nume_deno(2,:)
   stop
  end if
 end if
end subroutine reduce_frac

! simplify/reduce a combinatorial number C_{n}^{p}, e.g.
! C_{6}^{3} = (6*5*4)/(3*2*1) = 5*4
! C_{7}^{4} = C_{7}^{3} = (7*6*5)/(3*2*1) = 7*5
subroutine simplify_comb_cnp(n, p, p_min, int_prod)
 implicit none
 integer :: i
 integer, intent(in) :: n, p, p_min
 ! p_min = min(p, n-p)
 integer, intent(out) :: int_prod(p_min)
 integer, allocatable :: nume_deno(:,:)

 if(p==0 .or. p==n) then
  int_prod(1) = 1
  return
 else if(p==1 .or. p==n-1) then
  int_prod(1) = n
  return
 else if(p > n) then
  write(6,'(/,A)') 'ERROR in subroutine reduce_comb_cnp: p>n. Input values are &
                   &nonsense.'
  write(6,'(A,2I8)') 'n, p=', n, p
  stop
 end if

 allocate(nume_deno(2,p_min))

 do i = 1, p_min, 1
  nume_deno(1,i) = n - i + 1
  nume_deno(2,i) = p_min - i + 1
 end do ! for i

 call reduce_frac(p_min, nume_deno, .true.)
 int_prod = nume_deno(1,:)
 deallocate(nume_deno)
end subroutine simplify_comb_cnp

! combine two integer arrays with elements<2 discarded
subroutine combine_int_prod(n1, int_prod1, n2, int_prod2, n3, int_prod3)
 implicit none
 integer :: i, j, k
 integer, intent(in) :: n1, n2, n3
 integer, intent(in) :: int_prod1(n1), int_prod2(n2)
 integer, intent(out) :: int_prod3(n3)

 k = 0

 do i = 1, n1, 1
  j = int_prod1(i)
  if(j > 1) then
   k = k + 1
   int_prod3(k) = j
  end if
 end do ! for i

 do i = 1, n2, 1
  j = int_prod2(i)
  if(j > 1) then
   k = k + 1
   int_prod3(k) = j
  end if
 end do ! for i
end subroutine combine_int_prod

! check whether (C_{n}^{p1})*(C_{n}^{p2}) > (C_{m}^{q1})*(C_{m}^{q2})
function compare_cnp_cmq_prod(n, p1, p2, m, q1, q2) result(larger)
 implicit none
 integer :: i, k, k1, k2, p1_min, p2_min, q1_min, q2_min
 integer, intent(in) :: n, p1, p2, m, q1, q2
!f2py intent(in) :: n, p1, p2, m, q1, q2
 integer, allocatable :: int_prod1(:),int_prod2(:), int_prod3(:),int_prod4(:), &
  nume_deno(:,:)
 real(kind=8) :: res
 logical :: larger

 larger = .false.

 p1_min = MAX(1, MIN(p1, n-p1))
 p2_min = MAX(1, MIN(p2, n-p2))
 allocate(int_prod1(p1_min), int_prod2(p2_min))
 call simplify_comb_cnp(n, p1, p1_min, int_prod1)
 call simplify_comb_cnp(n, p2, p2_min, int_prod2)
 k1 = COUNT(int_prod1>1) + COUNT(int_prod2>1)
 allocate(int_prod3(k1))
 call combine_int_prod(p1_min,int_prod1, p2_min,int_prod2, k1,int_prod3)
 deallocate(int_prod1, int_prod2)

 q1_min = MAX(1, MIN(q1, m-q1))
 q2_min = MAX(1, MIN(q2, m-q2))
 allocate(int_prod1(q1_min), int_prod2(q2_min))
 call simplify_comb_cnp(m, q1, q1_min, int_prod1)
 call simplify_comb_cnp(m, q2, q2_min, int_prod2)
 k2 = COUNT(int_prod1>1) + COUNT(int_prod2>1)
 allocate(int_prod4(k2))
 call combine_int_prod(q1_min,int_prod1, q2_min,int_prod2, k2,int_prod4)
 deallocate(int_prod1, int_prod2)

 k = MAX(k1, k2)
 allocate(nume_deno(2,k), source=1)
 nume_deno(1,1:k1) = int_prod3
 nume_deno(2,1:k2) = int_prod4
 deallocate(int_prod3, int_prod4)

 call reduce_frac(k, nume_deno, .false.)
 ! PRODUCT(nume_deno(1,:)) is dangerous since it may exceed the range of
 ! integer(kind=4). So double precision division one by one is used below.
 res = 1d0
 do i = 1, k, 1
  k1 = nume_deno(1,i); k2 = nume_deno(2,i)
  if(k1 == 1) then
   if(k2 /= 1) res = res/DBLE(k2)
  else
   if(k2 == 1) then
    res = res*DBLE(k1)
   else
    res = res*DBLE(k1)/DBLE(k2)
   end if
  end if
 end do ! for i

 deallocate(nume_deno)
 if(res-1d0 > 1d-3) larger = .true.
end function compare_cnp_cmq_prod

! compare which active space size is larger by active orbitals, electrons and
! spin multiplicity
function compare_as_size(nacto1,nacte1,mult1, nacto2,nacte2,mult2) result(larger)
 implicit none
 integer :: nopen1, nacta1, nactb1, nopen2, nacta2, nactb2
 integer, intent(in) :: nacto1,nacte1,mult1, nacto2,nacte2,mult2
!f2py intent(in) :: nacto1,nacte1,mult1, nacto2,nacte2,mult2
 logical :: larger
 logical, external :: compare_cnp_cmq_prod

 if(nacto1==nacte1 .and. nacto2==nacte2 .and. mult1==mult2 .and. &
    nacto1>nacto2) then
  larger = .true.
  return
 end if

 nopen1 = mult1 - 1
 nactb1 = (nacte1 - nopen1)/2
 nacta1 = nacte1 - nactb1

 nopen2 = mult2 - 1
 nactb2 = (nacte2 - nopen2)/2
 nacta2 = nacte2 - nactb2

 larger = compare_cnp_cmq_prod(nacto1,nacta1,nactb1, nacto2,nacta2,nactb2)
end function compare_as_size

! sort an integer array by ascending/descending order
subroutine sort_int_array(n, a, ascending, idx)
 implicit none
 integer :: i, j, k, m
 integer, intent(in) :: n
 integer, intent(inout) :: a(n)
 integer, intent(out) :: idx(n)
 logical, intent(in) :: ascending

 forall(i = 1:n) idx(i) = i
 if(n == 1) return

 if(ascending) then
  do i = 1, n-1, 1
   k = a(i)
   do j = i+1, n, 1
    if(k > a(j)) then
     a(i) = a(j); m = a(j); a(j) = k; k = m
     m = idx(i); idx(i) = idx(j); idx(j) = m
    end if
   end do ! for j
  end do ! for i

 else ! descending order
  do i = 1, n-1, 1
   k = a(i)
   do j = i+1, n, 1
    if(k < a(j)) then
     a(i) = a(j); m = a(j); a(j) = k; k = m
     m = idx(i); idx(i) = idx(j); idx(j) = m
    end if
   end do ! for j
  end do ! for i
 end if
end subroutine sort_int_array

! sort a double precision array
subroutine sort_dp_array(n, a, ascending, idx)
 implicit none
 integer :: i, j, k
 integer, intent(in) :: n
 real(kind=8) :: r
 real(kind=8), intent(inout) :: a(n)
 integer, intent(out) :: idx(n)
 logical, intent(in) :: ascending

 forall(i = 1:n) idx(i) = i
 if(n == 1) return

 if(ascending) then
  do i = 1, n-1, 1
   r = a(i)
   do j = i+1, n, 1
    if(r > a(j)) then
     r = a(j); a(j) = a(i); a(i) = r
     k = idx(i); idx(i) = idx(j); idx(j) = k
    end if
   end do ! for j
  end do ! for i
 else ! descending order
  do i = 1, n-1, 1
   r = a(i)
   do j = i+1, n, 1
    if(r < a(j)) then
     r = a(j); a(j) = a(i); a(i) = r
     k = idx(i); idx(i) = idx(j); idx(j) = k
    end if
   end do ! for j
  end do ! for i
 end if
end subroutine sort_dp_array

! sort a set of MOs by the given eigenvalues (i.e. orbital energies or
! occupation numbers)
subroutine sort_mo_by_ev(nbf, nmo, mo, ev, new_mo, new_ev)
 implicit none
 integer :: i, j
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 real(kind=8) :: tmp_ev
 real(kind=8), intent(in) :: mo(nbf,nmo), ev(nmo)
!f2py intent(in) :: mo, ev
!f2py depend(nbf,nmo) :: mo
!f2py depend(nmo) :: ev
 real(kind=8), intent(out) :: new_mo(nbf,nmo), new_ev(nmo)
!f2py intent(out) :: new_mo, new_ev
!f2py depend(nbf,nmo) :: new_mo
!f2py depend(nmo) :: new_ev
 real(kind=8), allocatable :: tmp_mo(:)

 new_mo = mo; new_ev = ev
 allocate(tmp_mo(nbf))

 do i = 1, nmo-1, 1
  tmp_ev = new_ev(i)
  do j = i+1, nmo, 1
   if(new_ev(j) < tmp_ev) then
    new_ev(i) = new_ev(j); new_ev(j) = tmp_ev; tmp_ev = new_ev(i)
    tmp_mo = new_mo(:,i)
    new_mo(:,i) = new_mo(:,j)
    new_mo(:,j) = tmp_mo
   end if
  end do ! for j
 end do ! for i

 deallocate(tmp_mo)
end subroutine sort_mo_by_ev

! check whether a (double) complex matrix is hermitian
subroutine check_hermitian(n, a)
 implicit none
 integer, intent(in) :: n
!f2py intent(in) :: n
 complex(kind=8), intent(in) :: a(n,n)
!f2py intent(in) :: a
!f2py depend(n) :: a
 real(kind=8), allocatable :: b_real(:,:), b_imag(:,:)
 complex(kind=8), allocatable :: b(:,:)

 allocate(b(n,n))
 b = CONJG(TRANSPOSE(a)) - a
 allocate(b_real(n,n), source=DABS(REAL(b)))
 allocate(b_imag(n,n), source=DABS(AIMAG(b)))
 deallocate(b)
 write(6,'(A,F20.8)') 'Ave abs REAL:', SUM(b_real)/DBLE(n*n)
 write(6,'(A,F20.8)') 'Ave abs IMAG:', SUM(b_imag)/DBLE(n*n)
 deallocate(b_real, b_imag)
end subroutine check_hermitian

! check whether a (double) complex matrix is symmetric
subroutine check_symm_cmplx(n, a)
 implicit none
 integer, intent(in) :: n
!f2py intent(in) :: n
 complex(kind=8), intent(in) :: a(n,n)
!f2py intent(in) :: a
!f2py depend(n) :: a
 real(kind=8), allocatable :: b_real(:,:), b_imag(:,:)
 complex(kind=8), allocatable :: b(:,:)

 allocate(b(n,n))
 b = TRANSPOSE(a) - a
 allocate(b_real(n,n), source=DABS(REAL(b)))
 allocate(b_imag(n,n), source=DABS(AIMAG(b)))
 deallocate(b)
 write(6,'(A,F20.8)') 'Ave abs REAL:', SUM(b_real)/DBLE(n*n)
 write(6,'(A,F20.8)') 'Ave abs IMAG:', SUM(b_imag)/DBLE(n*n)
 deallocate(b_real, b_imag)
end subroutine check_symm_cmplx

! get upper triangle index pairs (j>=i), similar to numpy.triu_indices
subroutine get_triu_idx(n, map)
 implicit none
 integer :: i, j, k
 integer, intent(in) :: n
!f2py intent(in) :: n
 integer, intent(out) :: map(2,n*(n+1)/2)
!f2py intent(out) :: map
!f2py depend(n) :: map

 if(n < 1) then
  write(6,'(/,A)') 'ERROR in subroutine get_triu_idx: n<1 not allowed.'
  write(6,'(A,I0)') 'n=', n
  stop
 end if

 if(n < 100) then
  forall(i=1:n, j=1:n, j>=i) map(:,(2*n-i)*(i-1)/2+j) = [i,j]
 else
!$omp parallel do schedule(dynamic) default(shared) private(i,j,k)
  do i = 1, n, 1
   k = (2*n-i)*(i-1)/2
   do j = i, n, 1
    map(:,k+j) = [i,j]
   end do ! for j
  end do ! for i
!$omp end parallel do
 end if
end subroutine get_triu_idx

! get upper triangle index pairs (j>i), similar to numpy.triu_indices
subroutine get_triu_idx1(n, map)
 implicit none
 integer :: i, j, k
 integer, intent(in) :: n
!f2py intent(in) :: n
 integer, intent(out) :: map(2,n*(n-1)/2)
!f2py intent(out) :: map
!f2py depend(n) :: map

 if(n < 2) then
  write(6,'(/,A)') 'ERROR in subroutine get_triu_idx1: n<2 not allowed.'
  write(6,'(A,I0)') 'n=', n
  stop
 end if

 if(n < 100) then
  forall(i=1:n-1, j=1:n, j>i) map(:,(2*n-i)*(i-1)/2+j-i) = [i,j]
 else
!$omp parallel do schedule(dynamic) default(shared) private(i,j,k)
  do i = 1, n-1, 1
   k = (2*n-i)*(i-1)/2 - i
   do j = i+1, n, 1
    map(:,k+j) = [i,j]
   end do ! for j
  end do ! for i
!$omp end parallel do
 end if
end subroutine get_triu_idx1

!! generate the round robin ordering (DOI: 10.1109/EMPDP.1995.389182)
!subroutine init_round_robin_idx(n, np, map)
! implicit none
! integer :: i, j, k
! integer, intent(in) :: n, np
! integer, intent(out) :: map(2,np,2*np-1)
!
! forall(i = 1:np) map(:,i,1) = [2*i, 2*i-1]
! if(MOD(n,2) == 1) then
!  map(1,np,:) = 0
! else
!  map(1,np,:) = n
! end if
! k = 2*np - 1
!
! do i = 2, k, 1
!  do j = 1, np-2, 1
!   map(1,j,i) = map(1,j+1,i-1)
!  end do ! for j
!
!  map(1,np-1,i) = map(2,np,i-1)
!  map(2,1,i) = map(1,1,i-1)
!
!  do j = 2, np, 1
!   map(2,j,i) = map(2,j-1,i-1)
!  end do ! for j
! end do ! for i
!end subroutine init_round_robin_idx

! Generate the ring Jacobi ordering (DOI: 10.1109/EMPDP.1995.389182). This
! map array leads to more effective Jacobi rotations according to jxzou's test.
subroutine init_ring_jacobi_idx(n, np, map)
 implicit none
 integer :: i, j, k, m, dnp
 integer, intent(in) :: n, np
 integer, intent(out) :: map(2,np,4*np-2)

 dnp = 2*np

!$omp parallel do schedule(static) default(shared) private(i,j,k)
 do i = 1, np, 1
  j = 2*i
  k = dnp - j + 1
  map(:,i,1) = [k, k+1]
  map(:,i,dnp) = [j, j-1]
 end do ! for i
!$omp end parallel do

 if(MOD(n,2) == 1) then
  map(2,1,1) = 0
  map(1,np,dnp) = 0
 end if
 k = 2*np - 1

!$omp parallel sections private(i,m)
!$omp section
 do i = 1, k-1, 1
  map(:,:,i+1) = map(:,:,i)
  m = (i+1)/2
  map(:,m,i+1) = [map(2,m,i), map(1,m,i)]
  map(2,:,i+1) = CSHIFT(map(2,:,i+1), -1)
 end do ! for i
!$omp section
 do i = dnp, 2*k-1, 1
  map(:,:,i+1) = map(:,:,i)
  m = np - (i-dnp)/2
  map(:,m,i+1) = [map(2,m,i), map(1,m,i)]
  map(1,:,i+1) = CSHIFT(map(1,:,i+1), 1)
 end do ! for i
!$omp end parallel sections
end subroutine init_ring_jacobi_idx

! Diagonalize a real symmetric matrix and get all eigenvalues and eigenvectors.
! A = Ua(U^T). Eigenvectors U will be stored in the square matrix a, and eigenvalues
! in w() are in ascending order w(1)<=w(2)<=...
subroutine diag_get_e_and_vec(n, a, w)
 implicit none
 integer :: i, lwork, liwork
 integer, intent(in) :: n
 integer, allocatable :: iwork(:)
 real(kind=8), intent(inout) :: a(n,n)
 real(kind=8), intent(out) :: w(n)
 real(kind=8), allocatable :: work(:)

 if(n == 1) then
  w(1) = a(1,1)
  a(1,1) = 1d0
  return
 end if

 w = 0d0; lwork = -1; liwork = -1
 allocate(work(1), iwork(1))
 call dsyevd('V', 'U', n, a, n, w, work, lwork, iwork, liwork, i)
 lwork = CEILING(work(1))
 liwork = iwork(1)
 deallocate(work, iwork)
 if(lwork<1 .or. liwork<1) then
  write(6,'(/,A)') 'ERROR in subroutine diag_get_e_and_vec: lwork or liwork is &
                   &less than 1.'
  write(6,'(2(A,I0))') 'lwork=', lwork, ', liwork=', liwork
  stop
 end if
 allocate(work(lwork), iwork(liwork))
 call dsyevd('V', 'U', n, a, n, w, work, lwork, iwork, liwork, i)

 deallocate(work, iwork)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine diag_get_e_and_vec: info/=0 in dsyevd.'
  write(6,'(2(A,I0))') 'n=', n, ', info=', i
  stop
 end if
end subroutine diag_get_e_and_vec

! It seems that dsyevr is slower than dsyevd, according to my experience in
! 15315*15315. Eigenvectors U will be stored in the square matrix a, and
! eigenvalues in w() are in ascending order w(1)<=w(2)<=...
subroutine diag_get_e_and_vec2(n, a, w)
 implicit none
 integer :: i, m, lwork, liwork
 integer, intent(in) :: n
 integer, allocatable :: iwork(:), isuppz(:)
 real(kind=8), intent(inout) :: a(n,n)
 real(kind=8), intent(out) :: w(n)
 real(kind=8), allocatable :: U(:,:), work(:)

 if(n == 1) then
  w(1) = a(1,1)
  a(1,1) = 1d0
  return
 end if

 w = 0d0; lwork = -1; liwork = -1
 allocate(U(n,n), isuppz(2*n), work(1), iwork(1))
 call dsyevr('V', 'A', 'L', n, a, n, 0d0, 0d0, 0, 0, 1d-8, m, w, U, n, &
             isuppz, work, lwork, iwork, liwork, i)
 lwork = CEILING(work(1))
 liwork = iwork(1)
 deallocate(work, iwork)
 allocate(work(lwork), iwork(liwork))
 call dsyevr('V', 'A', 'L', n, a, n, 0d0, 0d0, 0, 0, 1d-8, m, w, U, n, &
             isuppz, work, lwork, iwork, liwork, i)

 deallocate(isuppz, work, iwork)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine diag_get_e_and_vec2: info/=0 in dsyevr.'
  write(6,'(2(A,I0))') 'n=', n, ', info=', i
  stop
 end if

 a = U
 deallocate(U)
end subroutine diag_get_e_and_vec2

! reverse eigenvalues and eigenvectors which were obtained from a previous
! matrix diagonalization
subroutine reverse_e_and_vec(n, u, w)
 implicit none
 integer :: i
 integer, intent(in) :: n
 real(kind=8), intent(inout) :: u(n,n), w(n)
 real(kind=8), allocatable :: y(:,:), z(:)

 allocate(z(n), source=w)
!$omp parallel do schedule(dynamic) default(shared) private(i)
 do i = 1, n, 1
  w(i) = z(n-i+1)
 end do ! for i
!$omp end parallel do
 deallocate(z)

 allocate(y(n,n), source=u)
!$omp parallel do schedule(dynamic) default(shared) private(i)
 do i = 1, n, 1
  u(:,i) = y(:,n-i+1)
 end do ! for i
!$omp end parallel do
 deallocate(y)
end subroutine reverse_e_and_vec

! find the number of independent MOs from the AO overlap integral matrix
subroutine get_nmo_from_ao_ovlp(nbf, ovlp, nmo)
 implicit none
 integer :: i
 integer, intent(in) :: nbf
!f2py intent(in) :: nbf
 integer, intent(out) :: nmo
!f2py intent(out) :: nmo
 real(kind=8), intent(in) :: ovlp(nbf,nbf)
!f2py intent(in) :: ovlp
!f2py depend(nbf) :: ovlp
 real(kind=8), parameter :: thres = 1d-6
 real(kind=8), allocatable :: S(:,:), w(:)

 allocate(S(nbf,nbf), source=ovlp)
 allocate(w(nbf))
 call diag_get_e_and_vec(nbf, S, w) ! w(1)<=w(2)<=...
 deallocate(S)

 nmo = 0
 do i = 1, nbf, 1
  if(w(i) > thres) exit
  nmo = nmo + 1
 end do ! for i
 deallocate(w)

 nmo = nbf - nmo
end subroutine get_nmo_from_ao_ovlp

! calculate Us(U^T), s is an array and will be enlarged into a diagonal matrix
subroutine calc_usut(n, s, U, usut)
 implicit none
 integer :: i
 integer, intent(in) :: n
 real(kind=8), intent(in) :: s(n), U(n,n)
 real(kind=8), intent(out) :: usut(n,n)
 real(kind=8), allocatable :: e(:,:), Ue(:,:)

 allocate(e(n,n), source=0d0)
 forall(i = 1:n) e(i,i) = s(i)
 allocate(Ue(n,n), source=0d0)
 call dsymm('R', 'L', n, n, 1d0, e, n, U, n, 0d0, Ue, n)
 deallocate(e)
 usut = 0d0
 call dgemm('N', 'T', n, n, n, 1d0, Ue, n, U, n, 0d0, usut, n)
 deallocate(Ue)
end subroutine calc_usut

! calculate only diagonal elements of Us(U^T)
subroutine calc_usut_diag_elem(n, s, u, d)
 implicit none
 integer :: i, j, k, m, n2
 integer, intent(in) :: n
 real(kind=8), intent(in) :: s(n), u(n,n)
 real(kind=8), intent(out) :: d(n)
 real(kind=8), allocatable :: u2(:,:)

 allocate(u2(n,n))
 n2 = n*n
 ! using one compound loop is faster than double loops (j,i)

 !$omp parallel do schedule(static) default(private) shared(n,n2,u,u2)
 do k = 1, n2, 1
  m = k/n
  if(k == n*m) then
   i = m
  else
   i = m + 1
  end if
  j = k - (i-1)*n
  u2(j,i) = u(j,i)*u(j,i)
 end do ! for k
 !$omp end parallel do

 d = 0d0
 call dgemv('N', n, n, 1d0, u2, n, s, 1, 0d0, d, 1)
 deallocate(u2)
end subroutine calc_usut_diag_elem

! calculate the difference matrix of two sets of coordinates
subroutine calc_coor_diff_mat(n, coor1, coor2, mat)
 implicit none
 integer :: i, j, k, m, n2
 integer, intent(in) :: n
!f2py intent(in) :: n
 real(kind=8) :: v(3)
 real(kind=8), intent(in) :: coor1(3,n), coor2(3,n)
!f2py intent(in) :: coor1, coor2
!f2py depend(n) :: coor1, coor2
 real(kind=8), intent(out) :: mat(n,n)
!f2py intent(out) :: mat
!f2py depend(n) :: mat

 n2 = n*n

 !$omp parallel do schedule(static) default(shared) private(i,j,k,m,v)
 do k = 1, n2, 1
  m = k/n
  if(k == n*m) then
   i = m
  else
   i = m + 1
  end if
  j = k - (i-1)*n
  v = coor1(:,j) - coor2(:,i)
  mat(j,i) = DSQRT(DOT_PRODUCT(v,v))
 end do ! for k
 !$omp end parallel do
end subroutine calc_coor_diff_mat

! solve the A^1/2 and A^(-1/2) for a real symmetric matrix A
! Note: the input matrix A must be symmetric
subroutine mat_dsqrt(n, a, calc_n_sqrt_a, sqrt_a, n_sqrt_a)
 implicit none
 integer :: i
 integer, intent(in) :: n
!f2py intent(in) :: n
 real(kind=8), parameter :: lin_dep = 1d-6
 ! 1D-6 is the default threshold of linear dependence in Gaussian and GAMESS
 ! But in PySCF and ORCA, one needs to manually adjust the threshold if linear
 ! dependence occurs.
 real(kind=8), intent(in) :: a(n,n)
!f2py intent(in) :: a
!f2py depend(n) :: a
 real(kind=8), intent(out) :: sqrt_a(n,n), n_sqrt_a(n,n) ! A^1/2, A(-1/2)
!f2py intent(out) :: sqrt_a, n_sqrt_a
!f2py depend(n) :: sqrt_a, n_sqrt_a
 real(kind=8), allocatable :: U(:,:), e(:)
 ! U: copy of a; eigenvectors from diagonalizing a
 ! e: eigenvalues
 logical, intent(in) :: calc_n_sqrt_a ! whether to calculate A^(-1/2)

 allocate(U(n,n), source=a)
 allocate(e(n))
 call diag_get_e_and_vec(n, U, e) ! e(1)<=e(2)<=...

 if(e(1) < -lin_dep) then
  deallocate(U, e)
  write(6,'(/,A)') 'ERROR in subroutine mat_dsqrt: too negative eigenvalue.'
  write(6,'(A,F16.9)') 'e(1)=', e(1)
  stop
 end if

 if(calc_n_sqrt_a .and. ANY(e<lin_dep)) then
  write(6,'(/,A)') 'ERROR in subroutine mat_dsqrt: linear dependency detected i&
                   &n matrix A.'
  write(6,'(A)') 'A^(-1/2) cannot be calculated.'
  write(6,'(A,I0)') 'n=', n
  write(6,'(A)') 'e='
  write(6,'(5(1X,ES15.8))') e
  stop
 end if

 ! Some tiny negative values are allowed, and they will be set to zero.
 do i = 1, n, 1
  if(e(i)>-lin_dep .and. e(i)<0d0) then
   e(i) = 0d0
  else
   e(i) = DSQRT(e(i))
  end if
 end do ! for i

 call calc_usut(n, e, U, sqrt_a)
 if(calc_n_sqrt_a) then
  e = 1d0/e
  call calc_usut(n, e, U, n_sqrt_a)
 end if

 deallocate(U, e)
end subroutine mat_dsqrt

! find the absolute value of the determinant of a square matrix by LU decomposition
! A=PLU, where
! P is a permutation matrix whose determinant is 1 or -1
! L is lower triangular with unit diagonal elements, so its determinant is 1
! U is upper triangular, so calculate the product of all diagonal elements
function abs_det(n, a) result(res)
 implicit none
 integer :: i
 integer, intent(in) :: n
 integer, allocatable :: ipiv(:)
 real(kind=8) :: res
 real(kind=8), intent(in) :: a(n,n)
 real(kind=8), allocatable :: a_copy(:,:), diag(:)

 res = 0d0
 allocate(ipiv(n), a_copy(n,n))
 a_copy = a
 call dgetrf(n, n, a_copy, n, ipiv, i) ! L/U is stored in a_copy
 deallocate(ipiv)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in function det: info/=0 in subroutine dgetrf.'
  write(6,'(A,I0)') 'info=', i
  stop
 end if

 allocate(diag(n))
 forall(i = 1:n) diag(i) = a_copy(i,i)
 deallocate(a_copy)

 res = PRODUCT(diag)
 deallocate(diag)
 if(res < 0d0) res = -res
end function abs_det

! find the inverse matrix of a square matrix A by solving systems of linear
! equations (the user should make sure that A is reversible)
subroutine inverse(n, a, inv_a)
 implicit none
 integer :: i
 integer, intent(in) :: n
 real(kind=8), intent(in) :: a(n,n)
 real(kind=8), intent(out) :: inv_a(n,n)
 real(kind=8), allocatable :: identity(:,:)

 allocate(identity(n,n), source=0d0)
 forall(i = 1:n) identity(i,i) = 1d0

 call solve_multi_lin_eqs(n, n, a, n, identity, inv_a)
 deallocate(identity)
end subroutine inverse

! find the inverse matrix of a square matrix A by diagonalization (the user
! should make sure that A is reversible)
subroutine inverse2(n, a, inv_a)
 implicit none
 integer :: i
 integer, intent(in) :: n
 real(kind=8), parameter :: thres = 1d-9
 real(kind=8), intent(in) :: a(n,n)
 real(kind=8), intent(out) :: inv_a(n,n)
 real(kind=8), allocatable :: w(:), u(:,:), ev(:,:), u_ev(:,:)

 inv_a = 0d0
 allocate(w(n))
 allocate(u(n,n), source=a)
 call diag_get_e_and_vec(n, u, w)

 if(ANY(DABS(w) < thres)) then
  write(6,'(/,A)') 'ERROR in subroutine inverse2: some eigenvalues are very clo&
                   &se to zero.'
  write(6,'(A)') 'Failed to calculate the inverse.'
  write(6,'(A,I0)') 'n=', n
  write(6,'(A)') 'w='
  write(6,'(5(1X,ES15.8))') w
  stop
 end if

 allocate(ev(n,n), source=0d0)
 forall(i = 1:n) ev(i,i) = 1d0/w(i)
 deallocate(w)
 allocate(u_ev(n,n), source=0d0)
 call dsymm('R', 'L', n, n, 1d0, ev, n, u, n, 0d0, u_ev, n)
 deallocate(ev)
 call dgemm('N', 'T', n, n, n, 1d0, u_ev, n, u, n, 0d0, inv_a, n)
 deallocate(u_ev, u)
end subroutine inverse2

subroutine newton_inv(n, a, inv_a)
 implicit none
 integer :: i
 integer, intent(in) :: n
 integer, parameter :: max_it = 999
 real(kind=8), intent(in) :: a(n,n)
 real(kind=8), intent(out) :: inv_a(n,n)
 real(kind=8), allocatable :: old_inv(:,:)
 real(kind=8), parameter :: thres = 1d-6

 allocate(old_inv(n,n), source=0d0)
 forall(i = 1:n) old_inv(i,i) = 1d0/a(i,i)

 do i = 1, max_it, 1
  call calc_sps(n, a, old_inv, inv_a)
  inv_a = 2d0*old_inv - inv_a
  if(SUM(DABS(old_inv - inv_a))/DBLE(n*n) < thres) exit
  old_inv = inv_a
 end do ! for i

 deallocate(old_inv)
 if(i-1 == max_it) then
  write(6,'(/,A)') 'ERROR in subroutine newton_inv: failed to converge.'
  stop
 else
  write(6,'(A,I3)') 'n_iter=', i
 end if
end subroutine newton_inv

! Solve systems of linear equations Ax=b with only 1 right-hand side, i.e.
! x and b are both vectors, not matrices.
subroutine solve_lin_eqs(a1, a2, a, b, x)
 implicit none
 integer, intent(in) :: a1, a2
 real(kind=8), intent(in) :: a(a1,a2), b(a1)
 real(kind=8), intent(out) :: x(a2)
 real(kind=8), allocatable :: bp(:,:), xp(:,:)

 allocate(bp(a1,1), xp(a2,1))
 bp(:,1) = b
 call solve_multi_lin_eqs(a1, a2, a, 1, bp, xp)
 x = xp(:,1)
 deallocate(xp, bp)
end subroutine solve_lin_eqs

! solving systems of linear equations with multiple right-hand sides
! Ax = b, x can be with multiple right-hand sides
subroutine solve_multi_lin_eqs(a1, a2, a, a3, b, x)
 implicit none
 integer :: i
 integer, intent(in) :: a1, a2, a3
 integer, allocatable :: ipiv(:)
 real(kind=8), intent(in) :: a(a1,a2), b(a1,a3)
 real(kind=8), intent(out) :: x(a2,a3)
 real(kind=8), allocatable :: a_copy(:,:), b_copy(:,:)

 x = 0d0
 allocate(a_copy(a1,a2), source=a)
 allocate(b_copy(a1,a3), source=b)
 allocate(ipiv(min(a1,a2)), source=0)

 call dgetrf(a1, a2, a_copy, a1, ipiv, i)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine solve_multi_lin_eqs: MKL subroutine dge&
                   &trf failed.'
  write(6,'(A,I0)') 'info=', i
  stop
 end if

 call dgetrs('N', a2, a3, a_copy, a1, ipiv, b_copy, a1, i)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine solve_multi_lin_eqs: MKL subroutine dge&
                   &trs failed.'
  write(6,'(A,I0)') 'info=', i
  stop
 end if
 deallocate(ipiv, a_copy)

 x = b_copy(1:a2, 1:a3)
 deallocate(b_copy)
end subroutine solve_multi_lin_eqs

! Cayley transformation K = (I-U)(I+U)^(-1), U = (I-K)(I+K)^(-1)
! For example, solve K using systems of linear equations (I+U)K = I-U. This
!  subroutine can be used for either K->U or U->K.
subroutine cayley_trans(n, u, k)
 implicit none
 integer, intent(in) :: n
 real(kind=8), intent(in) :: u(n,n)
 real(kind=8), intent(out) :: k(n,n)
 real(kind=8), allocatable :: one(:,:)

 allocate(one(n,n))
 call init_identity_mat(n, one)
 call solve_multi_lin_eqs(n, n, one+u, n, one-u, k)
 deallocate(one)
end subroutine cayley_trans

! perform SVD on a matrix
subroutine do_svd(m, n, a, u, vt, s)
 implicit none
 integer :: i, lwork
 integer, intent(in) :: m, n
 real(kind=8), intent(in) :: a(m,n)
 real(kind=8), intent(out) :: u(m,m), vt(n,n), s(m)
 real(kind=8), allocatable :: work(:), a_copy(:,:)

 u = 0d0; vt = 0d0; s = 0d0
 allocate(a_copy(m,n), source=a)

 lwork = -1
 allocate(work(1), source=0d0)
 call dgesvd('A', 'A', m, n, a_copy, m, s, u, m, vt, n, work, lwork, i)

 lwork = CEILING(work(1))
 deallocate(work)
 allocate(work(lwork),source=0d0)
 call dgesvd('A', 'A', m, n, a_copy, m, s, u, m, vt, n, work, lwork, i)
 deallocate(work, a_copy)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_svd: info/=0! Please check why.'
  write(6,'(A,I0)') 'info=', i
  stop
 end if
end subroutine do_svd

! perform SVD on the square matrix a, return U(V^T) and singular values
subroutine do_svd_get_uvt_s(n, a, uvt, s)
 implicit none
 integer, intent(in) :: n
 real(kind=8), intent(in) :: a(n,n)
 real(kind=8), intent(out) :: uvt(n,n), s(n)
 real(kind=8), allocatable :: u(:,:), vt(:,:)

 allocate(u(n,n), vt(n,n))
 call do_svd(n, n, a, u, vt, s)

 uvt = 0d0
 call dgemm('N', 'N', n, n, n, 1d0, u, n, vt, n, 0d0, uvt, n)
 deallocate(u, vt)
end subroutine do_svd_get_uvt_s

! calculate (C^T)SC, S is real symmetric, C is a vector
function calc_CiTSCi(nbf, C, S) result(res)
 implicit none
 integer, intent(in) :: nbf
 real(kind=8) :: ddot, res
 real(kind=8), intent(in) :: C(nbf), S(nbf,nbf)
 real(kind=8), allocatable :: SC(:)

 allocate(SC(nbf), source=0d0)
 call dsymv('U', nbf, 1d0, S, nbf, C, 1, 0d0, SC, 1)
 res = ddot(nbf, C, 1, SC, 1)
 deallocate(SC)
end function calc_CiTSCi

! symmetrize a double precision matrix
subroutine symmetrize_dmat(n, a)
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 real(kind=8), intent(inout) :: a(n,n)

 if(n == 1) return

!$omp parallel do schedule(dynamic) default(shared) private(i,j)
 do i = 1, n, 1
  do j = 1, i-1, 1
   a(i,j) = a(j,i)
  end do ! for j
 end do ! for i
!$omp end parallel do
end subroutine symmetrize_dmat

! symmetrize a double precision complex matrix
subroutine symmetrize_zmat(n, a)
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 complex(kind=8), intent(inout) :: a(n,n)

 if(n == 1) return

!$omp parallel do schedule(dynamic) default(shared) private(i,j)
 do i = 1, n, 1
  do j = 1, i-1, 1
   a(i,j) = a(j,i)
  end do ! for j
 end do ! for i
!$omp end parallel do
end subroutine symmetrize_zmat

! symmetrize the double precision complex matrices mo_zdip(i,:,:)
subroutine symmetrize_mo_zdip(nmo, mo_zdip)
 implicit none
 integer, intent(in) :: nmo
!f2py intent(in) :: nmo
 complex(kind=8), intent(inout) :: mo_zdip(3,nmo,nmo)
!f2py intent(in,out) :: mo_zdip
!f2py depend(nmo) :: mo_zdip

 call symmetrize_zmat(nmo, mo_zdip(1,:,:))
 call symmetrize_zmat(nmo, mo_zdip(2,:,:))
 call symmetrize_zmat(nmo, mo_zdip(3,:,:))
end subroutine symmetrize_mo_zdip

! calculate C(C^T)
subroutine calc_cct(nbf, nif, occ_mo, cct)
 implicit none
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: occ_mo(nbf,nif)
!f2py intent(in) :: occ_mo
!f2py depend(nbf,nif) :: occ_mo
 real(kind=8), intent(out) :: cct(nbf,nbf)
!f2py intent(out) :: cct
!f2py depend(nbf) :: cct

 cct = 0d0
 call dsyrk('U', 'N', nbf, nif, 1d0, occ_mo, nbf, 0d0, cct, nbf)
 call symmetrize_dmat(nbf, cct)
end subroutine calc_cct

! normalize an MO
subroutine normalize_mo(nbf, ao_ovlp, mo)
 implicit none
 integer, intent(in) :: nbf
 real(kind=8) :: a_2
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf)
 real(kind=8), intent(inout) :: mo(nbf)
 real(kind=8), external :: calc_CiTSCi

 a_2 = calc_CiTSCi(nbf, mo, ao_ovlp)
 mo = mo/DSQRT(a_2)
end subroutine normalize_mo

! normalize a set of MOs
subroutine normalize_mos(nbf, nif, ao_ovlp, mo)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nif
 real(kind=8) :: norm, ddot
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf)
 real(kind=8), intent(inout) :: mo(nbf,nif)
 real(kind=8), allocatable :: sc(:,:)

 allocate(sc(nbf,nif), source=0d0)
 call dsymm('L', 'L', nbf, nif, 1d0, ao_ovlp, nbf, mo, nbf, 0d0, sc, nbf)

 !$omp parallel do schedule(static) default(private) shared(nbf,nif,mo,sc)
 do i = 1, nif, 1
  norm = ddot(nbf, mo(:,i), 1, sc(:,i), 1)
  mo(:,i) = mo(:,i)/DSQRT(norm)
 end do ! for i
 !$omp end parallel do

 deallocate(sc)
end subroutine normalize_mos

! calculate/compute normalized PAOs
subroutine get_normalized_pao(nbf, nif, ao_ovlp, occ_mo, pao)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf), occ_mo(nbf,nif)
 real(kind=8), intent(out) :: pao(nbf,nbf)
 real(kind=8), allocatable :: cct(:,:)

 allocate(cct(nbf,nbf))
 call calc_cct(nbf, nif, occ_mo, cct)

 pao = 0d0
 forall(i = 1:nbf) pao(i,i) = 1d0
 call dsymm('L','L', nbf,nbf, -1d0,cct,nbf, ao_ovlp,nbf, 1d0,pao,nbf)
 deallocate(cct)

 call normalize_mos(nbf, nbf, ao_ovlp, pao) ! normalize each PAO
end subroutine get_normalized_pao

! calculate (A^T)BC
subroutine calc_atbc(n1, n2, n3, n4, a, b, c, atbc)
 implicit none
 integer, intent(in) :: n1, n2, n3, n4
!f2py intent(in) :: n1, n2, n3, n4
 real(kind=8), intent(in) :: a(n2,n1), b(n2,n3), c(n3,n4)
!f2py intent(in) :: a, b, c
!f2py depend(n1,n2) :: a
!f2py depend(n2,n3) :: b
!f2py depend(n3,n4) :: c
 real(kind=8), intent(out) :: atbc(n1,n4)
!f2py intent(out) :: atbc
!f2py depend(n1,n4) :: atbc
 real(kind=8), allocatable :: bc(:,:)

 allocate(bc(n2,n4), source=0d0)
 call dgemm('N', 'N', n2, n4, n3, 1d0, b, n2, c, n3, 0d0, bc, n2)
 atbc = 0d0
 call dgemm('T', 'N', n1, n4, n2, 1d0, a, n2, bc, n2, 0d0, atbc, n1)
 deallocate(bc)
end subroutine calc_atbc

! calculate diagonal elements of (C^T)SC, where S is real symmetric
subroutine calc_ctsc_diag(nbf, nif, s, c, diag)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8) :: ddot
 real(kind=8), intent(in) :: s(nbf,nbf), c(nbf,nif)
!f2py intent(in) :: s, c
!f2py depend(nbf) :: s
!f2py depend(nbf,nif) :: c
 real(kind=8), intent(out) :: diag(nif)
!f2py intent(out) :: diag
!f2py depend(nif) :: diag
 real(kind=8), allocatable :: sc(:,:)

 allocate(sc(nbf,nif), source=0d0)
 call dsymm('L', 'L', nbf, nif, 1d0, s, nbf, c, nbf, 0d0, sc, nbf)

 !$omp parallel do schedule(static) default(private) shared(nbf,nif,c,sc,diag)
 do i = 1, nif, 1
  diag(i) = ddot(nbf, c(:,i), 1, sc(:,i), 1)
 end do ! for i
 !$omp end parallel do

 deallocate(sc)
end subroutine calc_ctsc_diag

! calculate diagonal elements of each component of (C^T)DC, where D(i,:,:),
! i=1,2,3 is real symmetric
subroutine calc_ctdc_diag(nbf, nif, d, c, diag)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: d(3,nbf,nbf), c(nbf,nif)
!f2py intent(in) :: d, c
!f2py depend(nbf) :: d
!f2py depend(nbf,nif) :: c
 real(kind=8), intent(out) :: diag(3,nif)
!f2py intent(out) :: diag
!f2py depend(nif) :: diag
 real(kind=8), allocatable :: diag0(:), s(:,:)

 allocate(diag0(nbf), s(nbf,nbf))

 do i = 1, 3
  s = d(i,:,:)
  call calc_ctsc_diag(nbf, nif, s, c, diag0)
  diag(i,:) = diag0
 end do ! for i

 deallocate(diag0, s)
end subroutine calc_ctdc_diag

! calculate (C^T)SC, S must be real symmetric since dsymm is called
! C: nbf*nif, S: nbf*nbf
subroutine calc_CTSC(nbf, nif, C, S, CTSC)
 implicit none
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: C(nbf,nif), S(nbf,nbf)
!f2py intent(in) :: C, S
!f2py depend(nbf,nif) :: C
!f2py depend(nbf) :: S
 real(kind=8), intent(out) :: CTSC(nif,nif)
!f2py intent(out) :: CTSC
!f2py depend(nif) :: CTSC
 real(kind=8), allocatable :: SC(:,:)

 CTSC = 0d0
 allocate(SC(nbf,nif), source=0d0)
 call dsymm('L', 'L', nbf, nif, 1d0, S, nbf, C, nbf, 0d0, SC, nbf)
 call dgemm('T', 'N', nif, nif, nbf, 1d0, C, nbf, SC, nbf, 0d0, CTSC, nif)
 deallocate(SC)
end subroutine calc_CTSC

! calculate Cn(C^T), n must be real symmetric since dsymm is called
! C: nbf*nmo, n: nmo*nmo
subroutine calc_cnct(nbf, nmo, c, n, cnct)
 implicit none
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 real(kind=8), intent(in) :: c(nbf,nmo), n(nmo,nmo)
!f2py intent(in) :: c, n
!f2py depend(nbf,nmo) :: c
!f2py depend(nmo) :: n
 real(kind=8), intent(out) :: cnct(nbf,nbf)
!f2py intent(out) :: cnct
!f2py depend(nbf) :: cnct
 real(kind=8), allocatable :: cn(:,:)

 cnct = 0d0
 allocate(cn(nbf,nmo), source=0d0)
 call dsymm('R', 'L', nbf, nmo, 1d0, n, nmo, c, nbf, 0d0, cn, nbf)
 call dgemm('N', 'T', nbf, nbf, nmo, 1d0, cn, nbf, c, nbf, 0d0, cnct, nbf)
 deallocate(cn)
end subroutine calc_cnct

! Calculate (C^T)SC, where C is (double precision) real and S is (double)
! complex symmetric.
subroutine calc_ct_zs_c(nbf, nif, c, s, ctsc)
 implicit none
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: c(nbf,nif)
!f2py intent(in) :: c
!f2py depend(nbf,nif) :: c
 complex(kind=8), parameter :: zero=(0d0,0d0), one=(1d0,0d0)
 complex(kind=8), intent(in) :: s(nbf,nbf)
!f2py intent(in) :: s
!f2py depend(nbf) :: s
 complex(kind=8), intent(out) :: ctsc(nif,nif)
!f2py intent(out) :: ctsc
!f2py depend(nif) :: ctsc
 complex(kind=8), allocatable :: zc(:,:), sc(:,:)

 ctsc = zero
 ! Transform real matrix C to complex matrix. kind=8 is necessary in CMPLX()
 allocate(zc(nbf,nif), source=CMPLX(c, kind=8))
 allocate(sc(nbf,nif), source=zero)
 call zsymm('L', 'L', nbf, nif, one, s, nbf, zc, nbf, zero, sc, nbf)
 call zgemm('T', 'N', nif, nif, nbf, one, zc, nbf, sc, nbf, zero, ctsc, nif)
 deallocate(sc, zc)
end subroutine calc_ct_zs_c

! Update upper triangular part of i,j columns of a dipole integral matrix.
! Note: i<j is required.
subroutine update_up_tri_ij_dip(nmo, i, j, cos_a, sin_a, mo_dip)
 implicit none
 integer :: k
 integer, intent(in) :: nmo, i, j
 real(kind=8) :: tmp_dip(3)
 real(kind=8), intent(in) :: cos_a, sin_a
 real(kind=8), intent(inout) :: mo_dip(3,nmo,nmo)

 !$omp parallel sections private(k,tmp_dip)
 !$omp section
 do k = 1, i-1, 1
  tmp_dip = mo_dip(:,k,i)
  mo_dip(:,k,i) = cos_a*tmp_dip + sin_a*mo_dip(:,k,j)
  mo_dip(:,k,j) = cos_a*mo_dip(:,k,j) - sin_a*tmp_dip
 end do ! for k
 !$omp section
 do k = i+1, j-1, 1
  tmp_dip = mo_dip(:,i,k)
  mo_dip(:,i,k) = cos_a*tmp_dip + sin_a*mo_dip(:,k,j)
  mo_dip(:,k,j) = cos_a*mo_dip(:,k,j) - sin_a*tmp_dip
 end do ! for k
 !$omp section
 do k = j+1, nmo, 1
  tmp_dip = mo_dip(:,i,k)
  mo_dip(:,i,k) = cos_a*tmp_dip + sin_a*mo_dip(:,j,k)
  mo_dip(:,j,k) = cos_a*mo_dip(:,j,k) - sin_a*tmp_dip
 end do ! for k
 !$omp end parallel sections
end subroutine update_up_tri_ij_dip

! Update lower triangular part of i,j columns of a complex dipole integral matrix.
! Note: i<j is required.
subroutine update_up_tri_ij_zdip(nmo, i, j, cos_a, sin_a, mo_dip)
 implicit none
 integer :: k
 integer, intent(in) :: nmo, i, j
 real(kind=8), intent(in) :: cos_a, sin_a
 complex(kind=8) :: tmp_dip(3)
 complex(kind=8), intent(inout) :: mo_dip(3,nmo,nmo)

 !$omp parallel sections private(k,tmp_dip)
 !$omp section
 do k = 1, i-1, 1
  tmp_dip = mo_dip(:,k,i)
  mo_dip(:,k,i) = cos_a*tmp_dip + sin_a*mo_dip(:,k,j)
  mo_dip(:,k,j) = cos_a*mo_dip(:,k,j) - sin_a*tmp_dip
 end do ! for k
 !$omp section
 do k = i+1, j-1, 1
  tmp_dip = mo_dip(:,i,k)
  mo_dip(:,i,k) = cos_a*tmp_dip + sin_a*mo_dip(:,k,j)
  mo_dip(:,k,j) = cos_a*mo_dip(:,k,j) - sin_a*tmp_dip
 end do ! for k
 !$omp section
 do k = j+1, nmo, 1
  tmp_dip = mo_dip(:,i,k)
  mo_dip(:,i,k) = cos_a*tmp_dip + sin_a*mo_dip(:,j,k)
  mo_dip(:,j,k) = cos_a*mo_dip(:,j,k) - sin_a*tmp_dip
 end do ! for k
 !$omp end parallel sections
end subroutine update_up_tri_ij_zdip

! Update upper triangular part of i,j columns of a gross matrix.
! Note: i<j is required.
subroutine update_up_tri_ij_gross(natom, nmo, i, j, cos_a, sin_a, gross)
 implicit none
 integer :: k
 integer, intent(in) :: natom, nmo, i, j
 real(kind=8), allocatable :: tmp_g(:)
 real(kind=8), intent(in) :: cos_a, sin_a
 real(kind=8), intent(inout) :: gross(natom,nmo,nmo)

 allocate(tmp_g(natom))

!$omp parallel sections private(k,tmp_g)
!$omp section
 do k = 1, i-1, 1
  tmp_g = gross(:,k,i)
  gross(:,k,i) = cos_a*tmp_g + sin_a*gross(:,k,j)
  gross(:,k,j) = cos_a*gross(:,k,j) - sin_a*tmp_g
 end do ! for k
!$omp section
 do k = i+1, j-1, 1
  tmp_g = gross(:,i,k)
  gross(:,i,k) = cos_a*tmp_g + sin_a*gross(:,k,j)
  gross(:,k,j) = cos_a*gross(:,k,j) - sin_a*tmp_g
 end do ! for k
!$omp section
 do k = j+1, nmo, 1
  tmp_g = gross(:,i,k)
  gross(:,i,k) = cos_a*tmp_g + sin_a*gross(:,j,k)
  gross(:,j,k) = cos_a*gross(:,j,k) - sin_a*tmp_g
 end do ! for k
!$omp end parallel sections

 deallocate(tmp_g)
end subroutine update_up_tri_ij_gross

! transform AO dipole integrals into MO dipole integrals
subroutine ao2mo_dip(nbf, nmo, mo, ao_dip, mo_dip)
 implicit none
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 real(kind=8), intent(in) :: mo(nbf,nmo), ao_dip(3,nbf,nbf)
!f2py intent(in) :: mo, ao_dip
!f2py depend(nbf,nmo) :: mo
!f2py depend(nbf) :: ao_dip
 real(kind=8), intent(out) :: mo_dip(3,nmo,nmo)
!f2py intent(out) :: mo_dip
!f2py depend(nmo) :: mo_dip

 call calc_CTSC(nbf, nmo, mo, ao_dip(1,:,:), mo_dip(1,:,:))
 call calc_CTSC(nbf, nmo, mo, ao_dip(2,:,:), mo_dip(2,:,:))
 call calc_CTSC(nbf, nmo, mo, ao_dip(3,:,:), mo_dip(3,:,:))
end subroutine ao2mo_dip

! Transform complex AO dipole integrals into complex MO dipole integrals,
! where the MO coefficients are real.
! Note: ao_dip must be double complex and symmetric.
subroutine ao2mo_zdip(nbf, nmo, mo, ao_dip, mo_dip)
 implicit none
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 real(kind=8), intent(in) :: mo(nbf,nmo)
!f2py intent(in) :: mo
!f2py depend(nbf,nmo) :: mo
 complex(kind=8), intent(in) :: ao_dip(3,nbf,nbf)
!f2py intent(in) :: ao_dip
!f2py depend(nbf) :: ao_dip
 complex(kind=8), intent(out) :: mo_dip(3,nmo,nmo)
!f2py intent(out) :: mo_dip
!f2py depend(nmo) :: mo_dip

 call calc_ct_zs_c(nbf, nmo, mo, ao_dip(1,:,:), mo_dip(1,:,:))
 call calc_ct_zs_c(nbf, nmo, mo, ao_dip(2,:,:), mo_dip(2,:,:))
 call calc_ct_zs_c(nbf, nmo, mo, ao_dip(3,:,:), mo_dip(3,:,:))
end subroutine ao2mo_zdip

! TODO: change to mo_dip(nmo,nmo,3)
subroutine update_mo_zdip_by_u(nmo, u, mo_zdip)
 implicit none
 integer, intent(in) :: nmo
!f2py intent(in) :: nmo
 real(kind=8), intent(in) :: u(nmo,nmo)
!f2py intent(in) :: u
!f2py depend(nmo) :: u
 complex(kind=8), intent(inout) :: mo_zdip(3,nmo,nmo)
!f2py intent(in,out) :: mo_zdip
!f2py depend(nmo) :: mo_zdip
 complex(kind=8), allocatable :: r(:,:)

 allocate(r(nmo,nmo))
 call calc_ct_zs_c(nmo, nmo, u, mo_zdip(1,:,:), r)
 mo_zdip(1,:,:) = r
 call calc_ct_zs_c(nmo, nmo, u, mo_zdip(2,:,:), r)
 mo_zdip(2,:,:) = r
 call calc_ct_zs_c(nmo, nmo, u, mo_zdip(3,:,:), r)
 mo_zdip(3,:,:) = r
 deallocate(r)
end subroutine update_mo_zdip_by_u

! get the value of the Boys function sum_i <phi_i|r|phi_i>^2
subroutine get_fboys(nmo, mo_dip, fboys)
 implicit none
 integer :: i
 integer, intent(in) :: nmo
!f2py intent(in) :: nmo
 real(kind=8), intent(in) :: mo_dip(3,nmo,nmo)
!f2py intent(in) :: mo_dip
!f2py depend(nmo) :: mo_dip
 real(kind=8), intent(out) :: fboys
!f2py intent(out) :: fboys

 fboys = 0d0
!$omp parallel do schedule(static) default(private) shared(nmo,mo_dip) &
!$omp reduction(+:fboys)
 do i = 1, nmo, 1
  fboys = fboys + DOT_PRODUCT(mo_dip(:,i,i), mo_dip(:,i,i))
 end do ! for i
!$omp end parallel do
end subroutine get_fboys

! get the value of the Berry function sum_i <phi_i|r|phi_i>^2
subroutine get_fberry(nmo, mo_zdip, fberry)
 implicit none
 integer :: i
 integer, intent(in) :: nmo
!f2py intent(in) :: nmo
 complex(kind=8), intent(in) :: mo_zdip(3,nmo,nmo)
!f2py intent(in) :: mo_zdip
!f2py depend(nmo) :: mo_zdip
 real(kind=8), intent(out) :: fberry
!f2py intent(out) :: fberry
 real(kind=8), external :: cmplx_v3_square

 fberry = 0d0
!$omp parallel do schedule(static) default(private) shared(nmo,mo_zdip) &
!$omp reduction(+:fberry)
 do i = 1, nmo, 1
  fberry = fberry + cmplx_v3_square(mo_zdip(:,i,i))
 end do ! for i
!$omp end parallel do
end subroutine get_fberry

! calculate (C^T)S(C'), S must be real symmetric since dsymm is called
! C: nbf*nif  S: nbf*nbf, C': nbf*nif
subroutine calc_CTSCp(nbf, nif, C, S, Cp, CTSCp)
 implicit none
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: C(nbf,nif), S(nbf,nbf), Cp(nbf,nif)
!f2py intent(in) :: C, S, Cp
!f2py depend(nbf,nif) :: C, Cp
!f2py depend(nbf) :: S
 real(kind=8), intent(out) :: CTSCp(nif,nif)
!f2py intent(out) :: CTSCp
!f2py depend(nif) :: CTSCp
 real(kind=8), allocatable :: SCp(:,:)

 CTSCp = 0d0
 allocate(SCp(nbf,nif), source=0d0)
 call dsymm('L', 'L', nbf, nif, 1d0, S, nbf, Cp, nbf, 0d0, SCp, nbf)
 call dgemm('T', 'N', nif, nif, nbf, 1d0, C, nbf, SCp, nbf, 0d0, CTSCp, nif)
 deallocate(SCp)
end subroutine calc_CTSCp

subroutine calc_CTSCp2(nbf, nif1, nif2, C, S, Cp, CTSCp)
 implicit none
 integer, intent(in) :: nbf, nif1, nif2
 real(kind=8), intent(in) :: C(nbf,nif1), S(nbf,nbf), Cp(nbf,nif2)
 real(kind=8), intent(out) :: CTSCp(nif1,nif2)
 real(kind=8), allocatable :: SCp(:,:)

 CTSCp = 0d0
 allocate(SCp(nbf,nif2), source=0d0)
 call dsymm('L', 'L', nbf, nif2, 1d0, S, nbf, Cp, nbf, 0d0, SCp, nbf)
 call dgemm('T', 'N', nif1, nif2, nbf, 1d0, C, nbf, SCp, nbf, 0d0, CTSCp, nif1)
 deallocate(SCp)
end subroutine calc_CTSCp2

! calculate ((C_i)^T)S(C'_i), where C_i and C'_i are both vectors
function calc_CiTSCip(nbf1, nbf2, C, S, Cp) result(res)
 implicit none
 integer, intent(in) :: nbf1, nbf2
 real(kind=8), intent(in) :: C(nbf1), Cp(nbf2), S(nbf1,nbf2)
 real(kind=8), allocatable :: SCp(:)
 real(kind=8) :: ddot, res

 allocate(SCp(nbf1), source=0d0)
 call dgemv('N', nbf1, nbf2, 1d0, S, nbf1, Cp, 1, 0d0, SCp, 1)
 res = ddot(nbf1, C, 1, SCp, 1)
 deallocate(SCp)
end function calc_CiTSCip

! calculate CX(C^T), where X is a square matrix (symmetric is not required)
subroutine calc_CXCT(nbf, nmo, C, X, CXCT)
 implicit none
 integer, intent(in) :: nbf, nmo
 real(kind=8), intent(in) :: C(nbf,nmo), X(nmo,nmo)
 real(kind=8), intent(out) :: CXCT(nbf,nbf)
 real(kind=8), allocatable :: CX(:,:)

 CXCT = 0d0
 allocate(CX(nbf,nmo), source=0d0)
 call dgemm('N', 'N', nbf, nmo, nmo, 1d0, C, nbf, X, nmo, 0d0, CX, nbf)
 call dgemm('N', 'T', nbf, nbf, nmo, 1d0, CX, nbf, C, nbf, 0d0, CXCT, nbf)
 deallocate(CX)
end subroutine calc_CXCT

! calculate SPS, where S and P are both symmetric matrices
subroutine calc_sps(nbf, dm, s, sps)
 implicit none
 integer, intent(in) :: nbf
 real(kind=8), intent(in) :: dm(nbf,nbf), s(nbf,nbf)
 real(kind=8), intent(out) :: sps(nbf,nbf)
 real(kind=8), allocatable :: dm_s(:,:)

 sps = 0d0
 allocate(dm_s(nbf,nbf), source=0d0)
 call dsymm('R', 'L', nbf, nbf, 1d0, s, nbf, dm, nbf, 0d0, dm_s, nbf)
 call dsymm('L', 'L', nbf, nbf, 1d0, s, nbf, dm_s, nbf, 0d0, sps, nbf)
 deallocate(dm_s)
end subroutine calc_sps

! Calculate expectation values using the given unitary matrix and eigenvalues.
! E.g. LMO=CMO*U, where ev are orbital energies.
subroutine calc_expect_value(n, u, old_ev, new_ev)
 implicit none
 integer :: i
 integer, intent(in) :: n
!f2py intent(in) :: n
 real(kind=8), intent(in) :: u(n,n)
!f2py intent(in) :: u
!f2py depend(n) :: u
 real(kind=8), intent(in) :: old_ev(n)
!f2py intent(in) :: old_ev
!f2py depend(n) :: old_ev
 real(kind=8), intent(out) :: new_ev(n)
!f2py intent(out) :: new_ev
!f2py depend(n) :: new_ev
 real(kind=8), allocatable :: ev(:,:), utnu(:,:)

 allocate(ev(n,n), source=0d0)
 forall(i = 1:n) ev(i,i) = old_ev(i)
 allocate(utnu(n,n))
 call calc_CTSC(n, n, u, ev, utnu)
 deallocate(ev)
 forall(i = 1:n) new_ev(i) = utnu(i,i)
 deallocate(utnu)
end subroutine calc_expect_value

! Find a subset of MOs in mo1 to resemble mo2, which is a non-iterative method
! with nmo1>=nmo2. For example, mo2 contains RHF/STO-6G localized virtual MOs
! and mo1 is exactly the virtual space at cc-pVDZ. Now we want to find
! mo1(:,1:nmo2) which resembles mo2.
subroutine orb_resemble_ref1(nbf1, nmo1, mo1, nbf2, nmo2, mo2, cross_s, new_mo1)
 implicit none
 integer, intent(in) :: nbf1, nmo1, nbf2, nmo2
!f2py intent(in) :: nbf1, nmo1, nbf2, nmo2
 real(kind=8), intent(in) :: mo1(nbf1,nmo1), mo2(nbf2,nmo2), cross_s(nbf1,nbf2)
!f2py intent(in) :: mo1, mo2, cross_s
!f2py depend(nbf1,nmo1) :: mo1
!f2py depend(nbf2,nmo2) :: mo2
!f2py depend(nbf1,nbf2) :: cross_s
 real(kind=8), intent(out) :: new_mo1(nbf1,nmo1)
!f2py intent(out) :: new_mo1
!f2py depend(nbf1,nmo1) :: new_mo1
 real(kind=8), allocatable :: mo_ovlp(:,:), u(:,:), vt(:,:), rtmp(:,:), s(:)

 if(nmo1 < nmo2) then
  write(6,'(/,A)') 'ERROR in subroutine orb_resemble_ref1: this subroutine requ&
                   &ires nmo1>=nmo2.'
  write(6,'(2(A,I0))') 'But got nmo1=', nmo1, ', nmo2=', nmo2
  stop
 end if

 allocate(mo_ovlp(nmo1,nmo2))
 call calc_atbc(nmo1, nbf1, nbf2, nmo2, mo1, cross_s, mo2, mo_ovlp)
 allocate(u(nmo1,nmo1), vt(nmo2,nmo2), s(nmo1))
 call do_svd(nmo1, nmo2, mo_ovlp, u, vt, s)
 deallocate(mo_ovlp, s)
 new_mo1 = 0d0
 call dgemm('N','N', nbf1,nmo1,nmo1, 1d0,mo1,nbf1, u,nmo1, 0d0,new_mo1,nbf1)
 deallocate(u)

 allocate(rtmp(nbf1,nmo2), source=0d0)
 call dgemm('N','N', nbf1,nmo2,nmo2, 1d0,new_mo1(:,1:nmo2),nbf1, vt,nmo2, 0d0,&
            rtmp,nbf1)
 deallocate(vt)
 new_mo1(:,1:nmo2) = rtmp
 deallocate(rtmp)
end subroutine orb_resemble_ref1

! Update mo1 to resemble a set of MOs in mo2, which is a non-iterative method
! with nmo1<=nmo2. For example, mo2 contains orthogonal AOs (OAO) and mo1
! contains some MOs to be localized. Now we want to update mo1 which resembles
! nmo1 orbitals in mo2 as much as possible. mo1 and mo2 are assumed to share
! one AO overlap integral matrix.
!subroutine orb_resemble_ref3(nbf, nmo1, nmo2, mo1, mo2, ao_ovlp, new_mo1)
! implicit none
! integer, intent(in) :: nbf, nmo1, nmo2
!!f2py intent(in) :: nbf, nmo1, nmo2
! real(kind=8), intent(in) :: mo1(nbf,nmo1), mo2(nbf,nmo2), ao_ovlp(nbf,nbf)
!!f2py intent(in) :: mo1, mo2, ao_ovlp
!!f2py depend(nbf,nmo1) :: mo1
!!f2py depend(nbf,nmo2) :: mo2
!!f2py depend(nbf) :: ao_ovlp
! real(kind=8), intent(out) :: new_mo1(nbf,nmo1)
!!f2py intent(out) :: new_mo1
!!f2py depend(nbf,nmo1) :: new_mo1
!
! if(nmo1 > nmo2) then
!  write(6,'(/,A)') 'ERROR in subroutine orb_resemble_ref3: this subroutine requ&
!                   &ires nmo1<=nmo2.'
!  write(6,'(2(A,I0))') 'But got nmo1=', nmo1, ', nmo2=', nmo2
!  stop
! end if
!
! new_mo1 = mo1
!end subroutine orb_resemble_ref3

! compute the GAMMA_k matrix using the reference MOs and MOs of geometry/point k
! Note: here mo_ref and mo_k are both expressed at the orthogonal basis, i.e.
!  C' = (S^1/2)C
subroutine grassmann_C2GAMMA(nbf, nmo, mo_ref, mo_k)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nmo
 real(kind=8), intent(in) :: mo_ref(nbf,nmo)
 real(kind=8), intent(inout) :: mo_k(nbf,nmo) ! return G_k in mo_k
 real(kind=8), allocatable :: CTC(:,:), inv_CTC(:,:), u(:,:), vt(:,:), s(:), &
  s1(:,:), us(:,:), L_k(:,:)

 ! (mo_ref^T)mo_k
 allocate(CTC(nmo,nmo))
 call dgemm('T','N', nmo, nmo, nbf, 1d0, mo_ref, nbf, mo_k, nbf, 0d0, CTC, nmo)

 ! ((mo_ref^T)mo_k)^(-1)
 allocate(inv_CTC(nmo,nmo))
 call inverse(nmo, CTC, inv_CTC)
 ! Maybe using SVD to calculate the inverse would be more efficient
 deallocate(CTC)

 ! L_k = mo_k((mo_ref^T)mo_k)^(-1) - mo_ref
 allocate(L_k(nbf,nmo), source=mo_ref)
 call dgemm('N','N', nbf, nmo, nmo, 1d0, mo_k, nbf, inv_CTC, nmo, -1d0, L_k, nbf)
 deallocate(inv_CTC)

 allocate(u(nbf,nbf), vt(nmo,nmo), s(nbf))
 call do_svd(nbf, nmo, L_k, u, vt, s)
 deallocate(L_k)

 allocate(s1(nbf,nmo), source=0d0)
 forall(i = 1:nmo) s1(i,i) = DATAN(s(i))
 deallocate(s)

 allocate(us(nbf,nmo))
 call dgemm('N', 'N', nbf, nmo, nbf, 1d0, u, nbf, s1, nbf, 0d0, us, nbf)
 deallocate(u, s1)
 call dgemm('N', 'N', nbf, nmo, nmo, 1d0, us, nbf, vt, nmo, 0d0, mo_k, nbf)
 deallocate(vt, us)
end subroutine grassmann_C2GAMMA

! perform density matrix purification
subroutine purify_dm(nbf, S, P)
 implicit none
 integer :: i
 integer, intent(in) :: nbf
 integer, parameter :: max_it = 5
 real(kind=8) :: max_diff, mean_diff
 real(kind=8), intent(in) :: S(nbf,nbf)
 real(kind=8), intent(inout) :: P(nbf,nbf)
 real(kind=8), allocatable :: P0(:,:), SP(:,:), PSP(:,:), PSPSP(:,:)

 allocate(P0(nbf,nbf), PSP(nbf,nbf), SP(nbf,nbf), PSPSP(nbf,nbf))

 do i = 1, max_it, 1
  P0 = P
  call calc_sps(nbf, S, P, PSP)

  SP = 0d0
  call dsymm('L', 'L', nbf, nbf, 1d0, S, nbf, P, nbf, 0d0, SP, nbf)

  PSPSP = 0d0
  call dsymm('L', 'L', nbf, nbf, 1d0, PSP, nbf, SP, nbf, 0d0, PSPSP, nbf)

  P = 0.5d0*(3d0*PSP - PSPSP)

  P0 = DABS(P0 - P)
  max_diff = MAXVAL(P0)
  mean_diff = SUM(P0)/DBLE(nbf*nbf)
  write(6,'(2(A,F20.8))') 'max_v=', max_diff, ', mean_v=', mean_diff
  if(max_diff<1d-4 .and. mean_diff<1d-5) exit
 end do ! for i

 deallocate(P0, SP, PSP, PSPSP)
end subroutine purify_dm

! calculate density matrix using MO coefficients and occupation numbers
subroutine calc_dm_using_mo_and_on(nbf, nif, mo, noon, dm)
 implicit none
 integer :: u, v, k
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8) :: ddot
 real(kind=8), intent(in) :: mo(nbf,nif), noon(nif)
!f2py intent(in) :: mo, noon
!f2py depend(nif,nbf) :: mo
!f2py depend(nif) :: noon
 real(kind=8), intent(out) :: dm(nbf,nbf)
!f2py intent(out) :: dm
!f2py depend(nbf) :: dm
 real(kind=8), allocatable :: r(:)

 ! There is no need to initialize values in dm since each element in this
 ! array will be assigned a value below.
 allocate(r(nif))

!$omp parallel do schedule(dynamic) default(private) shared(nbf,nif,noon,mo,dm)
 do u = 1, nbf, 1
  forall(k = 1:nif) r(k) = noon(k)*mo(u,k)
  do v = 1, u, 1
   dm(v,u) = ddot(nif, r, 1, mo(v,:), 1)
  end do ! for v
 end do ! for u
!$omp end parallel do

 deallocate(r)
 call symmetrize_dmat(nbf, dm)
end subroutine calc_dm_using_mo_and_on

! get a random integer
! TODO: better random integer generator for Windows.
subroutine get_a_random_int(i)
 implicit none
 integer :: ierr, n, clock, fid
 integer, intent(out) :: i
 integer, allocatable :: seed(:)
 real(kind=4) :: r
 character(len=12), parameter :: fname = '/dev/urandom'

 clock = 0
 call RANDOM_SEED(size=n)
 allocate(seed(n), source=0)

#ifdef _WIN32
 call SYSTEM_CLOCK(count=clock)
 seed = clock
#else
 open(newunit=fid,file=fname,access='stream',form='unformatted',status='old',&
      iostat=ierr)
 if(ierr == 0) then
  read(fid) seed
  close(fid)
 else
  close(fid)
  call SYSTEM_CLOCK(count=seed(1))
 end if
#endif

 call RANDOM_SEED(put=seed)
 call RANDOM_NUMBER(r)
 deallocate(seed)

 i = CEILING(r*1e6)
end subroutine get_a_random_int

! Compute X=Us^(-1/2)(U^T) or Us^(-1/2) from AO overlap integral matrix S
subroutine solve_x_from_ao_ovlp(nbf, nif, S, X)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: S(nbf,nbf)
 real(kind=8), intent(out) :: X(nbf,nif)
 real(kind=8), allocatable :: U(:,:), w(:)

 if(nbf == nif) then
  ! no linear dependency, use symmetry orthogonalization
  allocate(U(nbf,nbf))
  call mat_dsqrt(nbf, S, .true., U, X)
  deallocate(U)
 else if(nbf > nif) then
  ! linear dependent, use canonical orthogonalization
  allocate(w(nbf))
  allocate(U(nbf,nbf), source=S)
  call diag_get_e_and_vec(nbf, U, w) ! ascending order
  call reverse_e_and_vec(nbf, U, w)  ! descending order
!$omp parallel do schedule(dynamic) default(shared) private(i)
  do i = 1, nif, 1
   X(:,i) = U(:,i)/DSQRT(w(i))
  end do ! for i
!$omp end parallel do
  deallocate(U, w)
 else
  write(6,'(/,A)') 'ERROR in subroutine solve_x_from_ao_ovlp: nbf>nif. Impossible.'
  write(6,'(2(A,I0))') 'nbf=', nbf, ', nif=', nif
  stop
 end if
end subroutine solve_x_from_ao_ovlp

! solve AO-based overlap matrix (S) from condition (C^T)SC=I
! Note: this subroutine only applies to nbf=nif, i.e. no linear dependence
subroutine solve_ovlp_from_ctsc(nbf, C, S)
 implicit none
 integer :: i
 integer, intent(in) :: nbf
 real(kind=8), intent(in) :: C(nbf,nbf)
 real(kind=8), intent(out) :: S(nbf,nbf)
 real(kind=8), allocatable :: SC(:,:)

 S = 0d0   ! this initialization is important
 forall(i = 1:nbf) S(i,i) = 1d0 ! unit matrix I

 allocate(SC(nbf,nbf))
 call solve_multi_lin_eqs(nbf, nbf, TRANSPOSE(C), nbf, S, SC)
 ! SC = X -> (C^T)S = X^T
 call solve_multi_lin_eqs(nbf, nbf, TRANSPOSE(C), nbf, TRANSPOSE(SC), S)
end subroutine solve_ovlp_from_ctsc

! solve AO-based overlap matrix (S) by calculating (C(C^T))^(-1)
! This subroutine has the same result as subroutine solve_ovlp_from_ctsc.
subroutine solve_ovlp_from_cct(nbf, C, S)
 implicit none
 integer, intent(in) :: nbf
 real(kind=8), intent(in) :: C(nbf,nbf)
 real(kind=8), intent(out) :: S(nbf,nbf)
 real(kind=8), allocatable :: cct(:,:)

 ! calculate C(C^T)
 allocate(cct(nbf,nbf))
 call calc_cct(nbf, nbf, C, cct)

 ! calculate (C(C^T))^(-1)
 call inverse(nbf, cct, S)
 deallocate(cct)
end subroutine solve_ovlp_from_cct

! solve AO-based Fock matrix (F) from condition (C^T)FC=E
subroutine solve_fock_from_ctfc(nbf, nif, C, E, F)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: C(nbf,nif), E(nif)
!f2py intent(in) :: C, E
!f2py depend(nbf,nif) :: C
!f2py depend(nif) :: E
 real(kind=8), intent(out) :: F(nbf,nbf)
!f2py intent(out) :: F
!f2py depend(nbf) :: F
 real(kind=8), allocatable :: FC(:,:), E1(:,:)

 if(nbf > nif) then
  write(6,'(/,A)') 'Warning from subroutine solve_fock_from_ctfc: nbf > nif. Ba&
                   &sis set linear'
  write(6,'(A)') 'dependency detected. The AO Fock matrix obtained by this subr&
                 &outine may be'
  write(6,'(A)') 'nonsense. You need to check.'
 end if

 allocate(E1(nif,nif), source=0d0)
 forall(i = 1:nif) E1(i,i) = E(i) ! diagonal matrix
 allocate(FC(nbf,nif))
 call solve_multi_lin_eqs(nif, nbf, TRANSPOSE(C), nif, E1, FC)
 deallocate(E1)
 ! FC = X -> (C^T)(F^T) = X^T, (F^T) = F
 call solve_multi_lin_eqs(nif, nbf, TRANSPOSE(C), nif, TRANSPOSE(FC), F)
end subroutine solve_fock_from_ctfc

! construct partial/all virtual orbitals using the PAO (projected atomic orbitals)
subroutine construct_vir(nbf, nif, idx, coeff, ovlp, new_mo)
 implicit none
 integer :: i, j, nvir
 integer, intent(in) :: nbf, nif, idx
!f2py intent(in) :: nbf, nif, idx
 ! nbf: the number of basis functions
 ! nif: the number of MOs
 ! idx: the beginning index (Fortran convention) of the virtual MOs
 real(kind=8), intent(in) :: coeff(nbf,nif), ovlp(nbf,nbf)
!f2py intent(in) :: coeff, ovlp
!f2py depend(nbf) :: ovlp
!f2py depend(nbf,nif) :: coeff
 real(kind=8), intent(out) :: new_mo(nbf,nif)
!f2py intent(out) :: new_mo
!f2py depend(nbf,nif) :: new_mo
 real(kind=8), allocatable :: p(:,:), v(:,:), s1(:,:), ev(:), x(:,:)
 ! V: projected atomic orbitals (PAO)
 ! P: density matrix of atomic basis functions, sigma_i(Cui*Cvi)
 ! Note that the index idx can be larger than ndb+1, in which case we only
 !  construct part of virtual orbitals.

 if(idx == nif+1) then
  write(6,'(/,A)') 'Warning in subroutine construct_vir: idx=nif+1 found.'
  write(6,'(A)') 'No need to construct virtual MOs.'
  new_mo = coeff
  return
 end if
 new_mo(:,1:idx-1) = coeff(:,1:idx-1) ! occupied space

 ! Step 1: P = sigma_i(Cui*Cvi)
 allocate(p(nbf,nbf), source=0d0)
 call dgemm('N','T', nbf,nbf,idx-1, 1d0,coeff(1:nbf,1:idx-1),nbf, &
            coeff(1:nbf,1:idx-1),nbf, 0d0,p,nbf)

 ! Step 2: V = 1 - PS
 allocate(v(nbf, nbf), source=0d0)
 forall(i = 1:nbf) v(i,i) = 1d0
 call dsymm('R', 'L', nbf, nbf, -1d0, ovlp, nbf, p, nbf, 1d0, v, nbf)
 deallocate(p)

 ! Step 3: S1 = (VT)SV
 allocate(s1(nbf,nbf))
 call calc_CTSC(nbf, nbf, v, ovlp, s1)

 ! Step 4: diagonalize S1 (note that S1 is symmetric) and get X
 allocate(ev(nbf))
 call diag_get_e_and_vec(nbf, s1, ev)
 nvir = nif - idx + 1
 forall(i = nbf-nvir+1:nbf) ev(i) = 1d0/DSQRT(ev(i))
 allocate(x(nbf,nvir), source=0d0)
 forall(i=1:nbf, j=1:nvir) x(i,j) = s1(i,nbf-nvir+j)*ev(nbf-nvir+j)
 deallocate(ev, s1)

 ! Step 5: get new virtual MO coefficients
 call dgemm('N','N', nbf,nvir,nbf, 1d0,v,nbf, x,nbf, 0d0,new_mo(:,idx:nif),nbf)
 deallocate(x, v)

 ! Step 6: check orthonormality
 allocate(x(nif,nif))
 call calc_CTSC(nbf, nif, new_mo, ovlp, x)

 forall(i = 1:nif) x(i,i) = x(i,i) - 1d0
 x = DABS(x)
 write(6,'(/,A)') 'The orthonormality of MOs after PAO construction:'
 write(6,'(A,F16.10)') 'maxv=', MAXVAL(x)
 write(6,'(A,F16.10)') 'abs_mean=', SUM(x)/DBLE(nif*nif)
 deallocate(x)
end subroutine construct_vir

subroutine get_occ_from_na_nb(nif, na, nb, occ)
 implicit none
 integer, intent(in) :: nif, na, nb
!f2py intent(in) :: nif, na, nb
 real(kind=8), intent(out) :: occ(nif)
!f2py intent(out) :: occ
!f2py depend(nif) :: occ

 occ(1:nb) = 2d0
 if(na > nb) occ(nb+1:na) = 1d0
 occ(na+1:) = 0d0
end subroutine get_occ_from_na_nb

subroutine get_occ_from_na_nb2(nif, na, nb, occ)
 implicit none
 integer, intent(in) :: nif, na, nb
!f2py intent(in) :: nif, na, nb
 real(kind=8), intent(out) :: occ(2,nif)
!f2py intent(out) :: occ
!f2py depend(nif) :: occ

 occ = 0d0
 occ(1,1:na) = 1d0
 occ(2,1:nb) = 1d0
end subroutine get_occ_from_na_nb2

! canonicalize MOs (digonalize the Fock matrix to get canonical MOs)
subroutine canonicalize_mo(nbf, nmo, f, old_mo, new_mo)
 implicit none
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 real(kind=8), intent(in) :: f(nbf,nbf), old_mo(nbf,nmo)
!f2py intent(in) :: f, old_mo
!f2py depend(nbf) :: f
!f2py depend(nbf,nmo) :: old_mo
 real(kind=8), intent(out) :: new_mo(nbf,nmo)
!f2py intent(out) :: new_mo
!f2py depend(nbf,nmo) :: new_mo
 real(kind=8), allocatable :: CTFC(:,:), orb_e(:)

 allocate(CTFC(nmo,nmo))
 call calc_CTSC(nbf, nmo, old_mo, f, CTFC)
 allocate(orb_e(nmo))
 call diag_get_e_and_vec(nmo, CTFC, orb_e)
 deallocate(orb_e)
 new_mo = 0d0
 call dgemm('N','N', nbf,nmo,nmo, 1d0,old_mo,nbf, CTFC,nmo, 0d0,new_mo,nbf)
 deallocate(CTFC)
end subroutine canonicalize_mo

! Calculate diagonal elements of Mulliken/Lowdin gross population matrix.
! TODO: if one wants to frequently call this subroutine to calculate the Lowdin
!  population, the rootS is supposed to be given as an optional input, so that
!  rootS would not be calculated multiple times.
subroutine calc_diag_gross_pop(natom, nbf, nif, bfirst, ao_ovlp, mo, popm, gross)
 implicit none
 integer :: i, j, i1, i2, i3
 integer, intent(in) :: natom, nbf, nif
 integer, intent(in) ::  bfirst(natom+1)
 real(kind=8) :: ddot
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf), mo(nbf,nif)
 real(kind=8), intent(out) :: gross(natom,nif)
 real(kind=8), allocatable :: rootS(:,:), SC(:,:)
 character(len=*), intent(in) :: popm

 select case(TRIM(popm))
 case('mulliken')
  allocate(SC(nbf,nif), source=0d0)
  call dsymm('L', 'L', nbf, nif, 1d0, ao_ovlp, nbf, mo, nbf, 0d0, SC, nbf)
!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(nif,natom,bfirst,mo,SC,gross)
  do i = 1, nif, 1
   do j = 1, natom, 1
    i1 = bfirst(j); i2 = bfirst(j+1) - 1
    i3 = i2 - i1 + 1
    gross(j,i) = ddot(i3, mo(i1:i2,i), 1, SC(i1:i2,i), 1)
   end do ! for j
  end do ! for i
!$omp end parallel do

 case('lowdin')
  allocate(rootS(nbf,nbf), SC(nbf,nbf))
  call mat_dsqrt(nbf, ao_ovlp, .false., rootS, SC)
  deallocate(SC)
  ! Note that Lowdin populations do not require S^(-1/2), but only S^1/2
  allocate(SC(nbf,nif), source=0d0) ! use SC to store (S^1/2)C
  call dsymm('L', 'L', nbf, nif, 1d0, rootS, nbf, mo, nbf, 0d0, SC, nbf)
  deallocate(rootS)
!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(nif,natom,bfirst,SC,gross)
  do i = 1, nif, 1
   do j = 1, natom, 1
    i1 = bfirst(j); i2 = bfirst(j+1) - 1
    i3 = i2 - i1 + 1
    gross(j,i) = ddot(i3, SC(i1:i2,i), 1, SC(i1:i2,i), 1)
   end do ! for j
  end do ! for i
!$omp end parallel do

 case default
  write(6,'(/,A)') 'ERROR in subroutine calc_diag_gross_pop: wrong population m&
                   &ethod provided.'
  write(6,'(A)') "Only 'mulliken' or 'lowdin' supported. But input popm="//popm
  stop
 end select

 deallocate(SC)
end subroutine calc_diag_gross_pop

! Calculate Mulliken/Lowdin gross population matrix.
! TODO: if one wants to frequently call this subroutine to calculate the Lowdin
!  population, the rootS is supposed to be given as an optional input, so that
!  rootS would not be calculated multiple times.
subroutine calc_gross_pop(natom, nbf, nif, bfirst, ao_ovlp, mo, popm, gross)
 implicit none
 integer :: i, j, k, m, i1, i2, i3, np
 integer, intent(in) :: natom, nbf, nif
 integer, intent(in) ::  bfirst(natom+1)
 integer, allocatable :: map(:,:)
 real(kind=8) :: ddot
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf), mo(nbf,nif)
 real(kind=8), intent(out) :: gross(natom,nif,nif)
 real(kind=8), allocatable :: rootS(:,:), SC(:,:)
 character(len=*), intent(in) :: popm

 np = nif*(nif+1)/2
 allocate(map(2,np))
 call get_triu_idx(nif, map)

 select case(TRIM(popm))
 case('mulliken')
  allocate(SC(nbf,nif), source=0d0)
  call dsymm('L', 'L', nbf, nif, 1d0, ao_ovlp, nbf, mo, nbf, 0d0, SC, nbf)

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(np,natom,map,bfirst,mo,SC,gross)
  do m = 1, np, 1
   i = map(1,m); j = map(2,m)
   do k = 1, natom, 1
    i1 = bfirst(k); i2 = bfirst(k+1) - 1
    i3 = i2 - i1 + 1
    gross(k,j,i) = 0.5d0*(ddot(i3, mo(i1:i2,i), 1, SC(i1:i2,j), 1) + &
                          ddot(i3, mo(i1:i2,j), 1, SC(i1:i2,i), 1) )
    gross(k,i,j) = gross(k,j,i)
   end do ! for k
  end do ! for m
!$omp end parallel do

 case('lowdin')
  allocate(rootS(nbf,nbf), SC(nbf,nbf))
  call mat_dsqrt(nbf, ao_ovlp, .false., rootS, SC)
  deallocate(SC)
  ! Note that Lowdin populations do not require S^(-1/2), but only S^1/2
  allocate(SC(nbf,nif), source=0d0) ! use SC to store (S^1/2)C
  call dsymm('L', 'L', nbf, nif, 1d0, rootS, nbf, mo, nbf, 0d0, SC, nbf)
  deallocate(rootS)

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(np,natom,map,bfirst,SC,gross)
  do m = 1, np, 1
   i = map(1,m); j = map(2,m)
   do k = 1, natom, 1
    i1 = bfirst(k); i2 = bfirst(k+1) - 1
    i3 = i2 - i1 + 1
    gross(k,j,i) = ddot(i3, SC(i1:i2,i), 1, SC(i1:i2,j), 1)
    gross(k,i,j) = gross(k,j,i)
   end do ! for k
  end do ! for m
!$omp end parallel do

 case default
  write(6,'(/,A)') 'ERROR in subroutine calc_gross_pop: wrong population method&
                   & provided.'
  write(6,'(A)') "Only 'mulliken' or 'lowdin' supported. But input popm="//popm
  stop
 end select

 deallocate(SC, map)
end subroutine calc_gross_pop

! Find the centers of each MO, using the population matrix. The population method
!  is determined when generating the pop array, so we do not need to know the
! population method in this subroutine.
! Note: the maximum number of centers for each MO is 4.
subroutine get_mo_center_from_pop(natom, nmo, pop, mo_center)
 implicit none
 integer :: i, j, k, m, ak(1)
 integer, intent(in) :: natom, nmo
 integer, intent(out) :: mo_center(0:4,nmo)
 real(kind=8) :: r
 real(kind=8), parameter :: diff = 0.15d0, pop_thres = 0.7d0
 ! diff: difference between the largest and the 2nd largest component
 real(kind=8), intent(in) :: pop(natom,nmo)

 mo_center = 0

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(natom,nmo,pop,mo_center)
 do i = 1, nmo, 1
  ! the largest component on an atom of an orbital
  ak = MAXLOC(pop(:,i)); k = ak(1); r = pop(k,i)
  mo_center(0,i) = 1; mo_center(1,i) = k; m = 1
  ! if this is lone pair, no need to check the 2nd largest component
  if(r > pop_thres) cycle

  ! find the 2nd largest component and so on
  do j = 1, natom, 1
   if(j == k) cycle
   if(r - pop(j,i) < diff) then
    m = m + 1
    if(m > 4) then
     write(6,'(/,A)') 'ERROR in subroutine get_mo_center_from_pop: ncenter>4. M&
                      &Os are too'
     write(6,'(A,2I7)') 'delocalized. natom, nmo=', natom, nmo
     stop
    end if
    mo_center(m,i) = j
   end if
  end do ! for j

  mo_center(0,i) = m
 end do ! for i
!$omp end parallel do

end subroutine get_mo_center_from_pop

! calculate the distance matrix from a set of coordinates
subroutine calc_dis_mat_from_coor(natom, coor, dis)
 implicit none
 integer :: i, j, k, m, n
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 integer, allocatable :: map(:,:)
 real(kind=8) :: r(3)
 real(kind=8), intent(in) :: coor(3,natom)
!f2py intent(in) :: coor
!f2py depend(natom) :: coor
 real(kind=8), intent(out) :: dis(natom,natom)
!f2py intent(out) :: dis
!f2py depend(natom) :: dis

 n = natom
 forall(i = 1:n) dis(i,i) = 0d0
 if(n == 1) return

 k = n*(n-1)/2
 allocate(map(2,k))
 call get_triu_idx1(n, map)

!$omp parallel do schedule(dynamic) default(private) shared(k,map,coor,dis)
 do m = 1, k, 1
  i = map(1,m); j = map(2,m)
  r = coor(:,i) - coor(:,j)
  dis(j,i) = DSQRT(DOT_PRODUCT(r,r))
  dis(i,j) = dis(j,i)
 end do ! for m
!$omp end parallel do

 deallocate(map)
end subroutine calc_dis_mat_from_coor

! calculate the distance matrix from a set of coordinates in a cell
subroutine calc_dis_mat_from_coor_pbc(natom, cell, coor, dis)
 implicit none
 integer :: i, j, k, m, n, np
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 integer, allocatable :: map(:,:)
 real(kind=8) :: tmp_dis, min_dis, r1(3), r2(3)
 real(kind=8), intent(in) :: cell(3,3), coor(3,natom)
!f2py intent(in) :: cell, coor
!f2py depend(natom) :: coor
 real(kind=8), intent(out) :: dis(natom,natom)
!f2py intent(out) :: dis
!f2py depend(natom) :: dis
 real(kind=8), allocatable :: coor1(:,:,:)

 n = natom
 forall(i = 1:n) dis(i,i) = 0d0
 if(n == 1) return
 allocate(coor1(3,2:27,natom))

!$omp parallel do schedule(static) default(shared) private(i)
 do i = 1, natom, 1
  coor1(:, 2,i) = coor(:,i) + cell(:,1)
  coor1(:, 3,i) = coor(:,i) + cell(:,2)
  coor1(:, 4,i) = coor(:,i) + cell(:,3)
  coor1(:, 5,i) = coor(:,i) - cell(:,1)
  coor1(:, 6,i) = coor(:,i) - cell(:,2)
  coor1(:, 7,i) = coor(:,i) - cell(:,3)
  coor1(:, 8,i) = coor(:,i) + cell(:,1) + cell(:,2)
  coor1(:, 9,i) = coor(:,i) + cell(:,1) - cell(:,2)
  coor1(:,10,i) = coor(:,i) - cell(:,1) + cell(:,2)
  coor1(:,11,i) = coor(:,i) - cell(:,1) - cell(:,2)
  coor1(:,12,i) = coor(:,i) + cell(:,1) + cell(:,3)
  coor1(:,13,i) = coor(:,i) + cell(:,1) - cell(:,3)
  coor1(:,14,i) = coor(:,i) - cell(:,1) + cell(:,3)
  coor1(:,15,i) = coor(:,i) - cell(:,1) - cell(:,3)
  coor1(:,16,i) = coor(:,i) + cell(:,2) + cell(:,3)
  coor1(:,17,i) = coor(:,i) + cell(:,2) - cell(:,3)
  coor1(:,18,i) = coor(:,i) - cell(:,2) + cell(:,3)
  coor1(:,19,i) = coor(:,i) - cell(:,2) - cell(:,3)
  coor1(:,20,i) = coor(:,i) + cell(:,1) + cell(:,2) + cell(:,3)
  coor1(:,21,i) = coor(:,i) + cell(:,1) + cell(:,2) - cell(:,3)
  coor1(:,22,i) = coor(:,i) + cell(:,1) - cell(:,2) + cell(:,3)
  coor1(:,23,i) = coor(:,i) - cell(:,1) + cell(:,2) + cell(:,3)
  coor1(:,24,i) = coor(:,i) + cell(:,1) - cell(:,2) - cell(:,3)
  coor1(:,25,i) = coor(:,i) - cell(:,1) + cell(:,2) - cell(:,3)
  coor1(:,26,i) = coor(:,i) - cell(:,1) - cell(:,2) + cell(:,3)
  coor1(:,27,i) = coor(:,i) - cell(:,1) - cell(:,2) - cell(:,3)
 end do ! for i
!$omp end parallel do

 np = n*(n-1)/2
 allocate(map(2,np))
 call get_triu_idx1(n, map)

!$omp parallel do schedule(static) default(private) &
!$omp shared(np,map,coor,coor1,dis)
 do m = 1, np, 1
  i = map(1,m); j = map(2,m)
  r1 = coor(:,i)
  r2 = r1 - coor(:,j)
  min_dis = DSQRT(DOT_PRODUCT(r2,r2))
  do k = 2, 27
   r2 = r1 - coor1(:,k,j)
   tmp_dis = DSQRT(DOT_PRODUCT(r2,r2))
   if(tmp_dis < min_dis) min_dis = tmp_dis
  end do ! for k
  dis(j,i) = min_dis
  dis(i,j) = min_dis
 end do ! for m
!$omp end parallel do

 deallocate(coor1, map)
end subroutine calc_dis_mat_from_coor_pbc

! check whether the input matrix is a diagonal matrix
function check_diagonal_mat(n, a) result(diagonal)
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 real(kind=8) :: r1, r2, s, maxv
 real(kind=8), intent(in) :: a(n,n)
 real(kind=8), parameter :: thres = 1d-5
 real(kind=8), allocatable :: abs_a(:,:)
 logical :: diagonal

 if(n == 1) then
  diagonal = .true.
  return
 end if

 allocate(abs_a(n,n), source=DABS(a))
 diagonal = .false.; s = 0d0; maxv = abs_a(2,1)

 do i = 1, n-1, 1
  do j = i+1, n, 1
   r1 = abs_a(j,i); r2 = abs_a(i,j)
   s = s + r1 + r2
   if(r1 > maxv) maxv = r1
   if(r2 > maxv) maxv = r2
  end do ! for j
 end do ! or i

 deallocate(abs_a)
 if(maxv<thres .and. s/DBLE(n*(n-1))<thres) diagonal = .true.
end function check_diagonal_mat

! compute the unitary matrix U between two sets of MOs
! --------------------------------------------------
!  nbf: the number of basis functions
!  nmo: the number of MOs
!  coeff: old MO Coefficients
!  lo_coeff: new MO Coefficients
!  u: the unitary(orthogonal) matrix to be computed
!  lo_coeff1 = coeff*U
! --------------------------------------------------
subroutine get_u(nbf, nmo, coeff, lo_coeff, u)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 integer, allocatable :: ipiv(:)
 real(kind=8), intent(in) :: coeff(nbf,nmo), lo_coeff(nbf,nmo)
!f2py intent(in) :: coeff, lo_coeff
!f2py depend(nbf,nmo) :: coeff, lo_coeff
 real(kind=8), intent(out) :: u(nmo,nmo)
!f2py intent(out) :: u
!f2py depend(nmo) :: u
 real(kind=8), allocatable :: coeff1(:,:), lo_coeff1(:,:)

 u = 0d0
 allocate(coeff1(nbf,nmo), source=coeff)
 allocate(lo_coeff1(nbf,nmo), source=lo_coeff)
 allocate(ipiv(min(nbf,nmo)), source=0)

 call dgetrf(nbf, nmo, coeff1, nbf, ipiv, i)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine get_u: dgetrf info/=0.'
  write(6,'(A,I0)') 'info=', i
  stop
 end if

 call dgetrs('N', nmo, nmo, coeff1, nbf, ipiv, lo_coeff1, nbf, i)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine get_u: dgetrs info/=0.'
  write(6,'(A,I0)') 'info=', i
  stop
 end if

 deallocate(ipiv, coeff1)
 u = lo_coeff1(1:nmo,1:nmo)
 deallocate(lo_coeff1)
end subroutine get_u

! perform QR factorization on a matrix A
subroutine qr_fac(m, n, A, Q, R)
 implicit none
 integer :: i, j, lwork
 integer, intent(in) :: m, n
 real(kind=8), intent(in) :: A(m,n)
 real(kind=8), intent(out) :: Q(m,n), R(n,n)
 real(kind=8), allocatable :: tau(:), work(:)

 if(m < n) then
  write(6,'(/,A)') 'ERROR in subroutine qr_fac: m < n. You can use A^T as the i&
                   &nput.'
  write(6,'(2(A,I0))') 'm=', m, ', n=', n
  stop
 end if

 Q = A; R = 0d0
 allocate(tau(n), source=0d0)
 lwork = -1
 allocate(work(1))
 call dgeqrf(m, n, Q, m, tau, work, lwork, i)
 lwork = CEILING(work(1))
 deallocate(work)
 allocate(work(lwork))
 call dgeqrf(m, n, Q, m, tau, work, lwork, i)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine qr_fac: failed to call dgeqrf.'
  write(6,'(A,I0)') 'info=', i
  stop
 end if

 deallocate(work)
!$omp parallel do schedule(dynamic) default(shared) private(i,j)
 do i = 1, n, 1
  do j = 1, i, 1
   R(j,i) = Q(j,i)
  end do ! for j
 end do ! for i
!$omp end parallel do

 ! calculate lwork again since a larger integer may be required here
 lwork = -1
 allocate(work(1))
 call dorgqr(m, n, n, Q, m, tau, work, lwork, i)
 lwork = CEILING(work(1))
 deallocate(work)
 allocate(work(lwork))
 call dorgqr(m, n, n, Q, m, tau, work, lwork, i)
 deallocate(tau, work)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine qr_fac: failed to call dorgqr.'
  write(6,'(A,I0)') 'info=', i
  stop
 end if
end subroutine qr_fac

! convert spin multiplicity to spin square <S^2>
pure function mult2ssquare(mult) result(ss)
 implicit none
 integer, intent(in) :: mult
 real(kind=8) :: ss

 ss = 0.25d0*DBLE(mult*mult - 1)
end function mult2ssquare

! `findloc` is in Fortran 2008, which cannot be used in gfortran 4.8.5, so use
! this function instead.
function find_1st_loc(n, a, i) result(k)
 implicit none
 integer :: k
 integer, intent(in) :: n, i
 integer, intent(in) :: a(n)

 do k = 1, n, 1
  if(a(k) == i) return
 end do ! for k

 if(k == n+1) k = 0
end function find_1st_loc

! Detect the number of MOs with all zero coefficients in the array mo.
! Note: detection begins from the last MO, i.e. mo(:,nif).
subroutine detect_zero_mo(nbf, nif, mo, nzero)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nif
 integer, intent(out) :: nzero
 real(kind=8), parameter :: thres = 1d-7
 real(kind=8), intent(in) :: mo(nbf,nif)

 nzero = 0
 do i = nif, 1, -1
  if(SUM(DABS(mo(:,i))) < thres) then
   nzero = nzero + 1
  else
   exit
  end if
 end do ! for i
end subroutine detect_zero_mo

! Move doubly occupied MOs (whose 1e expectation values are degenerate with
!  that of any bonding orbitals) below/close to the bonding orbitals. These
!  special doubly occupied MOs can be viewed as extra/new bonding orbitals.
!  Their corresponding antibonding orbitals may be found by calling other sub-
!  routines later.
! nb: the number of beta electrons
! npair: the number of bonding orbitals already exists
subroutine mv_dege_docc_below_bo(nbf, nb, npair, ev, mo, new_ev, new_mo)
 implicit none
 integer :: i, j, k, ndocc
 integer, intent(in) :: nbf, nb, npair
!f2py intent(in) :: nbf, nb, npair
 real(kind=8) :: r
 real(kind=8), parameter :: thres = 1d-3 ! |<i|F|i> - <j|F|j>|
 real(kind=8), intent(in) :: ev(nb), mo(nbf,nb)
!f2py intent(in) :: ev, mo
!f2py depend(nb) :: ev
!f2py depend(nbf,nb) :: mo
 real(kind=8), intent(out) :: new_ev(nb), new_mo(nbf,nb)
!f2py intent(out) :: new_ev, new_mo
!f2py depend(nb) :: new_ev
!f2py depend(nbf,nb) :: new_mo
 real(kind=8), allocatable :: mo_i(:)
 logical :: degenerate

 new_ev = ev; new_mo = mo
 ndocc = nb - npair; k = ndocc
 allocate(mo_i(nbf))

 do i = ndocc, 1, -1
  degenerate = .false.
  r = new_ev(i)

  do j = ndocc+1, nb, 1
   if(DABS(r - new_ev(j)) < thres) then
    degenerate = .true.
    exit
   end if
  end do ! for j

  if(degenerate) then
   if(i < k) then
    mo_i = new_mo(:,i)
    do j = i, k-1, 1
     new_ev(j) = new_ev(j+1)
     new_mo(:,j) = new_mo(:,j+1)
    end do ! for j
    new_ev(k) = r
    new_mo(:,k) = mo_i
   end if
   k = k - 1 ! remember to update k
  end if
 end do ! for i

 deallocate(mo_i)
end subroutine mv_dege_docc_below_bo

! generate natural orbitals from provided density matrix and overlap matrix
subroutine gen_no_from_dm_and_ao_ovlp(nbf, nif, P, ao_ovlp, noon, new_coeff)
 implicit none
 integer :: i, j, lwork, liwork
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 integer, allocatable :: isuppz(:), iwork(:)
 real(kind=8), intent(in) :: P(nbf,nbf), ao_ovlp(nbf,nbf)
!f2py intent(in) :: P, ao_ovlp
!f2py depend(nbf) :: P, ao_ovlp
 real(kind=8), intent(out) :: noon(nif), new_coeff(nbf,nif)
!f2py intent(out) :: noon, new_coeff
!f2py depend(nif) :: noon
!f2py depend(nbf,nif) :: new_coeff
 real(kind=8), allocatable :: S(:,:), sqrt_S(:,:), n_sqrt_S(:,:)
 real(kind=8), allocatable :: e(:), U(:,:), work(:)

 noon = 0d0; new_coeff = 0d0 ! initialization
 allocate(S(nbf,nbf), source=ao_ovlp)
 allocate(sqrt_S(nbf,nbf), n_sqrt_S(nbf,nbf))
 call mat_dsqrt(nbf, S, .true., sqrt_S, n_sqrt_S) ! solve S^1/2 and S^-1/2
 call calc_sps(nbf, P, sqrt_S, S) ! use S to store (S^1/2)P(S^1/2)
 deallocate(sqrt_S)

 lwork = -1; liwork = -1
 allocate(work(1), iwork(1), isuppz(2*nbf), e(nbf), U(nbf,nbf))
 call dsyevr('V', 'A',  'L', nbf, S, nbf, 0d0, 0d0, 0, 0, 1d-8, i, e, U, nbf, &
             isuppz, work, lwork, iwork, liwork, j)
 lwork = CEILING(work(1))
 liwork = iwork(1)
 deallocate(work, iwork)
 allocate(work(lwork), iwork(liwork))
 call dsyevr('V', 'A',  'L', nbf, S, nbf, 0d0, 0d0, 0, 0, 1d-8, i, e, U, nbf, &
             isuppz, work, lwork, iwork, liwork, j)
 deallocate(isuppz, work, iwork, S)
 ! eigenvalues in array e are in ascending order

 forall(i = 1:nif, e(nbf-i+1)>0d0) noon(i) = e(nbf-i+1)
 deallocate(e)

 call dgemm('N', 'N', nbf, nif, nbf, 1d0, n_sqrt_S, nbf, U(:,nbf-nif+1:nbf), &
            nbf, 0d0, new_coeff, nbf)
 deallocate(n_sqrt_S, U)

 ! reverse the order of MOs
 allocate(U(nbf,nif))
 forall(i = 1:nif) U(:,i) = new_coeff(:,nif-i+1)
 new_coeff = U
 deallocate(U)
end subroutine gen_no_from_dm_and_ao_ovlp

! solve the system of linear equations of Cayley K matrix via Ax=b
subroutine solve_cayley_k_diis(nmo,ndiis,ijmap,binfile, k_diis, k_old, cayley_k)
 implicit none
 integer :: i, j, k, n, npair, nfile, fid
 integer, intent(in) :: nmo, ndiis
 integer, intent(in) :: ijmap(2,nmo*(nmo-1)/2)
 real(kind=8), intent(in) :: k_diis(ndiis+1,ndiis+1)
 real(kind=8), intent(inout) :: k_old(nmo,nmo), cayley_k(nmo,nmo)
 real(kind=8) :: r1, r2
 real(kind=8), allocatable :: x(:), bv(:), older(:)
 ! b has been used by B, so use bv instead
 character(len=240), intent(in) :: binfile(2*ndiis-1)

 npair = nmo*(nmo-1)/2; nfile = 2*ndiis - 1
 allocate(bv(ndiis+1), source=0d0)
 bv(ndiis+1) = 1d0
 allocate(x(ndiis+1))
 call solve_lin_eqs(ndiis+1, ndiis+1, k_diis, bv, x)
 deallocate(bv)

 k_old = cayley_k ! remember to update k_old
 cayley_k = x(ndiis)*cayley_k
 allocate(older(npair))

 do k = ndiis+1, nfile, 1
  open(newunit=fid,file=TRIM(binfile(k)),status='old',form='unformatted')
  read(unit=fid) older
  close(fid)
  r1 = x(nfile+1-k)
!$omp parallel do schedule(dynamic) default(shared) private(i,j,n,r2)
  do n = 1, npair, 1
   i = ijmap(1,n); j = ijmap(2,n)
   r2 = cayley_k(j,i) + r1*older(n)
   cayley_k(j,i) = r2; cayley_k(i,j) = -r2 ! anti-symm
  end do ! for n
!$omp end parallel do
 end do ! for k

 deallocate(older)
end subroutine solve_cayley_k_diis

subroutine get_mo_and_u_from_cayley_k_diis(nbf, nmo, ndiis, ijmap, binfile, &
                                       k_diis, u, mo, k_old, cayley_k, u_tmp)
 implicit none
 integer, intent(in) :: nbf, nmo, ndiis
 integer, intent(in) :: ijmap(2,nmo*(nmo-1)/2)
 real(kind=8), intent(in) :: k_diis(ndiis+1,ndiis+1)
 real(kind=8), intent(inout) :: u(nmo,nmo), mo(nbf,nmo), k_old(nmo,nmo), &
  cayley_k(nmo,nmo)
 real(kind=8), intent(out) :: u_tmp(nmo,nmo)
 real(kind=8), allocatable :: u_bak(:,:), mo_bak(:,:)
 character(len=240), intent(in) :: binfile(2*ndiis-1)

 call solve_cayley_k_diis(nmo, ndiis, ijmap, binfile, k_diis, k_old, cayley_k)
 allocate(u_bak(nmo,nmo), source=u)
 call cayley_trans(nmo, cayley_k, u)
 u_tmp = 0d0
 call dgemm('T','N', nmo, nmo, nmo, 1d0, u_bak, nmo, u, nmo, 0d0, u_tmp, nmo)
 deallocate(u_bak)
 allocate(mo_bak(nbf,nmo), source=mo)
 mo = 0d0
 call dgemm('N','N', nbf, nmo, nmo, 1d0, mo_bak, nbf, u_tmp, nmo, 0d0, mo, nbf)
 deallocate(mo_bak)
end subroutine get_mo_and_u_from_cayley_k_diis

!subroutine merge_two_sets_of_t1(nocc1,nvir1,t1_1, nocc2,nvir2,t1_2, t1)
! implicit none
! integer :: nocc, nvir
! integer, intent(in) :: nocc1,nvir1, nocc2,nvir2
!!f2py intent(in) :: nocc1,nvir1, nocc2,nvir2
! real(kind=8), intent(in) :: t1_1(nocc1,nvir1), t1_2(nocc2,nvir2)
!!f2py intent(in) :: t1_1, t1_2
!!f2py depend(nocc1,nvir1) :: t1_1
!!f2py depend(nocc2,nvir2) :: t1_2
! real(kind=8), intent(out) :: t1(nocc1+nocc2,nvir1+nvir2)
!!f2py intent(out) :: t1
!!f2py depend(nocc1,nocc2,nvir1,nvir2) :: t1
!
! t1 = 0d0
! nocc = nocc1 + nocc2
! nvir = nvir1 + nvir2
! t1(1:nocc1,nvir2+1:nvir) = t1_1
! t1(nocc1+1:nocc,1:nvir2) = t1_2
!end subroutine merge_two_sets_of_t1
!
!subroutine merge_two_sets_of_t2(nocc1,nvir1,t2_1, nocc2,nvir2,t2_2, t2)
! implicit none
! integer :: nocc, nvir
! integer, intent(in) :: nocc1,nvir1, nocc2,nvir2
!!f2py intent(in) :: nocc1,nvir1, nocc2,nvir2
! real(kind=8), intent(in) :: t2_1(nocc1,nocc1,nvir1,nvir1)
!!f2py intent(in) :: t2_1
!!f2py depend(nocc1,nvir1) :: t2_1
! real(kind=8), intent(in) :: t2_2(nocc2,nocc2,nvir2,nvir2)
!!f2py intent(in) :: t2_2
!!f2py depend(nocc2,nvir2) :: t2_2
! real(kind=8), intent(out) :: t2(nocc1+nocc2,nocc1+nocc2,nvir1+nvir2,nvir1+nvir2)
!!f2py intent(out) :: t2
!!f2py depend(nocc1,nocc2,nvir1,nvir2) :: t2
!
! t2 = 0d0
! nocc = nocc1 + nocc2
! nvir = nvir1 + nvir2
! t2(1:nocc1,1:nocc1,nvir2+1:nvir,nvir2+1:nvir) = t2_1
! t2(nocc1+1:nocc,nocc1+1:nocc,1:nvir2,1:nvir2) = t2_2
!end subroutine merge_two_sets_of_t2
!
!! Rotate old t1 and t2 amplitudes according to old MOs and new MOs.
!! Note: in this subroutine, the old MOs and new MOs share the same geometry and
!!  basis set, so they share one AO overlap integral matrix
!subroutine rotate_t1_t2_amp(nbf, nmo, mo0, mo1, nocc, t1, t2, ovlp)
! implicit none
! integer :: nvir
! integer, intent(in) :: nbf, nmo, nocc
!!f2py intent(in) :: nbf, nmo, nocc
! real(kind=8), intent(in) :: mo0(nbf,nmo), mo1(nbf,nmo)
!!f2py intent(in) :: mo0, mo1
!!f2py depend(nbf,nmo) :: mo0, mo1
! real(kind=8), intent(inout) :: t1(nocc,nmo-nocc), t2(nocc,nocc,nmo-nocc,nmo-nocc)
!!f2py intent(in,out) :: t1, t2
!!f2py depend(nocc,nmo) :: t1, t2
! real(kind=8), intent(in) :: ovlp(nbf,nbf)
!!f2py intent(in) :: ovlp
!!f2py depend(nbf) :: ovlp
! real(kind=8), allocatable :: uvt0(:,:),uvt1(:,:),s0(:),s1(:),CTSCp(:,:),r(:,:)
!
! nvir = nmo - nocc
!
! ! perform SVD on the overlap of occ MOs
! allocate(CTSCp(nocc,nocc))
! call calc_CTSCp(nbf, nocc, mo0(:,1:nocc), ovlp, mo1(:,1:nocc), CTSCp)
! allocate(uvt0(nocc,nocc), s0(nocc))
! call do_svd_get_uvt_s(nocc, CTSCp, uvt0, s0)
! write(6,'(A)') 's0='
! write(6,'(5(1X,ES15.8))') s0
! deallocate(CTSCp, s0)
!
! ! perform SVD on the overlap of vir MOs
! allocate(CTSCp(nvir,nvir))
! call calc_CTSCp(nbf, nvir, mo0(:,nocc+1:nmo), ovlp, mo1(:,nocc+1:nmo), CTSCp)
! allocate(uvt1(nvir,nvir), s1(nvir))
! call do_svd_get_uvt_s(nvir, CTSCp, uvt1, s1)
! write(6,'(A)') 's1='
! write(6,'(5(1X,ES15.8))') s1
! deallocate(CTSCp, s1)
!
! ! get new t1, O(N^3)
! allocate(r(nocc,nvir), source=0d0)
! call dgemm('N','N', nocc,nvir,nvir, 1d0,t1,nocc, uvt1,nvir, 0d0,r,nocc)
! t1 = 0d0
! call dgemm('T','N', nocc,nvir,nocc, 1d0,uvt0,nocc, r,nocc, 0d0,t1,nocc)
! deallocate(r)
!
! ! get new t2, O(N^5), like ao2mo
! call update_t2_using_p_occ_p_vir(nocc, nvir, uvt0, uvt1, t2)
! deallocate(uvt0, uvt1)
!end subroutine rotate_t1_t2_amp
!
!subroutine update_t2_using_p_occ_p_vir(nocc, nvir, p_occ, p_vir, t2)
! implicit none
! integer :: i, j, a, b
! integer, intent(in) :: nocc, nvir
! real(kind=8), intent(in) :: p_occ(nocc,nocc), p_vir(nvir,nvir)
! real(kind=8), intent(inout) :: t2(nocc,nocc,nvir,nvir)
! real(kind=8), allocatable :: r(:,:), s(:,:), r1(:,:,:,:), r2(:,:,:,:)
!
! allocate(r(nocc,nocc), s(nocc,nocc), r1(nocc,nocc,nvir,nvir))
! do b = 1, nvir, 1
!  do a = 1, nvir, 1
!   r = t2(:,:,a,b); s = 0d0
!   call dgemm('T','N', nocc,nocc,nocc, 1d0,p_occ,nocc, r,nocc, 0d0,s,nocc)
!  end do ! for a
!  r1(:,:,a,b) = TRANSPOSE(s)
! end do ! for b
!
! allocate(r2(nocc,nocc,nvir,nvir))
! do b = 1, nvir, 1
!  do a = 1, nvir, 1
!   r = r1(:,:,a,b); s = 0d0
!   call dgemm('T','N', nocc,nocc,nocc, 1d0,p_occ,nocc, r,nocc, 0d0,s,nocc)
!  end do ! for a
!  r2(a,b,:,:) = s
! end do ! for b
! deallocate(r, s, r1)
!
! allocate(r(nvir,nvir), s(nvir,nvir), r1(nvir,nvir,nocc,nocc))
! do i = 1, nocc, 1
!  do j = 1, nocc, 1
!   r = r2(:,:,j,i); s = 0d0
!   call dgemm('T','N', nvir,nvir,nvir, 1d0,p_vir,nvir, r,nvir, 0d0,s,nvir)
!  end do ! for j
!  r1(:,:,j,i) = TRANSPOSE(s)
! end do ! for i
! deallocate(r2)
!
! allocate(r2(nvir,nvir,nocc,nocc))
! do i = 1, nocc, 1
!  do j = 1, nocc, 1
!   r = r1(:,:,j,i); s = 0d0
!   call dgemm('T','N', nvir,nvir,nvir, 1d0,p_vir,nvir, r,nvir, 0d0,s,nvir)
!  end do ! for j
!  r2(:,:,j,i) = TRANSPOSE(s)
! end do ! for i
! deallocate(r, s, r1)
!
! forall(i=1:nocc,j=1:nocc,a=1:nvir,b=1:nvir) t2(i,j,a,b) = r2(a,b,j,i)
! deallocate(r2)
!end subroutine update_t2_using_p_occ_p_vir

