! written by jxzou at 20220214: move math library wrappers into this file

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
    new_ev(i) = new_ev(j)
    new_ev(j) = tmp_ev
    tmp_ev = new_ev(i)
    tmp_mo = new_mo(:,i)
    new_mo(:,i) = new_mo(:,j)
    new_mo(:,j) = tmp_mo
   end if
  end do ! for j
 end do ! for i

 deallocate(tmp_mo)
end subroutine sort_mo_by_ev

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

 w = 0d0
 lwork = -1; liwork = -1
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
! 15315*15315.
subroutine diag_get_e_and_vec2(n, a, w)
 implicit none
 integer :: i, m, lwork, liwork
 integer, intent(in) :: n
 integer, allocatable :: iwork(:), isuppz(:)
 real(kind=8), intent(inout) :: a(n,n)
 real(kind=8), intent(out) :: w(n)
 real(kind=8), allocatable :: U(:,:), work(:)

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

 allocate(S(nbf,nbf), w(nbf))
 S = ovlp
 call diag_get_e_and_vec2(nbf, S, w)
 deallocate(S)

 nmo = 0
 do i = 1, nbf, 1
  if(w(i) > thres) exit
  nmo = nmo + 1
 end do ! for i
 deallocate(w)

 nmo = nbf - nmo
end subroutine get_nmo_from_ao_ovlp

! solve the A^1/2 and A^(-1/2) for a real symmetric matrix A
! Note: the input matrix A must be symmetric
subroutine mat_dsqrt(n, a0, sqrt_a, n_sqrt_a)
 implicit none
 integer :: i, m, lwork, liwork
 integer, intent(in) :: n
 integer, allocatable :: iwork(:), isuppz(:)
 real(kind=8), parameter :: lin_dep = 1d-6
 ! 1D-6 is the default threshold of linear dependence in Gaussian and GAMESS
 ! But in PySCF, one needs to manually adjust the threshold if linear dependence occurs
 real(kind=8), intent(in) :: a0(n,n)
 real(kind=8), intent(out) :: sqrt_a(n,n), n_sqrt_a(n,n)
 ! sqrt_a: A^1/2
 ! n_sqrt_a: A(-1/2)
 real(kind=8), allocatable :: a(:,:), e(:), U(:,:), Ue(:,:), work(:), e1(:,:)
 ! a: copy of a0
 ! e: eigenvalues; e1: all 0 except diagonal e(i)
 ! U: eigenvectors
 ! Ue: U*e

 allocate(a(n,n), source=a0)
 lwork = -1
 liwork = -1
 allocate(e(n), U(n,n), isuppz(2*n), work(1), iwork(1))
 call dsyevr('V', 'A', 'L', n, a, n, 0d0, 0d0, 0, 0, 1d-8, m, e, U, n, &
             isuppz, work, lwork, iwork, liwork, i)
 lwork = CEILING(work(1))
 liwork = iwork(1)
 deallocate(work, iwork)
 allocate(work(lwork), iwork(liwork))
 call dsyevr('V', 'A', 'L', n, a, n, 0d0, 0d0, 0, 0, 1d-8, m, e, U, n, &
             isuppz, work, lwork, iwork, liwork, i)

 deallocate(a, work, iwork, isuppz)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine mat_dsqrt: diagonalization failed.'
  write(6,'(A,I0)') 'i=', i
  stop
 end if

 if(e(1) < -1d-6) then
  write(6,'(/,A)') 'ERROR in subroutine mat_dsqrt: too negative eigenvalue.'
  write(6,'(A,F16.9)') 'e(1)=', e(1)
  stop
 end if

 allocate(e1(n,n), source=0d0)
 allocate(Ue(n,n), source=0d0)
 sqrt_a = 0d0
 forall(i=1:n, e(i)>0d0) e1(i,i) = DSQRT(e(i))
 call dsymm('R', 'L', n, n, 1d0, e1, n, U, n, 0d0, Ue, n)
 call dgemm('N', 'T', n, n, n, 1d0, Ue, n, U, n, 0d0, sqrt_a, n)

 e1 = 0d0
 n_sqrt_a = 0d0
 forall(i=1:n, e(i)>=lin_dep) e1(i,i) = 1d0/DSQRT(e(i))
 call dsymm('R', 'L', n, n, 1d0, e1, n, U, n, 0d0, Ue, n)
 call dgemm('N', 'T', n, n, n, 1d0, Ue, n, U, n, 0d0, n_sqrt_a, n)

 deallocate(e, e1, U, Ue)
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
  call calc_SPS(n, a, old_inv, inv_a)
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

! calculate (C^T)SC, S must be real symmetric since dsymm is called
! C: nbf*nif  S: nbf*nbf
subroutine calc_CTSC(nbf, nif, C, S, CTSC)
 implicit none
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: C(nbf,nif), S(nbf,nbf)
 real(kind=8), intent(out) :: CTSC(nif,nif)
 real(kind=8), allocatable :: SC(:,:)

 CTSC = 0d0
 allocate(SC(nbf,nif), source=0d0)
 call dsymm('L', 'L', nbf, nif, 1d0, S, nbf, C, nbf, 0d0, SC, nbf)
 call dgemm('T', 'N', nif, nif, nbf, 1d0, C, nbf, SC, nbf, 0d0, CTSC, nif)
 deallocate(SC)
end subroutine calc_CTSC

! calculate (C^T)S(C'), S must be real symmetric since dsymm is called
! C: nbf*nif  S: nbf*nbf, C': nbf*nif
subroutine calc_CTSCp(nbf, nif, C, S, Cp, CTSCp)
 implicit none
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: C(nbf,nif), S(nbf,nbf), Cp(nbf,nif)
 real(kind=8), intent(out) :: CTSCp(nif,nif)
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
subroutine calc_SPS(nbf, P, S, SPS)
 implicit none
 integer, intent(in) :: nbf
 real(kind=8), intent(in) :: P(nbf,nbf), S(nbf,nbf)
 real(kind=8), intent(out) :: SPS(nbf,nbf)
 real(kind=8), allocatable :: PS(:,:)

 SPS = 0d0
 allocate(PS(nbf,nbf), source=0d0)
 call dsymm('R', 'L', nbf, nbf, 1d0, S, nbf, P, nbf, 0d0, PS, nbf)
 call dsymm('L', 'L', nbf, nbf, 1d0, S, nbf, PS, nbf, 0d0, SPS, nbf)
 deallocate(PS)
end subroutine calc_SPS

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
  call calc_SPS(nbf, S, P, PSP)

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
 forall(u=1:nbf, v=1:nbf, v<u) dm(u,v) = dm(v,u)
end subroutine calc_dm_using_mo_and_on

! get a random integer
subroutine get_a_random_int(i)
 implicit none
 integer :: n, clock
 integer, intent(out) :: i
 integer, allocatable :: seed(:)
 real(kind=4) :: r

 call RANDOM_SEED(size=n)
 allocate(seed(n))
 call SYSTEM_CLOCK(count=clock)
 seed = clock
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
  call mat_dsqrt(nbf, S, U, X)
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

! solver AO-based overlap matrix (S) from condition (C^T)SC=I
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

! solver AO-based overlap matrix (S) by calculating (C(C^T))^(-1)
! This subroutine has the same result as subroutine solve_ovlp_from_ctsc.
subroutine solve_ovlp_from_cct(nbf, C, S)
 implicit none
 integer, intent(in) :: nbf
 real(kind=8), intent(in) :: C(nbf,nbf)
 real(kind=8), intent(out) :: S(nbf,nbf)
 real(kind=8), allocatable :: CCT(:,:)

 ! calculate C(C^T)
 allocate(CCT(nbf,nbf), source=0d0)
 call dgemm('N','T', nbf,nbf,nbf, 1d0,C,nbf, C,nbf, 0d0,CCT,nbf)

 ! calculate (C(C^T))^(-1)
 call inverse(nbf, CCT, S)
 deallocate(CCT)
end subroutine solve_ovlp_from_cct

! solver AO-based Fock matrix (F) from condition (C^T)FC=E
subroutine solve_fock_from_ctfc(nbf, nif, C, E, F)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: C(nbf,nif), E(nif)
 real(kind=8), intent(out) :: F(nbf,nbf)
 real(kind=8), allocatable :: FC(:,:), E1(:,:)

 allocate(E1(nif,nif), source=0d0)
 forall(i = 1:nif) E1(i,i) = E(i) ! diagonal matrix
 allocate(FC(nbf,nif))
 call solve_multi_lin_eqs(nif, nbf, TRANSPOSE(C), nif, E1, FC)
 deallocate(E1)
 ! FC = X -> (C^T)(F^T) = X^T, (F^T) = F
 call solve_multi_lin_eqs(nif, nbf, TRANSPOSE(C), nif, TRANSPOSE(FC), F)
end subroutine solve_fock_from_ctfc

! symmetrize a double precision matrix
subroutine symmetrize_dmat(n, a)
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 real(kind=8), intent(inout) :: a(n,n)

 if(n == 1) return
 forall(i=1:n-1, j=1:n, j>i) a(i,j) = a(j,i)
end subroutine symmetrize_dmat

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

! canonicalize MOs
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

! calculate the distance matrix from a set of coordinates
subroutine cal_dis_mat_from_coor(natom, coor, dis)
 implicit none
 integer :: i, j, k, m, n
 integer, intent(in) :: natom
 integer, allocatable :: map(:,:)
 real(kind=8) :: r(3)
 real(kind=8), intent(in) :: coor(3,natom)
 real(kind=8), intent(out) :: dis(natom,natom)

 n = natom
 forall(i = 1:n) dis(i,i) = 0d0
 k = n*(n-1)/2
 allocate(map(2,k))
 forall(i=1:n-1, j=1:n, j>i) map(:,(2*n-i)*(i-1)/2+j-i) = [i,j]

!$omp parallel do schedule(dynamic) default(private) shared(k,map,coor,dis)
 do m = 1, k, 1
  i = map(1,m); j = map(2,m)
  r = coor(:,i) - coor(:,j)
  dis(j,i) = DSQRT(DOT_PRODUCT(r,r))
  dis(i,j) = dis(j,i)
 end do ! for m
!$omp end parallel do

 deallocate(map)
end subroutine cal_dis_mat_from_coor

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
 integer, allocatable :: ipiv(:)
 real(kind=8), intent(in) :: coeff(nbf,nmo), lo_coeff(nbf,nmo)
 real(kind=8), intent(out) :: u(nmo,nmo)
 real(kind=8), allocatable :: coeff1(:,:), lo_coeff1(:,:)

 u = 0d0
 allocate(coeff1(nbf,nmo), lo_coeff1(nbf,nmo))
 coeff1 = coeff
 lo_coeff1 = lo_coeff
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
 lwork = work(1)
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
 lwork = work(1)
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
!! Rotate old t1 and t2 amplitudes accoording to old MOs and new MOs.
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
!!f2py intent(inout) :: t1, t2
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

