! written by jxzou at 20220214: move math library wrappers into this file

! diagonalize a real symmetric matrix and get all eigenvalues and eigenvectors.
! eigenvalues in w() are in ascending order w(1)<=w(2)<=...
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
  return
 end if
 w = 0d0

 ! ?syevd: Computes all eigenvalues and (optionally) all eigenvectors of a real
 ! symmetric matrix using divide and conquer algorithm.
 ! Syntax FORTRAN 77:
 ! call dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
 lwork = -1; liwork = -1
 allocate(work(1), iwork(1))
 call dsyevd('V', 'U', n, a, n, w, work, lwork, iwork, liwork, i)
 lwork = CEILING(work(1))
 liwork = iwork(1)
 deallocate(work, iwork)
 allocate(work(lwork), iwork(liwork))
 call dsyevd('V', 'U', n, a, n, w, work, lwork, iwork, liwork, i)

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine diag_get_e_and_vec: info/=0 in dsyevd.'
  write(6,'(2(A,I0))') 'n=', n, ', info=', i
  stop
 end if
 deallocate(work,iwork)
end subroutine diag_get_e_and_vec

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
 ! call dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z,
 ! ldz, isuppz, work, lwork, iwork, liwork, info)
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
  write(6,'(A)') 'ERROR in subroutine mat_dsqrt: diagonalization failed.'
  write(6,'(A,I0)') 'i=', i
  stop
 end if

 if(e(1) < -1d-6) then
  write(6,'(A)') 'ERROR in subroutine mat_dsqrt: too negative eigenvalue.'
  write(6,'(A,F16.9)') 'e(1)=', e(1)
  stop
 end if

 allocate(e1(n,n), source=0d0)
 allocate(Ue(n,n), source=0d0)
 sqrt_a = 0d0
 forall(i=1:n, e(i)>0d0) e1(i,i) = DSQRT(e(i))
 ! call dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 call dsymm('R', 'L', n, n, 1d0, e1, n, U, n, 0d0, Ue, n)
 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 call dgemm('N', 'T', n, n, n, 1d0, Ue, n, U, n, 0d0, sqrt_a, n)

 e1 = 0d0
 n_sqrt_a = 0d0
 forall(i=1:n, e(i)>=lin_dep) e1(i,i) = 1d0/DSQRT(e(i))
 ! call dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 call dsymm('R', 'L', n, n, 1d0, e1, n, U, n, 0d0, Ue, n)
 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
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

! solve the inverse matrix of a square matrix A (the user should make sure that
! A is reversible)
subroutine inverse(n, a, inv)
 implicit none
 integer :: i
 integer, intent(in) :: n
 real(kind=8), intent(in) :: a(n,n)
 real(kind=8), intent(out) :: inv(n,n)
 real(kind=8), allocatable :: identity(:,:)

 allocate(identity(n,n), source=0d0)
 forall(i = 1:n) identity(i,i) = 1d0
 call solve_multi_lin_eqs(n, n, a, n, identity, inv)
 deallocate(identity)
end subroutine inverse

! solving systems of linear equations with multiple right-hand sides
! Ax = b, x can be with multiple right-hand sides
subroutine solve_multi_lin_eqs(a1, a2, a, a3, b, x)
 implicit none
 integer :: i
 integer, intent(in) :: a1, a2, a3
!f2py intent(in) :: a1, a2, a3
 integer, allocatable :: ipiv(:)
 real(kind=8), intent(in) :: a(a1,a2), b(a1,a3)
!f2py intent(in) :: a, b
!f2py depend(a1,a2) :: a
!f2py depend(a1,a3) :: b
 real(kind=8), intent(out) :: x(a2,a3)
!f2py intent(out) :: x
!f2py depend(a2,a3) :: x
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

 x = b_copy(1:a2,1:a3)
 deallocate(b_copy)
end subroutine solve_multi_lin_eqs

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
 call dsymm('L', 'U', nbf, nif, 1d0, S, nbf, C, nbf, 0d0, SC, nbf)
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
 call dsymm('L', 'U', nbf, nif, 1d0, S, nbf, Cp, nbf, 0d0, SCp, nbf)
 call dgemm('T', 'N', nif, nif, nbf, 1d0, C, nbf, SCp, nbf, 0d0, CTSCp, nif)
 deallocate(SCp)
end subroutine calc_CTSCp

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

