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
 integer, parameter :: iout = 6
 integer, intent(in) :: n
 integer, allocatable :: iwork(:), isuppz(:)

 real(kind=8), parameter :: lin_dep = 1d-6
 ! 1.0D-6 is the default threshold of linear dependence in Gaussian and GAMESS
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
  write(iout,'(A)') 'ERROR in subroutine mat_dsqrt: diagonalization failed.'
  write(iout,'(A,I0)') 'i=', i
  stop
 end if

 if(e(1) < -1d-6) then
  write(iout,'(A)') 'ERROR in subroutine mat_dsqrt: too negative eigenvalue.'
  write(iout,'(A,F16.9)') 'e(1)=', e(1)
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

! solving systems of linear equations with multiple right-hand sides
! Ax = b, x can be with multiple right-hand sides
subroutine solve_multi_lin_eqs(a1, a2, a, a3, b, x)
 implicit none
 integer :: i
 integer, intent(in) :: a1, a2, a3
 integer, parameter :: iout = 6
 integer, allocatable :: ipiv(:)
 real(kind=8), intent(in) :: a(a1,a2), b(a1,a3)
 real(kind=8), intent(out) :: x(a2,a3)
 real(kind=8), allocatable :: a_copy(:,:), b_copy(:,:)

 ! mkl Syntax FORTRAN 77:
 ! ?getrf: Computes the LU factorization of a general m-by-n matrix
 ! call dgetrf(m, n, a, lda, ipiv, info)

 ! ?getrs: Solves a system of linear equations with an LU-factored square
 !  matrix, with multiple right-hand sides.
 ! call dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
 allocate(a_copy(a1,a2), source=a)
 allocate(b_copy(a1,a3), source=b)
 allocate(ipiv(min(a1,a2)), source=0)

 call dgetrf(a1, a2, a_copy, a1, ipiv, i)
 if(i /= 0) then
  write(iout,'(/,A)') 'ERROR in subroutine solve_multi_lin_eqs: MKL subroutine&
                     & dgetrf failed.'
  write(iout,'(A,I0)') 'info=', i
  stop
 end if

 call dgetrs('N', a2, a3, a_copy, a1, ipiv, b_copy, a1, i)
 deallocate(ipiv, a_copy)

 x = b_copy(1:a2,1:a3)
 deallocate(b_copy)
 return
end subroutine solve_multi_lin_eqs

