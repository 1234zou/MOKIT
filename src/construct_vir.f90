! written by jxzou at 20180406: construct part of or all virtual orbitals using the PAO (projected atomic orbitals)

! updated by jxzou at 20191215: simplify code

! This subroutine is designed to be imported as a module by Python.
!
! For GNU compiler, use
! --------------------------------------------------------------
!  f2py -m construct_vir -c construct_vir.f90 --link-lapack_opt
! --------------------------------------------------------------
!
! For Intel compiler, use
! -----------------------------------------------------------------------------------------------------
!  f2py -m construct_vir -c construct_vir.f90 --link-lapack_opt --fcompiler=intelem --compiler=intelem
! -----------------------------------------------------------------------------------------------------
!
!  to compile this file (a construct_vir.so file will be generated). Then in Python
!  you can import the construct_vir module.

subroutine construct_vir(nbf, nif, idx, coeff, ovlp, new_coeff)
 implicit none
 integer, parameter :: iout = 6
 integer :: i, j, nvir
 integer :: nbf, nif, idx
!f2py intent(in) :: nbf, nif, idx
 ! nbf: the number of basis functions
 ! nif: the number of MOs
 ! idx: the beginning index (Fortran convention) of the virtual MOs

 real(kind=8) :: coeff(nbf,nif), ovlp(nbf,nbf), new_coeff(nbf,nif)
!f2py depend(nbf,nif) coeff, new_coeff
!f2py depend(nbf) ovlp
!f2py intent(in) :: ovlp
!f2py intent(in,copy) :: coeff
!f2py intent(out) :: new_coeff

 real(kind=8), allocatable :: p(:,:), v(:,:), sv(:,:)
 real(kind=8), allocatable :: s1(:,:), ev(:), x(:,:)
 ! V: projected atomic orbitals (PAO)
 ! P: density matrix of atomic basis functions, sigma_i(Cui*Cvi)
 ! Note that the index i can be larger than nocc, in which case we only
 !  construct part of virtual orbitals.

 ! Step 1: P = sigma_i(Cui*Cvi)
 ! ?gemm: Computes a matrix-matrix product with general matrices
 ! Syntax FORTRAN 77:
 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 allocate(p(nbf,nbf), source=0.0d0)
 call dgemm('N', 'T', nbf, nbf, idx-1, 1.0d0, coeff(1:nbf,1:idx-1), nbf, coeff(1:nbf,1:idx-1), nbf, 0.0d0, p, nbf)

 ! Step 2: V = 1 - PS
 ! ?symm: Computes a matrix-matrix product where one input matrix is symmetric
 ! Syntax FORTRAN 77:
 ! call dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 allocate(v(nbf, nbf), source=0.0d0)
 forall(i = 1:nbf) v(i,i) = 1.0d0
 call dsymm('R', 'U', nbf, nbf, -1.0d0, ovlp, nbf, p, nbf, 1.0d0, v, nbf)
 deallocate(p)

 ! Step 3: S1 = (VT)SV
 allocate(sv(nbf,nbf), source=0.0d0)
 allocate(s1(nbf,nbf), source=0.0d0)
 call dsymm('L', 'U', nbf, nbf, 1.0d0, ovlp, nbf, v, nbf, 0.0d0, sv, nbf)
 call dgemm('T', 'N', nbf, nbf, nbf, 1.0d0, v, nbf, sv, nbf, 0.0d0, s1, nbf)
 deallocate(sv)

 ! Step 4: diagonalize S1 (note that S1 is symmetric) and get X
 allocate(ev(nbf))
 call diag_get_e_and_vec(nbf, s1, ev)
 nvir = nif - idx + 1
 forall(i = nbf-nvir+1:nbf) ev(i) = 1.0d0/DSQRT(ev(i))
 allocate(x(nbf,nvir), source=0.0d0)
 forall(i=1:nbf, j=1:nvir) x(i,j) = s1(i,nbf-nvir+j)*ev(nbf-nvir+j)
 deallocate(ev, s1)

 ! Step 5: get new virtual MO coefficients
 call dgemm('N', 'N', nbf, nvir, nbf, 1.0d0, v, nbf, x, nbf, 0.0d0, coeff(:,idx:nif), nbf)
 deallocate(x, v)
 new_coeff = coeff

 ! Step 6: check orthonormality
 allocate(s1(nbf,nif), source=0.0d0)
 allocate(x(nif,nif), source=0.0d0)
 call dsymm('L', 'U', nbf, nif, 1.0d0, ovlp, nbf, coeff, nbf, 0.0d0, s1, nbf)
 call dgemm('T', 'N', nif, nif, nbf, 1.0d0, coeff, nbf, s1, nbf, 0.0d0, x, nif)
 deallocate(s1)
 forall(i = 1:nif) x(i,i) = x(i,i) - 1.0d0
 x = DABS(x)
 write(iout,'(/,A)') 'The orthonormality of Alpha MO after PAO construction:'
 write(iout,'(A,F16.10)') 'maxv=', MAXVAL(x)
 write(iout,'(A,F16.10)') 'abs_mean=', SUM(x)/DBLE(nif*nif)
 deallocate(x)

 return
end subroutine construct_vir

! subroutine diag_get_e_and_vec: diagonalize a real symmetric matrix and get all eigenvalues and eigenvectors
subroutine diag_get_e_and_vec(n, a, w)
 implicit none
 integer :: info, lwork, liwork
 integer, intent(in) :: n
 integer, allocatable :: iwork(:)
 real(kind=8), intent(inout) :: a(n,n)
 real(kind=8), intent(out) :: w(n)
 real(kind=8), allocatable :: work(:)

 ! ?syevd: Computes all eigenvalues and (optionally) all eigenvectors of a real
 ! symmetric matrix using divide and conquer algorithm.
 ! Syntax FORTRAN 77:
 ! call dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
 lwork = -1
 liwork = -1
 allocate(work(1), iwork(1))
 call dsyevd('V', 'U', n, a, n, w, work, lwork, iwork, liwork, info)
 lwork = CEILING(work(1))
 liwork = iwork(1)
 deallocate(work, iwork)
 allocate(work(lwork), iwork(liwork))
 call dsyevd('V', 'U', n, a, n, w, work, lwork, iwork, liwork, info)
 deallocate(work,iwork)
 return
end subroutine diag_get_e_and_vec

