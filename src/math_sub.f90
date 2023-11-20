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

! calculate density matrix using MO coefficients and occupation numbers
subroutine calc_dm_using_mo_and_on(nbf, nif, mo, noon, dm)
 implicit none
 integer :: i, j, k
 integer, intent(in) :: nbf, nif
 real(kind=8), parameter :: thres = 1d-8
 real(kind=8), intent(in) :: mo(nbf,nif), noon(nif)
 real(kind=8), intent(out) :: dm(nbf,nbf)

 dm = 0d0 ! initialization

 do i = 1, nbf, 1
  do j = 1, i, 1
   do k = 1, nif, 1
    if(DABS(noon(k)) < thres) cycle
    dm(j,i) = dm(j,i) + noon(k)*mo(j,k)*mo(i,k)
   end do ! for k
  end do ! for j
 end do ! for i
 ! Note that dm(i,j) is not assigned
end subroutine calc_dm_using_mo_and_on

! get a random integer
subroutine get_a_random_int(i)
 implicit none
 integer :: n, clock
 integer, intent(out) :: i
 integer, allocatable :: seed(:)
 real(kind=4) :: r

 call random_seed(size=n)
 allocate(seed(n))
 call system_clock(count=clock)
 seed = clock
 call random_seed(put=seed)
 call random_number(r)
 deallocate(seed)

 i = CEILING(r*1e6)
end subroutine get_a_random_int

subroutine merge_two_sets_of_t1(nocc1,nvir1,t1_1, nocc2,nvir2,t1_2, t1)
 implicit none
 integer :: nocc, nvir
 integer, intent(in) :: nocc1,nvir1, nocc2,nvir2
!f2py intent(in) :: nocc1,nvir1, nocc2,nvir2
 real(kind=8), intent(in) :: t1_1(nocc1,nvir1), t1_2(nocc2,nvir2)
!f2py intent(in) :: t1_1, t1_2
!f2py depend(nocc1,nvir1) :: t1_1
!f2py depend(nocc2,nvir2) :: t1_2
 real(kind=8), intent(out) :: t1(nocc1+nocc2,nvir1+nvir2)
!f2py intent(out) :: t1
!f2py depend(nocc1,nocc2,nvir1,nvir2) :: t1

 t1 = 0d0
 nocc = nocc1 + nocc2
 nvir = nvir1 + nvir2
 t1(1:nocc1,nvir2+1:nvir) = t1_1
 t1(nocc1+1:nocc,1:nvir2) = t1_2
end subroutine merge_two_sets_of_t1

subroutine merge_two_sets_of_t2(nocc1,nvir1,t2_1, nocc2,nvir2,t2_2, t2)
 implicit none
 integer :: nocc, nvir
 integer, intent(in) :: nocc1,nvir1, nocc2,nvir2
!f2py intent(in) :: nocc1,nvir1, nocc2,nvir2
 real(kind=8), intent(in) :: t2_1(nocc1,nocc1,nvir1,nvir1)
!f2py intent(in) :: t2_1
!f2py depend(nocc1,nvir1) :: t2_1
 real(kind=8), intent(in) :: t2_2(nocc2,nocc2,nvir2,nvir2)
!f2py intent(in) :: t2_2
!f2py depend(nocc2,nvir2) :: t2_2
 real(kind=8), intent(out) :: t2(nocc1+nocc2,nocc1+nocc2,nvir1+nvir2,nvir1+nvir2)
!f2py intent(out) :: t2
!f2py depend(nocc1,nocc2,nvir1,nvir2) :: t2

 t2 = 0d0
 nocc = nocc1 + nocc2
 nvir = nvir1 + nvir2
 t2(1:nocc1,1:nocc1,nvir2+1:nvir,nvir2+1:nvir) = t2_1
 t2(nocc1+1:nocc,nocc1+1:nocc,1:nvir2,1:nvir2) = t2_2
end subroutine merge_two_sets_of_t2

! Rotate old t1 and t2 amplitudes accoording to old MOs and new MOs.
! Note: in this subroutine, the old MOs and new MOs share the same geometry and
!  basis set, so they share one AO overlap integral matrix
subroutine rotate_t1_t2_amp(nbf, nmo, mo0, mo1, nocc, t1, t2, ovlp)
 implicit none
 integer :: nvir
 integer, intent(in) :: nbf, nmo, nocc
!f2py intent(in) :: nbf, nmo, nocc
 real(kind=8), intent(in) :: mo0(nbf,nmo), mo1(nbf,nmo)
!f2py intent(in) :: mo0, mo1
!f2py depend(nbf,nmo) :: mo0, mo1
 real(kind=8), intent(inout) :: t1(nocc,nmo-nocc), t2(nocc,nocc,nmo-nocc,nmo-nocc)
!f2py intent(inout) :: t1, t2
!f2py depend(nocc,nmo) :: t1, t2
 real(kind=8), intent(in) :: ovlp(nbf,nbf)
!f2py intent(in) :: ovlp
!f2py depend(nbf) :: ovlp
 real(kind=8), allocatable :: uvt0(:,:),uvt1(:,:),s0(:),s1(:),CTSCp(:,:),r(:,:)

 nvir = nmo - nocc

 ! perform SVD on the overlap of occ MOs
 allocate(CTSCp(nocc,nocc))
 call calc_CTSCp(nbf, nocc, mo0(:,1:nocc), ovlp, mo1(:,1:nocc), CTSCp)
 allocate(uvt0(nocc,nocc), s0(nocc))
 call do_svd_get_uvt_s(nocc, CTSCp, uvt0, s0)
 write(6,'(A)') 's0='
 write(6,'(5(1X,ES15.8))') s0
 deallocate(CTSCp, s0)

 ! perform SVD on the overlap of vir MOs
 allocate(CTSCp(nvir,nvir))
 call calc_CTSCp(nbf, nvir, mo0(:,nocc+1:nmo), ovlp, mo1(:,nocc+1:nmo), CTSCp)
 allocate(uvt1(nvir,nvir), s1(nvir))
 call do_svd_get_uvt_s(nvir, CTSCp, uvt1, s1)
 write(6,'(A)') 's1='
 write(6,'(5(1X,ES15.8))') s1
 deallocate(CTSCp, s1)

 ! get new t1, O(N^3)
 allocate(r(nocc,nvir), source=0d0)
 call dgemm('N','N', nocc,nvir,nvir, 1d0,t1,nocc, uvt1,nvir, 0d0,r,nocc)
 t1 = 0d0
 call dgemm('T','N', nocc,nvir,nocc, 1d0,uvt0,nocc, r,nocc, 0d0,t1,nocc)
 deallocate(r)

 ! get new t2, O(N^5), like ao2mo
 call update_t2_using_p_occ_p_vir(nocc, nvir, uvt0, uvt1, t2)
 deallocate(uvt0, uvt1)
end subroutine rotate_t1_t2_amp

subroutine update_t2_using_p_occ_p_vir(nocc, nvir, p_occ, p_vir, t2)
 implicit none
 integer :: i, j, a, b
 integer, intent(in) :: nocc, nvir
 real(kind=8), intent(in) :: p_occ(nocc,nocc), p_vir(nvir,nvir)
 real(kind=8), intent(inout) :: t2(nocc,nocc,nvir,nvir)
 real(kind=8), allocatable :: r(:,:), s(:,:), r1(:,:,:,:), r2(:,:,:,:)

 allocate(r(nocc,nocc), s(nocc,nocc), r1(nocc,nocc,nvir,nvir))
 do b = 1, nvir, 1
  do a = 1, nvir, 1
   r = t2(:,:,a,b); s = 0d0
   call dgemm('T','N', nocc,nocc,nocc, 1d0,p_occ,nocc, r,nocc, 0d0,s,nocc)
  end do ! for a
  r1(:,:,a,b) = TRANSPOSE(s)
 end do ! for b

 allocate(r2(nocc,nocc,nvir,nvir))
 do b = 1, nvir, 1
  do a = 1, nvir, 1
   r = r1(:,:,a,b); s = 0d0
   call dgemm('T','N', nocc,nocc,nocc, 1d0,p_occ,nocc, r,nocc, 0d0,s,nocc)
  end do ! for a
  r2(a,b,:,:) = s
 end do ! for b
 deallocate(r, s, r1)

 allocate(r(nvir,nvir), s(nvir,nvir), r1(nvir,nvir,nocc,nocc))
 do i = 1, nocc, 1
  do j = 1, nocc, 1
   r = r2(:,:,j,i); s = 0d0
   call dgemm('T','N', nvir,nvir,nvir, 1d0,p_vir,nvir, r,nvir, 0d0,s,nvir)
  end do ! for j
  r1(:,:,j,i) = TRANSPOSE(s)
 end do ! for i
 deallocate(r2)

 allocate(r2(nvir,nvir,nocc,nocc))
 do i = 1, nocc, 1
  do j = 1, nocc, 1
   r = r1(:,:,j,i); s = 0d0
   call dgemm('T','N', nvir,nvir,nvir, 1d0,p_vir,nvir, r,nvir, 0d0,s,nvir)
  end do ! for j
  r2(:,:,j,i) = TRANSPOSE(s)
 end do ! for i
 deallocate(r, s, r1)

 forall(i=1:nocc,j=1:nocc,a=1:nvir,b=1:nvir) t2(i,j,a,b) = r2(a,b,j,i)
 deallocate(r2)
end subroutine update_t2_using_p_occ_p_vir

