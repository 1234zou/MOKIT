! written by jxzou at 20180410
! modified by jxzou at 20190910:
!  1) optimize the code
!  2) add subroutine mo_svd_qcmo, diagonalize SVOs into QCMOs (quasi-canonical MOs)

!  For GNU compiler, use
! ------------------------------------------------
!  f2py -m mo_svd -c mo_svd.f90 --link-lapack_opt
! ------------------------------------------------
!  For INTEL compiler, use
! ---------------------------------------------------------------------------------------
!  f2py -m mo_svd -c mo_svd.f90 --link-lapack_opt --fcompiler=intelem --compiler=intelem
! ---------------------------------------------------------------------------------------

module mo_ovlp_and_svd
 implicit none
 integer, parameter :: iout = 6
contains

 ! Compute the overlap of two sets of MOs
 ! Note that these two sets MOs can have different number of basis functions and number of MOs
 subroutine get_mo_basis_ovlp2(nbf1, nmo1, nbf2, nmo2, coeff1, coeff2, ao_ovlp, mo_ovlp)
  implicit none
  integer, intent(in) :: nbf1, nmo1, nbf2, nmo2
  real(kind=8), intent(in) :: coeff1(nbf1,nmo1), coeff2(nbf2,nmo2)
  real(kind=8), intent(in) :: ao_ovlp(nbf1,nbf2)
  real(kind=8), intent(out) :: mo_ovlp(nmo1,nmo2)
  real(kind=8), allocatable :: temp(:,:)

  mo_ovlp = 0d0
  allocate(temp(nbf1,nmo2), source=0d0)
  ! ?gemm: Computes a matrix-matrix product with general matrices
  ! Syntax FORTRAN 77:
  ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  call dgemm('N', 'N', nbf1, nmo2, nbf2, 1d0, ao_ovlp, nbf1, coeff2, nbf2, 0d0, temp, nbf1)
  call dgemm('T', 'N', nmo1, nmo2, nbf1, 1d0, coeff1, nbf1, temp, nbf1, 0d0, mo_ovlp, nmo1)

  deallocate(temp)
  return
 end subroutine get_mo_basis_ovlp2

 ! Perform SVD on two sets of MOs and rotate them using obtained unitary matrice U and V_T
 ! Note that these two sets MOs can have different number of basis functions and number of MOs
 subroutine svd_and_rotate2(nbf1, nmo1, nbf2, nmo2, coeff1, coeff2, mo_ovlp, sv, reverse, mo_e)
  implicit none
  integer :: i
  integer, intent(in) :: nbf1, nmo1, nbf2, nmo2
  real(kind=8), intent(inout) :: coeff1(nbf1,nmo1), coeff2(nbf2,nmo2)
  real(kind=8), intent(inout), optional :: mo_e(nmo1)
  real(kind=8), intent(in) :: mo_ovlp(nmo1,nmo2)
  real(kind=8), intent(out) :: sv(nmo1)
  real(kind=8), allocatable :: sv2(:)
  real(kind=8), allocatable :: u(:,:), vt(:,:), u2(:,:), vt2(:,:)
  logical, intent(in) :: reverse

  sv = 0d0
  allocate(u(nmo1,nmo1), source=0d0)
  allocate(vt(nmo2,nmo2), source=0d0)

  ! perform SVD on mo_ovlp
  call svd_on_ovlp(nmo1, nmo2, mo_ovlp, u, vt, sv)

  if(reverse) then
   allocate(sv2(nmo1), source=0d0)
   forall(i = 1:nmo1)
    sv2(i) = sv(nmo1-i+1)
   end forall
   sv = sv2
   deallocate(sv2)

   allocate(u2(nmo1,nmo1), source=0d0)
   allocate(vt2(nmo2,nmo2), source=0d0)
   forall(i = 1:nmo1)
    u2(:,i) = u(:,nmo1-i+1)
   end forall
   u = u2

   forall(i = 1:nmo2)
    vt2(i,:) = vt(nmo2-i+1,:)
   end forall
   vt = vt2

   deallocate(u2, vt2)
  end if

  ! rotate the original orbitals
  allocate(u2(nbf1,nmo1), source=0d0)
  allocate(vt2(nbf2,nmo2), source=0d0)
  ! ?gemm: Computes a matrix-matrix product with general matrices
  ! Syntax FORTRAN 77:
  ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  call dgemm('N', 'N', nbf1, nmo1, nmo1, 1d0, coeff1, nbf1, u, nmo1, 0d0, u2, nbf1)
  call dgemm('N', 'T', nbf2, nmo2, nmo2, 1d0, coeff2, nbf2, vt, nmo2, 0d0, vt2, nbf2)
  coeff1 = u2
  coeff2 = vt2
  deallocate(u2, vt, vt2)

  if(present(mo_e)) then
   ! e' = (U^T)eU
   allocate(vt(nmo1,nmo1), source=0d0)
   allocate(u2(nmo1,nmo1), source=0d0)
   forall(i = 1:nmo1) vt(i,i) = mo_e(i)
   allocate(vt2(nmo1,nmo1), source=0d0)
   call dgemm('T', 'N', nmo1, nmo1, nmo1, 1d0, u, nmo1, vt, nmo1, 0d0, vt2, nmo1)
   call dgemm('N', 'N', nmo1, nmo1, nmo1, 1d0, vt2, nmo1, u, nmo1, 0d0, vt, nmo1)
   u2 = vt
   deallocate(vt2)

   ! diagonalize (1:nmo2,1:nmo2) part of symmetric matrix e', obtaining a new U matrix
   allocate(vt2(nmo2,nmo2), source=0d0)
   vt2 = vt(1:nmo2,1:nmo2)
   deallocate(vt)
   allocate(sv2(nmo2), source=0d0)
   call diag_get_e_and_vec(nmo2, vt2, sv2)
   mo_e = 0d0
   mo_e(1:nmo2) = sv2
   deallocate(sv2)
   ! Note that orthonormal eigenvectors are now stored in matrix vt2

   ! C' = CU done in (1:nmo2,1:nmo2) block
   allocate(vt(nbf1,nmo2))
   call dgemm('N', 'N', nbf1, nmo2, nmo2, 1d0, coeff1(:,1:nmo2), nbf1, vt2, nmo2, 0d0, vt, nbf1)
   coeff1(:,1:nmo2) = vt
   deallocate(vt, vt2)

   ! diagonalize the left of symmetric matrix e', obtaining a new U matrix
   i = nmo1 - nmo2
   allocate(vt2(i,i), source=0d0)
   vt2 = u2(nmo2+1:nmo1,nmo2+1:nmo1)
   deallocate(u2)
   allocate(sv2(i), source=0d0)
   call diag_get_e_and_vec(i, vt2, sv2)
   mo_e(nmo2+1:nmo1) = sv2
   deallocate(sv2)
   ! Note that orthonormal eigenvectors are now stored in matrix vt2

   ! C' = CU done in (nmo2+1:nmo1,nmo2+1:nmo1) block
   allocate(vt(nbf1,i))
   call dgemm('N', 'N', nbf1, i, i, 1d0, coeff1(:,nmo2+1:nmo1), nbf1, vt2, i, 0d0, vt, nbf1)
   coeff1(:,nmo2+1:nmo1) = vt
   deallocate(vt, vt2)

!   allocate(vt2(nmo1,nmo1), source=0d0)
!   vt2 = vt(1:nmo1,1:nmo1)
!   deallocate(vt)
!   allocate(sv2(nmo1), source=0d0)
!   call diag_get_e_and_vec(nmo1, vt2, sv2)
!   mo_e = 0d0
!   mo_e(1:nmo1) = sv2
!   deallocate(sv2)
!   ! Note that orthonormal eigenvectors are now stored in matrix vt2
!
!   ! C' = CU done in (1:nmo2,1:nmo2) block
!   allocate(vt(nbf1,nmo1))
!   call dgemm('N', 'N', nbf1, nmo1, nmo1, 1d0, coeff1(:,1:nmo1), nbf1, vt2, nmo1, 0d0, vt, nbf1)
!   coeff1(:,1:nmo1) = vt
!   deallocate(vt2)
  end if

  deallocate(u)
  return
 end subroutine svd_and_rotate2

 ! perform SVD on overlap matrix s
 subroutine svd_on_ovlp(m, n, a, u, vt, s)
  implicit none
  integer :: lwork, info
  integer, intent(in) :: m, n
  integer, parameter :: iout = 6
  real(kind=8), intent(in) :: a(m,n)
  real(kind=8), intent(out) :: u(m,m), vt(n,n), s(m)
  real(kind=8), allocatable :: work(:), a_copy(:,:)
 
  info = 0
  u = 0d0; vt = 0d0; s = 0d0
  allocate(a_copy(m,n), source=a)
 
  ! ?gesvd: Computes the singular value decomposition of a general rectangular matrix
  ! Syntax FORTRAN 77:
  ! call dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
  lwork = -1
  allocate(work(1), source=0d0)
  call dgesvd('A', 'A', m, n, a_copy, m, s, u, m, vt, n, work, lwork, info)
 
  lwork = CEILING(work(1))
  deallocate(work)
  allocate(work(lwork),source=0d0)
  call dgesvd('A', 'A', m, n, a_copy, m, s, u, m, vt, n, work, lwork, info)
 
  if(info /= 0) then
   write(iout,'(A)') 'ERROR in subroutine svd_on_ovlp: info/=0! Please check why.'
   write(iout,'(A5,I0)') 'info=', info
   stop
  end if
 
  deallocate(work, a_copy)
  return
 end subroutine svd_on_ovlp

end module mo_ovlp_and_svd

! perform SVD on two sets of MOs, and get new MOs
subroutine mo_svd(nbf1, nmo1, nbf2, nmo2, coeff1, coeff2, ao_ovlp, reverse)
 use mo_ovlp_and_svd
 implicit none
 integer :: i
 integer :: nbf1, nmo1, nbf2, nmo2
!f2py intent(in) :: nbf1, nmo1, nbf2, nmo2

 real(kind=8) :: coeff1(nbf1,nmo1), coeff2(nbf2,nmo2)
!f2py intent(inout) :: coeff1, coeff2
!f2py depend(nbf1,nmo1) :: coeff1
!f2py depend(nbf2,nmo2) :: coeff2

 real(kind=8) :: ao_ovlp(nbf1,nbf2)
!f2py intent(in) :: ao_ovlp
!f2py depend(nbf1,nbf2) :: ao_ovlp

 logical :: reverse
!f2py intent(in) :: reverse

 real(kind=8), allocatable :: mo_ovlp(:,:), sv(:)

 ! compute MO basis overlap
 allocate(mo_ovlp(nmo1,nmo2), source=0d0)
 call get_mo_basis_ovlp2(nbf1, nmo1, nbf2, nmo2, coeff1, coeff2, ao_ovlp, mo_ovlp)

 ! perform SVD and get new MO
 allocate(sv(nmo1), source=0d0)
 call svd_and_rotate2(nbf1, nmo1, nbf2, nmo2, coeff1, coeff2, mo_ovlp, sv, reverse)

 write(iout,'(/,A)') 'Singular values:'
 write(iout,'(5(1X,ES15.8))') (sv(i),i=1,nmo1)
 deallocate(mo_ovlp, sv)
 return
end subroutine mo_svd

! 1) project a set of MOs of a small basis set onto the (partially) occupied MOs of
!    a large basis set, leading to virtual orbitals of the small basis set
! 2) project virtual MOs of the large basis set onto those orbitals of the small
!    basis set
subroutine proj_occ_get_act_vir(nbf1, nmo1, nbf2, na_np, S2, cross_S, coeff1, coeff)
 use mo_ovlp_and_svd, only: iout, get_mo_basis_ovlp2, svd_and_rotate2
 implicit none
 integer :: i, j, nmo2, nvir1, nvir2
 integer :: nbf1, nmo1, nbf2, na_np
!f2py intent(in) nbf1, nmo1, nbf2, na_np
 real(kind=8) :: S2(nbf2,nbf2), cross_S(nbf1,nbf2)
!f2py intent(in,copy) S2
!f2py intent(in) cross_S
!f2py depend(nbf2) S2
!f2py depend(nbf1,nbf2) cross_S
 real(kind=8) :: coeff1(nbf1,nmo1), coeff(nbf1,nmo1)
!f2py intent(in,copy) coeff1
!f2py intent(out) :: coeff
!f2py depend(nbf1,nmo1) coeff, coeff1
 real(kind=8), allocatable :: coeff2(:,:), mo_ovlp(:,:), sv(:), sv2(:,:)
 
! na_np: the number of alpha orbitals + unoccupied localized UNOs in coeff1
! coeff1: MOs of large basis set
! coeff2: MOs of small basis set
 coeff = 0d0

 ! U^(T)SU = s, X = Us^(-1/2), a orthonormal set of MOs of small basis set
 allocate(sv(nbf2))
 call diag_get_e_and_vec(nbf2, S2, sv)
 nmo2 = COUNT(sv > 1d-6)
 if(na_np >= nmo2) then
  write(iout,'(/,A)') 'ERROR in subroutine proj_occ_get_act_vir: na_np>=nmo2.'
  write(iout,'(A,5I5)') 'nbf1,nmo1,nbf2,nmo2,na_np=',nbf1,nmo1,nbf2,nmo2,na_np
  stop
 end if
 if(nmo2 < nbf2) write(iout,'(A)') 'Warning from subroutine proj_occ_get_act_&
                                   &vir: nmo2 < nbf2, rare case.'
 allocate(coeff2(nbf2,nmo2), sv2(nmo2,nmo2))
 sv2 = 0d0
 j = nbf2 - nmo2
 forall(i = 1:nmo2) sv2(i,i) = 1d0/DSQRT(sv(j+i))
 call dgemm('N', 'N', nbf2, nmo2, nmo2, 1d0, S2(j+1:nbf2,j+1:nbf2), nbf2, sv2,&
            nmo2, 0d0, coeff2, nbf2)
 deallocate(sv, sv2)

 allocate(mo_ovlp(nmo2,na_np))
 call get_mo_basis_ovlp2(nbf2, nmo2, nbf1, na_np, coeff2, coeff1(:,1:na_np), &
                         TRANSPOSE(cross_S), mo_ovlp)

 coeff(:,1:na_np) = coeff1(:,1:na_np)
 allocate(sv(nmo2))
 call svd_and_rotate2(nbf2,nmo2, nbf1,na_np, coeff2,coeff,mo_ovlp, sv, .false.)
 write(iout,'(/,A)') 'Singular values of projecting small basis MOs onto large&
                    & basis (partially) occupied MOs:'
 write(iout,'(5(1X,ES15.8))') (sv(i),i=1,nmo2)
 deallocate(sv, mo_ovlp)

 nvir1 = nmo1 - na_np
 nvir2 = nmo2 - na_np
 allocate(mo_ovlp(nvir1,nvir2))
 call get_mo_basis_ovlp2(nbf1, nvir1, nbf2, nvir2, coeff1(:,na_np+1:nmo1), &
                         coeff2(:,na_np+1:nmo2), cross_S, mo_ovlp)
 allocate(sv(nvir1))
 call svd_and_rotate2(nbf1, nvir1, nbf2, nvir2, coeff1(:,na_np+1:nmo1), &
                      coeff2(:,na_np+1:nmo2), mo_ovlp, sv, .false.)
 write(iout,'(/,A)') 'Singular values of projecting large basis set virtual MOs&
                    & onto small basis virtual MOs:'
 write(iout,'(5(1X,ES15.8))') (sv(i),i=1,nvir1)
 deallocate(sv, mo_ovlp)

 coeff = coeff1
 return
end subroutine proj_occ_get_act_vir

! perform SVD on two sets of MOs, get new MOs, and diagonalize part of MOs to
!  obtain quasi-canonical MOs (with corresponding energies)
subroutine mo_svd_qcmo(nbf1, nmo1, nbf2, nmo2, coeff1, coeff2, ao_ovlp, mo_e)
 use mo_ovlp_and_svd
 implicit none
 integer :: i
 integer :: nbf1, nmo1, nbf2, nmo2
!f2py intent(in) :: nbf1, nmo1, nbf2, nmo2

 real(kind=8) :: coeff1(nbf1,nmo1), coeff2(nbf2,nmo2), mo_e(nmo1)
!f2py intent(inout) :: coeff1, coeff2, mo_e
!f2py depend(nbf1,nmo1) :: coeff1
!f2py depend(nbf2,nmo2) :: coeff2
!f2py depend(nmo1) :: mo_e

 real(kind=8) :: ao_ovlp(nbf1,nbf2)
!f2py intent(in) :: ao_ovlp
!f2py depend(nbf1,nbf2) :: ao_ovlp

 real(kind=8), allocatable :: mo_ovlp(:,:), sv(:)

 ! compute MO basis overlap
 allocate(mo_ovlp(nmo1,nmo2), source=0d0)
 call get_mo_basis_ovlp2(nbf1, nmo1, nbf2, nmo2, coeff1, coeff2, ao_ovlp, mo_ovlp)

 ! perform SVD and get new MO
 allocate(sv(nmo1), source=0d0)
 call svd_and_rotate2(nbf1, nmo1, nbf2, nmo2, coeff1, coeff2, mo_ovlp, sv, .false., mo_e)

 write(iout,'(/,A)') 'Singular values:'
 write(iout,'(5(1X,ES15.8))') (sv(i),i=1,nmo1)

 deallocate(mo_ovlp, sv)
 return
end subroutine mo_svd_qcmo

