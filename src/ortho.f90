! written by jxzou at 20191029
! updated by jxzou at 20200403: move check_orthonormal into this file

! Currently this file contains
! subroutine check_orthonormal: check whether a given set of MOs are orthonormal
! subroutine can_ortho: generate canonical orthonormalized atomic orbitals
! subroutine sym_ortho: generate symmetric orthonormalized atomic orbitals

! TO DO: Schmidt orthogonalization

! check whether a given set of MOs are orthonormal
subroutine check_orthonormal(nbf, nif, coeff, S)
 implicit none
 integer :: i, j, i0, j0
 integer :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8) :: coeff(nbf,nif), S(nbf,nbf)
!f2py intent(in) :: coeff, S
!f2py depend(nbf,nif) :: coeff
!f2py depend(nbf) :: S
 real(kind=8) :: maxv
 real(kind=8), allocatable :: C_T_S_C(:,:)

 allocate(C_T_S_C(nif,nif))
 call calc_CTSC(nbf, nif, coeff, S, C_T_S_C)

 forall(i=1:nif) C_T_S_C(i,i) = C_T_S_C(i,i) - 1d0
 forall(i=1:nif,j=1:nif) C_T_S_C(i,j) = DABS(C_T_S_C(i,j))

 i0 = 1; j0 = 1
 maxv = C_T_S_C(1,1)

 do i = 1, nif, 1
  do j = 1, nif, 1
   if(C_T_S_C(j,i) > maxv) then
    maxv = C_T_S_C(j,i)
    i0 = i; j0 = j
   end if
  end do ! for j
 end do ! for i

 write(6,'(/,2(A,I4,1X),A5,ES15.8)') 'Orthonormality check: j=',j0,'i=',i0,'maxv=',maxv
 deallocate(C_T_S_C)
end subroutine check_orthonormal

! check whether a given set of complex MOs are orthonormal
! Note: 
! 1) this subroutine is specially designed for complex GHF MOs in PySCF, where 
!  the basis oerder of a MO is (ar, ai)1, (ar, ai)2, ..., (br, bi)1, (br, bi)2
! 2) do not use this subroutine for complex GHF MOs in Gaussian, where the
!  basis oerder of a MO is (ar, ai, br, bi)1, (ar, ai, br, bi)2, ...
subroutine check_cghf_orthonormal(nbf, nif, coeff, S)
 implicit none
 integer :: i, j, i0, j0
 integer :: nbf, nif
!f2py intent(in) :: nbf, nif
 complex(kind=8), intent(in) :: coeff(nbf,nif)
!f2py intent(in) :: coeff
!f2py depend(nbf,nif) :: coeff
 real(kind=8), intent(in) :: S(nbf/2,nbf/2)
!f2py intent(in) :: S
!f2py depend(nbf) :: S
 real(kind=8) :: rtmp, maxv
 real(kind=8), allocatable :: C_r(:,:), C_i(:,:), A(:,:), A_b(:,:)
 real(kind=8), allocatable :: B(:,:), B_2(:,:)

 allocate(C_r(nbf,nif), source=REAL(coeff))  ! real part
 allocate(C_i(nbf,nif), source=AIMAG(coeff)) ! imaginary part

 allocate(A(nif,nif))
 call calc_CTSCp(nbf/2, nif, C_r(1:nbf/2,:), S, C_i(1:nbf/2,:), A)

 allocate(A_b(nif,nif))
 call calc_CTSCp(nbf/2, nif, C_r(nbf/2+1:nbf,:), S, C_i(nbf/2+1:nbf,:), A_b)

 A = A + A_b
 deallocate(A_b)
 A = A - TRANSPOSE(A)

 rtmp = SUM(DABS(A))/DBLE(nif*nif)
 deallocate(A)
 write(6,'(/,A,ES15.8)') 'SUM(DABS(A))/nif^2=', rtmp
 if(rtmp > 1d-7) then
  write(6,'(/,A)') 'Warning in subroutine check_cghf_orthonormal: rtmp>1D-7.'
  write(6,'(A)') 'SUM(DABS(A))/nif^2 should be zero. Orthonormality not good.'
 end if

 allocate(B(nif,nif))
 call calc_CTSC(nbf/2, nif, C_r(1:nbf/2,:), S, B)

 allocate(B_2(nif,nif))
 call calc_CTSC(nbf/2, nif, C_r(nbf/2+1:nbf,:), S, B_2)
 deallocate(C_r)
 B = B + B_2

 call calc_CTSC(nbf/2, nif, C_i(1:nbf/2,:), S, B_2)
 B = B + B_2

 call calc_CTSC(nbf/2, nif, C_i(nbf/2+1:nbf,:), S, B_2)
 deallocate(C_i)
 B = B + B_2
 deallocate(B_2)

 forall(i = 1:nif) B(i,i) = B(i,i)- 1d0
 B = DABS(B)

 do i = 1, nif, 1
  do j = 1, nif, 1
   if(B(j,i) > maxv) then
    maxv = B(j,i)
    i0 = i; j0 = j
   end if
  end do ! for j
 end do ! for i

 write(6,'(/,2(A,I4,1X),A5,ES15.8)') 'Orthonormality check: j=',j0,'i=',i0,'maxv=',maxv
 deallocate(B)
end subroutine check_cghf_orthonormal

! orthonormalize a set of non-orthogonal MOs
! Note: this subroutine can only be used for no linear dependence
subroutine orthonormalize_orb(nbf, nif, ao_ovlp, old_mo, new_mo)
 implicit none
 integer :: i, j
 integer :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8) :: ao_ovlp(nbf,nbf), old_mo(nbf,nif), new_mo(nbf,nif)
!f2py depend(nbf) :: ao_ovlp
!f2py intent(in) :: ao_ovlp
!f2py depend(nbf,nif) :: mo_coeff, new_mo
!f2py intent(in) :: mo_coeff
!f2py intent(out) :: new_mo
 real(kind=8), allocatable :: X(:,:), ev(:)
 real(kind=8), allocatable :: Sp(:,:) ! S', S prime

 new_mo = 0d0
 allocate(Sp(nif,nif))
 call calc_CTSC(nbf, nif, old_mo, ao_ovlp, Sp) !(C^T)SC = S'
 ! S' is not I, so C is not orthonormal

 allocate(ev(nif), source=0d0)
 call diag_get_e_and_vec(nif, Sp, ev) ! S' = U(s')(U^T), U stored in Sp

 if(ANY(ev < 0d0)) then
  write(6,'(A)') 'ERROR in subroutine orthonormalize_orb: some eigenvalues are&
                & negative. ev='
  write(6,'(5(1X,ES15.8))') ev
  stop
 end if

 forall(i = 1:nif) ev(i) = 1d0/DSQRT(ev(i))
 allocate(X(nif,nif))   ! X = U((s')^(-1/2))
 forall(i=1:nif, j=1:nif) X(i,j) = Sp(i,j)*ev(j)
 deallocate(Sp, ev)

 ! C' = CX
 call dgemm('N','N', nbf,nif,nif, 1d0,old_mo,nbf, X,nif, 0d0,new_mo,nbf)
 deallocate(X)
end subroutine orthonormalize_orb

! generate canonical orthonormalized atomic orbitals
subroutine can_ortho(nbf, nif, ao_ovlp, mo_coeff)
 implicit none
 integer :: i, nif0
 integer :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), parameter :: thresh = 1d-6
 real(kind=8) :: ao_ovlp(nbf,nbf), mo_coeff(nbf,nif)
!f2py depend(nbf) :: ao_ovlp
!f2py intent(in,copy) :: ao_ovlp
!f2py depend(nbf,nif) :: mo_coeff
!f2py intent(out) :: mo_coeff
 real(kind=8), allocatable :: U(:,:), s(:)

 mo_coeff = 0d0
 allocate(s(nbf), source=0d0)
 call diag_get_e_and_vec(nbf, ao_ovlp, s) ! S = UsU^T, U stored in ao_ovlp

 nif0 = COUNT(s > thresh)
 if(nif0 /= nif) then
  write(6,'(A)') 'ERROR in subroutine can_ortho: nif /= nif0. This is because&
                   & the linear dependence threshold here and outside subroutine&
                   & is not consistent.'
  write(6,'(A,ES15.8)') 'Default threshold here:', thresh
  stop
 end if

 ! reverse the eigenvalues
 allocate(U(nbf,1), source=0d0)
 forall(i = 1:nbf) U(i,1) = s(nbf-i+1)
 s = U(:,1)
 deallocate(U)
 ! now s1 > s2 > s3 > ...

 ! reverse the eigenvectors, according to the descending order of eigenvalues
 allocate(U(nbf, nbf), source=0d0)
 forall(i = 1:nbf) U(:,i) = ao_ovlp(:,nbf-i+1)

 ! compute s^(-1/2) (only the first nif ones), stored as diagonal in ao_ovlp
 ao_ovlp = 0d0
 forall(i = 1:nif) ao_ovlp(i,i) = 1d0/DSQRT(s(i))
 deallocate(s)
 ! ao_ovlp now is a nif*nif diagonal matrix

 ! compute Us^(-1/2), where s^(-1/2) is symmetric (in fact, diagonal)
 call dsymm('R', 'U', nbf, nif, 1d0, ao_ovlp, nif, U, nbf, 0d0, mo_coeff, nbf)

 deallocate(U)
end subroutine can_ortho

! generate symmetric orthonormalized atomic orbitals
subroutine sym_ortho(nbf, ao_ovlp, mo_coeff)
 implicit none
 integer :: i
 integer :: nbf
!f2py intent(in) :: nbf
 real(kind=8), parameter :: thresh = 1d-6
 real(kind=8) :: ao_ovlp(nbf,nbf), mo_coeff(nbf,nbf)
!f2py depend(nbf) :: ao_ovlp, mo_coeff
!f2py intent(in,copy) :: ao_ovlp
!f2py intent(out) :: mo_coeff
 real(kind=8), allocatable :: X(:,:), U(:,:), s(:)

 mo_coeff = 0d0
 allocate(s(nbf), source=0d0)
 call diag_get_e_and_vec(nbf, ao_ovlp, s) ! S = UsU^T, U stored in ao_ovlp

 if(ANY(s<thresh)) then
  write(6,'(A)') 'Warning: linear dependence detected in symmetric orthogonalization!'
  write(6,'(A,ES15.8)') 'Default threshold here:', thresh
 end if

 ! reverse the eigenvalues
 allocate(X(nbf,1), source=0d0)
 forall(i = 1:nbf) X(i,1) = s(nbf-i+1)
 s = X(:,1)
 deallocate(X)
 ! now s1 > s2 > s3 > ...

 ! reverse the eigenvectors, according to the descending order of eigenvalues
 allocate(U(nbf, nbf), source=0d0)
 forall(i = 1:nbf) U(:,i) = ao_ovlp(:,nbf-i+1)

 ! compute s^(-1/2), stored as diagonal in ao_ovlp
 ao_ovlp = 0d0
 forall(i = 1:nbf) ao_ovlp(i,i) = 1d0/DSQRT(s(i))
 deallocate(s)

 ! compute Us^(-1/2), where s^(-1/2) is symmetric (in fact, diagonal)
 allocate(X(nbf,nbf), source=0d0)
 call dsymm('R', 'U', nbf, nbf, 1d0, ao_ovlp, nbf, U, nbf, 0d0, X, nbf)

 ! compute X1*U^T, where X1 is Us^(-1/2)
 call dgemm('N', 'T', nbf, nbf, nbf, 1d0, X, nbf, U, nbf, 0d0, mo_coeff, nbf)

 deallocate(X, U)
end subroutine sym_ortho

