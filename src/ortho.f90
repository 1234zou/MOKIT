! written by jxzou at 20191029
! updated by jxzou at 20200403: move check_orthonormal into this file
! TODO: Schmidt orthogonalization

! check whether a given set of MOs are orthonormal
subroutine check_orthonormal(nbf, nif, coeff, S)
 implicit none
 integer :: i, j, k2(2)
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: coeff(nbf,nif), S(nbf,nbf)
!f2py intent(in) :: coeff, S
!f2py depend(nbf,nif) :: coeff
!f2py depend(nbf) :: S
 real(kind=8) :: maxv
 real(kind=8), allocatable :: ctsc(:,:)

 allocate(ctsc(nif,nif))
 call calc_CTSC(nbf, nif, coeff, S, ctsc)

 do i = 1, nif, 1
  ctsc(i,i) = ctsc(i,i) - 1d0
 end do ! for i

 ctsc = DABS(ctsc)
 k2 = MAXLOC(ctsc)
 j = k2(1); i = k2(2)
 maxv = ctsc(j,i)
 deallocate(ctsc)

 write(6,'(/,2(A,I0),A,ES15.8)') 'Orthonormality check: j=', j, ', i=', i,&
                                 ', maxv=', maxv
 if(maxv > 1d-2) then
  write(6,'(/,A)') 'Warning: severe non-orthonormal problem!'
  write(6,'(2(A,I0))') 'nbf=', nbf, ', nif=', nif
 end if
end subroutine check_orthonormal

! check whether a given set of complex MOs are orthonormal
! Note: 
! 1) this subroutine is specially designed for complex GHF MOs in PySCF, where 
!  the basis oerder of a MO is (ar, ai)1, (ar, ai)2, ..., (br, bi)1, (br, bi)2
! 2) do not use this subroutine for complex GHF MOs in Gaussian, where the
!  basis oerder of a MO is (ar, ai, br, bi)1, (ar, ai, br, bi)2, ...
subroutine check_cghf_orthonormal(nbf, nif, coeff, S)
 implicit none
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 complex(kind=8), intent(in) :: coeff(nbf,nif)
!f2py intent(in) :: coeff
!f2py depend(nbf,nif) :: coeff
 real(kind=8), intent(in) :: S(nbf/2,nbf/2)
!f2py intent(in) :: S
!f2py depend(nbf) :: S
 real(kind=8) :: rtmp
 real(kind=8), allocatable :: C_r(:,:), C_i(:,:), A(:,:), A_b(:,:)

 allocate(C_r(nbf,nif), source=REAL(coeff))  ! real part
 allocate(C_i(nbf,nif), source=AIMAG(coeff)) ! imaginary part

 allocate(A(nif,nif))
 call calc_CTSCp(nbf/2, nif, C_r(1:nbf/2,:), S, C_i(1:nbf/2,:), A)

 allocate(A_b(nif,nif))
 call calc_CTSCp(nbf/2, nif, C_r(nbf/2+1:nbf,:), S, C_i(nbf/2+1:nbf,:), A_b)

 A = DABS(A) + DABS(A_b)
 deallocate(A_b)

 rtmp = SUM(A)/DBLE(nif*nif)
 deallocate(A)
 write(6,'(/,A,ES15.8)') 'SUM(DABS(A))/nif^2=', rtmp
 if(rtmp > 1d-6) then
  write(6,'(/,A)') 'Warning in subroutine check_cghf_orthonormal: rtmp>1d-6.'
  write(6,'(A)') 'SUM(DABS(A))/nif^2 should be zero. Orthogonality not good.'
 end if

 call check_orthonormal(nbf/2, nif, C_r(1:nbf/2,:), S)
 call check_orthonormal(nbf/2, nif, C_r(nbf/2+1:nbf,:), S)
 deallocate(C_r)

 call check_orthonormal(nbf/2, nif, C_i(1:nbf/2,:), S)
 call check_orthonormal(nbf/2, nif, C_i(nbf/2+1:nbf,:), S)
 deallocate(C_i)
end subroutine check_cghf_orthonormal

! Perform canonical/symmetric orthonormalization on a set of non-orthogonal MOs.
! This subroutine is different from subroutines can_ortho and sym_ortho below.
subroutine orthonormalize_orb(sym_ortho, prt_warn, nbf, nmo, ao_ovlp, old_mo, &
                              new_mo)
 implicit none
 integer :: i, j, nif
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 real(kind=8), parameter :: thres = 1d-6
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf), old_mo(nbf,nmo)
!f2py intent(in) :: ao_ovlp, old_mo
!f2py depend(nbf) :: ao_ovlp
!f2py depend(nbf,nmo) :: old_mo
 real(kind=8), intent(out) :: new_mo(nbf,nmo)
!f2py intent(out) :: new_mo
!f2py depend(nbf,nmo) :: new_mo
 real(kind=8), allocatable :: X(:,:), Y(:,:), ev(:)
 real(kind=8), allocatable :: Sp(:,:) ! S', S prime
 logical, intent(in) :: sym_ortho, prt_warn
!f2py intent(in) :: sym_ortho, prt_warn
 ! sym_ortho: True/False for symmetric/canonical orthonormalization

 if(nmo > nbf) then
  write(6,'(/,A)') 'ERROR in subroutine orthonormalize_orb: nmo>nbf. Impossible.'
  write(6,'(A,L1)') 'sym_ortho=', sym_ortho
  write(6,'(2(A,I0))') 'nbf=', nbf, ', nmo=', nmo
  stop
 end if

 if(nmo == 1) then ! only one MO, normalize and return
  new_mo = old_mo
  call normalize_mo(nbf, ao_ovlp, new_mo(:,1))
  return
 end if

 new_mo = 0d0 ! initialization
 allocate(Sp(nmo,nmo))
 call calc_CTSC(nbf, nmo, old_mo, ao_ovlp, Sp) ! (C^T)SC = S'
 ! S' is not I, so C is not orthonormal

 allocate(ev(nmo), source=0d0)
 call diag_get_e_and_vec(nmo, Sp, ev) ! S' = U(s')(U^T), U stored in Sp
 call reverse_e_and_vec(nmo, Sp, ev)  ! in descending order

 nif = COUNT(ev > thres)
 if(nmo > nif) then
  if(prt_warn) then
   write(6,'(/,A)') REPEAT('-',79)
   write(6,'(A)') 'Warning from subroutine orthonormalize_orb: the input integer&
                  & nmo is larger'
   write(6,'(A)') 'than independent MOs nif. Independent MOs will be stored in n&
                  &ew_mo(:,1:nif),'
   write(6,'(A)') 'while new_mo(:,nif+1:nmo) will be set to zero. This is not an&
                  & error, but just'
   write(6,'(A)') 'to remind you that you are supposed to know what you are calc&
                  &ulating.'
   write(6,'(A)') REPEAT('-',79)
   write(6,'(3(A,I0))') 'nbf=', nbf, ', nif=', nif, ', nmo=', nmo
  end if
  if(sym_ortho) then
   write(6,'(/,A)') 'ERROR in subroutine orthonormalize_orb: symmetric orthogon&
                    &alization cannot'
   write(6,'(A)') 'be used when there is linear dependency among these MOs. You&
                  & can use canonical'
   write(6,'(A)') 'orthogonalization instead. Eigenvalues='
   write(6,'(5(1X,ES15.8))') ev
   deallocate(Sp, ev)
   stop
  end if
 end if

 forall(i = 1:nif) ev(i) = 1d0/DSQRT(ev(i))
 allocate(X(nmo,nmo), source=0d0) ! X = U((s')^(-1/2))
 forall(i=1:nmo, j=1:nif) X(i,j) = Sp(i,j)*ev(j)
 deallocate(ev)

 if(sym_ortho) then ! X = X(U^T)
  allocate(Y(nmo,nmo), source=0d0)
  call dgemm('N','T', nmo,nmo,nmo, 1d0,X,nmo, Sp,nmo, 0d0,Y,nmo)
  X = Y
  deallocate(Y)
 end if
 deallocate(Sp)

 ! C' = CX
 call dgemm('N','N', nbf,nmo,nmo, 1d0,old_mo,nbf, X,nmo, 0d0,new_mo,nbf)
 deallocate(X)
end subroutine orthonormalize_orb

! Generate orthonormalized atomic orbitals (OAO) from the given AO overlap
! integral matrix. If there exists linear dependency in ao_ovlp, the output
! nif will be smaller than nbf, and mo(:,1:nif) are orthonormalized MOs.
subroutine gen_ortho_mo(nbf, ao_ovlp, nif, mo)
 implicit none
 integer :: i, k
 integer, intent(in) :: nbf
!f2py intent(in) :: nbf
 integer, intent(out) :: nif
!f2py intent(out) :: nif
 real(kind=8), parameter :: thresh = 1d-6 ! commonly used threshold
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf)
!f2py intent(in) :: ao_ovlp
!f2py depend(nbf) :: ao_ovlp
 real(kind=8), intent(out) :: mo(nbf,nbf)
!f2py intent(out) :: mo
!f2py depend(nbf) :: mo
 real(kind=8), allocatable :: diag(:,:), U(:,:), s(:)

 allocate(U(nbf,nbf), source=ao_ovlp)
 allocate(s(nbf))
 call diag_get_e_and_vec(nbf, U, s)
 ! s(1) <= s(2) <= s(3) <= ...
 nif = COUNT(s > thresh)

 if(nbf > nif) then ! linear dep
  ! compute s^(-1/2), stored as diagonal elements
  allocate(diag(nif,nif), source=0d0)
  k = nbf - nif
  forall(i = 1:nif) diag(i,i) = 1d0/DSQRT(s(k+i))
  deallocate(s)
  mo = 0d0
  ! compute Us^(-1/2)
  call dsymm('R','L',nbf,nif,1d0,diag,nif,U(:,k+1:nbf),nbf,0d0,mo(:,1:nif),nbf)
  deallocate(diag, U)
 else if(nbf == nif) then ! no linear dep
  s = 1d0/DSQRT(s)
  call calc_usut(nbf, s, U, mo)
  deallocate(s, U)
 else
  write(6,'(/,A)') 'ERROR in subroutine gen_ortho_mo: nbf<nif. Impossible.'
  write(6,'(A,2I8)') 'nbf, nif=', nbf, nif
  stop
 end if
end subroutine gen_ortho_mo

