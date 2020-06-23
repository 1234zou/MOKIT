! written by jxzou at 20200404: do population analysis

! This subroutine is designed to be called as a module by Python.
!  For INTEL compiler, use
! ---------------------------------------------------------------------------------
!  f2py -m pop -c pop.f90 --link-lapack_opt --fcompiler=intelem --compiler=intelem
! ---------------------------------------------------------------------------------
!  For GNU compiler, use
! ------------------------------------------
!  f2py -m pop -c pop.f90 --link-lapack_opt
! ------------------------------------------

!  to compile this file (a pop.so file will be generated). Then in
!  Python you can import the pop module.

! perform Mulliken population analysis
subroutine mulliken(nshl, shl2atm, ang, ibas, cart, nbf, P, S, natom, eff_nuc)
 implicit none
 integer :: i, j
 integer, parameter :: iout = 6

 integer :: nshl, nbf, natom
!f2py intent(in) :: nshl, nbf, natom

 integer :: shl2atm(nshl), ang(nshl), ibas(nshl), eff_nuc(natom)
!f2py intent(in) :: shl2atm, ibas, eff_nuc
!f2py intent(in,copy) :: ang
!f2py depend(nshl) :: shl2atm, ang, ibas
!f2py depend(natom) :: eff_nuc

 integer, allocatable :: bfirst(:) ! size natom
 ! bfirst: the beginning index of basis func. of each atom

 real(kind=8) :: P(nbf,nbf), S(nbf,nbf)
!f2py intent(in) :: P, S
!f2py depend(nbf) :: P, S

 real(kind=8) :: rtmp, ddot
 real(kind=8), allocatable :: gross(:) ! size natom

 logical :: cart
!f2py intent(in) :: cart

 if(cart) then
  forall(i = 1:nshl) ang(i) = (ang(i)+1)*(ang(i)+2)/2
 else
  forall(i = 1:nshl) ang(i) = 2*ang(i) + 1
 end if

 j = DOT_PRODUCT(ang, ibas)
 if(j /= nbf) then
  write(iout,'(A)') 'ERROR in subroutine mulliken_pop: number of basis function is&
                   & inconsistent between j and nbf.'
  write(iout,'(2(A,I0))') 'j=', j, ', nbf=', nbf
  stop
 end if

 allocate(bfirst(natom+1), source=0)
 bfirst(1) = 1

 do i = 1, nshl, 1
  j = shl2atm(i) + 2
  bfirst(j) = bfirst(j) + ang(i)*ibas(i)
 end do ! for i

 do i = 2, natom+1, 1
  bfirst(i) = bfirst(i) + bfirst(i-1)
 end do ! for i

 allocate(gross(natom), source=0.0d0)

 do i = 1, natom, 1
  rtmp = 0.0d0

  do j = bfirst(i), bfirst(i+1)-1, 1
   rtmp = rtmp + ddot(nbf, P(j,:), 1, S(:,j), 1)
  end do ! for j

  gross(i) = rtmp
 end do ! for i

 deallocate(bfirst)
 gross = eff_nuc - gross

 write(iout,'(/,A)') 'Mulliken population:'

 do i = 1, natom, 1
  write(iout,'(I5,1X,F11.6)') i, gross(i)
 end do ! for i

 deallocate(gross)
 return
end subroutine mulliken

! perform lowdin/lowedin population analysis
subroutine lowdin(nshl, shl2atm, ang, ibas, cart, nbf, P, S, natom, eff_nuc)
 implicit none
 integer :: i, j, m, lwork, liwork
 integer, parameter :: iout = 6

 integer :: nshl, nbf, natom
!f2py intent(in) :: nshl, nbf, natom

 integer :: shl2atm(nshl), ang(nshl), ibas(nshl), eff_nuc(natom)
!f2py intent(in) :: shl2atm, ibas, eff_nuc
!f2py intent(in,copy) :: ang
!f2py depend(nshl) :: shl2atm, ang, ibas
!f2py depend(natom) :: eff_nuc

 integer, allocatable :: bfirst(:) ! size natom
 ! bfirst: the beginning index of basis func. of each atom

 integer, allocatable :: iwork(:), isuppz(:)

 real(kind=8) :: P(nbf,nbf), S(nbf,nbf)
!f2py intent(in) :: P
!f2py intent(in,copy) :: S
!f2py depend(nbf) :: P, S

 real(kind=8) :: rtmp, ddot
 real(kind=8), allocatable :: gross(:) ! size natom
 real(kind=8), allocatable :: work(:), e(:), ev(:,:), sqrt_S_P(:,:)
 ! e: eigenvalues
 ! ev: eigenvectors
 ! sqrt_S_P: S^(1/2)*P

 logical :: cart
!f2py intent(in) :: cart

 if(cart) then
  forall(i = 1:nshl) ang(i) = (ang(i)+1)*(ang(i)+2)/2
 else
  forall(i = 1:nshl) ang(i) = 2*ang(i) + 1
 end if

 j = DOT_PRODUCT(ang, ibas)
 if(j /= nbf) then
  write(iout,'(A)') 'ERROR in subroutine mulliken_pop: number of basis function is&
                   & inconsistent between j and nbf.'
  write(iout,'(2(A,I0))') 'j=', j, ', nbf=', nbf
  stop
 end if

 allocate(bfirst(natom+1), source=0)
 bfirst(1) = 1

 do i = 1, nshl, 1
  j = shl2atm(i) + 2
  bfirst(j) = bfirst(j) + ang(i)*ibas(i)
 end do ! for i

 do i = 2, natom+1, 1
  bfirst(i) = bfirst(i) + bfirst(i-1)
 end do ! for i

 ! call dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z,
 ! ldz, isuppz, work, lwork, iwork, liwork, info)
 allocate(e(nbf), ev(nbf, nbf), isuppz(2*nbf))
 lwork = -1
 liwork = -1
 allocate(work(1), iwork(1))
 call dsyevr('V', 'A', 'L', nbf, S, nbf, 0.0d0, 0.0d0, 0, 0, 1.0d-6, m, e, ev,&
             nbf, isuppz, work, lwork, iwork, liwork, i)
 lwork = CEILING(work(1))
 liwork = iwork(1)
 deallocate(work, iwork)
 allocate(work(lwork), iwork(liwork))
 call dsyevr('V', 'A', 'L', nbf, S, nbf, 0.0d0, 0.0d0, 0, 0, 1.0d-6, m, e, ev,&
             nbf, isuppz, work, lwork, iwork, liwork, i)
 deallocate(isuppz, work, iwork)

 S = 0.0d0
 forall(i = 1:nbf) S(i,i) = DSQRT(DABS(e(i)))
 deallocate(e)
 allocate(sqrt_S_P(nbf,nbf))
 ! call dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 call dsymm('R', 'U', nbf, nbf, 1.0d0, S, nbf, ev, nbf, 0.0d0, sqrt_S_P, nbf)
 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 call dgemm('N', 'T', nbf, nbf, nbf, 1.0d0, sqrt_S_P, nbf, ev, nbf, 0.0d0, S, nbf)
 deallocate(ev, sqrt_S_P)

 allocate(gross(natom), source=0.0d0)
 allocate(e(nbf))

 do i = 1, natom, 1
  rtmp = 0.0d0

  do j = bfirst(i), bfirst(i+1)-1, 1
   e = S(j,:)

   do m = 1, nbf, 1
    rtmp = rtmp + ddot(nbf,e,1,P(:,m),1)*S(m,j)
   end do ! for m

  end do ! for j

  gross(i) = rtmp
 end do ! for i

 deallocate(bfirst, e)
 gross = eff_nuc - gross

 write(iout,'(/,A)') 'Lowdin population:'

 do i = 1, natom, 1
  write(iout,'(I5,1X,F11.6)') i, gross(i)
 end do ! for i

 deallocate(gross)
 return
end subroutine lowdin

