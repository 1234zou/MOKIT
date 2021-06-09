! written by jxzou at 20190723
! updated by jxzou at 20200411: add Pipek-Mezey orbital localization (DOI: 10.1063/1.456588)
! updated by jxzou at 20200413: add Cholesky decomposition LMOs (DOI: 10.1063/1.2360264)

! Note: before PySCF-1.6.4, its dumped .molden file is wrong when using Cartesian functions.

!  For GNU compiler, use
! ----------------------------------------
!  f2py -m lo -c lo.f90 --link-lapack_opt
! ----------------------------------------
!  For INTEL compiler, use
! -------------------------------------------------------------------------------
!  f2py -m lo -c lo.f90 --link-lapack_opt --fcompiler=intelem --compiler=intelem
! -------------------------------------------------------------------------------

! generate natural orbitals from density matrix and overlap matrix
subroutine no(nbf, nif, P, S, noon, new_coeff)
 implicit none
 integer :: i, j, lwork, liwork
 integer, parameter :: iout = 6
 integer :: nbf, nif
!f2py intent(in) :: nbf, nif
 ! nbf: the number of atomic basis functions
 ! nif: the number of independent basis functions, i.e., the number of MOs
 ! Note: nif<nbf when linear dependence occurs (be cautious)

 integer, allocatable :: isuppz(:), iwork(:)

 real(kind=8) :: new_coeff(nbf,nif), P(nbf,nbf), S(nbf,nbf)
!f2py intent(in) :: P
!f2py intent(in,copy) :: S
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nif) :: new_coeff
!f2py depend(nbf) :: P, S

 real(kind=8) :: noon(nif)
!f2py intent(out) :: noon
!f2py depend(nif) :: noon

 real(kind=8), allocatable :: sqrt_S(:,:), n_sqrt_S(:,:), PS12(:,:)
 real(kind=8), allocatable :: e(:), U(:,:), work(:)

 allocate(sqrt_S(nbf,nbf), n_sqrt_S(nbf,nbf))
 call mat_dsqrt(nbf, S, sqrt_S, n_sqrt_S) ! solve S^1/2 and S^-1/2

 allocate(PS12(nbf,nbf), source=0d0)
 ! call dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 call dsymm('R', 'L', nbf, nbf, 1d0, sqrt_S, nbf, P, nbf, 0d0, PS12, nbf)
 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 call dgemm('N', 'N', nbf, nbf, nbf, 1d0, sqrt_S, nbf, PS12, nbf, 0d0, S, nbf)
 ! use S to store (S^1/2)P(S^1/2)

 deallocate(PS12, sqrt_S)

 ! call dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z,
 ! ldz, isuppz, work, lwork, iwork, liwork, info)
 lwork = -1
 liwork = -1
 allocate(work(1), iwork(1), isuppz(2*nbf), e(nbf), U(nbf,nbf))
 call dsyevr('V', 'A',  'L', nbf, S, nbf, 0d0, 0d0, 0, 0, 1d-8, i, e, &
             U, nbf, isuppz, work, lwork, iwork, liwork, j)
 lwork = CEILING(work(1))
 liwork = iwork(1)
 deallocate(work, iwork)
 allocate(work(lwork), iwork(liwork))
 call dsyevr('V', 'A',  'L', nbf, S, nbf, 0d0, 0d0, 0, 0, 1d-8, i, e, &
             U, nbf, isuppz, work, lwork, iwork, liwork, j)
 deallocate(isuppz, work, iwork)
 ! eigenvalues in array e are in ascending order

 noon = 0d0
 forall(i = 1:nif, e(nbf-i+1)>0d0) noon(i) = e(nbf-i+1)
 write(iout,'(A)') 'Natural Orbital Occupancy Numbers (NOON):'
 write(iout,'(5(1X,ES15.8))') (noon(i),i=1,nif)
 deallocate(e)

 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 call dgemm('N', 'N', nbf, nif, nbf, 1d0, n_sqrt_S, nbf, U(:,nbf-nif+1:nbf), nbf, 0d0, new_coeff, nbf)
 deallocate(n_sqrt_S, U)

 ! reverse the order of MOs
 allocate(U(nbf,nif))
 forall(i = 1:nif) U(:,i) = new_coeff(:,nif-i+1)
 new_coeff = U
 deallocate(U)
 return
end subroutine no

! compute MO-based density matrix
! Note: the result of this subroutine is identical to that of subroutine
!       solve_ON_matrix, but using a different formula
subroutine get_mo_based_dm(nbf, nif, coeff, S, P, dm)
 implicit none
 integer :: i, j
 integer :: nbf, nif
!f2py intent(in) :: nbf, nif
 integer, parameter :: iout = 6
 real(kind=8) :: coeff(nbf,nif), S(nbf,nbf), P(nbf,nbf), dm(nif,nif)
!f2py intent(in) :: coeff, S, P
!f2py intent(out) :: dm
!f2py depend(nbf,nif) :: coeff
!f2py depend(nbf) :: S, P
!f2py depend(nif) :: dm
 real(kind=8), allocatable :: SC(:,:), PSC(:,:)

 allocate(SC(nbf,nif), source=0d0)
 ! call dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 call dsymm('L', 'L', nbf, nif, 1d0, S, nbf, coeff, nbf, 0d0, SC, nbf)

 allocate(PSC(nbf,nif), source=0d0)
 call dsymm('L', 'L', nbf, nif, 1d0, P, nbf, SC, nbf, 0d0, PSC, nbf)

 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 call dgemm('T', 'N', nif, nif, nbf, 1d0, SC, nbf, PSC, nbf, 0d0, dm, nif)

 write(iout,'(A)') 'MO-based density matrix (Final one electron symbolic density matrix):'
 do i = 1, nif, 1
  do j = i, nif, 1
   write(iout,'(2I4,F15.8)') j, i, dm(j,i)
  end do ! for j
 end do ! for i

 deallocate(SC, PSC)
 return
end subroutine get_mo_based_dm

! use Cholesky factorization/decomposition of the density matrix to generate LMOs
subroutine cholesky(nbf, nif, coeff, new_coeff)
 implicit none
 integer :: i, j, t0, t1, time
 integer, parameter :: iout = 6
 integer :: nbf, nif
!f2py intent(in) :: nbf, nif

 real(kind=8) :: coeff(nbf,nif), new_coeff(nbf,nif)
!f2py intent(in) :: coeff
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nif) :: coeff, new_coeff

 real(kind=8) :: ddot, rtmp1, rtmp2
 real(kind=8), allocatable :: P(:,:)

 t0 = time()
 new_coeff = coeff
 if(nif == 1) then
  write(iout,'(A)') 'Warning in subroutine cholesky: only 1 orbital. Nothing to do.'
  stop
 end if

 allocate(P(nbf,nbf), source=0d0)
 call dgemm('N', 'T', nbf, nbf, nif, 1d0, coeff, nbf, coeff, nbf, 0d0, P, nbf)

 do i = 1, nif, 1
  rtmp1 = P(i,i) - ddot(i-1, new_coeff(i,1:i-1), 1, new_coeff(i,1:i-1), 1)

  if(rtmp1 < 1.0d-12) then
   write(iout,'(A)') 'ERROR in subroutine cholesky: density matrix not positive semidefinite.'
   write(iout,'(A)') 'rtmp1 < 1.0d-12.'
   stop
  end if
  rtmp1 = DSQRT(rtmp1)

  do j = 1, nbf, 1
   rtmp2 = ddot(i-1, new_coeff(j,1:i-1), 1, new_coeff(i,1:i-1), 1)
   new_coeff(j,i) = (P(j,i)-rtmp2)/rtmp1
  end do ! for j
 end do ! for i

 deallocate(P)
 t1 = time()
 write(iout,'(A,I0)') 'Decomposition time(sec):', t1-t0
 return
end subroutine cholesky

! perform Boys orbital localization (Jacobian 2*2 rotations) on a set of MOs
subroutine boys(nbf, nif, coeff, mo_dipole, new_coeff)
 implicit none
 integer :: t0, t1, time
 integer :: nbf, nif
!f2py intent(in) :: nbf, nif
 integer, parameter :: iout = 6

 real(kind=8) :: coeff(nbf,nif), new_coeff(nbf,nif)
 real(kind=8) :: mo_dipole(3,nif,nif)
!f2py intent(in) :: coeff
!f2py intent(in,copy) :: mo_dipole
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nif) :: coeff, new_coeff
!f2py depend(nif) :: mo_dipole
 ! coeff: all orbitals
 ! new_coeff: updated coeff
 ! mo_dipole: the whole MO basis dipole integrals matrix

 t0 = time()
 write(iout,'(/,A)') 'Boys orbital localization begins:'
 new_coeff = coeff

 if(nif == 1) then
  write(iout,'(A)') 'Warning in subroutine boys: only 1 orbital. No rotation.'
  return
 end if

 call serial2by2(nbf, nif, new_coeff, 3, mo_dipole)

 t1 = time()
 write(iout,'(A,I0)') 'Localization time(sec):', t1-t0
 return
end subroutine boys

! perform Pipek-Mezey orbital localization (Jacobian 2*2 rotations) on a set of MOs
subroutine pm(nshl, shl2atm, ang, ibas, cart, nbf, nif, coeff, S, pop, new_coeff)
 implicit none
 integer :: i, j, k, natom, t0, t1, time, i1, i2, i3
 integer :: lwork, liwork
 integer, parameter :: iout = 6
 integer :: nshl, nbf, nif
!f2py intent(in) :: nshl, nbf, nif

 integer :: shl2atm(nshl), ang(nshl), ibas(nshl)
!f2py intent(in) :: shl2atm, ibas
!f2py intent(in,copy) :: ang
!f2py depend(nshl) :: shl2atm, ang, ibas

 integer, allocatable :: bfirst(:) ! size natom
 ! bfirst: the beginning index of basis func. of each atom
 integer, allocatable :: iwork(:), isuppz(:)

 character(len=*), intent(in) :: pop ! 'mulliken' or 'lowdin'

 real(kind=8) :: coeff(nbf,nif), S(nbf,nbf), new_coeff(nbf,nif)
!f2py intent(in) :: coeff, S
!f2py intent(out) :: new_coeff
!f2py depend(nbf) :: S
!f2py depend(nbf,nif) :: coeff, new_coeff

 real(kind=8) :: ddot, rtmp
 real(kind=8), allocatable :: gross(:,:,:), SC(:,:), rootS(:,:)
 ! gross0: old gross population of an orthonormalized MO, (natom,nif,nif)
 ! SC: MATMUL(S,C)
 ! rootS: S^-1/2
 ! e: eigenvalues of S after diagonalization
 ! ev: eigenvectors
 real(kind=8), allocatable :: e(:), ev(:,:), work(:)

 logical :: cart
!f2py intent(in) :: cart

 t0 = time()
 write(iout,'(/,A)') 'PM orbital localization begins: using '//TRIM(pop)//' population'
 new_coeff = coeff

 if(nif == 1) then
  write(iout,'(A)') 'Warning in subroutine pm: only 1 orbital. No rotation.'
  return
 end if

 if(ANY(ang<0)) then
  write(iout,'(A)') 'ERROR in subroutine pm: there exists ang(i)<0.'
  stop
 end if

 if(cart) then
  forall(i = 1:nshl) ang(i) = (ang(i)+1)*(ang(i)+2)/2
 else
  forall(i = 1:nshl) ang(i) = 2*ang(i) + 1
 end if

 j = DOT_PRODUCT(ang, ibas)
 if(j /= nbf) then
  write(iout,'(A)') 'ERROR in subroutine pm: number of basis function is&
                   & inconsistent between j and nbf.'
  write(iout,'(2(A,I0))') 'j=', j, ', nbf=', nbf
  stop
 end if

 natom = shl2atm(nshl) + 1
 write(iout,'(A,I0)') 'natom=', natom
 allocate(bfirst(natom+1), source=0)
 bfirst(1) = 1
 do i = 1, nshl, 1
  j = shl2atm(i) + 2
  bfirst(j) = bfirst(j) + ang(i)*ibas(i)
 end do ! for i

 do i = 2, natom+1, 1
  bfirst(i) = bfirst(i) + bfirst(i-1)
 end do ! for i

 if(pop == 'mulliken') then
  allocate(SC(nbf,nif), source=0d0)
  call dsymm('L', 'L', nbf, nif, 1d0, S, nbf, coeff, nbf, 0d0, SC, nbf)
  allocate(gross(natom,nif,nif), source=0d0)
 
  do i = 1, nif, 1
   do j = i, nif, 1
    do k = 1, natom, 1
     i1 = bfirst(k); i2 = bfirst(k+1)-1; i3 = bfirst(k+1)-bfirst(k)
     rtmp = ddot(i3, coeff(i1:i2,i), 1, SC(i1:i2,j), 1) + &
            ddot(i3, coeff(i1:i2,j), 1, SC(i1:i2,i), 1)
 
     rtmp = 0.5d0*rtmp
     gross(k,j,i) = rtmp
     gross(k,i,j) = rtmp
    end do ! for k
   end do ! for j
  end do ! for i

 else if(pop == 'lowdin') then
  allocate(e(nbf), ev(nbf, nbf), isuppz(2*nbf))
  allocate(SC(nbf,nbf), source=S) ! use SC to temporarily store S

  ! call dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z,
  ! ldz, isuppz, work, lwork, iwork, liwork, info)
  lwork = -1
  liwork = -1
  allocate(work(1), iwork(1))
  call dsyevr('V', 'A', 'L', nbf, SC, nbf, 0d0, 0d0, 0, 0, 1d-8, k, e, &
              ev, nbf, isuppz, work, lwork, iwork, liwork, i)
  lwork = CEILING(work(1))
  liwork = iwork(1)
  deallocate(work, iwork)
  allocate(work(lwork), iwork(liwork))
  call dsyevr('V', 'A', 'L', nbf, SC, nbf, 0d0, 0d0, 0, 0, 1d-8, k, e, &
              ev, nbf, isuppz, work, lwork, iwork, liwork, i)
  deallocate(SC, isuppz, work, iwork)

  allocate(rootS(nbf,nbf), source=0d0)
  forall(i=1:nbf) rootS(i,i) = DSQRT(DABS(e(i)))
  deallocate(e)
  allocate(SC(nbf,nbf)) ! use SC to temporarily store US^1/2
  call dsymm('R', 'L', nbf, nbf, 1d0, rootS, nbf, ev, nbf, 0d0, SC, nbf)
  call dgemm('N', 'T', nbf, nbf, nbf, 1d0, SC, nbf, ev, nbf, 0d0, rootS, nbf)
  deallocate(ev, SC)
  allocate(SC(nbf,nif), source=0d0) ! use SC to temporarily store (S^-1/2)C
  call dsymm('L', 'L', nbf, nif, 1d0, rootS, nbf, coeff, nbf, 0d0, SC, nbf)
  deallocate(rootS)

  allocate(gross(natom,nif,nif), source=0d0)

  do i = 1, nif, 1
   do j = i, nif, 1
    do k = 1, natom, 1
     i1 = bfirst(k); i2 = bfirst(k+1)-1; i3 = bfirst(k+1)-bfirst(k)
     rtmp = ddot(i3, SC(i1:i2,i), 1, SC(i1:i2,j), 1)
     gross(k,j,i) = rtmp
     gross(k,i,j) = rtmp
    end do ! for k
   end do ! for j
  end do ! for i

 else
  write(iout,'(A)') 'ERROR in subroutine pm: wrong population method provided.'
  write(iout,'(A)') "Only 'mulliken' or 'lowdin' supported. But input pop="//pop
  stop
 end if

 deallocate(bfirst, SC)

 call serial2by2(nbf, nif, new_coeff, natom, gross)
 deallocate(gross)

 t1 = time()
 write(iout,'(A,I0)') 'Localization time(sec):', t1-t0
 return
end subroutine pm

! perform serial 2-by-2 rotation on given MOs
subroutine serial2by2(nbf, nif, coeff, ncomp, mo_dipole)
 implicit none
 integer :: i, j, k, niter
 integer, intent(in) :: nbf, nif, ncomp
 integer, parameter :: iout = 6, niter_max = 9999
 ! niter_max: max number of iterations

 real(kind=8), parameter :: threshold1 = 1d-7, threshold2 = 1d-6
 ! threshold1: determine whether to rotate (and update MOs, dipole integrals)
 ! threshold2: determine whether rotation/localization converged
 ! if threshold2 is set to 1d-5, the localization is sometimes not thorough
 real(kind=8), parameter :: PI = 4d0*DATAN(1d0)

 real(kind=8), intent(inout) :: coeff(nbf,nif), mo_dipole(ncomp,nif,nif)
 ! for Boys, ncomp = 3
 ! for PM,   ncomp = natom, mo_dipole is actually population matrix

 real(kind=8) :: ddot, rtmp, sum_change
 real(kind=8) :: Aij, Bij, alpha, sin_4a, cos_a, sin_a
 real(kind=8) :: cc, ss, sin_2a, cos_2a
 real(kind=8), allocatable :: dipole(:,:,:), motmp(:,:), vtmp(:,:), vdiff(:)
 ! dipole: for Boys, store updated dipole integrals matrix
 !         for PM,   store updated population matrix
 ! motmp: store 2 MOs
 ! vtmp: array dipole(:,y,x)
 ! vdiff := vtmp(:,1) - vtmp(:,3)

 allocate(dipole(ncomp,nif,2), motmp(nbf,2), vtmp(ncomp,3), vdiff(ncomp))

 ! perform 2*2 rotation
 niter = 0
 do while(niter <= niter_max)
  sum_change = 0d0

  do i = 1, nif-1, 1
   do j = i+1, nif, 1
    vtmp(:,1) = mo_dipole(:,i,i)
    vtmp(:,2) = mo_dipole(:,j,i)
    vtmp(:,3) = mo_dipole(:,j,j)
    vdiff = vtmp(:,1) - vtmp(:,3)
    Aij = ddot(ncomp,vtmp(:,2),1,vtmp(:,2),1) - 0.25d0*ddot(ncomp,vdiff,1,vdiff,1)
    Bij = ddot(ncomp, vdiff, 1, vtmp(:,2), 1)
    rtmp = HYPOT(Aij, Bij)
    sin_4a = Bij/rtmp
    sin_4a = MAX(-1d0, MIN(sin_4a, 1d0)) ! in case of numerical error
    alpha = DASIN(sin_4a) ! [-PI/2,PI/2]
    if(Aij > 0d0) then
     alpha = PI - alpha
    else if(Aij<0d0 .and. Bij<0d0) then
     alpha = 2d0*PI + alpha
    end if
    alpha = 0.25d0*alpha
    ! if theta/alpha is very close to zero or PI/2, not to rotate
    if(alpha<threshold1 .or. 0.5d0*PI-alpha<threshold1) cycle

    cos_a = DCOS(alpha)
    sin_a = DSIN(alpha)
    sum_change = sum_change + Aij + rtmp

    ! update two orbitals
    motmp(:,1) = coeff(:,i)
    motmp(:,2) = coeff(:,j)
    coeff(:,i) = cos_a*motmp(:,1) + sin_a*motmp(:,2)
    coeff(:,j) = cos_a*motmp(:,2) - sin_a*motmp(:,1)

    ! update corresponding dipole integrals, only indices in range to be updated
    cc = cos_a*cos_a
    ss = sin_a*sin_a
    sin_2a = 2d0*sin_a*cos_a
    cos_2a = cc - ss
    dipole(:,i,1) = cc*vtmp(:,1) + ss*vtmp(:,3) + sin_2a*vtmp(:,2)
    dipole(:,j,2) = ss*vtmp(:,1) + cc*vtmp(:,3) - sin_2a*vtmp(:,2)
    dipole(:,j,1) = cos_2a*vtmp(:,2) - 0.5d0*sin_2a*vdiff
    dipole(:,i,2) = dipole(:,j,1)

    do k = 1, nif, 1
     if(k==i .or. k==j) cycle
     dipole(:,k,1) = cos_a*mo_dipole(:,k,i) + sin_a*mo_dipole(:,k,j)
     dipole(:,k,2) = cos_a*mo_dipole(:,k,j) - sin_a*mo_dipole(:,k,i)
    end do ! for k

    mo_dipole(:,:,i) = dipole(:,:,1)
    mo_dipole(:,:,j) = dipole(:,:,2)
    mo_dipole(:,i,:) = dipole(:,:,1)
    mo_dipole(:,j,:) = dipole(:,:,2)
   end do ! for j
  end do ! for i

  niter = niter + 1
  write(iout,'(A,I5,A,F15.7)') 'niter=', niter, ', sum_change=', sum_change
  if(sum_change < threshold2) exit
 end do ! for while

 deallocate(motmp,vdiff,vtmp,dipole)

 if(niter <= niter_max) then
  write(iout,'(A)') 'Localization converged successfully.'
 else
  write(iout,'(A)') 'ERROR in subroutine serial2by2: niter_max exceeded.'
  stop
 end if
 return
end subroutine serial2by2

! get the value of the modified Boys function
subroutine get_mboys(nif, ncore, npair, nopen, mo_dipole)
 implicit none
 integer i, j, nocc
 integer, parameter :: iout = 6
 integer nif, ncore, npair, nopen
!f2py intent(in) :: nif, nocc, npair, nopen
 real(kind=8) :: mo_dipole(3,nif,nif)
!f2py intent(in) :: mo_dipole
!f2py depend(nif) mo_dipole
 real(kind=8) :: ddot, fBoys, temp_dipole(3)

 nocc = ncore + npair + nopen

 fBoys = 0d0
 j = ncore + npair
 do i = 1, j, 1
  temp_dipole =  mo_dipole(1:3,i,i)
  fBoys = fBoys + ddot(3, temp_dipole, 1, temp_dipole, 1)
 end do
 fBoys = DSQRT(fBoys/DBLE(npair))
 write(iout,'(A,F13.6)') 'In occ, Modified f(Boys)=', fBoys
 fBoys = 0d0
 j = nocc + npair
 do i = nocc+1, j, 1
  temp_dipole =  mo_dipole(1:3,i,i)
  fBoys = fBoys + ddot(3, temp_dipole, 1, temp_dipole, 1)
 end do
 fBoys = DSQRT(fBoys/DBLE(npair))
 write(iout,'(A,F13.6)') 'In vir, Modified f(Boys)=', fBoys
 return
end subroutine get_mboys

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
 return
end subroutine mat_dsqrt

! perform immediate Boys localization by diagonalizing the DxDx+DyDy+DzDz
subroutine boys_diag(nbf, nmo, mo_coeff, mo_dipole, new_coeff)
 implicit none
 integer :: i, lwork, liwork
 integer :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 integer, parameter :: iout = 6
 integer, allocatable :: iwork(:)
 real(kind=8) :: mo_coeff(nbf,nmo), mo_dipole(3,nmo,nmo), new_coeff(nbf,nmo)
!f2py intent(in) :: mo_coeff, mo_dipole
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nmo) :: mo_coeff, new_coeff
!f2py depend(nmo) :: mo_dipole
 real(kind=8), allocatable :: f(:,:)
 real(kind=8), allocatable :: w(:), work(:)

 allocate(f(nmo,nmo), source=0d0)
 ! call dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 call dsymm('L','L',nmo,nmo, 1d0,mo_dipole(1,:,:),nmo, mo_dipole(1,:,:),nmo, 0d0,f,nmo)
 call dsymm('L','L',nmo,nmo, 1d0,mo_dipole(2,:,:),nmo, mo_dipole(2,:,:),nmo, 1d0,f,nmo)
 call dsymm('L','L',nmo,nmo, 1d0,mo_dipole(3,:,:),nmo, mo_dipole(3,:,:),nmo, 1d0,f,nmo)

 ! call dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
 allocate(work(1), iwork(1))
 lwork = -1
 liwork = -1
 call dsyevd('V', 'U', nmo, f, nmo, w, work, lwork, iwork, liwork, i)
 lwork = CEILING(work(1))
 liwork = iwork(1)
 deallocate(work, iwork)
 allocate(w(nmo), work(lwork), iwork(liwork))
 call dsyevd('V', 'U', nmo, f, nmo, w, work, lwork, iwork, liwork, i)
 deallocate(w, work, iwork)

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine boys_diag: info /= 0 in dsyevd.'
  write(iout,'(A,I0)') 'info = ', i
  stop
 end if

 call dgemm('N','N',nbf,nmo,nmo, 1d0,mo_coeff,nbf, f,nmo, 0d0,new_coeff,nbf)
 deallocate(f)
 return
end subroutine boys_diag

subroutine solve_boys_lamda_matrix(nbf, nmo, coeff, lo_coeff, mo_dipole)
 implicit none
 integer :: i, j
 integer :: nbf, nmo
 integer, parameter :: iout = 6
!f2py intent(in) :: nbf, nmo
 real(kind=8) :: coeff(nbf,nmo), lo_coeff(nbf,nmo), mo_dipole(3,nmo,nmo)
!f2py intent(in) :: coeff, lo_coeff, mo_dipole
!f2py depend(nbf,nmo) :: coeff, lo_coeff
!f2py depend(nmo) ::  mo_dipole
 real(kind=8), allocatable :: f(:,:), U(:,:), fU(:,:), lamda(:,:)

 allocate(f(nmo,nmo), source=0d0)
 ! call dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 call dsymm('L','L',nmo,nmo, 1d0,mo_dipole(1,:,:),nmo, mo_dipole(1,:,:),nmo, 0d0,f,nmo)
 call dsymm('L','L',nmo,nmo, 1d0,mo_dipole(2,:,:),nmo, mo_dipole(2,:,:),nmo, 1d0,f,nmo)
 call dsymm('L','L',nmo,nmo, 1d0,mo_dipole(3,:,:),nmo, mo_dipole(3,:,:),nmo, 1d0,f,nmo)

 allocate(U(nmo,nmo), source=0d0)
 call get_u(nbf, nmo, coeff, lo_coeff, U) 
 allocate(fU(nmo,nmo), source=0d0)
 call dsymm('L','L',nmo,nmo, 1d0,f,nmo, U,nmo, 0d0,fU,nmo)
 deallocate(f)
 allocate(lamda(nmo,nmo), source=0d0)
 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 call dgemm('T', 'N', nmo,nmo,nmo, -4d0, U,nmo, fU,nmo, 0d0, lamda,nmo)
 deallocate(fU, U)

 do i = 1, nmo, 1
  write(iout,'(20F14.4)') (lamda(j,i),j=1,i)
 end do ! for i
 deallocate(lamda)
 return
end subroutine solve_boys_lamda_matrix

subroutine get_u(nbf, nmo, coeff, lo_coeff, u)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nmo
 integer, allocatable :: ipiv(:)
 real(kind=8), intent(in) :: coeff(nbf,nmo), lo_coeff(nbf,nmo)
 real(kind=8), intent(out) :: u(nmo,nmo)
 real(kind=8), allocatable :: coeff1(:,:), lo_coeff1(:,:)

 ! mkl Syntax FORTRAN 77:
 ! ?getrf: Computes the LU factorization of a general m-by-n matrix
 ! call dgetrf(m, n, a, lda, ipiv, info)

 ! ?getrs: Solves a system of linear equations with an LU-factored square
 !  matrix, with multiple right-hand sides.
 ! call dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)

 ! ?gemm: Computes a matrix-matrix product with general matrix
 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

 ! find the unitary (orthogonal) matrix between coeff1 and lo_coeff1,
 !  where coeff1*U = lo_coeff1
 allocate(coeff1(nbf,nmo), lo_coeff1(nbf,nmo))
 coeff1 = coeff
 lo_coeff1 = lo_coeff
 allocate(ipiv(min(nbf,nmo)), source=0)

 call dgetrf(nbf, nmo, coeff1, nbf, ipiv, i)
 call dgetrs('N', nmo, nmo, coeff1, nbf, ipiv, lo_coeff1, nbf, i)
 deallocate(ipiv, coeff1)

 u = lo_coeff1(1:nmo,1:nmo)
 deallocate(lo_coeff1)
 return
end subroutine get_u

