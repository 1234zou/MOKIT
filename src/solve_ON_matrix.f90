! written by jxzou at 20200502: compute the (non-diagonal) occupation number
!  matrix of a set of MOs

! Note:
! (1) the target MOs are hold in .fch(k) file fname1
! (2) the Natural Orbitals (NOs) are stored in .fch(k) file fname2
! (3) the NOON must also be stored in .fch(k) file fname2 ('Alpha Orbital Energies' section)
! (4) the digonal elements of occupation number matrix will be written into
!     the 'Alpha Orbital Energies' section of file fname1
! (5) the whole non-diagonal occupation number matrix will be output into a file

program main
 implicit none
 integer :: i, idx(2)
! idx(1)/idx(2): the start/final NO index
 character(len=8) :: buf
 character(len=240) :: fname1, fname2

 i = iargc()
 if(.not. (i==2 .or. i==4)) then
  write(6,'(/,A)') ' ERROR in subroutine solve_ON_matrix: wrong command line arguments.'
  write(6,'(A)') ' Format: solve_ON_matrix MO_file NO_file [idx1] [idx2]'
  write(6,'(/,A)') ' Example 1(Gaussian): solve_ON_matrix mo.fchk no.fchk'
  write(6,'(/,A)') ' Example 2(Gaussian): solve_ON_matrix mo.fchk no.fchk 19 24'
  write(6,'(/,A)') ' Example 3(OpenMolcas): solve_ON_matrix mo.UnaOrb no.RasOrb.1'
  write(6,'(/,A,/)') ' Example 4(OpenMolcas): solve_ON_matrix mo.UnaOrb no.RasOrb.1 19 24'
  stop
 end if

 idx = 0
 call getarg(1, fname1)
 call getarg(2, fname2)
 if(i == 4) then
  call getarg(3, buf)
  read(buf,*) idx(1)
  call getarg(4, buf)
  read(buf,*) idx(2)
 end if

 call solve_ON_matrix(fname1, fname2, idx)
end program main

subroutine solve_ON_matrix(fname1, fname2, idx)
 implicit none
 integer :: i, j, fid, idx1, idx2, nif, nbf, nmo
! nif: total number of MOs
! nbf: total number of basis functions
! nmo: the number of NOs, <= nif, e.g. in CAS(6,6), nmo=6 is enough
 integer, intent(in) :: idx(2)
 character(len=240) :: fname
 character(len=240), intent(in) :: fname1, fname2
 real(kind=8), allocatable :: mo(:,:), no(:,:), noon(:)
 real(kind=8), allocatable :: U(:,:), n(:,:), nU(:,:)
 logical :: gau

 gau = .false.
 i = index(fname1,'.fch', back=.true.)
 if(i /= 0) gau = .true. ! two Gaussian .fch(k) files

 if(gau) then
  call read_nbf_and_nif_from_fch(fname1, nbf, nif)
  allocate(mo(nbf,nif), no(nbf,nif), noon(nif))

  call read_mo_from_fch(fname1, nbf, nif, 'a', mo)
  call read_mo_from_fch(fname2, nbf, nif, 'a', no)
  call read_eigenvalues_from_fch(fname2, nif, 'a', noon)

 else ! two OpenMolcas Orb files
  call read_nbf_and_nif_from_orb(fname1, nbf, nif)
  allocate(mo(nbf,nif), no(nbf,nif), noon(nif))

  call read_mo_from_orb(fname1, nbf, nif, 'a', mo)
  call read_mo_from_orb(fname2, nbf, nif, 'a', no)
  call read_on_from_orb(fname2, nif, 'a', noon)
 end if

 if(idx(1) == 0) then
  idx1 = 1; idx2 = nif; nmo = nif
 else
  idx1 = idx(1); idx2 = idx(2)
  nmo = idx2 - idx1 + 1
  if(nmo < 2) then
   write(6,'(A)') 'ERROR in subroutine solve_ON_matrix: >=2 NOs required.'
   stop
  end if
 end if

 ! solve mo = no*U
 allocate(U(nmo,nmo), source=0d0)
 call get_u(nbf, nmo, no(:,idx1:idx2), mo(:,idx1:idx2), U)
 deallocate(no)

 ! check whether matrix U is unitary
 call check_unitary(nmo, U)

 allocate(n(nmo,nmo), source=0d0)
 forall(i=1:nmo) n(i,i) = noon(idx1+i-1)

 ! D = (U^T)nU
 allocate(nU(nmo,nmo), source=0d0)
 ! call dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 call dsymm('L', 'L', nmo, nmo, 1d0, n, nmo, U, nmo, 0d0, nU, nmo)
 n = 0d0
 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 call dgemm('T', 'N', nmo, nmo, nmo, 1d0, U, nmo, nU, nmo, 0d0, n, nmo)
 deallocate(U, nU)

 ! update the noon matrix
 forall(i=1:nmo) noon(idx1+i-1) = n(i,i)

 if(gau) then
  j = index(fname1,'.fch',back=.true.)
 else
  j = index(fname1,'.',back=.true.)
 end if
 fname = fname1(1:j-1)//'_D.txt'

 open(newunit=fid,file=TRIM(fname),status='replace')
 write(fid,'(A)') 'Occupation matrix:'
 do i = 1, nmo, 1
  do j = i, nmo, 1
   write(fid,'(2(I4,1X),F9.5)') j, i, n(j,i)
  end do ! for j
 end do ! for i
 close(fid)
 deallocate(n)

 if(gau) then ! Gaussian .fch(k) case
  call write_eigenvalues_to_fch(fname1, nif, 'a', noon, .true.)
 else ! Molcas/OpenMolcas case
  call write_on_to_orb(fname1, nif, 'a', noon, .true.)
 end if

 deallocate(mo, noon)
end subroutine solve_ON_matrix

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
end subroutine get_u

! check whether matrix U is unitary
subroutine check_unitary(nmo, U)
 implicit none
 integer :: i, j
 integer, intent(in) :: nmo
 real(kind=8), parameter :: thresh = 1.3d-6
 real(kind=8), intent(in) :: U(nmo,nmo)
 real(kind=8), allocatable :: UUT(:,:), UTU(:,:)

 allocate(UUT(nmo,nmo), source=0d0)
 call dgemm('N', 'T', nmo, nmo, nmo, 1d0, U, nmo, U, nmo, 0d0, UUT, nmo)
 forall(i = 1:nmo) UUT(i,i) = UUT(i,i) - 1d0

 do i = 1, nmo, 1
  do j = i, nmo, 1
   if(DABS(UUT(j,i)) > thresh) then
    write(6,'(A)') 'ERROR in subroutine check_unitary: input U is not unitary.'
    write(6,'(A)') 'Possible reasons: wrong indices or wrong MOs provided.'
    write(6,'(2(A,I0),A,F15.8)') 'j=', j, ', i=', i, ', UUT(j,i)=', UUT(j,i)
    stop
   end if
  end do ! for j
 end do ! for i

 deallocate(UUT)
 allocate(UTU(nmo,nmo), source=0d0)
 call dgemm('T', 'N', nmo, nmo, nmo, 1d0, U, nmo, U, nmo, 0d0, UTU, nmo)
 forall(i = 1:nmo) UTU(i,i) = UTU(i,i) - 1d0

 do i = 1, nmo, 1
  do j = i, nmo, 1
   if(DABS(UTU(j,i)) > thresh) then
    write(6,'(A)') 'ERROR in subroutine check_unitary: input U is not unitary.'
    write(6,'(A)') 'Possible reasons: wrong indices or wrong MOs provided.'
    write(6,'(2(A,I0),A,F15.8)') 'j=', j, ', i=', i, ', UTU(j,i)=', UTU(j,i)
    stop
   end if
  end do ! for j
 end do ! for i

 deallocate(UTU)
end subroutine check_unitary

