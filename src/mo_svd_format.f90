! written by jxzou at 20200607
! This version of mo_svd read overlap, dipole matrix, etc from formatted files
!  (like the Gaussian .log/.out file with iop(3/33=1), OpenMolcas output file
!  with 'print' keyword in &SEWARD) and performs SVD on the overlap of two sets of MOs.

! Another version file mo_svd.f90 reads overlap matrix from memory when running PySCF,
! thus that file is compiled by f2py wrapper. While this file is compiled by gfortran/ifort.

! Three arguments are: orbital file1, orbital file2, and output file containing
! AO-basis overlap. If the first argument ends with .fch, this will be identified
! as the Gaussian case. Otherwise it will be identified as the OpenMolcas case.

program main
 implicit none
 integer :: i, idx1, idx2
 character(len=10) :: str
 character(len=240) :: fname1, fname2, ovlp_file

 i = iargc()
 if(i /= 5) then
  write(6,'(/,A)') ' ERROR in subroutine mo_svd: wrong command line arguments!'
  write(6,'(/,A)') ' Example 1(Gaussian): mo_svd a.fch b.fch a.log idx1 idx2'
  write(6,'(/,A,/)') ' Example 2(OpenMolcas): mo_svd a.INPORB b.RasOrb a.out idx1 idx2'
  stop
 end if

 fname1 = ' '
 fname2 = ' '
 ovlp_file = ' '
 call getarg(1, fname1)
 call getarg(2, fname2)
 call getarg(3, ovlp_file)
 call getarg(4, str)
 read(str,*) idx1
 call getarg(5, str)
 read(str,*) idx2

 call mo_svd(fname1, fname2, ovlp_file, idx1, idx2)
end program main

subroutine mo_svd(fname1, fname2, ovlp_file, idx1, idx2)
 implicit none
 integer :: i, j, nbf, nif, nmo
 integer, intent(in) :: idx1, idx2
 character(len=240), intent(in) :: fname1, fname2, ovlp_file
 real(kind=8), allocatable :: S(:,:), mo_ovlp(:,:), SC(:,:)
 real(kind=8), allocatable :: coeff1(:,:), coeff2(:,:)
 real(kind=8), allocatable :: u(:,:), vt(:,:), ev(:)
 logical :: gau

 i = index(fname1, '.txt', back=.true.)
 if(i /= 0) then ! Gaussian case
  gau = .true.
  call read_nbf_and_nif_from_fch(fname1, nbf, nif)
 else ! i == 0, OpenMolcas case
  gau = .false.
  call read_nbf_and_nif_from_orb(fname1, nbf, nif)
 end if

 nmo = idx2 - idx1 + 1
 allocate(coeff1(nbf,nif), source=0d0)
 allocate(coeff2(nbf,nif), source=0d0)
 allocate(S(nbf,nbf), source=0d0)

 if(gau) then ! Gaussian case
  call read_mo_from_fch(fname1, nbf, nif, 'a', coeff1)
  call read_mo_from_fch(fname2, nbf, nif, 'a', coeff2)
!  call read_mo_from_chk_txt(fname1, nbf, nif, 'a', coeff1)
!  call read_mo_from_chk_txt(fname2, nbf, nif, 'a', coeff2)
  call read_int1e_from_gau_log(ovlp_file, 1, nbf, S)
 else         ! OpenMolcas case
  call read_mo_from_orb(fname1, nbf, nif, 'a', coeff1)
  call read_mo_from_orb(fname2, nbf, nif, 'a', coeff2)
  call read_ovlp_from_molcas_out(ovlp_file, nbf, S)
 end if

 ! check orthonormality
 allocate(mo_ovlp(nif,nif), source=0d0)
 mo_ovlp = MATMUL(TRANSPOSE(coeff1), MATMUL(S,coeff1))
 write(*,'(A)') 'orthonormality of file1:'
 do i = 1, nif, 1
  do j = i, nif, 1
   if(j==i) mo_ovlp(j,i) = mo_ovlp(j,i) - 1d0
   if(DABS(mo_ovlp(j,i)) > 1d-6) write(6,'(2I6,F15.8)') j, i, mo_ovlp(j,i)
  end do ! for j
 end do ! for i
 mo_ovlp = MATMUL(TRANSPOSE(coeff2), MATMUL(S,coeff2))
 write(*,'(A)') 'orthonormality of file2:'
 do i = 1, nif, 1
  do j = i, nif, 1
   if(j==i) mo_ovlp(j,i) = mo_ovlp(j,i) - 1d0
   if(DABS(mo_ovlp(j,i)) > 1d-6) write(6,'(2I6,F15.8)') j, i, mo_ovlp(j,i)
  end do ! for j
 end do ! for i
 deallocate(mo_ovlp)
 ! check done

 ! compute the MO-basis overlap matrix (C1^T)SC2
 allocate(SC(nbf,nmo), source=0d0)
 call dsymm('L', 'U', nbf, nmo, 1d0, S, nbf, coeff2(:,idx1:idx2), nbf, 0d0, SC, nbf)
 deallocate(S, coeff2)
 allocate(mo_ovlp(nmo,nmo), source=0d0)
 call dgemm('T', 'N', nmo, nmo, nbf, 1d0, coeff1(:,idx1:idx2), nbf, SC, nbf, 0d0, mo_ovlp, nmo)
 deallocate(SC, coeff1)

 ! perform SVD on the MO-basis overlap matrix
 allocate(u(nmo,nmo), vt(nmo,nmo), ev(nmo))
 call svd_on_ovlp(nmo, nmo, mo_ovlp, u, vt, ev)
 deallocate(mo_ovlp, u, vt)

 write(6,'(/,A)') 'SVD analysis of two sets of MOs:'
 write(6,'(A,ES15.8)') 'The smallest singular value:', MINVAL(ev)
 i = COUNT(ev < 0.1d0)
 write(6,'(A,I0)') 'Number of singular values< 0.1: ', i
 i = COUNT(ev < 0.01d0)
 write(6,'(A,I0)') 'Number of singular values<0.01: ', i

 write(6,'(A)') 'All singular values:'
 write(6,'(5(1X,ES15.8))') (ev(i),i=1,nmo)
 deallocate(ev)
end subroutine mo_svd

! perform SVD on a given overlap matrix
subroutine svd_on_ovlp(m, n, a, u, vt, s)
 implicit none
 integer :: i, lwork
 integer, intent(in) :: m,n
 real(kind=8), intent(in) :: a(m,n)
 real(kind=8), intent(out) :: u(m,m), vt(n,n), s(m)
 real(kind=8), allocatable :: work(:)
 real(kind=8) :: a_copy(m,n)

 a_copy = a
 u = 0d0; vt = 0d0; s = 0d0

 ! ?gesvd: Computes the singular value decomposition of a general rectangular matrix
 ! Syntax FORTRAN 77:
 ! call dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
 lwork = -1
 allocate(work(1))
 work = 0
 call dgesvd('A', 'A', m, n, a_copy, m, s, u, m, vt, n, work, lwork, i)
 lwork = CEILING(work(1))
 deallocate(work)
 allocate(work(lwork))
 call dgesvd('A', 'A', m, n, a_copy, m, s, u, m, vt, n, work, lwork, i)

 deallocate(work)
 if(i /= 0) then
  write(6,'(A)') 'ERROR: info /= 0 in subroutine svd_on_ovlp. Please check why.'
  write(6,'(A5,I0)') 'info=',i
  stop
 end if
end subroutine svd_on_ovlp

