! written by jxzou at 20180331
! generate UHF natural orbtials (UNO) from UHF canonical orbitals

! updated by jxzou at 20180420: to support nalpha > nbeta
! updated by jxzou at 20180520: modify the intent of array alpha_coeff from (inout) to (in,copy)
! updated by jxzou at 20180825: fix the bug when only 1 pair
! updated by jxzou at 20191215: delete SVD generation of virtual/inactive MOs (PAO outside recommanded)
! updated by jxzou at 20200426: add NOON output

! This subroutine is designed to be imported as a module in Python.
!
! For GNU compiler, use
! --------------------------------------------
!  f2py -m uno -c pyuno.f90 --link-lapack_opt
! --------------------------------------------
!
! For Intel compiler, use
! -----------------------------------------------------------------------------------
!  f2py -m uno -c pyuno.f90 --link-lapack_opt --fcompiler=intelem --compiler=intelem
! -----------------------------------------------------------------------------------
!
!  to compile this file (a uno.so file will be generated). Then in Python
!  you can import the uno module.

subroutine uno(nbf, nif, nalpha, nbeta, alpha_coeff, beta_coeff, ao_ovlp, idx, noon, uno_coeff)
 implicit none
 integer :: i, outid
 integer :: ndb, nact, nact0, nopen, nocc
 ! ndb: the number of doubly occupied MOs
 ! nact: the number of active orbitals
 ! nact = nact0 + nopen
 ! nocc = nalpha + nact0

 integer :: nbf, nif, nalpha, nbeta
!f2py intent(in) :: nbf, nif, nalpha, nbeta
 ! nbf: the number of basis functions
 ! nif: the number of independent functions, i.e., the number of MOs
 ! nalpha: the number of alpha electrons
 ! nbeta: the number of beta electrons

 integer :: idx
!f2py intent(out) :: idx

 integer, allocatable :: idx1(:), idx2(:)

 real(kind=8), parameter :: on_criteria = 0.99999d0
 real(kind=8) :: maxv, abs_mean
 real(kind=8) :: ao_ovlp(nbf,nbf)
!f2py intent(in) :: ao_ovlp
!f2py depend(nbf) :: ao_ovlp

 real(kind=8) :: alpha_coeff(nbf,nif), beta_coeff(nbf,nif)
!f2py intent(in,copy) :: alpha_coeff
!f2py intent(in) :: beta_coeff
!f2py depend(nbf,nif) :: alpha_coeff, beta_coeff

 real(kind=8) :: uno_coeff(nbf,nif), noon(nif)
!f2py intent(out) :: uno_coeff, noon
!f2py depend(nbf,nif) :: uno_coeff
!f2py depend(nif) :: noon

 real(kind=8), allocatable :: mo_basis_ovlp(:,:)
 real(kind=8), allocatable :: alpha_occ(:,:), beta_occ(:,:)
 real(kind=8), allocatable :: sv_occ(:)
 character(len=7), parameter :: outname = 'uno.out'

 nopen = nalpha - nbeta
 if(nopen < 0) then
  write(6,'(A)') 'ERROR in subroutine uno: nalpha < nbeta.'
  stop
 end if

 open(newunit=outid,file=outname,status='replace')
 write(outid,'(A,F10.6)') 'ON_criteria=', on_criteria
 write(outid,'(A6,I6)') 'nopen=', nopen

 ! check the orthonormality of initial Alpha and Beta MO, respectively
 allocate(mo_basis_ovlp(nif,nif))
 call get_mo_basis_ovlp(nif, nif, nbf, alpha_coeff, alpha_coeff, ao_ovlp, mo_basis_ovlp)
 call check_unity(nif, mo_basis_ovlp, maxv, abs_mean)
 write(outid,'(/,A)') 'The orthonormality of initial Alpha MO:'
 write(outid,'(A,F16.10)') 'maxv=', maxv
 write(outid,'(A,F16.10)') 'abs_mean=', abs_mean

 call get_mo_basis_ovlp(nif, nif, nbf, beta_coeff, beta_coeff, ao_ovlp, mo_basis_ovlp)
 call check_unity(nif, mo_basis_ovlp, maxv, abs_mean)
 deallocate(mo_basis_ovlp)
 write(outid,'(/,A)') 'The orthonormality of initial Beta MO:'
 write(outid,'(A,F16.10)') 'maxv=', maxv
 write(outid,'(A,F16.10)') 'abs_mean=', abs_mean
 ! check orthonormality done

 ! allocate arrays for alpha and beta occupied orbitals, prepare for SVD
 allocate(alpha_occ(nbf,nalpha), source=alpha_coeff(:,1:nalpha))
 allocate(beta_occ(nbf,nbeta), source=beta_coeff(:,1:nbeta))

 ! calculate the overlap between alpha and beta occupied spatial orbitals
 allocate(mo_basis_ovlp(nalpha,nbeta))
 call get_mo_basis_ovlp(nalpha, nbeta, nbf, alpha_occ, beta_occ, ao_ovlp, mo_basis_ovlp)
 ! calculate done

 ! do SVD on the alpha_beta_ovlp of occupied spatial orbitals
 allocate(sv_occ(nalpha))
 call svd_and_rotate(nalpha, nbeta, nbf, alpha_occ, beta_occ, mo_basis_ovlp, sv_occ, .false.)
 deallocate(mo_basis_ovlp)
 ! SVD done in occ space

 nact = COUNT(sv_occ < on_criteria)
 ndb = nalpha - nact
 nact0 = nact - nopen
 nocc = nalpha + nact0
 idx = nocc + 1
 write(outid,'(/,A6,I6)') 'ndb  =', ndb
 write(outid,'(A6,I6)')   'nact =', nact
 write(outid,'(A6,I6)')   'nact0=', nact0
 write(outid,'(A6,I6)')   'idx  =', idx

 ! generate NOON (Natural Orbital Occupation Number)
 noon = 0.0d0
 forall(i = 1:nalpha) noon(i) = 1.0d0 + sv_occ(i)
 forall(i = 1:nact0) noon(nalpha+i) = 2.0d0 - noon(nbeta-i+1)
 deallocate(sv_occ)

 ! copy the doubly occupied MOs
 alpha_coeff(:,1:ndb) = alpha_occ(:,1:ndb)

 ! copy the singly occupied MO
 if(nopen > 0) alpha_coeff(:,nbeta+1:nalpha) = alpha_occ(:,nbeta+1:nalpha)

 ! transform the corresponding orbitals to UNOs in occ space
 allocate(idx1(nact0), idx2(nact0))
 forall(i = 1:nact0)
  idx1(i) = ndb + i
  idx2(i) = nocc + 1 - i
 end forall
 forall(i = 1:nact0)
  alpha_coeff(:,idx1(i)) = (alpha_occ(:,idx1(i)) + beta_occ(:,idx1(i)))/DSQRT(2.0d0*noon(idx1(i)))
  alpha_coeff(:,idx2(i)) = (alpha_occ(:,idx1(i)) - beta_occ(:,idx1(i)))/DSQRT(2.0d0*noon(idx2(i)))
 end forall
 deallocate(idx1, idx2, alpha_occ, beta_occ)
 ! done transform in occ space

 ! Set virtual MOs to be zero. They may be calculated in following PAO constructions
 if(nocc < nif) alpha_coeff(:,nocc+1:nif) = 0.0d0

 ! check the orthonormality of final Alpha MO
 allocate(mo_basis_ovlp(nocc,nocc))
 call get_mo_basis_ovlp(nocc, nocc, nbf, alpha_coeff(:,1:nocc), alpha_coeff(:,1:nocc), ao_ovlp, mo_basis_ovlp)
 call check_unity(nocc, mo_basis_ovlp, maxv, abs_mean)
 deallocate(mo_basis_ovlp)
 write(outid,'(/,A)') 'The orthonormality of final Alpha MO:'
 write(outid,'(A,F16.10)') 'maxv=', maxv
 write(outid,'(A,F16.10)') 'abs_mean=', abs_mean

 ! print NOON
 write(outid,'(/,A)') "NOON: (can be copied into 'Alpha Orbital Energies' section of .fchk)"
 write(outid,'(5(1X,ES15.8))') (noon(i),i=1,nif)
 close(outid)

 ! copy the UNO coefficients
 uno_coeff = alpha_coeff
 return
end subroutine uno

subroutine get_mo_basis_ovlp(na, nb, nbf, c_alpha, c_beta, ao_ovlp, mo_basis_ovlp)
 implicit none
 integer, intent(in) :: na, nb, nbf
 real(kind=8), intent(in) :: c_alpha(nbf,na), c_beta(nbf,nb)
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf)
 real(kind=8), intent(out) :: mo_basis_ovlp(na,nb)
 real(kind=8) :: s_c_beta(nbf, nb)

 s_c_beta = 0.0d0
 mo_basis_ovlp = 0.0d0

 ! ?symm: Computes a matrix-matrix product where one input matrix is symmetric
 ! ?gemm: Computes a matrix-matrix product with general matrices
 ! Syntax FORTRAN 77:
 ! call dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 call dsymm('L', 'U', nbf, nb, 1.0d0, ao_ovlp, nbf, c_beta, nbf, 0.0d0, s_c_beta, nbf)
 call dgemm('T', 'N', na, nb, nbf, 1.0d0, c_alpha, nbf, s_c_beta, nbf, 0.0d0, mo_basis_ovlp, na)
 return
end subroutine get_mo_basis_ovlp

subroutine svd_on_ovlp(m, n, a, u, vt, s)
 implicit none
 integer :: lwork, info
 integer, intent(in) :: m,n
 real(kind=8), intent(in) :: a(m,n)
 real(kind=8), intent(out) :: u(m,m), vt(n,n), s(m)
 real(kind=8), allocatable :: work(:)
 real(kind=8) :: a_copy(m,n)

 info = 0
 a_copy = a
 u = 0.0d0; vt = 0.0d0; s = 0.0d0

 ! ?gesvd: Computes the singular value decomposition of a general rectangular matrix
 ! Syntax FORTRAN 77:
 ! call dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
 lwork = -1
 allocate(work(1))
 work = 0
 call dgesvd('A', 'A', m, n, a_copy, m, s, u, m, vt, n, work, lwork, info)
 lwork = CEILING(work(1))
 deallocate(work)
 allocate(work(lwork))
 call dgesvd('A', 'A', m, n, a_copy, m, s, u, m, vt, n, work, lwork, info)

 deallocate(work)
 if(info /= 0) then
  write(6,'(A)') 'ERROR: info /= 0 in subroutine svd_on_ovlp. Please check why.'
  write(6,'(A5,I8)') 'info=',info
 end if
 return
end subroutine svd_on_ovlp

subroutine svd_and_rotate(nalpha,nbeta,nbf,c_alpha,c_beta,alpha_beta_ovlp,sv,reverse)
 implicit none
 integer :: i
 integer, intent(in) :: nalpha, nbeta, nbf
 real(kind=8), intent(inout) :: c_alpha(nbf,nalpha), c_beta(nbf,nbeta)
 real(kind=8), intent(in) :: alpha_beta_ovlp(nalpha,nbeta)
 real(kind=8), intent(out) :: sv(nalpha)
 real(kind=8) :: u(nalpha,nalpha), vt(nbeta,nbeta)
 real(kind=8), allocatable :: sv2(:), u2(:,:), vt2(:,:)
 logical, intent(in) :: reverse

 sv = 0.0d0
 u = 0.0d0
 vt = 0.0d0

 ! perform SVD on alpha_beta_ovlp of occupied or virtual orbitals
 call svd_on_ovlp(nalpha, nbeta, alpha_beta_ovlp, u, vt, sv)

 ! rotate occupied or virtual orbitals
 if(reverse) then
  allocate(sv2(nalpha))
  forall(i = 1:nalpha) sv2(i) = sv(nalpha-i+1)
  sv = sv2
  deallocate(sv2)
  allocate(u2(nalpha,nalpha), vt2(nbeta,nbeta))
  forall(i = 1:nalpha) u2(:,i) = u(:,nalpha-i+1)
  u = u2
  forall(i = 1:nbeta) vt2(i,:) = vt(nbeta-i+1,:)
  vt = vt2
  deallocate(u2, vt2)
 end if

 ! ?gemm: Computes a matrix-matrix product with general matrices
 ! Syntax FORTRAN 77:
 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 allocate(u2(nbf,nalpha), source=0.0d0)
 allocate(vt2(nbf,nbeta), source=0.0d0)
 call dgemm('N', 'N', nbf, nalpha, nalpha, 1.0d0, c_alpha, nbf, u, nalpha, 0.0d0, u2, nbf)
 call dgemm('N', 'T', nbf, nbeta, nbeta, 1.0d0, c_beta, nbf, vt, nbeta, 0.0d0, vt2, nbf)

 c_alpha = u2
 c_beta = vt2
 deallocate(u2, vt2)
 return
end subroutine svd_and_rotate

subroutine check_unity(n, a, maxv, abs_mean)
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 real(kind=8), intent(in) :: a(n,n)
 real(kind=8), intent(out) :: maxv, abs_mean
 real(kind=8) :: b(n,n)

 b = a
 forall(i = 1:n) b(i,i) = b(i,i) - 1.0d0
 forall(i=1:n, j=1:n) b(i,j) = DABS(b(i,j))

 maxv = MAXVAL(b)
 abs_mean = SUM(b)/DBLE(n*n)
 return
end subroutine check_unity

