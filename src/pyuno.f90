! written by jxzou at 20180331
! generate UHF natural orbtials (UNO) from UHF canonical orbitals

! updated by jxzou at 20180420: to support nalpha > nbeta
! updated by jxzou at 20180520: modify the intent of array alpha_coeff from (inout) to (in,copy)
! updated by jxzou at 20180825: fix the bug when only 1 pair
! updated by jxzou at 20191215: delete SVD generation of virtual/inactive MOs (PAO outside recommanded)
! updated by jxzou at 20200426: add NOON output
! updated by jxzou at 20210518: add an intent(in) parameter ON_thres
! updated by jxzou at 20220711: change ON_thres to uno_thres

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

subroutine uno(nbf, nif, nalpha, nbeta, alpha_coeff, beta_coeff, ao_ovlp, &
               uno_thres, idx, noon, uno_coeff)
 implicit none
 integer :: i, outid
 integer :: ndb, nact, nact0, nopen, nocc
 ! ndb: the number of doubly occupied MOs
 ! nact: the number of active occupied orbitals
 ! nact = nact0 + nopen
 ! nocc = nalpha + nact0

 integer :: nbf, nif, nalpha, nbeta
!f2py intent(in) :: nbf, nif, nalpha, nbeta
 ! nbf: the number of basis functions
 ! nif: the number of independent functions, i.e., the number of MOs
 ! nalpha: the number of alpha electrons
 ! nbeta: the number of beta electrons

 ! array idx is in Fortran convention, i.e., begins from 1, not 0
 integer :: idx(3) ! 1: start of pair, 2: start of virtual, 3: nopen
!f2py intent(out) :: idx
 integer, allocatable :: idx1(:), idx2(:)

 real(kind=8) :: maxv, abs_mean
 real(kind=8), parameter :: ON_criteria = 1d-5
 real(kind=8) :: uno_thres ! UNO occupation number threshold
!f2py intent(in) :: uno_thres

 real(kind=8) :: ao_ovlp(nbf,nbf)
!f2py intent(in) :: ao_ovlp
!f2py depend(nbf) :: ao_ovlp

 real(kind=8) :: alpha_coeff(nbf,nif), beta_coeff(nbf,nif)
!f2py intent(in,copy) :: alpha_coeff
!f2py intent(in) :: beta_coeff
!f2py depend(nbf,nif) :: alpha_coeff, beta_coeff

 real(kind=8) :: noon(nif), uno_coeff(nbf,nif)
!f2py intent(out) :: noon, uno_coeff
!f2py depend(nif) :: noon
!f2py depend(nbf,nif) :: uno_coeff
 real(kind=8), allocatable :: mo_basis_ovlp(:,:), alpha_occ(:,:), beta_occ(:,:)
 real(kind=8), allocatable :: sv_occ0(:), sv_occ(:)
 character(len=7), parameter :: outname = 'uno.out'
 character(len=62), parameter :: on_warn1 = 'Warning in subroutine uno: uno_th&
                                            &res deviates from ON_criteria.'
 character(len=35), parameter :: on_warn2 = 'You better know what you are doing.'

 noon = 0d0
 nopen = nalpha - nbeta
 if(nopen < 0) then
  write(6,'(A)') 'ERROR in subroutine uno: nalpha < nbeta.'
  write(6,'(2(A,I0))') 'nalpha=', nalpha, ', nbeta=', nbeta
  stop
 end if

 open(newunit=outid,file=outname,status='replace')
 write(outid,'(A4,I5)') 'nbf=', nbf
 write(outid,'(A4,I5)') 'nif=', nif
 write(outid,'(A,F11.7)') 'ON_criteria=', ON_criteria
 write(outid,'(A,F11.7)') 'uno_thres=', uno_thres
 if(DABS(uno_thres - ON_criteria) > 1d-5) then
  write(outid,'(/,A)') on_warn1
  write(outid,'(A)') on_warn2
  write(6,'(/,A)') on_warn1
  write(6,'(A)') on_warn2
 end if

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

 if(nbeta == 0) then ! no beta electrons, return
  forall(i = 1:nopen) noon(i) = 1d0
  uno_coeff = alpha_coeff
  idx = [1, nopen+1, nopen]
  write(outid,'(/,A6,I5)') 'ndb  =', 0
  write(outid,'(A6,I5)')   'nact =', nalpha
  write(outid,'(A6,I5)')   'nact0=', 0
  write(outid,'(A6,3I5)')  'idx  =', idx
  close(outid)
  return
 end if
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
 write(6,'(/,A)') 'Singular values from SVD of Alpha/Beta MOs:'
 write(6,'(5(1X,ES15.8))') (sv_occ(i), i=1,nalpha)
 ! SVD done in occ space

 allocate(sv_occ0(nalpha), source=sv_occ)
 nact = COUNT(sv_occ < 1d0-ON_criteria)
 ndb = nalpha - nact
 nact0 = nact - nopen
 nocc = nalpha + nact0
 idx = [ndb+1, nocc+1, nopen]
 write(outid,'(/,A6,I5)') 'ndb  =', ndb
 write(outid,'(A6,I5)')   'nact =', nact
 write(outid,'(A6,I5)')   'nact0=', nact0
 write(outid,'(A6,3I5)')  'idx  =', idx

 ! generate NOON (Natural Orbital Occupation Number)
 forall(i = 1:nalpha) noon(i) = 1d0 + sv_occ(i)
 forall(i = 1:nact0) noon(nalpha+i) = 2d0 - noon(nbeta-i+1)
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
  alpha_coeff(:,idx1(i)) = (alpha_occ(:,idx1(i)) + beta_occ(:,idx1(i)))/DSQRT(2d0*noon(idx1(i)))
  alpha_coeff(:,idx2(i)) = (alpha_occ(:,idx1(i)) - beta_occ(:,idx1(i)))/DSQRT(2d0*noon(idx2(i)))
 end forall
 deallocate(idx1, idx2, alpha_occ, beta_occ)
 ! done transform in occ space

 ! Set virtual MOs to be zero. They may be calculated in following PAO constructions
 if(nocc < nif) alpha_coeff(:,nocc+1:nif) = 0d0

 ! check the orthonormality of final Alpha MO
 allocate(mo_basis_ovlp(nocc,nocc))
 call get_mo_basis_ovlp(nocc, nocc, nbf, alpha_coeff(:,1:nocc), alpha_coeff(:,1:nocc), ao_ovlp, mo_basis_ovlp)
 call check_unity(nocc, mo_basis_ovlp, maxv, abs_mean)
 deallocate(mo_basis_ovlp)
 write(outid,'(/,A)') 'The orthonormality of final Alpha MO:'
 write(outid,'(A,F16.10)') 'maxv=', maxv
 write(outid,'(A,F16.10)') 'abs_mean=', abs_mean

 ! copy the UNO coefficients
 uno_coeff = alpha_coeff

 ! now let's update ndb, nact, nact0, idx according to input uno_thres
 if(DABS(uno_thres - ON_criteria) > 1d-5) then
  nact = COUNT(sv_occ0 < 1d0-uno_thres)
  ndb = nalpha - nact
  nact0 = nact - nopen
  nocc = nalpha + nact0
  idx = [ndb+1, nocc+1, nopen]
 end if

 write(outid,'(/,A6,I5)') 'ndb  =', ndb
 write(outid,'(A6,I5)')   'nact =', nact
 write(outid,'(A6,I5)')   'nact0=', nact0
 write(outid,'(A6,3I5)')  'idx  =', idx
 close(outid)
 deallocate(sv_occ0)
end subroutine uno

subroutine get_mo_basis_ovlp(na, nb, nbf, c_alpha, c_beta, ao_ovlp, mo_basis_ovlp)
 implicit none
 integer, intent(in) :: na, nb, nbf
 real(kind=8), intent(in) :: c_alpha(nbf,na), c_beta(nbf,nb)
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf)
 real(kind=8), intent(out) :: mo_basis_ovlp(na,nb)
 real(kind=8) :: s_c_beta(nbf, nb)

 s_c_beta = 0d0
 mo_basis_ovlp = 0d0

 ! ?symm: Computes a matrix-matrix product where one input matrix is symmetric
 ! ?gemm: Computes a matrix-matrix product with general matrices
 ! Syntax FORTRAN 77:
 ! call dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 call dsymm('L', 'U', nbf, nb, 1d0, ao_ovlp, nbf, c_beta, nbf, 0d0, s_c_beta, nbf)
 call dgemm('T', 'N', na, nb, nbf, 1d0, c_alpha, nbf, s_c_beta, nbf, 0d0, mo_basis_ovlp, na)
end subroutine get_mo_basis_ovlp

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

 sv = 0d0
 u = 0d0
 vt = 0d0

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
 allocate(u2(nbf,nalpha), source=0d0)
 allocate(vt2(nbf,nbeta), source=0d0)
 call dgemm('N', 'N', nbf, nalpha, nalpha, 1d0, c_alpha, nbf, u, nalpha, 0d0, u2, nbf)
 call dgemm('N', 'T', nbf, nbeta, nbeta, 1d0, c_beta, nbf, vt, nbeta, 0d0, vt2, nbf)

 c_alpha = u2
 c_beta = vt2
 deallocate(u2, vt2)
end subroutine svd_and_rotate

subroutine check_unity(n, a, maxv, abs_mean)
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 real(kind=8), intent(in) :: a(n,n)
 real(kind=8), intent(out) :: maxv, abs_mean
 real(kind=8) :: b(n,n)

 b = a
 forall(i = 1:n) b(i,i) = b(i,i) - 1d0
 forall(i=1:n, j=1:n) b(i,j) = DABS(b(i,j))

 maxv = MAXVAL(b)
 abs_mean = SUM(b)/DBLE(n*n)
end subroutine check_unity

! calculated the SVD singular values of overlap of two sets of MOs
! Note: the input coeff1 and coeff2 must have the same dimension (nbf,nif)
subroutine svd_of_two_mo(nbf, nif, coeff1, coeff2, S)
 implicit none
 integer :: i
 integer :: nbf, nif
!f2py intent(in) :: nbf, nif

 real(kind=8), intent(in) :: coeff1(nbf,nif), coeff2(nbf,nif), S(nbf,nbf)
!f2py intent(in) :: coeff1, coeff2, S
!f2py depend(nbf,nif) :: coeff1, coeff2
!f2py depend(nbf) :: S

 real(kind=8), allocatable :: mo_ovlp(:,:), u(:,:), vt(:,:), ev(:)

 allocate(mo_ovlp(nif,nif), source=0d0)
 call get_mo_basis_ovlp(nif, nif, nbf, coeff1, coeff2, S, mo_ovlp)

 allocate(u(nif,nif), vt(nif,nif), ev(nif))
 call svd_on_ovlp(nif, nif, mo_ovlp, u, vt, ev)
 deallocate(mo_ovlp, u, vt)

 write(6,'(/,A)') 'SVD analysis of two sets of MOs:'
 write(6,'(A,ES15.8)') 'The smallest singular value:', MINVAL(ev)
 i = COUNT(ev < 1d-1)
 write(6,'(A,I0)') 'Number of singular values< 0.1: ', i
 i = COUNT(ev < 1d-2)
 write(6,'(A,I0)') 'Number of singular values<0.01: ', i

 write(6,'(A)') 'All singular values:'
 write(6,'(5(1X,ES15.8))') (ev(i),i=1,nif)
 deallocate(ev)
end subroutine svd_of_two_mo

