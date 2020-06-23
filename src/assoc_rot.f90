! written by jxzou at 20180616

! This subroutine/program will first compute the U matrix between two sets of MO (coeff1 and lo_coeff1)
!  then apply U onto coeff2 and the result is stored in new_coeff2

!  For GNU compiler, use
! ------------------------------------------------------
!  f2py -m assoc_rot -c assoc_rot.f90 --link-lapack_opt
! ------------------------------------------------------
!  For INTEL compiler, use
! ---------------------------------------------------------------------------------------------
!  f2py -m assoc_rot -c assoc_rot.f90 --link-lapack_opt --fcompiler=intelem --compiler=intelem
! ---------------------------------------------------------------------------------------------

subroutine assoc_rot(nbf, nmo, coeff1, lo_coeff1, coeff2, new_coeff2)
 implicit none
 integer i, info
 integer nbf, nmo
!f2py intent(in) :: nbf, nmo
 integer, allocatable :: ipiv(:)
 real(kind=8) coeff1(nbf,nmo), lo_coeff1(nbf,nmo), coeff2(nbf,nmo)
 real(kind=8) new_coeff2(nbf,nmo)
!f2py intent(in,copy) :: coeff1, lo_coeff1
!f2py intent(in) :: coeff2
!f2py intent(out) :: new_coeff2
!f2py depend(nbf,nmo) coeff1, lo_coeff1, coeff2, new_coeff2

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
 allocate(ipiv(min(nbf,nmo)))
 ipiv = 0
 call dgetrf(nbf, nmo, coeff1, nbf, ipiv, info)
 call dgetrs('N', nmo, nmo, coeff1, nbf, ipiv, lo_coeff1, nbf, info)
 deallocate(ipiv)

 ! reverse the coeff2
 forall(i = 1:nmo)
  coeff1(:,i) = coeff2(:,nmo-i+1)
 end forall
 coeff2 = coeff1
 ! rotate the coeff2
 new_coeff2 = 0.0d0
 call dgemm('N', 'N', nbf, nmo, nmo, 1.0d0, coeff2, nbf, lo_coeff1, nbf, 0.0d0, new_coeff2, nbf)
 ! reverse the coeff2 again
 forall(i = 1:nmo)
  coeff1(:,i) = new_coeff2(:,nmo-i+1)
 end forall
 new_coeff2 = coeff1
 return
end subroutine assoc_rot

