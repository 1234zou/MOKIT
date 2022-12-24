! written by jxzou at 20181012
! Note:
! 1) This is different from subroutine assoc_rot;
! 2) input integers 'rot1, rot2, ref1, ref2' are in Python/C convention
!    i.e., strart with [0], and [rot1,rot2), [ref1,ref2)

! update 20181020: support the number of virtual orbitals larger than that of occupied ones
! update 20181022: support either rotating virtual orbitals or occupied ones

! This subroutine/program will perform serial 2*2 rotations on active (un)occupied MOs by
!  maximizing the sum of orbital transition dipoles between an active occupied MO and
!  an active unoccupied MO.

!  For GNU compiler, use
! ------------------------------------
!  f2py -m assoc_loc -c assoc_loc.f90
! ------------------------------------
!  For INTEL compiler, use
! ---------------------------------------------------------------------------
!  f2py -m assoc_loc -c assoc_loc.f90 --fcompiler=intelem --compiler=intelem
! ---------------------------------------------------------------------------

! perform associated separated localization on a set of MOs
subroutine assoc_loc(nbf, nif, ref1, ref2, rot1, rot2, coeff, mo_dipole, new_coeff)
 implicit none
 integer i, j, k, niter
 integer nrot, nref
 integer nbf, nif, rot1, rot2, ref1, ref2
!f2py intent(in) :: nbf, nif, rot1, rot2, ref1, ref2
 ! nbf: the number of atomic basis functions
 ! nif: the number of independent basis functions, i.e., the number of MOs
 ! rot1: the begin index of orbitals to be rotated
 ! rot2: the end index of orbitals to be rotated
 ! ref1: the begin index of reference orbitals
 ! ref2: the end index of reference orbitals
 integer, parameter :: niter_max = 10000
 real(kind=8) coeff(nbf,nif), new_coeff(nbf,nif)
 real(kind=8) mo_dipole(3,nif,nif)
!f2py intent(in) :: coeff
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nif) coeff, new_coeff
!f2py depend(nif) mo_dipole
 ! coeff: all orbitals
 ! new_coeff: updated coeff
 ! mo_dipole: the whole MO basis dipole integrals matrix
 real(kind=8), parameter :: threshold1 = 1.0d-8, threshold2 = 1.0d-5
 ! threshold1: threshold to decide whether to rotate (and update MOs, dipole integrals)
 ! threshold2: threshold to decide if rotation/localization converged
 real(kind=8) rtmp, decrease
 real(kind=8) vtmp(3,4), motmp(nbf,2), diptmp(ref2-ref1,2)
 real(kind=8) Aij, Bij, cos_a, cos_theta, sin_theta
 real(kind=8), allocatable :: dipole(:,:,:)
 !                                  3,ref,rot
 ! dipole: the MO dipole integrals will be used, part of the matrix mo_dipole

 nrot = rot2 - rot1
 nref = ref2 - ref1

 if(nref > nrot) then
  write(6,'(A)') 'ERROR in subroutine assoc_loc: the number of reference orbitals is &
  &larger than that of rotated orbitals. Not allowed.'
  write(6,'(4(A,I3))') 'ref1=',ref1,', ref2=',ref2,', rot1=',rot1,', rot2=',rot2
  return
 end if

 allocate(dipole(3,nref,nrot))
 dipole = 0.0d0

 new_coeff = coeff
 do i = 1, nrot, 1
  do j = 1, nref, 1
   dipole(:,j,i) = mo_dipole(:,ref2-j+1,rot1+i)
  end do
 end do
 ! Note that Fortran conventional is considered in 'ref2-j+1, rot1+i'

 ! perform 2*2 rotation
 niter = 0
 do while(niter <= niter_max)
  decrease = 0.0d0
  do i = 1, nref, 1
   do j = i+1, nrot, 1
    if(j > nref) then
     vtmp(:,1) = dipole(:,i,i)
     vtmp(:,2) = dipole(:,i,j)
     Aij = DOT_PRODUCT(vtmp(:,1),vtmp(:,1)) - DOT_PRODUCT(vtmp(:,2),vtmp(:,2))
     Bij = DOT_PRODUCT(vtmp(:,1),vtmp(:,2))
    else
     vtmp(:,1) = dipole(:,i,i)
     vtmp(:,2) = dipole(:,j,j)
     vtmp(:,3) = dipole(:,i,j)
     vtmp(:,4) = dipole(:,j,i)
     Aij = DOT_PRODUCT(vtmp(:,1),vtmp(:,1)) + DOT_PRODUCT(vtmp(:,2),vtmp(:,2)) - &
         & DOT_PRODUCT(vtmp(:,3),vtmp(:,3)) - DOT_PRODUCT(vtmp(:,4),vtmp(:,4))
     Bij = DOT_PRODUCT(vtmp(:,1),vtmp(:,3)) - DOT_PRODUCT(vtmp(:,2),vtmp(:,4))
    end if
    Aij = 0.5d0*Aij
    rtmp = DSQRT(Aij*Aij + Bij*Bij)
    cos_a = Aij/rtmp
    cos_theta = DSQRT(0.5d0*(1.0d0+cos_a))
    if(DABS(1.0d0-cos_theta) < threshold1) cycle
    ! if theta is very close to zero, not to rotate
    decrease = decrease + Aij - rtmp
    rtmp = DSQRT(0.5d0*(1.0d0-cos_a))
    if(Bij > 0.0d0) then
     sin_theta = rtmp
    else if(Bij < 0.0d0) then
     sin_theta = -rtmp
    else
     cos_theta = 1.0d0
     sin_theta = 0.0d0
    end if
    !write(6,*) 'Aij=', Aij, 'Bij=', Bij
    !write(6,*) 'sin_theta=', sin_theta, 'cos_theta=', cos_theta
    ! update two orbitals
    motmp(:,1) = new_coeff(:,rot1+i)
    motmp(:,2) = new_coeff(:,rot1+j)
    new_coeff(:,rot1+i) = cos_theta*motmp(:,1) + sin_theta*motmp(:,2)
    new_coeff(:,rot1+j) = cos_theta*motmp(:,2) - sin_theta*motmp(:,1)
    ! update corresponding dipole integrals
    do k = 1, 3
     diptmp(:,1) = dipole(k,:,i)
     diptmp(:,2) = dipole(k,:,j)
     dipole(k,:,i) = cos_theta*diptmp(:,1) + sin_theta*diptmp(:,2)
     dipole(k,:,j) = cos_theta*diptmp(:,2) - sin_theta*diptmp(:,1)
    end do
   end do
  end do
  niter = niter + 1
  write(6,'(A,I5,A,F13.6)') 'niter=', niter, ', decrease=', decrease
  if(-decrease < threshold2) exit
  if(nref == 1) exit
 end do

 if(niter <= niter_max) then
  write(6,'(A)') 'Associated localization converged successfully.'
 else
  write(6,'(A,I5)') 'niter_max=', niter_max
  write(6,'(A)') 'Associated localization fails to converge.'
 end if

 deallocate(dipole)
end subroutine assoc_loc

