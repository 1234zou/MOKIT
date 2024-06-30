! written by jxzou at 20180616

module polyn_info
 implicit none
 real(kind=8), parameter :: zero = 1d-9
 real(kind=8), parameter :: alpha = 0.01d0
 real(kind=8), parameter :: threshold1 = 1d-9, threshold2 = 1d-6
! threshold1: threshold to decide whether to rotate (and update MOs, dipole integrals)
! threshold2: threshold to decide if rotation/localization converged
end module polyn_info

module ao2mo_arrays
 implicit none
 integer :: nbf, ndb, nvir, nuv, nls, nls2, nia, nsia
 integer, allocatable :: idxuv(:,:), idxls(:,:), idxls2(:,:), idxia(:,:), &
  idxsia(:,:)
 real(kind=8), parameter :: thres = 1d-12
 real(kind=8), allocatable :: x1(:,:,:,:), x2(:,:,:), x3(:,:), y1(:,:,:,:), &
  y2(:,:,:), y3(:,:)
 logical, allocatable :: small(:,:), small_dm(:,:)
end module ao2mo_arrays

! This subroutine/program will first compute the U matrix between two sets of MO
! (coeff1 and lo_coeff1), then apply U onto coeff2 and the result is stored in
! new_coeff2
subroutine assoc_rot(nbf, nmo, coeff1, lo_coeff1, coeff2, new_coeff2)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 integer, allocatable :: ipiv(:)
 real(kind=8) :: coeff1(nbf,nmo), lo_coeff1(nbf,nmo), coeff2(nbf,nmo)
!f2py intent(in,copy) :: coeff1, lo_coeff1
!f2py intent(in) :: coeff2
!f2py depend(nbf,nmo) :: coeff1, lo_coeff1, coeff2
 real(kind=8), intent(out) :: new_coeff2(nbf,nmo)
!f2py intent(out) :: new_coeff2
!f2py depend(nbf,nmo) :: new_coeff2

 ! find the unitary (orthogonal) matrix between coeff1 and lo_coeff1,
 !  where coeff1*U = lo_coeff1
 allocate(ipiv(min(nbf,nmo)), source=0)
 call dgetrf(nbf, nmo, coeff1, nbf, ipiv, i)
 call dgetrs('N', nmo, nmo, coeff1, nbf, ipiv, lo_coeff1, nbf, i)
 deallocate(ipiv)

 ! reverse the coeff2
 forall(i = 1:nmo) coeff1(:,i) = coeff2(:,nmo-i+1)
 coeff2 = coeff1
 ! rotate the coeff2
 new_coeff2 = 0d0
 call dgemm('N', 'N', nbf, nmo, nmo, 1d0, coeff2, nbf, lo_coeff1, nbf, 0d0, &
            new_coeff2, nbf)
 ! reverse the coeff2 again
 forall(i = 1:nmo) coeff1(:,i) = new_coeff2(:,nmo-i+1)

 new_coeff2 = coeff1
end subroutine assoc_rot

! written by jxzou at 20181012
! update 20181020: support the number of virtual orbitals larger than that of occupied ones
! update 20181022: support either rotating virtual orbitals or occupied ones

! This subroutine/program will perform serial 2*2 rotations on active (un)occupied MOs by
!  maximizing the sum of orbital transition dipoles between an active occupied MO and
!  an active unoccupied MO.
! Note:
! 1) This is different from subroutine assoc_rot;
! 2) input integers 'rot1, rot2, ref1, ref2' are in Python/C convention
!    i.e., strart with [0], and [rot1,rot2), [ref1,ref2)
! 3) orbitals obtained from this subroutine is delocalized, so I rename this subroutine 
!    to assoc_loc2 and write a new subroutine assoc_loc
! perform associated separated localization on a set of MOs
subroutine assoc_loc2(nbf, nif, ref1, ref2, rot1, rot2, coeff, mo_dipole, &
                      new_coeff)
 use polyn_info, only: threshold1, threshold2
 implicit none
 integer :: i, j, k, niter, nrot, nref
 integer, intent(in) ::  nbf, nif, ref1, ref2, rot1, rot2
!f2py intent(in) :: nbf, nif, ref1, ref2, rot1, rot2
 integer, parameter :: niter_max = 10000
 real(kind=8), intent(in) :: coeff(nbf,nif) ! input orbitals
!f2py intent(in) :: coeff
!f2py depend(nbf,nif) :: coeff
 real(kind=8), intent(out) :: new_coeff(nbf,nif) ! updated orbitals
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nif) :: new_coeff
 real(kind=8), intent(in) :: mo_dipole(3,nif,nif) ! MO basis dipole integrals
!f2py intent(in) :: mo_dipole
!f2py depend(nif) mo_dipole
 real(kind=8) :: rtmp, increase, vtmp(3,4), motmp(nbf,2), diptmp(ref2-ref1,2)
 real(kind=8) :: Aij, Bij, cos_a, cos_theta, sin_theta
 real(kind=8), allocatable :: dipole(:,:,:) ! size (3,ref,rot)

 nrot = rot2 - rot1
 nref = ref2 - ref1

 if(nref > nrot) then
  write(6,'(/,A)') 'ERROR in subroutine assoc_loc2: the number of reference or&
                   &bitals is larger'
  write(6,'(A)') 'than that of rotated orbitals. Not allowed.'
  write(6,'(4(A,I3))') 'ref1=',ref1,', ref2=',ref2,', rot1=',rot1,', rot2=',rot2
  stop
 end if

 new_coeff = coeff
 allocate(dipole(3,nref,nrot))
 forall(i=1:nrot,j=1:nref) dipole(:,j,i) = mo_dipole(:,ref2-j+1,rot1+i)
 ! Note that Fortran conventional is considered in 'ref2-j+1, rot1+i'

 ! perform 2*2 rotation
 niter = 0
 do while(niter <= niter_max)
  increase = 0d0
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
    cos_theta = DSQRT(0.5d0*(1d0+cos_a))
    if(DABS(1d0-cos_theta) < threshold1) cycle
    ! if theta is very close to 0, not to rotate

    increase = increase + rtmp - Aij
    rtmp = DSQRT(0.5d0*(1d0-cos_a))
    if(Bij > 0d0) then
     sin_theta = -rtmp
    else if(Bij < 0d0) then
     sin_theta = rtmp
    else
     sin_theta = 0d0; cos_theta = 1d0
    end if

    ! update two orbitals
    motmp(:,1) = new_coeff(:,rot1+i)
    motmp(:,2) = new_coeff(:,rot1+j)
    new_coeff(:,rot1+i) = cos_theta*motmp(:,1) - sin_theta*motmp(:,2)
    new_coeff(:,rot1+j) = sin_theta*motmp(:,1) + cos_theta*motmp(:,2)

    ! update corresponding dipole integrals
    do k = 1, 3
     diptmp(:,1) = dipole(k,:,i)
     diptmp(:,2) = dipole(k,:,j)
     dipole(k,:,i) = cos_theta*diptmp(:,1) - sin_theta*diptmp(:,2)
     dipole(k,:,j) = sin_theta*diptmp(:,1) + cos_theta*diptmp(:,2)
    end do ! for k
   end do ! for j
  end do ! for i

  niter = niter + 1
  write(6,'(A,I5,A,F13.6)') 'niter=', niter, ', increase=', increase
  if(increase < threshold2) exit
  if(nref == 1) exit
 end do ! for while

 if(niter <= niter_max) then
  write(6,'(A)') 'Associated localization converged successfully.'
 else
  write(6,'(A,I5)') 'niter_max=', niter_max
  write(6,'(A)') 'Associated localization fails to converge.'
 end if

 deallocate(dipole)
end subroutine assoc_loc2

! perform associated separated localization on a set of MOs by maximizing the
!  target function
!   f = sum_i(<i_ref|r|i_rot>^2) + alpha*(sum_i(<i_rot|r|i_rot>)^2)
!  which means large orbital transition dipole moments meanwhile localized
! nbf: the number of atomic basis functions
! nif: the number of independent basis functions, i.e. the number of MOs
! ref1: the begin index of reference orbitals
! ref2: the end index of reference orbitals
! rot1: the begin index of orbitals to be rotated
! rot2: the end index of orbitals to be rotated
! coeff: all MO coefficients of a molecule
! new_coeff: all MO coefficients with coeff(:,rot1+1:rot2) updated
! mo_dipole: MO-based dipole integrals of all MOs
subroutine assoc_loc(nbf, nif, ref1, ref2, rot1, rot2, coeff, mo_dipole, &
                     new_coeff)
 use polyn_info, only: alpha, threshold1, threshold2
 implicit none
 integer :: i, j, k, niter, nrot, nref
 integer, intent(in) ::  nbf, nif, ref1, ref2, rot1, rot2
!f2py intent(in) :: nbf, nif, ref1, ref2, rot1, rot2
 integer, parameter :: niter_max = 9999
 real(kind=8), intent(in) :: coeff(nbf,nif) ! input orbitals
!f2py intent(in) :: coeff
!f2py depend(nbf,nif) :: coeff
 real(kind=8), intent(out) :: new_coeff(nbf,nif) ! updated orbitals
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nif) :: new_coeff
 real(kind=8), intent(in) :: mo_dipole(3,nif,nif) ! MO basis dipole integrals
!f2py intent(in) :: mo_dipole
!f2py depend(nif) :: mo_dipole
 real(kind=8) :: increase, vt(3,6), motmp(nbf,2), sin_2t, cos_2t
 real(kind=8) :: Aij, Bij, Cij, Dij, cos_theta, sin_theta, cc, ss
 real(kind=8), allocatable :: diptmp(:,:), r(:,:,:), d(:,:,:)
                                      ! (3,nref,nrot) (3,nrot,nrot)
 nrot = rot2 - rot1
 nref = ref2 - ref1
 if(nref > nrot) then
  write(6,'(/,A)') 'ERROR in subroutine assoc_loc: the number of reference or&
                   &bitals is larger'
  write(6,'(A)') 'than that of rotated orbitals. Not allowed.'
  write(6,'(4(A,I3))') 'ref1=',ref1,', ref2=',ref2,', rot1=',rot1,', rot2=',rot2
  stop
 end if

 new_coeff = coeff
 allocate(r(3,nref,nrot), d(3,nrot,nrot))
 d = mo_dipole(:,rot1+1:rot2,rot1+1:rot2)
 forall(i=1:nrot,j=1:nref) r(:,j,i) = mo_dipole(:,ref2-j+1,rot1+i)
 ! Note that Fortran conventional is considered in 'ref2-j+1, rot1+i'

 ! Perform 2*2 rotation. Note: r is not symmetric, d is symmetric
 niter = 0
 do while(niter <= niter_max)
  increase = 0d0
  do i = 1, nref, 1
   do j = i+1, nrot, 1
    if(j > nref) then
     vt(:,1) = r(:,i,i); vt(:,2) = r(:,i,j); vt(:,3) = d(:,i,i)
     vt(:,4) = d(:,j,i); vt(:,5) = d(:,j,j); vt(:,6) = vt(:,3)-vt(:,5)
     Aij = DOT_PRODUCT(vt(:,1),vt(:,1)) - DOT_PRODUCT(vt(:,2),vt(:,2)) + &
           alpha*(DOT_PRODUCT(vt(:,3),vt(:,3)) - DOT_PRODUCT(vt(:,5),vt(:,5)))
     Bij = DOT_PRODUCT(vt(:,1),vt(:,2)) + alpha*DOT_PRODUCT(vt(:,3)+vt(:,5),vt(:,4))
     Cij = 0.125d0*alpha*(DOT_PRODUCT(vt(:,6),vt(:,6)) - 4d0*DOT_PRODUCT(vt(:,4),vt(:,4)))
     Dij = 0.5d0*alpha*DOT_PRODUCT(vt(:,6),vt(:,4))
    else
     vt(:,1) = r(:,i,i); vt(:,2) = r(:,j,j); vt(:,3) = r(:,i,j)
     vt(:,4) = r(:,j,i); vt(:,5) = d(:,j,i); vt(:,6) = d(:,i,i)-d(:,j,j)
     Aij = DOT_PRODUCT(vt(:,1),vt(:,1)) + DOT_PRODUCT(vt(:,2),vt(:,2)) - &
           DOT_PRODUCT(vt(:,3),vt(:,3)) - DOT_PRODUCT(vt(:,4),vt(:,4))
     Bij = DOT_PRODUCT(vt(:,1),vt(:,3)) - DOT_PRODUCT(vt(:,2),vt(:,4))
     Cij = 0.25d0*alpha*(DOT_PRODUCT(vt(:,6),vt(:,6)) - 4d0*DOT_PRODUCT(vt(:,5),vt(:,5)))
     Dij = alpha*DOT_PRODUCT(vt(:,5),vt(:,6))
    end if

    Aij = 0.5d0*Aij
    call find_cos_quartic_poly_maximum(Aij,Bij,Cij,Dij,cos_theta,sin_theta, cc)
    if(cc < threshold1) cycle
    ! if the function change is too tiny, not to rotate, no matter what are cos_x
    ! or sin_x
    increase = increase + cc

    ! update two orbitals
    motmp(:,1) = new_coeff(:,rot1+i)
    motmp(:,2) = new_coeff(:,rot1+j)
    new_coeff(:,rot1+i) = cos_theta*motmp(:,1) - sin_theta*motmp(:,2)
    new_coeff(:,rot1+j) = sin_theta*motmp(:,1) + cos_theta*motmp(:,2)

    ! update dipole integrals array r
    allocate(diptmp(ref2-ref1,2), source=0d0)
    do k = 1, 3
     diptmp(:,1) = r(k,:,i)
     diptmp(:,2) = r(k,:,j)
     r(k,:,i) = cos_theta*diptmp(:,1) - sin_theta*diptmp(:,2)
     r(k,:,j) = sin_theta*diptmp(:,1) + cos_theta*diptmp(:,2)
    end do ! for k
    deallocate(diptmp)

    ! update dipole integrals array d
    cc = cos_theta*cos_theta
    ss = sin_theta*sin_theta
    sin_2t = 2d0*sin_theta*cos_theta
    cos_2t = cc - ss
    allocate(diptmp(3,3), source=0d0)
    diptmp(:,1) = d(:,i,i)
    diptmp(:,2) = d(:,j,i)
    diptmp(:,3) = d(:,j,j)
    d(:,i,i) = cc*diptmp(:,1) + ss*diptmp(:,3) - sin_2t*diptmp(:,2)
    d(:,j,j) = ss*diptmp(:,1) + cc*diptmp(:,3) + sin_2t*diptmp(:,2)
    d(:,j,i) = cos_2t*diptmp(:,2) + 0.5d0*sin_2t*vt(:,6)
    d(:,i,j) = d(:,j,i)
    do k = 1, nrot, 1
     if(k==i .or. k==j) cycle
     diptmp(:,1) = d(:,k,i)
     diptmp(:,2) = d(:,k,j)
     d(:,k,i) = cos_theta*diptmp(:,1) - sin_theta*diptmp(:,2)
     d(:,k,j) = sin_theta*diptmp(:,1) + cos_theta*diptmp(:,2)
    end do ! for k

    deallocate(diptmp)
   end do ! for j
  end do ! for i

  niter = niter + 1
  write(6,'(A,I4,A,F13.6)') 'niter=', niter, ', increase=', increase
  if(increase < threshold2) exit
  if(nref == 1) exit
 end do ! for while

 if(niter <= niter_max) then
  write(6,'(A)') 'Associated localization converged successfully.'
 else
  write(6,'(A,I5)') 'niter_max=', niter_max
  write(6,'(A)') 'Associated localization fails to converge.'
 end if

 deallocate(r, d)
end subroutine assoc_loc

! find maximum of f(x) = Acos2x - Bsin2x + Ccos4x - Dsin4x - A - C
subroutine find_cos_quartic_poly_maximum(a, b, c, d, cos_x, sin_x, y)
 implicit none
 integer :: i, nroot
 real(kind=8) :: a1, b1, c1, d1, e1, c2, s2, c4, s4, c22, y1, y2, root(4)
 real(kind=8), intent(in) :: a, b, c, d
 real(kind=8), intent(out) :: cos_x, sin_x, y
 logical :: valid_root

 cos_x = 1d0; sin_x = 0d0   ! initialization
 y = -(2d0*DABS(A) + DABS(B) + 2d0*DABS(C) + DABS(D))
 ! a negative value which is definitely < f(x)_min

 ! f'(x)=0 => -2Asin2x - 2Bcos2x - 4Csin4x - 4Dcos4x = 0
 !        i.e.  Asin2x + Bcos2x + 2Csin4x + 2Dcos4x = 0
 ! Let t = cos2x, we can obtain the folowing equation
 ! 16(C^2+D^2)t^4 + 8(AC+BD)t^3 + (A^2+B^2-16C^2-16D^2)t^2 - 4(2AC+BD)t +
 ! 4D^2-A^2 = 0
 ! Note here -1 <= t <= 1
 a1 = 16d0*(c*c + d*d)
 b1 = 8d0*(a*c + b*d)
 c1 = a*a + b*b - a1
 d1 = -4d0*(2d0*a*c + b*d)
 e1 = 4d0*d*d - a*a
 call general_quartic_solver(a1, b1, c1, d1, e1, nroot, root)
 ! The obtained roots are in ascending order

 ! This quartic equation may have some root beyond [-1,1], which should be discarded
 call discard_root(nroot, root)
 valid_root = .false.

 do i = 1, nroot, 1
  ! It is possible that c2 is slightly larger thant 1.0 due to numerical error,
  ! so we need MIN() and MAX() here.
  c2 = MAX(-1d0, MIN(root(i),1d0)) ! cos2x
  c22 = c2*c2                      ! (cos2x)^2
  c4 = 2d0*c22 - 1d0               ! cos4x
  s2 = DSQRT(1d0 - c22)            ! sin2x
  ! determine the sign of sin2x
  y1 = 2d0*d - 4d0*d*c22 - b*c2
  y2 = a + 4d0*c*c2
  if((y1<0d0 .and. y2>0d0) .or. (y1>0d0 .and. y2<0d0)) s2 = -s2
  ! done determine
  s4 = 2d0*c2*s2        ! sin4x
  y1 = a*c2 - b*s2 + c*c4 - d*s4 - a - c
  if(y1 < 0d0) cycle
  if(y1 > y) then
   valid_root = .true.
   y = y1
   cos_x = DSQRT(0.5d0*(1d0+c2))
   sin_x = DSQRT(0.5d0*(1d0-c2))
   if(s2 < 0d0) sin_x = -sin_x
  end if
 end do ! for i

 ! If there is no valid root, we must reset y=0 since it was set to a large
 ! negative value when initialization
 if(.not. valid_root) y = 0d0
end subroutine find_cos_quartic_poly_maximum

! written by jxzou 20170328, modified at 20170505, modified again at 20170901
! merged into assoc_rot.f90 at 20230415

! Find roots of Ax^4 + Bx^3 + Cx^2 + Dx + E = 0. Return roots in ascending order
subroutine general_quartic_solver(a, b, c, d, e, nroot, root)
 use polyn_info, only: zero
 implicit none
 integer :: i, j, nroot_c, nroot_q1, nroot_q2
 integer, intent(out) :: nroot
 real(kind=8), intent(in) :: a, b, c, d, e
 real(kind=8), intent(out) :: root(4)
 real(kind=8) :: b_a, c_a, d_a, e_a, b_c, c_c, d_c, root_c(3), y
 real(kind=8) :: a_p, b_p, c_p, b_q1, b_q2, c_q1, c_q2
 real(kind=8) :: root_q1(2), root_q2(2), temp_value(4), tmpv
 logical :: alive(4)

 nroot = 0; root = 0d0   ! initialization
 alive = .false.
 if(DABS(a) < zero) alive(1) = .true.

 if(alive(1)) then   ! a = 0, degrade to Bx^3 + Cx^2 + Dx + E = 0
  call general_cubic_solver(b, c, d, e, nroot, root(1:3))
  return
 end if

 if(DABS(b) < zero) alive(2) = .true.
 if(DABS(d) < zero) alive(3) = .true.
 if(alive(2) .and. alive(3)) then ! b = d = 0, Ax^4 + Cx^2 + E = 0
  call general_biquadratic_solver(a, c, e, nroot, root)
  return
 end if

 ! coefficient transformation and solve a cubic equation
 b_a = b/a ; c_a = c/a ; d_a = d/a ; e_a = e/a
 b_c = -c_a
 c_c = b_a*d_a - 4d0*e_a
 d_c = 4d0*c_a*e_a - d_a*d_a - b_a*b_a*e_a
 call general_cubic_solver(1d0, b_c, c_c, d_c, nroot_c, root_c)
 y = root_c(1)
 ! done transform and solve cubic

 ! The cubic equation above may have multiple solutions, but only one non-zero
 ! root is needed. We are looking for a proper transformation, not for real
 ! roots of the original quartic equation. Such transformation may be more than
 ! 1, but only one is needed.

 ! coefficient transformation and solve two quadratic equations
 a_p = 0.25d0*b_a*b_a - c_a + y
 b_p = 0.5d0*b_a*y - d_a
 alive = .false.
 if(DABS(a_p) < zero) alive(1) = .true.
 if(DABS(b_p) < zero) alive(2) = .true.

 if(alive(1) .and. alive(2)) then ! a_p = b_p = 0
  c_p = 0.25d0*y*y-e
  if(c_p < -zero) then
   nroot = 0
   return
  end if
  nroot = 4
  tmpv = b_a*b_a - 8d0*(y + 2d0*c_p)
  if(tmpv < 0d0) then
   root(1) = 0d0
  else
   root(1) = DSQRT(tmpv)
  end if
  root(2) = -root(1)
  tmpv = b_a*b_a - 8d0*(y - 2d0*c_p)
  if(tmpv < 0d0) then
   root(3) = 0d0
  else
   root(3) = DSQRT(tmpv)
  end if
  root(4) = -root(3)
  root = 0.25d0*(root - b_a)
 else
  ! a_p /= 0 or b_p /= 0
  temp_value(1) = DSQRT(DABS(a_p))
  temp_value(2) = DSQRT(DABS(0.25d0*y*y - e_a))
  if(b_p < 0d0) temp_value(2) = -temp_value(2)
  temp_value(3) = 0.5d0*b_a
  temp_value(4) = 0.5d0*y
  b_q1 = temp_value(3) - temp_value(1)
  b_q2 = temp_value(3) + temp_value(1)
  c_q1 = temp_value(4) - temp_value(2)
  c_q2 = temp_value(4) + temp_value(2)
  call general_quadratic_solver(1d0, b_q1, c_q1, nroot_q1, root_q1)
  call general_quadratic_solver(1d0, b_q2, c_q2, nroot_q2, root_q2)
  ! done transform and solve quadratic
  ! combine all roots
  nroot = nroot_q1 + nroot_q2
  if(nroot_q1 == 2) then
   root = [root_q1(1), root_q1(2), root_q2(1), root_q2(2)]
  else if(nroot_q1 == 1) then
   root(1:3) = [root_q1(1), root_q2(1), root_q2(2)]
  else
   root(1:2) = root_q2
  end if
 end if

 ! delete identical roots
 if(nroot > 1) then
  alive = .true.

  do i = 1, nroot-1, 1
   if(.not. alive(i)) cycle
   do j = i+1, nroot, 1
    if(.not. alive(j)) cycle
    if(DABS(root(j) - root(i)) < zero) then
     root(j) = 0d0
     alive(j) = .false.
    end if
   end do ! for j
  end do ! for i
  temp_value = 0d0
  j = 0
  do i = 1, nroot, 1
   if(.not. alive(i)) cycle
   j = j + 1
   temp_value(j) = root(i)
  end do ! for i
  root = temp_value
  nroot = j ! remember to update nroot
 end if

 call sort_darray(nroot, root)
end subroutine general_quartic_solver

! Find roots of Ax^3 + Bx^2 + Cx + D = 0. Return roots in ascending order
subroutine general_cubic_solver(a, b, c, d, nroot, root)
 use polyn_info, only: zero
 implicit none
 integer, intent(out) :: nroot
 real(kind=8), intent(in) :: a, b, c, d
 real(kind=8), intent(out) :: root(3)
 real(kind=8) :: p, q, a2, b2, delta, sqrt_delta
 real(kind=8) :: theta, cos_theta, temp_value, tmpv1, tmpv2
 real(kind=8), parameter :: PI = 4d0*DATAN(1d0)

 root = 0d0 ! initialization

 if(a<zero .and. a>-zero) then ! degrade to Bx^2 + Cx + D = 0
  call general_quadratic_solver(b, c, d, nroot, root(1:2))
  return
 end if

 ! coefficient transformation
 a2 = a**2
 b2 = b**2
 temp_value = b/(3d0*a)
 p = (3d0*a*c - b2)/(3d0*a2)
 q = (2d0*b2*b - 9d0*a*b*c + 27d0*a2*d)/(27d0*a2*a)
 ! transformation done

 delta = 0.25d0*(q*q) + (p**3)/27d0

 if(delta > zero) then
  nroot = 1
  sqrt_delta = DSQRT(delta)
  tmpv1 = - 0.5d0*q + sqrt_delta
  tmpv2 = - 0.5d0*q - sqrt_delta
  call dcurt(tmpv1)
  call dcurt(tmpv2)
  root(1) = tmpv1 + tmpv2 - temp_value
 else if(delta < -zero) then
  if(p > 0d0) p = 0d0
  nroot = 3
  cos_theta = -0.5d0*q*DSQRT(-27d0*p)/(p*p)
  cos_theta = MIN(MAX(-1d0,cos_theta),1d0) ! in case of value >1 or <-1
  theta = DACOS(cos_theta)/3
  q = 2d0*DSQRT(-p/3d0)
  cos_theta = 2d0*PI/3d0
  root(1) = q*DCOS(theta) - temp_value
  root(2) = q*DCOS(theta + cos_theta) - temp_value
  root(3) = q*DCOS(theta - cos_theta) - temp_value
  call sort_darray(3, root)
 else ! delta = 0
  nroot = 2
  tmpv1 = 0.5d0*q
  call dcurt(tmpv1)
  root(2) = tmpv1
  !root(2) = DSQRT(-p/3d0) ! this sometimes cause errors
  root(1) = -2d0*root(2) - temp_value
  root(2) = root(2) - temp_value
  ! delete identical roots
  if(DABS(root(1) - root(2)) <zero) then
   nroot = 1
   root(2) = 0d0
  else if(root(1) > root(2)) then
   tmpv1 = root(1); root(1) = root(2); root(2) = tmpv1
  end if
 end if
end subroutine general_cubic_solver

! Find roots of Ax^4 + Bx^2 + C = 0. Return roots in ascending order
subroutine general_biquadratic_solver(a, b, c, nroot, root)
 use polyn_info, only: zero
 implicit none
 integer, intent(out) :: nroot
 real(kind=8) :: r1, r2
 real(kind=8), intent(in) :: a, b, c
 real(kind=8), intent(out) :: root(4)

 nroot = 0; root = 0d0 ! initialization

 if(a<zero .and. a>-zero) then ! degrade to Bx^2 + C = 0
  call general_quadratic_solver(b, 0d0, c, nroot, root(1:2))
  return
 end if

 call general_quadratic_solver(a, b, c, nroot, root(1:2))

 if(nroot == 2) then
  if(root(1) > zero) then ! root(2)>root(1)>0
   nroot = 4
   r1 = DSQRT(root(1))
   r2 = DSQRT(root(2))
   root = [-r2, -r1, r1, r2]
  else if(root(1) < -zero) then ! root(1)<0
   if(root(2) > zero) then ! root(2)>0
    nroot = 2
    root(2) = DSQRT(root(2))
    root(1) = -root(2)
   else if(root(2) > -zero) then ! root(2)=0
    nroot = 1
   end if
  else ! root(1) = 0
   if(root(2) > zero) then ! root(2)>0
    nroot = 3
    r1 = DSQRT(root(2))
    root(1:3) = [-r1, 0d0, r1]
   else if(root(2) > -zero) then ! root(2)=0
    nroot = 1
   end if
  end if
 else if(nroot == 1) then
  if(root(1) > zero) then       ! root(1)>0
   nroot = 2
   root(2) = DSQRT(root(1))
   root(1) = -root(2)
  else if(root(1) < -zero) then ! root(1)<0
   nroot = 0; root = 0d0
  else
   nroot = 1                    ! root(1)=0
  end if
 end if
end subroutine general_biquadratic_solver

! Find roots of Ax^2 + Bx + C = 0. Return roots in ascending order
subroutine general_quadratic_solver(a, b, c, nroot, root)
 use polyn_info, only: zero
 implicit none
 integer, intent(out) :: nroot
 real(kind=8), intent(in) :: a, b, c
 real(kind=8), intent(out) :: root(2)
 real(kind=8) :: delta
 
 nroot = 0; root = 0d0   ! initialization

 if(a<zero .and. a>-zero) then ! a = 0
  if(b<zero .and. b>-zero) then ! b = 0
   if(c<zero .and. c>-zero) then ! c = 0
    nroot = 1
   else   ! c /= 0
    return
   end if
  else   ! b /= 0, Bx + C = 0
   nroot = 1
   root(1) = -c/b
  end if
 else   ! a /= 0, Ax^2 + Bx + C = 0
  delta = b*b - 4d0*a*c

  if(delta > zero) then       ! two different solutions
   nroot = 2
   delta = DSQRT(delta)
   if(a > 0d0) then
    root(1) = -0.5d0*(b + delta)/a
    root(2) = -0.5d0*(b - delta)/a
   else ! a < 0d0
    root(1) = -0.5d0*(b - delta)/a
    root(2) = -0.5d0*(b + delta)/a
   end if
  else if(delta > -zero) then ! two identical solutions
   nroot = 1
   root(1) = -0.5d0*b/a
  end if
  ! delta < -zero corresponds to no solution, i.e. nroot=0
 end if
end subroutine general_quadratic_solver

! double precision for cubic root
! In expression 'a**(1d0/3d0)', if a<0d0, it may cause errors in GNU compiler,
!  so calling subroutine dcurt is safer
subroutine dcurt(num)
 implicit none
 real(kind=8), intent(inout) :: num

 if(num > 0d0) then
  num = num**(1d0/3d0)
 else
  num = (-num)**(1d0/3d0)
  num = -num
 end if
end subroutine dcurt

! sort an double precision array
subroutine sort_darray(n, a)
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 real(kind=8) :: r
 real(kind=8), intent(inout) :: a(n)

 if(n < 2) return

 do i = 1, n-1, 1
  r = a(i)
  do j = i+1, n, 1
   if(a(j) < r) then
    r = a(j); a(j) = a(i); a(i) = r
   end if
  end do ! for j
 end do ! for i
end subroutine sort_darray

! discard elements which are beyond the range [-1,1] in a real(kind=8) array
subroutine discard_root(n, root)
 use polyn_info, only: zero
 implicit none
 integer :: i, j, k
 integer, intent(inout) :: n
 real(kind=8), intent(inout) :: root(n)
 real(kind=8), allocatable :: root1(:)
 logical, allocatable :: del(:)

 allocate(del(n))
 del = .false.
 ! use a threshold zero, not simply root(i)<-1d0, root(i)>1d0
 do i = 1, n, 1
  if(root(i)+1d0<-zero .or. root(i)-1d0>zero) del(i) = .true.
 end do ! for i

 k = COUNT(del .eqv. .true.)
 if(k == 0) then
  deallocate(del)
  return
 end if

 allocate(root1(n-k))
 j = 0
 do i = 1, n, 1
  if(.not. del(i)) then
   j = j + 1
   root1(j) = root(i)
  end if
 end do ! for i

 n = n - k
 root(1:n) = root1
 deallocate(del, root1)
end subroutine discard_root

! initialize compound-index arrays
subroutine init_idx_arrays()
 use ao2mo_arrays, only: nbf, ndb, nvir, nuv, nls, nls2, nia, nsia, idxuv, &
  idxls, idxls2, idxia, idxsia
 implicit none
 integer :: a, i, k, u, v, l, s

 nuv = nbf*nbf; nia = ndb*nvir; nsia = nbf*nia
 nls = nbf*(nbf+1)/2; nls2 = nbf*(nbf-1)/2
 allocate(idxuv(2,nuv),idxls(2,nls),idxls2(2,nls2),idxia(2,nia),idxsia(3,nsia))

!$omp parallel do schedule(static) default(private) shared(nbf,idxuv)
 do v = 1, nbf, 1
  k = (v-1)*nbf
  forall(u = 1:nbf) idxuv(:,k+u) = [u,v]
 end do ! for v
!$omp end parallel do

!$omp parallel do schedule(dynamic) default(private) shared(nbf,idxls)
 do s = 1, nbf, 1
  k = (s-1)*s/2
  forall(l = 1:s) idxls(:,k+l) = [l,s]
 end do ! for s
!$omp end parallel do

!$omp parallel do schedule(dynamic) default(private) shared(nbf,idxls2)
 do s = 2, nbf, 1
  k = (s-2)*(s-1)/2
  forall(l = 1:s-1) idxls2(:,k+l) = [l,s]
 end do ! for s
!$omp end parallel do

!$omp parallel do schedule(dynamic) default(shared) private(u,k,v)
 do u = ndb, 1, -1
  k = (ndb-u)*nvir
  forall(v = 1:nvir) idxia(:,k+v) = [v,u]
 end do ! for u
!$omp end parallel do

!$omp parallel do schedule(dynamic) default(shared) private(s,i,k,a)
 do s = 1, nbf, 1
  do i = 1, ndb, 1
   k = (s-1)*nia + (i-1)*nvir
   forall(a = 1:nvir) idxsia(:,k+a) = [a,i,s]
  end do ! for i
 end do ! for s
!$omp end parallel do
end subroutine init_idx_arrays

subroutine free_ao2mo_arrays()
 use ao2mo_arrays, only: x1, y1, idxuv, idxls, idxls2, idxia, idxsia, small, &
  small_dm
 implicit none

 deallocate(x1, y1, idxuv, idxls, idxls2, idxia, idxsia, small, small_dm)
end subroutine free_ao2mo_arrays

! build the G matrix, i.e. the two-electron part of the Fock matrix
subroutine build_ao_g(nbf, dm, eri, ao_g)
 use ao2mo_arrays, only: nuv, nls, nls2, thres, idxuv, idxls, idxls2, small, &
  small_dm
 implicit none
 integer :: i, j, u, v, l, s
 integer, intent(in) :: nbf
 real(kind=8) :: r1, r2, r3
 real(kind=8), intent(in) :: dm(nbf,nbf), eri(nbf,nbf,nbf,nbf)
 real(kind=8), intent(out) :: ao_g(nbf,nbf)

 if(allocated(small)) then
!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(nls,nuv,idxls,idxuv,small,small_dm,dm,eri,ao_g)
  do i = 1, nls, 1
   l = idxls(1,i); s = idxls(2,i)
   r1 = 0d0
   do j = 1, nuv, 1
    if(small(j,i)) cycle
    u = idxuv(1,j); v = idxuv(2,j)
    if(small_dm(u,v)) cycle
    r1 = r1 + dm(u,v)*(eri(u,v,l,s) - 0.5d0*eri(u,s,v,l))
   end do ! for j
   ao_g(l,s) = r1
  end do ! for i
!$omp end parallel do

 else ! small is not allocated
  allocate(small(nuv,nls))
  small = .false.

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(nls,nuv,idxls,idxuv,small,dm,eri,ao_g)
  do i = 1, nls, 1
   l = idxls(1,i); s = idxls(2,i)
   r1 = 0d0
   do j = 1, nuv, 1
    u = idxuv(1,j); v = idxuv(2,j)
    r2 = eri(u,v,l,s) - 0.5d0*eri(u,s,v,l)
    if(DABS(r2) < thres) then
     small(j,i) = .true.
    else
     r3 = dm(u,v)
     if(DABS(r3) > thres) r1 = r1 + r3*r2
    end if
   end do ! for j
   ao_g(l,s) = r1
  end do ! for i
!$omp end parallel do
 end if

!$omp parallel do schedule(dynamic) default(shared) private(i,l,s)
 do i = 1, nls2, 1
  l = idxls2(1,i); s = idxls2(2,i)
  ao_g(s,l) = ao_g(l,s)
 end do ! for i
!$omp end parallel do
end subroutine build_ao_g

! (iv|ls) = X_1(v,i,l,s) = sum_u (C_ui)*(uv|ls)
subroutine ao2mo_x1(nbf, ndb, mo, eri, x1)
 use ao2mo_arrays, only: nls, nls2, idxls, idxls2
 implicit none
 integer :: i, l, s
 integer, intent(in) :: nbf, ndb
 real(kind=8), intent(in) :: mo(nbf,ndb), eri(nbf,nbf,nbf,nbf)
 real(kind=8), intent(out) :: x1(nbf,ndb,nbf,nbf)
 real(kind=8), allocatable :: mat1(:,:), mat2(:,:)

 allocate(mat1(nbf,nbf), mat2(ndb,nbf))

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(nbf,ndb,nls,idxls,mo,eri,x1)
  do i = 1, nls, 1
   l = idxls(1,i); s = idxls(2,i)
   mat1 = eri(:,:,l,s); mat2 = 0d0
   call dgemm('T','N',ndb,nbf,nbf,1d0,mo,nbf,mat1,nbf,0d0,mat2,ndb)
   x1(:,:,l,s) = TRANSPOSE(mat2)
  end do ! for i
!$omp end parallel do

 deallocate(mat1, mat2)

!$omp parallel do schedule(dynamic) default(shared) private(i,l,s)
 do i = 1, nls2, 1
  l = idxls2(1,i); s = idxls2(2,i)
  x1(:,:,s,l) = x1(:,:,l,s)
 end do ! for i
!$omp end parallel do
end subroutine ao2mo_x1

subroutine update_x1_y1(cos_t, sin_t, i, a)
 use ao2mo_arrays, only: x1, y1
 implicit none
 integer, intent(in) :: i, a
 real(kind=8) :: tan_t, t21c
 real(kind=8), intent(in) :: cos_t, sin_t

 x1(:,i,:,:) = cos_t*x1(:,i,:,:) - sin_t*y1(:,a,:,:)
 tan_t = sin_t/cos_t
 t21c = (tan_t*tan_t + 1d0)*cos_t
 y1(:,a,:,:) = tan_t*x1(:,i,:,:) + t21c*y1(:,a,:,:)
end subroutine update_x1_y1

! (ii|ls) = X_2(l,i,s) = sum_v (C_vi)*X_1(v,i,l,s)
subroutine ao2mo_x2(nbf, nmo, mo, x1, x2)
 use ao2mo_arrays, only: nls, nls2, idxls, idxls2
 implicit none
 integer :: i, j, l, s
 integer, intent(in) :: nbf, nmo
 real(kind=8) :: ddot
 real(kind=8), intent(in) :: mo(nbf,nmo), x1(nbf,nmo,nbf,nbf)
 real(kind=8), intent(out) :: x2(nbf,nmo,nbf)

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(nls,nbf,nmo,idxls,mo,x1,x2)
 do i = 1, nls, 1
  l = idxls(1,i); s = idxls(2,i)
   do j = 1, nmo, 1
    x2(l,j,s) = ddot(nbf, mo(:,j), 1, x1(:,j,l,s), 1)
   end do ! for j
 end do ! for i
!$omp end parallel do

!$omp parallel do schedule(dynamic) default(shared) private(i,l,s)
 do i = 1, nls2, 1
  l = idxls2(1,i); s = idxls2(2,i)
  x2(s,:,l) = x2(l,:,s)
 end do ! for i
!$omp end parallel do
end subroutine ao2mo_x2

! (ii|is) = X_3(s,i) = sum_l (C_li)*X_2(l,i,s)
subroutine ao2mo_x3(nbf, nmo, mo, x2, x3)
 implicit none
 integer :: i, s
 integer, intent(in) :: nbf, nmo
 real(kind=8) :: ddot
 real(kind=8), intent(in) :: mo(nbf,nmo), x2(nbf,nmo,nbf)
 real(kind=8), intent(out) :: x3(nbf,nmo)

!$omp parallel do schedule(dynamic) default(private) shared(nbf,nmo,x2,x3,mo)
 do s = 1, nbf, 1
  do i = 1, nmo, 1
   x3(s,i) = ddot(nbf, mo(:,i), 1, x2(:,i,s), 1)
  end do ! for i
 end do ! for s
!$omp end parallel do
end subroutine ao2mo_x3

! (ii|ii) = X_4(i) = sum_s (C_si)*X_3(s,i)
subroutine ao2mo_x4(nbf, nmo, mo, x3, x4)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nmo
 real(kind=8) :: ddot
 real(kind=8), intent(in) :: mo(nbf,nmo), x3(nbf,nmo)
 real(kind=8), intent(out) :: x4(nmo)

!$omp parallel do schedule(dynamic) default(private) shared(nbf,nmo,x3,x4,mo)
 do i = 1, nmo, 1
  x4(i) = ddot(nbf, mo(:,i), 1, x3(:,i), 1)
 end do ! for i
!$omp end parallel do
end subroutine ao2mo_x4

! (ij|ls) = X_2'(l,j,i,s) = sum_v (C_vj)*X_1(v,i,l,s)
subroutine ao2mo_x2p(nbf, nmo, nmo1, mo, x1, x2p)
 use ao2mo_arrays, only: nls, nls2, idxls, idxls2
 implicit none
 integer :: i, l, s
 integer, intent(in) :: nbf, nmo, nmo1
 real(kind=8), intent(in) :: mo(nbf,nmo1), x1(nbf,nmo,nbf,nbf)
 real(kind=8), intent(out) :: x2p(nbf,nmo1,nmo,nbf)
 real(kind=8), allocatable :: mat1(:,:), mat2(:,:)

 allocate(mat1(nbf,nmo), mat2(nmo1,nmo))

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(nls,nbf,nmo,nmo1,idxls,mo,x1,x2p)
 do i = 1, nls, 1
  l = idxls(1,i); s = idxls(2,i)
  mat1 = x1(:,:,l,s); mat2 = 0d0
  call dgemm('T','N',nmo1,nmo,nbf,1d0,mo,nbf,mat1,nbf,0d0,mat2,nmo1)
  x2p(l,:,:,s) = mat2
 end do ! for i
!$omp end parallel do

 deallocate(mat1, mat2)

!$omp parallel do schedule(dynamic) default(shared) private(i,l,s)
 do i = 1, nls2, 1
  l = idxls2(1,i); s = idxls2(2,i)
  x2p(s,:,:,l) = x2p(l,:,:,s)
 end do ! for i
!$omp end parallel do
end subroutine ao2mo_x2p

! (ii|js) = X_3'(s,j,i) = sum_l (C_lj)*X_2(l,i,s)
subroutine ao2mo_x3p(nbf, nmo, nmo1, mo, x2, x3p)
 implicit none
 integer :: s
 integer, intent(in) :: nbf, nmo, nmo1
 real(kind=8), intent(in) :: mo(nbf,nmo1), x2(nbf,nmo,nbf)
 real(kind=8), intent(out) :: x3p(nbf,nmo1,nmo)
 real(kind=8), allocatable :: mat1(:,:), mat2(:,:)

 allocate(mat1(nbf,nmo), mat2(nmo1,nmo))

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(nbf,nmo,nmo1,x2,x3p,mo)
 do s = 1, nbf, 1
  mat1 = x2(:,:,s); mat2 = 0d0
  call dgemm('T','N',nmo1,nmo,nbf,1d0,mo,nbf,mat1,nbf,0d0,mat2,nmo1)
  x3p(s,:,:) = mat2
 end do ! for s
!$omp end parallel do

 deallocate(mat1, mat2)
end subroutine ao2mo_x3p

! (ia|is) = X_3''(s,a,i) = sum_l (C_li)*X_2'(l,a,i,s)
subroutine ao2mo_x3pp(nbf, ndb, nvir, mo, x2p, x3pp)
 use ao2mo_arrays, only: nsia, idxsia
 implicit none
 integer :: a, i, k, s
 integer, intent(in) :: nbf, ndb, nvir
 real(kind=8) :: ddot
 real(kind=8), intent(in) :: mo(nbf,ndb), x2p(nbf,nvir,ndb,nbf)
 real(kind=8), intent(out) :: x3pp(nbf,nvir,ndb)

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(nbf,nsia,idxsia,mo,x2p,x3pp)
 do k = 1, nsia, 1
  a = idxsia(1,k); i = idxsia(2,k); s = idxsia(3,k)
  x3pp(s,a,i) = ddot(nbf, mo(:,i), 1, x2p(:,a,i,s), 1)
 end do ! for k
!$omp end parallel do
end subroutine ao2mo_x3pp

! (ii|ij) = X_4'(j,i) = sum_s (C_sj)*X_3(s,i)
subroutine ao2mo_x4p(nbf, nmo, nmo1, mo, x3, x4p)
 implicit none
 integer, intent(in) :: nbf, nmo, nmo1
 real(kind=8), intent(in) :: mo(nbf,nmo1), x3(nbf,nmo)
 real(kind=8), intent(out) :: x4p(nmo1,nmo)

 x4p = 0d0
 call dgemm('T', 'N', nmo1, nmo, nbf, 1d0, mo, nbf, x3, nbf, 0d0, x4p, nmo1)
end subroutine ao2mo_x4p

! (ii|jj) = X_4''(j,i) = sum_s (C_sj)*X_3'(s,j,i)
subroutine ao2mo_x4pp(nbf, nmo, nmo1, mo, x3p, x4pp)
 implicit none
 integer :: i, j
 integer, intent(in) :: nbf, nmo, nmo1
 real(kind=8) :: ddot
 real(kind=8), intent(in) :: mo(nbf,nmo1), x3p(nbf,nmo1,nmo)
 real(kind=8), intent(out) :: x4pp(nmo1,nmo)

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(nbf,nmo,nmo1,x3p,x4pp,mo)
 do i = 1, nmo, 1
  do j = 1, nmo1, 1
   x4pp(j,i) = ddot(nbf, mo(:,j), 1, x3p(:,j,i), 1)
  end do ! for j
 end do ! for i
!$omp end parallel do
end subroutine ao2mo_x4pp

! (ij|ij) = X_4'''(j,i) = sum_s (C_sj)*X_3''(s,j,i)
subroutine ao2mo_x4ppp(nbf, nmo, nmo1, mo, x3pp, x4ppp)
 implicit none
 integer, intent(in) :: nbf, nmo, nmo1
 real(kind=8), intent(in) :: mo(nbf,nmo1), x3pp(nbf,nmo1,nmo)
 real(kind=8), intent(out) :: x4ppp(nmo1,nmo)

 call ao2mo_x4pp(nbf, nmo, nmo1, mo, x3pp, x4ppp)
end subroutine ao2mo_x4ppp

! calculate arrays g_iiia, g_iiaa, g_iaia, g_iaaa
subroutine ao2mo_ia1(nbf, ndb, nvir, mo, eri, g_iiii, g_iiia, g_iiaa, g_iaia, &
                     g_iaaa, g_aaaa)
 use ao2mo_arrays, only: x1, x2, x3, y1, y2, y3
 implicit none
 integer :: nif
 integer, intent(in) :: nbf, ndb, nvir
 real(kind=8), intent(in) :: mo(nbf,ndb+nvir), eri(nbf,nbf,nbf,nbf)
 real(kind=8), intent(out) :: g_iiii(ndb), g_iiia(nvir,ndb), g_iiaa(nvir,ndb), &
  g_iaia(nvir,ndb), g_iaaa(nvir,ndb), g_aaaa(nvir)
 real(kind=8), allocatable :: x2p(:,:,:,:),x3p(:,:,:),x3pp(:,:,:),g_iaaa_t(:,:)

 nif = ndb + nvir
 allocate(x1(nbf,ndb,nbf,nbf))
 call ao2mo_x1(nbf, ndb, mo(:,1:ndb), eri, x1)    ! (iv|ls)
 allocate(x2(nbf,ndb,nbf))
 call ao2mo_x2(nbf, ndb, mo(:,1:ndb), x1, x2)     ! (ii|ls)
 allocate(x3(nbf,ndb))
 call ao2mo_x3(nbf, ndb, mo(:,1:ndb), x2, x3)     ! (ii|is)
 call ao2mo_x4(nbf, ndb, mo(:,1:ndb), x3, g_iiii) ! (ii|ii)
 call ao2mo_x4p(nbf, ndb, nvir, mo(:,ndb+1:nif), x3, g_iiia) ! (ii|ia)
 deallocate(x3)
 allocate(x3p(nbf,nvir,ndb))
 call ao2mo_x3p(nbf, ndb, nvir, mo(:,ndb+1:nif), x2, x3p) ! (ii|as)
 deallocate(x2)
 call ao2mo_x4pp(nbf, ndb, nvir, mo(:,ndb+1:nif), x3p, g_iiaa) ! (ii|aa)
 deallocate(x3p)
 allocate(x2p(nbf,nvir,ndb,nbf))
 call ao2mo_x2p(nbf, ndb, nvir, mo(:,ndb+1:nif), x1, x2p) ! (ia|ls)
 allocate(x3pp(nbf,nvir,ndb))
 call ao2mo_x3pp(nbf, ndb, nvir, mo(:,1:ndb), x2p, x3pp) ! (ia|is)
 deallocate(x2p)
 call ao2mo_x4ppp(nbf, ndb, nvir, mo(:,ndb+1:nif), x3pp, g_iaia) ! (ia|ia)
 deallocate(x3pp)

 allocate(y1(nbf,nvir,nbf,nbf))
 call ao2mo_x1(nbf, nvir, mo(:,ndb+1:nif), eri, y1)    ! (av|ls)
 allocate(y2(nbf,nvir,nbf))
 call ao2mo_x2(nbf, nvir, mo(:,ndb+1:nif), y1, y2)     ! (aa|ls)
 allocate(y3(nbf,nvir))
 call ao2mo_x3(nbf, nvir, mo(:,ndb+1:nif), y2, y3)     ! (aa|as)
 deallocate(y2)
 call ao2mo_x4(nbf, nvir, mo(:,ndb+1:nif), y3, g_aaaa) ! (aa|aa)
 allocate(g_iaaa_t(ndb,nvir))
 call ao2mo_x4p(nbf, nvir, ndb, mo(:,1:ndb), y3, g_iaaa_t) ! (aa|ai)
 g_iaaa = TRANSPOSE(g_iaaa_t)
 deallocate(y3, g_iaaa_t)
end subroutine ao2mo_ia1

subroutine ao2mo_ia2(nbf, ndb, nvir, mo, g_iiia, g_iiaa, g_iaia, g_iaaa)
 use ao2mo_arrays, only: x1, x2, x3, y1, y2, y3
 implicit none
 integer :: nif
 integer, intent(in) :: nbf, ndb, nvir
 real(kind=8), intent(in) :: mo(nbf,ndb+nvir)
 real(kind=8), intent(out) :: g_iiia(nvir,ndb), g_iiaa(nvir,ndb), &
                              g_iaia(nvir,ndb), g_iaaa(nvir,ndb)
 real(kind=8), allocatable :: x2p(:,:,:,:),x3p(:,:,:),x3pp(:,:,:),g_iaaa_t(:,:)

 nif = ndb + nvir
 allocate(x2(nbf,ndb,nbf))
 call ao2mo_x2(nbf, ndb, mo(:,1:ndb), x1, x2) ! (ii|ls)
 allocate(x3(nbf,ndb))
 call ao2mo_x3(nbf, ndb, mo(:,1:ndb), x2, x3) ! (ii|is)
 call ao2mo_x4p(nbf, ndb, nvir, mo(:,ndb+1:nif), x3, g_iiia) ! (ii|ia)
 deallocate(x3)
 allocate(x3p(nbf,nvir,ndb))
 call ao2mo_x3p(nbf, ndb, nvir, mo(:,ndb+1:nif), x2, x3p)
 deallocate(x2)
 call ao2mo_x4pp(nbf, ndb, nvir, mo(:,ndb+1:nif), x3p, g_iiaa)
 deallocate(x3p)
 allocate(x2p(nbf,nvir,ndb,nbf))
 call ao2mo_x2p(nbf, ndb, nvir, mo(:,ndb+1:nif), x1, x2p)
 allocate(x3pp(nbf,nvir,ndb))
 call ao2mo_x3pp(nbf, ndb, nvir, mo(:,1:ndb), x2p, x3pp)
 deallocate(x2p)
 call ao2mo_x4ppp(nbf, ndb, nvir, mo(:,ndb+1:nif), x3pp, g_iaia)
 deallocate(x3pp)

 allocate(y2(nbf,nvir,nbf))
 call ao2mo_x2(nbf, nvir, mo(:,ndb+1:nif), y1, y2) ! (aa|ls)
 allocate(y3(nbf,nvir))
 call ao2mo_x3(nbf, nvir, mo(:,ndb+1:nif), y2, y3) ! (aa|as)
 deallocate(y2)
 allocate(g_iaaa_t(ndb,nvir))
 call ao2mo_x4p(nbf, nvir, ndb, mo(:,1:ndb), y3, g_iaaa_t) ! (aa|ai)
 g_iaaa = TRANSPOSE(g_iaaa_t)
 deallocate(y3, g_iaaa_t)
end subroutine ao2mo_ia2

! calculate the electronic energy using dm and Hcore+F
! Note: nuclear repulsion is not included in e.
subroutine calc_elec_c(nbf, dm, h_plus_f, e)
 implicit none
 integer :: u
 integer, intent(in) :: nbf
 real(kind=8) :: ddot
 real(kind=8), intent(in) :: dm(nbf,nbf), h_plus_f(nbf,nbf)
 real(kind=8), intent(out) :: e

 e = 0d0
!$omp parallel do schedule(dynamic) default(private) reduction(+:e) &
!$omp shared(nbf,dm,h_plus_f)
 do u = 1, nbf, 1
  e = e + ddot(nbf, dm(:,u), 1, h_plus_f(u,:), 1)
 end do ! for u
!$omp end parallel do
 e = 0.5d0*e
end subroutine calc_elec_c

! update only two elements: g_iiii(i) and g_aaaa(a)
subroutine update_g4i_g4a(cos_t, sin_t, v, g4i, g4a)
 implicit none
 integer :: i
 real(kind=8), intent(in) :: cos_t, sin_t, v(6)
 real(kind=8), intent(out) :: g4i, g4a
 real(kind=8), parameter :: cons(6) = [0.125d0,0.5d0,0.25d0,0.5d0,0.5d0,0.125d0]
 real(kind=8) :: ddot, cos2t, sin2t, cos4t, sin4t, cos2tp, sin2tp, cos4tp, r, &
  r1(6), r2(6)

 sin2t = 2d0*sin_t*cos_t; cos2t = 2d0*cos_t*cos_t - 1d0
 sin4t = 2d0*sin2t*cos2t; cos4t = 2d0*cos2t*cos2t - 1d0
 sin2tp = 2d0*sin2t; cos2tp = 4d0*cos2t
 cos4tp = cos4t + 3d0; r = 1d0 - cos4t

 r1 = [cos4tp+cos2tp,-sin4t-sin2tp,r,r,sin4t-sin2tp,cos4tp-cos2tp]
 r2 = [cos4tp-cos2tp,sin2tp-sin4t,r,r,sin4t+sin2tp,cos4tp+cos2tp]
 forall(i = 1:6)
  r1(i) = cons(i)*r1(i)
  r2(i) = cons(i)*r2(i)
 end forall
 g4i = ddot(6, r1, 1, v, 1)
 g4a = ddot(6, r2, 1, v, 1)
end subroutine update_g4i_g4a

! calculate the change of density matrix after the i-a rotation
subroutine calc_ddm(nbf, cos_t, sin_t, mo_i, mo_a, ddm)
 use ao2mo_arrays, only: nls2, thres, idxls2, small_dm
 implicit none
 integer :: i, u, v
 integer, intent(in) :: nbf
 real(kind=8) :: r
 real(kind=8), intent(in) :: cos_t, sin_t, mo_i(nbf), mo_a(nbf)
 real(kind=8), intent(out) :: ddm(nbf,nbf)

 small_dm = .false.; r = 2d0*sin_t

!$omp parallel do schedule(dynamic) default(shared) private(u,v)
 do u = 1, nbf, 1
  do v = 1, u, 1
   ddm(v,u) = r*( sin_t*(mo_a(u)*mo_a(v) - mo_i(u)*mo_i(v)) - &
                  cos_t*(mo_i(u)*mo_a(v) + mo_i(v)*mo_a(u)) )
   if(DABS(ddm(v,u)) < thres) small_dm(v,u) = .true.
  end do ! for v
 end do ! for u
!$omp end parallel do

!$omp parallel do schedule(dynamic) default(shared) private(i,v,u)
 do i = 1, nls2, 1
  v = idxls2(1,i); u = idxls2(2,i)
  ddm(u,v) = ddm(v,u)
  small_dm(u,v) = small_dm(v,u)
 end do ! for i
!$omp end parallel do
end subroutine calc_ddm

subroutine update_XA(ndb, nvir, i, a, g_iiii, g_aaaa, fii, faa, X, AA)
 use ao2mo_arrays, only: nia, idxia
 implicit none
 integer :: j, k, b
 integer, intent(in) :: ndb, nvir, i, a
 real(kind=8), intent(in) :: g_iiii(ndb), g_aaaa(nvir), fii(ndb,ndb), &
  faa(nvir,nvir)
 real(kind=8), intent(inout) :: X(nvir,ndb), AA(nvir,ndb)

!$omp parallel do schedule(dynamic) default(shared) private(j)
 do j = 1, ndb, 1
  X(a,j) = 0.5d0*(g_iiii(j) + g_aaaa(a))
 end do ! for j
!$omp end parallel do

!$omp parallel do schedule(dynamic) default(shared) private(b)
 do b = 1, nvir, 1
  X(b,i) = 0.5d0*(g_iiii(i) + g_aaaa(b))
 end do ! for b
!$omp end parallel do

!$omp parallel do schedule(dynamic) default(shared) private(k,b,j)
 do k = 1, nia, 1
  b = idxia(1,k); j = idxia(2,k)
  AA(b,j) = faa(b,b) - fii(j,j)
 end do ! for k
!$omp end parallel do
end subroutine update_XA

! RHF orbital optimization using the Jacobian 2-by-2 rotation method
! TODO: add DIIS
subroutine jacob22_rhf(ndb, nbf, nif, hcore, eri, mo, new_mo)
 use ao2mo_arrays, only: nbf0=>nbf, ndb0=>ndb, nvir0=>nvir, nia, idxia, small_dm
 implicit none
 integer :: i, k, a, nvir, niter
 integer(kind=4) :: hostnm
 integer, intent(in) :: ndb, nbf, nif
!f2py intent(in) :: ndb, nbf, nif
 integer, parameter :: niter_max = 999
 real(kind=8) :: thres, cos_t, sin_t, v(6), lower_e, inc1, inc2, elec_e
 real(kind=8), parameter :: thres1 = 1d-10, thres2 = 1d-8, thres3 = 1d-3 ! a.u.
 real(kind=8), intent(in) :: hcore(nbf,nbf), eri(nbf,nbf,nbf,nbf), mo(nbf,nif)
!f2py intent(in) :: hcore, eri, mo
!f2py depend(nbf) :: hcore, eri
!f2py depend(nbf,nif) :: mo
 real(kind=8), intent(out) :: new_mo(nbf,nif)
!f2py intent(out) :: new_mo
!f2py depend(nbf,nif) :: new_mo
 real(kind=8), allocatable :: dm(:,:), ddm(:,:), ao_g(:,:), dg(:,:), f(:,:), &
  fii(:,:), faa(:,:), fai(:,:), X(:,:), AA(:,:), BB(:,:), C(:,:), D(:,:)
 ! a is already used as an integer, here AA is used for the matrix A(a,i)
 ! ddm: P_new-P, the change of density matrix, whose elements are smaller
 ! dg: the change of matrix G
 real(kind=8), allocatable :: g_iiii(:), g_iiia(:,:), g_iiaa(:,:), g_iaia(:,:),&
  g_iaaa(:,:), g_aaaa(:), mo_i(:), mo_a(:)
 character(len=8) :: hostname
 character(len=24) :: data_string
 logical :: check_conv

 i = hostnm(hostname)
 call fdate(data_string)
 write(6,'(A)') 'HOST '//TRIM(hostname)//', '//TRIM(data_string)

 new_mo = mo
 nvir = nif - ndb
 write(6,'(2(A,I0))') 'ndb=', ndb, ', nvir=', nvir
 nbf0 = nbf; ndb0 = ndb; nvir0 = nvir
 call init_idx_arrays()

 ! calculate the density matrix
 allocate(dm(nbf,nbf), source=0d0)
 call dgemm('N','T',nbf,nbf,ndb,2d0,mo(:,1:ndb),nbf,mo(:,1:ndb),nbf,0d0,dm,nbf)

 write(6,'(/,A)') 'Building the AO Fock matrix...'
 ! calculate the G matrix
 allocate(ao_g(nbf,nbf))
 call build_ao_g(nbf, dm, eri, ao_g)

 ! calculate the AO-based Fock matrix
 allocate(f(nbf,nbf), source=hcore+ao_g)
 deallocate(ao_g)

 ! calculate the initial electronic energy
 call calc_elec_c(nbf, dm, hcore+f, elec_e)
 write(6,'(A,F21.9)') 'The initial electronic energy is: ', elec_e
 write(6,'(A)') 'Note: the nuclear repulsion energy is not included.'

 ! calculate the MO-based Fock matrix
 allocate(fii(ndb,ndb), faa(nvir,nvir), fai(nvir,ndb))
 call calc_CTSC(nbf, ndb, mo(:,1:ndb), f, fii)
 call calc_CTSC(nbf, nvir, mo(:,ndb+1:nif), f, faa)
 call calc_CTSCp2(nbf, nvir, ndb, mo(:,ndb+1:nif), f, mo(:,1:ndb), fai)

 ! AO -> MO integral transformation
 allocate(g_iiii(ndb), g_iiia(nvir,ndb), g_iiaa(nvir,ndb), g_iaia(nvir,ndb), &
          g_iaaa(nvir,ndb), g_aaaa(nvir))
 write(6,'(/,A)') 'AO->MO integral transformation...'
 call ao2mo_ia1(nbf, ndb, nvir, mo, eri, g_iiii, g_iiia, g_iiaa, g_iaia, &
                g_iaaa, g_aaaa)

 write(6,'(/,A)') 'Construct arrays X, A, B, C, D...'
 allocate(X(nvir,ndb), AA(nvir,ndb), BB(nvir,ndb), C(nvir,ndb), D(nvir,ndb))
 D = 0.5d0*(g_iaaa - g_iiia)
 BB = -2d0*(fai + D)

!$omp parallel do schedule(dynamic) default(shared) private(k,a,i)
 do k = 1, nia, 1
  a = idxia(1,k); i = idxia(2,k)
  X(a,i) = 0.5d0*(g_iiii(i) + g_aaaa(a))
  AA(a,i) = faa(a,a) - fii(i,i)
 end do ! for k
!$omp end parallel do

 AA = AA + g_iaia - 2d0*g_iiaa + X
 C = 0.25d0*(g_iiaa + 2d0*g_iaia - X)

 allocate(mo_i(nbf), mo_a(nbf), ddm(nbf,nbf), dg(nbf,nbf), small_dm(nbf,nbf))
 write(6,'(/,A)') 'Jacobian 2*2 orbital rotations...'

 niter = 0; inc1 = 0d0; thres = thres3; check_conv = .false.
 do while(niter <= niter_max)
  inc2 = 0d0

  do k = 1, nia, 1
   a = idxia(1,k); i = idxia(2,k)
   call find_cos_quartic_poly_maximum(AA(a,i),BB(a,i),C(a,i),D(a,i),cos_t,sin_t,lower_e)
   if(lower_e < thres) cycle
   inc2 = inc2 + lower_e
   ! update various matrices
   mo_i = new_mo(:,i); mo_a = new_mo(:,ndb+a)
   new_mo(:,i) = cos_t*mo_i - sin_t*mo_a
   new_mo(:,ndb+a) = sin_t*mo_i + cos_t*mo_a
   call calc_ddm(nbf, cos_t, sin_t, mo_i, mo_a, ddm)
   call build_ao_g(nbf, ddm, eri, dg) ! maybe time-consuming
   f = f + dg
   call calc_CTSC(nbf, ndb, new_mo(:,1:ndb), f, fii)
   call calc_CTSC(nbf, nvir, new_mo(:,ndb+1:nif), f, faa)
   call calc_CTSCp2(nbf, nvir, ndb, new_mo(:,ndb+1:nif), f, new_mo(:,1:ndb), fai)
   v = [g_iiii(i),g_iiia(a,i),g_iiaa(a,i),g_iaia(a,i),g_iaaa(a,i),g_aaaa(a)]
   call update_g4i_g4a(cos_t, sin_t, v, g_iiii(i), g_aaaa(a))
   call update_x1_y1(cos_t, sin_t, i, a)
   call ao2mo_ia2(nbf, ndb, nvir, new_mo, g_iiia, g_iiaa, g_iaia, g_iaaa)
   D = 0.5d0*(g_iaaa - g_iiia)
   BB = -2d0*(fai + D)
   call update_XA(ndb, nvir, i, a, g_iiii, g_aaaa, fii, faa, X, AA)
   AA = AA + g_iaia - 2d0*g_iiaa + X
   C = 0.25d0*(g_iiaa + 2d0*g_iaia - X)
  end do ! for k

  inc1 = inc1 + inc2
  niter = niter + 1
  write(6,'(A,I3,A,E7.1,A,F15.9)') 'niter=',niter,', thres=',thres,', E_lower=',inc2
  if(check_conv) then
   if(inc2 < thres2) exit
  else
   if(thres > 1.1d0*thres1) then
    thres = MIN(inc2, MAX(thres*0.1d0,inc2*1d-3))
    if(thres < thres1) then
     thres = thres1
     check_conv = .true.
    end if
   else
    check_conv = .true.
   end if
  end if
 end do ! for while

 deallocate(dm, ddm, f, fii, faa, fai, X, AA, BB, C, D, mo_i, mo_a, dg, g_iiii,&
            g_iiia, g_iiaa, g_iaia, g_iaaa, g_aaaa)
 call free_ao2mo_arrays()
 elec_e = elec_e - inc1
 write(6,'(A,F21.9)') 'The final electronic energy is: ', elec_e

 if(niter <= niter_max) then
  write(6,'(A)') 'RHF converged successfully.'
 else
  write(6,'(A,I5)') 'niter_max=', niter_max
  write(6,'(A)') 'RHF fails to converge.'
 end if
 call fdate(data_string)
 write(6,'(/,A)') TRIM(data_string)
end subroutine jacob22_rhf

