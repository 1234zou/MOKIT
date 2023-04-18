! written by jxzou at 20180616

module polyn_info
 implicit none
 real(kind=8), parameter :: zero = 1d-8
 real(kind=8), parameter :: alpha = -0.75d0
 real(kind=8), parameter :: threshold1 = 1d-8, threshold2 = 1d-5
! threshold1: threshold to decide whether to rotate (and update MOs, dipole integrals)
! threshold2: threshold to decide if rotation/localization converged
end module polyn_info

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
 real(kind=8) :: new_coeff2(nbf,nmo)
!f2py intent(in,copy) :: coeff1, lo_coeff1
!f2py intent(in) :: coeff2
!f2py intent(out) :: new_coeff2
!f2py depend(nbf,nmo) coeff1, lo_coeff1, coeff2, new_coeff2

 ! find the unitary (orthogonal) matrix between coeff1 and lo_coeff1,
 !  where coeff1*U = lo_coeff1
 allocate(ipiv(min(nbf,nmo)), source=0)
 call dgetrf(nbf, nmo, coeff1, nbf, ipiv, i)
 call dgetrs('N', nmo, nmo, coeff1, nbf, ipiv, lo_coeff1, nbf, i)
 deallocate(ipiv)

 ! reverse the coeff2
 forall(i = 1:nmo)
  coeff1(:,i) = coeff2(:,nmo-i+1)
 end forall
 coeff2 = coeff1
 ! rotate the coeff2
 new_coeff2 = 0d0
 call dgemm('N', 'N', nbf, nmo, nmo, 1d0, coeff2, nbf, lo_coeff1, nbf, 0d0, &
            new_coeff2, nbf)
 ! reverse the coeff2 again
 forall(i = 1:nmo)
  coeff1(:,i) = new_coeff2(:,nmo-i+1)
 end forall

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
 integer, intent(in) ::  nbf, nif, rot1, rot2, ref1, ref2
!f2py intent(in) :: nbf, nif, rot1, rot2, ref1, ref2
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
  return
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
! nif: the number of independent basis functions, i.e., the number of MOs
! rot1: the begin index of orbitals to be rotated
! rot2: the end index of orbitals to be rotated
! ref1: the begin index of reference orbitals
! ref2: the end index of reference orbitals
! coeff: all MO coefficients of a molecule
! new_coeff: all MO coefficients with coeff(:,rot1,rot2) updated
! mo_dipole: MO based dipole integrals of all MOs
subroutine assoc_loc(nbf, nif, ref1, ref2, rot1, rot2, coeff, mo_dipole, &
                     new_coeff)
 use polyn_info, only: alpha, threshold1, threshold2
 implicit none
 integer :: i, j, k, niter, nrot, nref
 integer, intent(in) ::  nbf, nif, rot1, rot2, ref1, ref2
!f2py intent(in) :: nbf, nif, rot1, rot2, ref1, ref2
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
 real(kind=8) :: increase, vt(3,6), motmp(nbf,2), sin_2t, cos_2t
 real(kind=8) :: Aij, Bij, Cij, Dij, cos_a, cos_theta, sin_theta, cc, ss
 real(kind=8), allocatable :: diptmp(:,:)
 real(kind=8), allocatable :: r(:,:,:) ! size (3,ref,rot)
 real(kind=8), allocatable :: d(:,:,:) ! size (3,rot,rot)

 nrot = rot2 - rot1
 nref = ref2 - ref1

 if(nref > nrot) then
  write(6,'(/,A)') 'ERROR in subroutine assoc_loc: the number of reference or&
                   &bitals is larger'
  write(6,'(A)') 'than that of rotated orbitals. Not allowed.'
  write(6,'(4(A,I3))') 'ref1=',ref1,', ref2=',ref2,', rot1=',rot1,', rot2=',rot2
  return
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
    call find_cos_quartic_poly_maximum(Aij, Bij, Cij, Dij, cos_theta, sin_theta, cc)
    increase = increase + cc
    if(DABS(1d0-cos_theta) < threshold1) cycle
    ! if theta is very close to 0, not to rotate

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
    d(:,j,i) = cos_2t*diptmp(:,2) - 0.5d0*sin_2t*vt(:,6)
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

 deallocate(r, d)
end subroutine assoc_loc

! find maximum of f(x) = Acos2x - Bsin2x + Ccos4x - Dsin4x - A - C
subroutine find_cos_quartic_poly_maximum(a, b, c, d, cos_x, sin_x, y)
 implicit none
 integer :: i, nroot
 real(kind=8) :: a1, b1, c1, d1, e1, c2, s2, c4, s4, c22, y1, y2, root(4)
 real(kind=8), intent(in) :: a, b, c, d
 real(kind=8), intent(out) :: cos_x, sin_x, y

 sin_x = 0d0; cos_x = 1d0; y = 0d0 ! initialization

 ! f'(x)=0 => -2Asin2x - 2Bcos2x - 4Csin4x - 4Dcos4x = 0
 !        i.e.  Asin2x + Bcos2x + 2Csin4x + 2Dcos4x = 0
 ! let t = cos2x, we can obtain the folowing equation
 ! 16(C^2+D^2)t^4 + 8(AC+BD)t^3 + (A^2+B^2-16C^2-16D^2)t^2 - 4(2AC+BD)t +
 ! 4D^2-A^2 = 0
 a1 = 16d0*(c*c + d*d)
 b1 = 8d0*(a*c + b*d)
 c1 = a*a + b*b - a1
 d1 = -4d0*(2d0*a*c + b*d)
 e1 = 4d0*d*d - a*a
 call general_quartic_solver(a1, b1, c1, d1, e1, nroot, root)

 do i = 1, nroot, 1
  c2 = root(i)          ! cos2x
  c22 = c2*c2           ! (cos2x)^2
  c4 = 2d0*c22 - 1d0    ! cos4x
  s2 = DSQRT(1d0 - c22) ! sin2x
  ! determine the sign of sin2x
  y1 = 2d0*d - 4d0*d*c22 - b*c2
  y2 = a + 4d0*c*c2
  if((y1<0d0 .and. y2>0d0) .or. (y1>0d0 .and. y2<0d0)) s2 = -s2
  ! done determine
  s4 = 2d0*c2*s2        ! sin4x
  y1 = a*c2 - b*s2 + c*c4 - d*s4 - a - c
  if(y1 < 0d0) cycle
  if(y1 > y) then
   y = y1
   cos_x = DSQRT(0.5d0*(1d0+c2))
   sin_x = DSQRT(0.5d0*(1d0-c2))
   if(s2 < 0d0) sin_x = -sin_x
  end if
 end do ! for i
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
 real(kind=8) :: b_a, c_a, d_a, e_a, a_c, b_c, c_c, d_c, root_c(3), y
 real(kind=8) :: a_p, b_p, c_p, a_q1, b_q1, b_q2, c_q1, c_q2
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
 a_c = 1d0
 b_c = -c_a
 c_c = b_a*d_a - 4d0*e_a
 d_c = 4d0*c_a*e_a - d_a*d_a - b_a*b_a*e_a
 call general_cubic_solver(a_c, b_c, c_c, d_c, nroot_c, root_c)
 y = root_c(1)
 ! done transform and solve cubic

 ! coefficient transformation and solve two quadratic equations
 a_p = 0.25d0*b_a*b_a - c_a + y
 b_p = 0.5d0*b_a*y - d_a
 alive = .false.
 if(a_p<zero .and. a_p>-zero) alive(1) = .true.
 if(b_p<zero .and. b_p>-zero) alive(2) = .true.

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
  ! a_p /= 0 and b_p /= 0
  a_q1 = 1d0
  temp_value(1) = DSQRT(DABS(a_p))
  temp_value(2) = 0.5d0*b_p*temp_value(1)/a_p
  temp_value(3) = 0.5d0*b_a
  temp_value(4) = 0.5d0*y
  b_q1 = temp_value(3) - temp_value(1)
  b_q2 = temp_value(3) + temp_value(1)
  c_q1 = temp_value(4) - temp_value(2)
  c_q2 = temp_value(4) + temp_value(2)
  call general_quadratic_solver(a_q1, b_q1, c_q1, nroot_q1, root_q1)
  call general_quadratic_solver(a_q1, b_q2, c_q2, nroot_q2, root_q2)
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
  cos_theta = min(max(-1d0,cos_theta),1d0) ! in case of value >1 or <-1
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
 real(kind=8) :: r1, r2, delta
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

