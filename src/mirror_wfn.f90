! written by jxzou at 20230516: generate the wave function of the mirror of a
! molecule

module bas_rot
 implicit none
 real(kind=8), parameter :: root3 = DSQRT(3d0), root5 = DSQRT(5d0), &
  root10 = DSQRT(10d0), root14 = DSQRT(14d0), root15 = DSQRT(15d0), &
  root35 = DSQRT(35d0), root70 = DSQRT(70d0), root105 = DSQRT(105d0)
 real(kind=8), parameter :: vec7f(3,7) = RESHAPE([1d0,0d0,0d0, 0d0,0d0,1d0, &
  1d0,root3,1d0, 1d0,0d0,1d0, 0d0,1d0,1d0, 1d0,1d0,0d0, root3,1d0,root5], [3,7])
 real(kind=8), parameter :: vec9g(3,9) = RESHAPE([1d0,0d0,0d0, 0d0,0d0,1d0, &
  1d0,root3,1d0, 1d0,0d0,1d0, 0d0,1d0,1d0, 1d0,1d0,0d0, root3,1d0,root5, &
  root5,1d0,root3, root5,root3,1d0], [3,9])
 real(kind=8), parameter :: vec11h(3,11) = RESHAPE([1d0,0d0,0d0, 0d0,0d0,1d0, &
  1d0,root3,1d0, 1d0,0d0,1d0, 0d0,1d0,1d0, 1d0,1d0,0d0, root3,1d0,root5, &
  root5,1d0,root3, root5,root3,1d0, root3,10d0,-root14, root3,-root14,10d0], &
  [3,11])
 real(kind=8), allocatable :: A_7f(:,:), A_9g(:,:), A_11h(:,:), invA_7f(:,:), &
  invA_9g(:,:), invA_11h(:,:)
contains

! Hard-coding of invA_7f, invA_9g, or invA_11h is so tedious that typos will
!  probably occur, so here we solve the three inverse matrices numerically.
! The invA_5d is hard-coded in the subroutine get_rot_mat_5d.
subroutine get_invA(n)
 implicit none
 integer :: i
 integer, intent(in) :: n

 select case(n)
 case(7)
  allocate(A_7f(7,7), invA_7f(7,7))
  forall(i = 1:7) A_7f(:,i) = get7f_vector(vec7f(:,i))
  call inverse(7, A_7f, invA_7f)
  deallocate(A_7f)
 case(9)
  allocate(A_9g(9,9), invA_9g(9,9))
  forall(i = 1:9) A_9g(:,i) = get9g_vector(vec9g(:,i))
  call inverse(9, A_9g, invA_9g)
  deallocate(A_9g)
 case(11)
  allocate(A_11h(11,11), invA_11h(11,11))
  forall(i = 1:11) A_11h(:,i) = get11h_vector(vec11h(:,i))
  call inverse(11, A_11h, invA_11h)
  deallocate(A_11h)
 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine get_invA: invalid n=', n
  stop
 end select
end subroutine get_invA

! Input the rotation matrix of Cartesian coordinates, rotate the pre-defined five
!  3D basis vectors, calculate values of spherical harmonic functions using
!  rotated basis vectors, and calculate the rotation matrix for 5D, i.e. MA^(-1).
! Ref1: http://www.iro.umontreal.ca/~derek/files/SHRot_mine.pdf
! Ref2: http://filmicworlds.com/blog/simple-and-fast-spherical-harmonic-rotation
! Ref3: https://zhuanlan.zhihu.com/p/51267461?utm_id=0
subroutine get_rot_mat_5d(rotation, rot5d)
 implicit none
 integer :: i
 real(kind=8) :: vecs(3,5)
 real(kind=8), intent(in) :: rotation(3,3)
 real(kind=8), intent(out) :: rot5d(5,5)
 real(kind=8), parameter :: r = 1d0/DSQRT(2d0)
 real(kind=8), parameter :: s = 4d0*DSQRT(DATAN(1d0)/15d0)
 real(kind=8), parameter :: vecs0(3,5) = &
  RESHAPE([1d0,0d0,0d0, 0d0,0d0,1d0, r,r,0d0, r,0d0,r, 0d0,r,r], [3,5])
 ! The inverse of square matrix A is hard-coded here
 real(kind=8), parameter :: invA(5,5) = RESHAPE([0d0,root3,0d0,0d0,0d0, &
  -1d0,-1d0,0d0,2d0,0d0, 1d0,0d0,0d0,0d0,2d0, 2d0,1d0,0d0,0d0,0d0, &
  0d0,1d0,2d0,0d0,0d0], [5,5])

 vecs = MATMUL(TRANSPOSE(rotation), vecs0)

 ! calculate values one by one column
 forall(i = 1:5) rot5d(:,i) = get5d_vector(vecs(:,i))

 rot5d = MATMUL(rot5d, s*invA)
end subroutine get_rot_mat_5d

subroutine get_rot_mat_6d(rotation, rot6d)
 implicit none
 real(kind=8), intent(in) :: rotation(3,3)
 real(kind=8), intent(out) :: rot6d(6,6)

 rot6d = 0d0
end subroutine get_rot_mat_6d

subroutine get_rot_mat_7f(rotation, rot7f)
 implicit none
 integer :: i
 real(kind=8) :: vecs(3,7)
 real(kind=8), intent(in) :: rotation(3,3)
 real(kind=8), intent(out) :: rot7f(7,7)

 vecs = MATMUL(TRANSPOSE(rotation), vec7f)
 forall(i = 1:7) rot7f(:,i) = get7f_vector(vecs(:,i))
 rot7f = MATMUL(rot7f, invA_7f)
 deallocate(invA_7f)
end subroutine get_rot_mat_7f

subroutine get_rot_mat_9g(rotation, rot9g)
 implicit none
 integer :: i
 real(kind=8) :: vecs(3,9)
 real(kind=8), intent(in) :: rotation(3,3)
 real(kind=8), intent(out) :: rot9g(9,9)

 vecs = MATMUL(TRANSPOSE(rotation), vec9g)
 forall(i = 1:9) rot9g(:,i) = get9g_vector(vecs(:,i))
 rot9g = MATMUL(rot9g, invA_9g)
 deallocate(invA_9g)
end subroutine get_rot_mat_9g

subroutine get_rot_mat_11h(rotation, rot11h)
 implicit none
 integer :: i
 real(kind=8) :: vecs(3,11)
 real(kind=8), intent(in) :: rotation(3,3)
 real(kind=8), intent(out) :: rot11h(11,11)

 vecs = MATMUL(TRANSPOSE(rotation), vec11h)
 forall(i = 1:11) rot11h(:,i) = get11h_vector(vecs(:,i))
 rot11h = MATMUL(rot11h, invA_11h)
 deallocate(invA_11h)
end subroutine get_rot_mat_11h

! input a 3D-vector (x,y,z), calculate the values of spherical harmonic 5D functions
pure function get5d_vector(v1) result(v2)
 implicit none
 real(kind=8) :: r2, v2(5)
 real(kind=8), intent(in) :: v1(3)

 r2 = DOT_PRODUCT(v1, v1)
 v2(1) = (3d0*v1(3)*v1(3) - r2)/DSQRT(3d0)
 v2(2) = 2d0*v1(1)*v1(3)
 v2(3) = 2d0*v1(2)*v1(3)
 v2(4) = v1(1)*v1(1) - v1(2)*v1(2)
 v2(5) = 2d0*v1(1)*v1(2)
 v2 = v2*0.25d0*DSQRT(3.75d0/DATAN(1d0))/r2
end function get5d_vector

pure function get7f_vector(v1) result(v2)
 implicit none
 real(kind=8) :: x2, y2, z2, r2, r3, fivez2, v2(7)
 real(kind=8), intent(in) :: v1(3)

 x2 = v1(1)*v1(1); y2 = v1(2)*v1(2); z2 = v1(3)*v1(3)
 r2 = x2 + y2 + z2
 r3 = r2*DSQRT(r2)
 fivez2 = 5d0*z2
 v2(1) = DSQRT(2d0)*(fivez2 - 3d0*r2)*v1(3)
 v2(2) = DSQRT(3d0)*(fivez2 - r2)*v1(1)
 v2(3) = DSQRT(3d0)*(fivez2 - r2)*v1(2)
 v2(4) = DSQRT(30d0)*(x2 - y2)*v1(3)
 v2(5) = 2d0*DSQRT(30d0)*v1(1)*v1(2)*v1(3)
 v2(6) = DSQRT(5d0)*v1(1)*(x2 - 3d0*y2)
 v2(7) = DSQRT(5d0)*v1(2)*(3d0*x2 - y2)
 v2 = v2*0.125d0*DSQRT(3.5d0/DATAN(1d0))/r3
end function get7f_vector

pure function get9g_vector(v1) result(v2)
 implicit none
 real(kind=8) :: x2, y2, z2, r2, r4, seven_z2, xy, x2_y2, v2(9)
 real(kind=8), intent(in) :: v1(3)

 x2 = v1(1)*v1(1); y2 = v1(2)*v1(2); z2 = v1(3)*v1(3)
 xy = v1(1)*v1(2)
 x2_y2 = x2 - y2
 r2 = x2 + y2 + z2
 r4 = r2*r2
 seven_z2 = 7d0*z2
 v2(1) = 5d0*(seven_z2 - 6d0*r2)*z2 + 3d0*r4
 v2(2) = 2d0*root10*(seven_z2 - 3d0*r2)*v1(1)*v1(3)
 v2(3) = 2d0*root10*(seven_z2 - 3d0*r2)*v1(2)*v1(3)
 v2(4) = 2d0*root5*(seven_z2 - r2)*x2_y2
 v2(5) = 4d0*root5*(seven_z2 - r2)*xy
 v2(6) = 2d0*root70*(x2 - 3d0*y2)*v1(1)*v1(3)
 v2(7) = 2d0*root70*(3d0*x2 - y2)*v1(2)*v1(3)
 v2(8) = root35*(x2_y2 + 2d0*xy)*(x2_y2 - 2d0*xy)
 v2(9) = 4d0*root35*x2_y2*xy
 v2 = v2*3d0/(r4*32d0*DSQRT(DATAN(1d0)))
end function get9g_vector

pure function get11h_vector(v1) result(v2)
 implicit none
 real(kind=8) :: x2, y2, z2, r, r2, r3, r5, xy, xyz, x2_y2, z_r, z2_r2, z4_r4, &
  t1, t2, t3, v2(11)
 real(kind=8), intent(in) :: v1(3)

 x2 = v1(1)*v1(1); y2 = v1(2)*v1(2); z2 = v1(3)*v1(3)
 r2 = x2 + y2 + z2
 r = DSQRT(r2); r3 = r2*r; r5 = r3*r2
 z_r = v1(3)/r
 z2_r2 = z_r*z_r
 z4_r4 = z2_r2*z2_r2
 xy = v1(1)*v1(2); xyz = xy*v1(3)
 x2_y2 = x2 - y2
 t1 = 21d0*z4_r4 - 14d0*z2_r2 + 1d0
 t2 = 3d0*z2_r2 - 1d0
 t3 = 9d0*z2_r2 - 1d0
 v2(1) = 2d0*(63d0*(z_r**5) - 70d0*(z_r**3) + 15d0*z_r)
 v2(2) = 2d0*root15*t1*v1(1)/r
 v2(3) = 2d0*root15*t1*v1(2)/r
 v2(4) = 4d0*root105*t2*x2_y2*v1(3)/r3
 v2(5) = 8d0*root105*t2*xyz/r3
 v2(6) = root70*t3*(x2 - 3d0*y2)*v1(1)/r3
 v2(7) = root70*t3*(3d0*x2 - y2)*v1(2)/r3
 v2(8) = 6d0*root35*(x2_y2 + 2d0*xy)*(x2_y2 - 2d0*xy)*v1(3)/r5
 v2(9) = 24d0*root35*x2_y2*xyz/r5
 v2(10) = 3d0*root14*(x2*x2 - 10d0*x2*y2 + 5d0*y2*y2)*v1(1)/r5
 v2(11) = 3d0*root14*(5d0*x2*x2 - 10d0*x2*y2 + y2*y2)*v1(2)/r5
 v2 = v2*DSQRT(2.75d0/DATAN(1d0))/32d0
end function get11h_vector
end module bas_rot

! Get the wave function of a rotated molecule. The original coordinates are
!  provided in fchname, while the new coordinates are provided in coor_file
!  (which is a .xyz or .gjf file).
! Note: the one-to-one correspondence of element symbols in two files will not
!  be checked, so it is assumed that you already ensure the correspondence.
! Assuming the input file is a.fch, the output file would be a_r.fch
subroutine rotate_atoms_wfn(fchname, coor_file)
 implicit none
 integer :: i, natom, natom1
 real(kind=8), allocatable :: coor(:,:)
 character(len=240) :: new_fch
 character(len=240), intent(in) :: fchname, coor_file
!f2py intent(in) :: fchname, coor_file

 i = index(fchname, '.fch', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine rotate_atoms_wfn: '.fch' suffix not fou&
                   &nd in"
  write(6,'(A)') 'filename '//TRIM(fchname)
  stop
 end if
 ! the filename to store rotated molecule and wfn
 new_fch = fchname(1:i-1)//'_r.fch'

 i = index(coor_file, '.xyz', back=.true.)
 if(i > 0) then
  call read_natom_from_xyz(coor_file, natom)
 else ! not .xyz, assume it be .gjf
  call read_natom_from_gjf(coor_file, natom)
 end if

 call read_natom_from_fch(fchname, natom1)
 if(natom1 /= natom) then
  write(6,'(/,A)') 'ERROR in subroutine rotate_atoms_wfn: the number of atoms i&
                   &s not equal to'
  write(6,'(A)') 'each other in two files.'
  write(6,'(A,I0)') 'fchname = '//TRIM(fchname)//', natom1 = ', natom1
  write(6,'(A,I0)') 'coor_file = '//TRIM(coor_file)//', natom = ', natom
  stop
 end if

 allocate(coor(3,natom))
 call read_coor_from_gjf_or_xyz(coor_file, natom, coor)

 call rotate_atoms_wfn2(fchname, natom, coor, new_fch)
 deallocate(coor)
end subroutine rotate_atoms_wfn

! Get the wave function of a rotated molecule. The original coordinates are
!  provided in fchname, while the new coordinates are provided in the array
!  coor. The rotated molecule and wfn are stored in new_fch.
! Note: the one-to-one correspondence of element symbols will not be checked,
!  so it is assumed that you already ensure the correspondence.
! Assuming the input file is a.fch, the output file would be a_r.fch
subroutine rotate_atoms_wfn2(fchname, natom, coor, new_fch)
 use fch_content, only: nbf, nif, is_uhf, ncontr, shell_type, alpha_coeff, &
  beta_coeff, tot_dm, spin_dm, check_uhf_in_fch, read_fch, coor0=>coor
 use bas_rot, only: get_invA, get_rot_mat_5d, get_rot_mat_6d, get_rot_mat_7f, &
  get_rot_mat_9g, get_rot_mat_11h
 implicit none
 integer :: i, j, nif1
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 integer, allocatable :: no_mark(:), mark(:,:)
 real(kind=8) :: rmsd_v, trans1(3), trans2(3), rotation(3,3), rot5d(5,5), &
  rot6d(6,6), rot7f(7,7), rot9g(9,9), rot11h(11,11)
 real(kind=8), intent(in) :: coor(3,natom)
!f2py intent(in) :: coor
!f2py depend(natom) :: coor
 real(kind=8), allocatable :: rtmp(:,:), coeff(:,:)
 character(len=240), intent(in) :: fchname, new_fch
!f2py intent(in) :: fchname, new_fch

 call check_uhf_in_fch(fchname, is_uhf)
 call read_fch(fchname, is_uhf)

 call rmsd(natom, coor0, coor, rmsd_v, trans1, trans2, rotation)
 ! coor0 will be translated to its own center and rotated after calling the
 ! subroutine rmsd

 if(rmsd_v > 2d-2) then
  write(6,'(/,A)') 'Warning in subroutine rotate_atoms_wfn2: RMSD value > 0.02.'
  write(6,'(A)') 'Anyway, the program will continue, but you should be aware of&
                 & what you are'
  write(6,'(A,I0)') 'doing. fchname='//TRIM(fchname)//', natom=', natom
 end if

 ! now translate coor0 to the frame of coor
 forall(i = 1:natom) coor0(:,i) = coor0(:,i) - trans2

 ! merge alpha_coeff and beta_coeff for easier manipulation
 nif1 = nif
 if(is_uhf) nif1 = 2*nif
 allocate(coeff(nbf,nif1))
 coeff(:,1:nif) = alpha_coeff
 deallocate(alpha_coeff)
 if(is_uhf) then
  coeff(:,nif+1:nif1) = beta_coeff
  deallocate(beta_coeff)
 end if

 ! get the basis function marks of 3P/5D/6D/7F/10F/9G/15G/11H/21H from the
 ! integer array shell_type
 allocate(no_mark(9), mark(ncontr,9))
 call read_pdfgh_mark_from_shltyp(ncontr, shell_type, no_mark, mark)

 if(ANY(shell_type > 1)) then
  write(6,'(/,A)') 'ERROR in subroutine rotate_atoms_wfn2: Cartesian-type basis&
                  & functions 6D'
  write(6,'(A)') '(or higher) are not supported currently.'
  stop
 end if

 ! the rotation of MO coefficients of 3P angular momentum are the same to that
 ! of Cartesian coordinates, i.e. left-multiplied by R^T
 allocate(rtmp(3,nif1), source=0d0)
 do i = 1, no_mark(1), 1
  j = mark(i,1)
  rtmp = coeff(j:j+2,:)
  call dgemm('T', 'N', 3, nif1, 3, 1d0,rotation,3, rtmp,3, 0d0,coeff(j:j+2,:),3)
 end do ! for i
 deallocate(rtmp)

 ! The rotation of MO coefficients of angular momentum higher than 3P are a little
 ! bit of tedious. There are two ways to calculate the rotated MO coefficients:
 ! (1) using the Wigner-D matrix; (2) do not use the Wigner-D matrix, but solve
 ! and save the inverse matrix of basis function values (which can be hard-coded).
 ! Here the second approach is adopted.

 ! 5D
 if(no_mark(2) > 0) then
  call get_rot_mat_5d(rotation, rot5d)
  allocate(rtmp(5,nif1), source=0d0)
  do i = 1, no_mark(2), 1
   j = mark(i,2)
   rtmp = coeff(j:j+4,:)
   call dgemm('N', 'N', 5, nif1, 5, 1d0,rot5d,5, rtmp,5, 0d0,coeff(j:j+4,:),5)
  end do ! for i
  deallocate(rtmp)
 end if

 ! 6D
 if(no_mark(3) > 0) then
  call get_rot_mat_6d(rotation, rot6d)
  allocate(rtmp(6,nif1), source=0d0)
  do i = 1, no_mark(3), 1
   j = mark(i,3)
   rtmp = coeff(j:j+5,:)
   call dgemm('N', 'N', 6, nif1, 6, 1d0,rot6d,6, rtmp,6, 0d0,coeff(j:j+5,:),6)
  end do ! for i
  deallocate(rtmp)
 end if

 ! 7F
 if(no_mark(4) > 0) then
  call get_invA(7)
  call get_rot_mat_7f(rotation, rot7f)
  allocate(rtmp(7,nif1), source=0d0)
  do i = 1, no_mark(4), 1
   j = mark(i,4)
   rtmp = coeff(j:j+6,:)
   call dgemm('N', 'N', 7, nif1, 7, 1d0,rot7f,7, rtmp,7, 0d0,coeff(j:j+6,:),7)
  end do ! for i
  deallocate(rtmp)
 end if

 ! 9G
 if(no_mark(6) > 0) then
  call get_invA(9)
  call get_rot_mat_9g(rotation, rot9g)
  allocate(rtmp(9,nif1), source=0d0)
  do i = 1, no_mark(6), 1
   j = mark(i,6)
   rtmp = coeff(j:j+8,:)
   call dgemm('N', 'N', 9, nif1, 9, 1d0,rot9g,9, rtmp,9, 0d0,coeff(j:j+8,:),9)
  end do ! for i
  deallocate(rtmp)
 end if

 ! 11H
 if(no_mark(8) > 0) then
  call get_invA(11)
  call get_rot_mat_11h(rotation, rot11h)
  allocate(rtmp(11,nif1), source=0d0)
  do i = 1, no_mark(8), 1
   j = mark(i,8)
   rtmp = coeff(j:j+10,:)
   call dgemm('N', 'N', 11, nif1, 11, 1d0,rot11h,11, rtmp,11, 0d0,coeff(j:j+10,:),11)
  end do ! for i
  deallocate(rtmp)
 end if

 deallocate(no_mark, mark)
 allocate(alpha_coeff(nbf,nif), source=coeff(:,1:nif))
 if(is_uhf) allocate(beta_coeff(nbf,nif), source=coeff(:,nif+1:nif1))
 deallocate(coeff)

 ! currently the tot_dm and spin_dm are not rotated, so set to 0
 tot_dm = 0d0
 if(is_uhf) spin_dm = 0d0

 call write_fch(new_fch)
end subroutine rotate_atoms_wfn2

! calculate the RMSD value of two sets of coordinates
subroutine rmsd(natom, coor1, coor2, rmsd_v, trans1, trans2, rotation)
 implicit none
 integer :: i, j, lwork
 integer, intent(in) :: natom
 real(kind=8), intent(inout) :: coor1(3,natom)
 real(kind=8), intent(in) :: coor2(3,natom)
 real(kind=8), intent(out) :: trans1(3), trans2(3), rotation(3,3)
 real(kind=8), intent(out) :: rmsd_v
 real(kind=8) :: d, H(3,3), u(3,3), vt(3,3), s(3), uvt(3,3), unity(3,3)
 real(kind=8), allocatable :: work(:), new_coor1(:,:), coor(:,:)

 u = 0d0; vt = 0d0; s = 0d0; uvt = 0d0; H = 0d0 ! initialization

 ! translate coor1 to its origin
 call move_coor_to_center(natom, coor1, trans1)

 ! make a copy of coor2 as coor. translate coor to its origin, and keep coor2 fixed
 allocate(coor(3,natom), source=coor2)
 call move_coor_to_center(natom, coor, trans2)

 call dgemm('N', 'T', 3, 3, natom, 1d0, coor1, 3, coor, 3, 0d0, H, 3)

 ! SVD on H
 lwork = 30
 allocate(work(lwork), source=0d0)
 ! call dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
 call dgesvd('A', 'A', 3, 3, H, 3, s, u, 3, vt, 3, work, lwork, i)
 deallocate(work)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine rmsd: MKL dgesvd error info /= 0.'
  write(6,'(A,I0)') 'info = ', i
  stop
 end if

 uvt = MATMUL(u, vt) ! U(V^T)

 ! d = |U(V^T)|
 d = uvt(1,1)*(uvt(2,2)*uvt(3,3) - uvt(2,3)*uvt(3,2)) - &
     uvt(1,2)*(uvt(2,1)*uvt(3,3) - uvt(2,3)*uvt(3,1)) + &
     uvt(1,3)*(uvt(2,1)*uvt(3,2) - uvt(2,2)*uvt(3,1))
 unity = 0d0
 forall(i = 1:3) unity(i,i) = 1d0
 if(d < 0d0) unity(3,3) = -1d0

 ! R = UD(V^T), using array uvt to hold UD
 uvt  = MATMUL(u, unity)
 rotation = MATMUL(uvt, vt)

 ! P' = PR
 allocate(new_coor1(3,natom), source=0d0)
 call dgemm('T', 'N', 3, natom, 3, 1d0, rotation,3, coor1,3, 0d0, new_coor1, 3)
 coor1 = new_coor1
 deallocate(new_coor1)

 coor = coor - coor1
 allocate(work(natom))
 forall(i = 1:natom) work(i) = DOT_PRODUCT(coor(:,i), coor(:,i))
 deallocate(coor)

 rmsd_v = DSQRT(SUM(work)/DBLE(natom))
 deallocate(work)
end subroutine rmsd

! calculate the geometry center (centeroid) of a set of coordinates,
! return the corresponding translational vector and translated coordinates
subroutine move_coor_to_center(natom, coor, trans)
 implicit none
 integer :: i
 integer, intent(in) :: natom
 real(kind=8), intent(inout) :: coor(3,natom)
 real(kind=8), intent(out) :: trans(3)

 do i = 1, 3
  trans(i) = -SUM(coor(i,:))
 end do ! for i

 trans = trans/DBLE(natom)
 forall(i = 1:natom) coor(:,i) = coor(:,i) + trans
end subroutine move_coor_to_center

! .chk -> .chk, a wrapper of {formchk, mirror_wfn, and unfchk}
subroutine mirror_c2c(chkname)
 use util_wrapper, only: formchk, unfchk
 implicit none
 integer :: i
 character(len=240) :: fchname1, fchname2
 character(len=240), intent(in) :: chkname
!f2py intent(in) :: chkname

 call formchk(chkname)
 i = index(chkname, '.chk', back=.true.)
 fchname1 = chkname(1:i-1)//'.fch'
 fchname2 = chkname(1:i-1)//'_m.fch'
 call mirror_wfn(fchname1)
 call unfchk(fchname2)
end subroutine mirror_c2c

! Get the wave function of the mirror of a molecule.
! The z-component of Cartesian coordinates are multiplied by -1.
! Some of the MO coefficients on basis functions are multiplied by -1
!  5D: D+1, D-1
!  7F: F0, F+2, F-2
!  9G: G+1, G-1, G+3, G-3
! 11H: H+2, H-2, H+4, H-4
subroutine mirror_wfn(fchname)
 use fch_content, only: nbf, nif, ncontr, alpha_coeff, beta_coeff, shell_type, &
  is_uhf, coor, tot_dm, spin_dm, check_uhf_in_fch, read_fch
 implicit none
 integer :: i, j, k, nif1
 integer, parameter :: a1(3) = [1,4,5]
 integer, parameter :: a2(4) = [3,6,9,10]
 integer, parameter :: a3(4) =[2,3,6,7]
 integer, parameter :: a4(6) = [2,4,6,8,11,13]
 integer, parameter :: a5(4) = [4,5,8,9]
 integer, parameter :: a6(9) = [1,3,5,8,10,12,14,17,19]
 character(len=240) :: new_fch
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 real(kind=8), allocatable :: coeff(:,:)
 logical, allocatable :: pos(:)
 logical :: uhf

 call require_file_exist(fchname)
 i = index(fchname, '.fch', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine mirror_wfn: '.fch' suffix not found in&
                   & file "//TRIM(fchname)
  stop
 end if
 new_fch = fchname(1:i-1)//'_m.fch'

 call check_uhf_in_fch(fchname, uhf)
 is_uhf = uhf ! be careful of this variable

 call read_fch(fchname, uhf)
 coor(3,:) = -coor(3,:) ! z*(-1)

 nif1 = nif
 if(uhf) nif1 = 2*nif
 allocate(coeff(nbf,nif1))
 if(uhf) then
  coeff(:,1:nif) = alpha_coeff
  coeff(:,nif+1:nif1) = beta_coeff
 else
  coeff = alpha_coeff
 end if

 allocate(pos(nbf)) ! record the positions to be positive or to be negative
 pos = .true.
 k = 0

 do i = 1, ncontr, 1
  select case(shell_type(i))
  case( 0)   ! S
   k = k + 1
  case( 1)   ! 3P
   k = k + 3
   pos(k) = .false.
  case(-1)   ! SP or L
   k = k + 4
   pos(k) = .false.
  case(-2)   ! 5D
   pos(k+2:k+3) = .false.
   k = k + 5
  case( 2)   ! 6D
   k = k + 6
   pos(k-1:k) = .false.
  case(-3)   ! 7F
   forall(j = 1:3) pos(k+a1(j)) = .false.
   k = k + 7
  case( 3)   ! 10F
   forall(j = 1:4) pos(k+a2(j)) = .false.
   k = k + 10
  case(-4)   ! 9G
   forall(j = 1:4) pos(k+a3(j)) = .false.
   k = k + 9
  case( 4)   ! 15G
   forall(j = 1:6) pos(k+a4(j)) = .false.
   k = k + 15
  case(-5)   ! 11H
   forall(j = 1:4) pos(k+a5(j)) = .false.
   k = k + 11
  case( 5)   ! 21H
   forall(j = 1:9) pos(k+a6(j)) = .false.
   k = k + 21
  case default
   write(6,'(/,A,I0)') 'ERROR in subroutine mirror_wfn: invalid shell_type(i)=',&
                        shell_type(i)
   stop
  end select
 end do ! for i

 if(k /= nbf) then
  write(6,'(/,A)') 'ERROR in subroutine mirror_wfn: internal inconsistency.'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  write(6,'(2(A,I0))') 'nbf=', nbf, ', k=', k
  stop
 end if

 ! adjust MO coefficients
 forall(i=1:nbf, (.not.pos(i))) coeff(i,:) = -coeff(i,:)
 if(uhf) then
  alpha_coeff = coeff(:,1:nif)
  beta_coeff = coeff(:,nif+1:nif1)
 else
  alpha_coeff = coeff
 end if
 deallocate(coeff)

 ! adjust Total SCF Density
 forall(i=1:nbf, j=1:nbf, (pos(i).neqv.pos(j))) tot_dm(j,i) = -tot_dm(j,i)

 ! adjust Spin SCF Density, if any
 if(allocated(spin_dm)) then
  forall(i=1:nbf, j=1:nbf, (pos(i).neqv.pos(j)))
   spin_dm(j,i) = -spin_dm(j,i)
  end forall
 end if
 deallocate(pos)

 ! generate a new .fch file
 call write_fch(new_fch)
end subroutine mirror_wfn

! generate geometries using linear interpolation of Cartesian coordinates
! fname1: .gjf file which includes the initial geometry
! fname2: .gjf file which includes the final geometry
! Note: the spatial orientations of the initial and the final geometry can be
!  different, since the RMSD will be called firstly to align two molecules. But
!  the one-to-one correspondence of atoms in two files must be ensured by the
!  user(e.g. C1-C1, H2-H2).
subroutine geom_lin_intrplt(gjfname1, gjfname2, n)
 implicit none
 integer :: i, k, natom0, natom, charge, mult
 integer, intent(in) :: n ! the number of geometries to be generated
!f2py intent(in) :: n
 integer, allocatable :: nuc(:)
 real(kind=8) :: r, rmsd_v, trans1(3), trans2(3), rotation(3,3)
 real(kind=8), allocatable :: coor1(:,:), coor2(:,:), coor(:,:)
 character(len=240) :: gjfname
 character(len=240), intent(in) :: gjfname1, gjfname2
!f2py intent(in) :: gjfname1, gjfname2
 character(len=240), allocatable :: elem(:)

 k = index(gjfname1, '.gjf', back=.true.)
 if(k == 0) then
  write(6,'(/,A)') "ERROR in subroutine geom_lin_intrplt: '.gjf' suffix not fou&
                   &nd in filename "//TRIM(gjfname1)
  stop
 end if
 call read_natom_from_gjf(gjfname1, natom0)
 call read_natom_from_gjf(gjfname2, natom)

 if(natom0 /= natom) then
  write(6,'(/,A)') 'ERROR in subroutine geom_lin_intrplt: the number of atoms a&
                   &re inconsistent'
  write(6,'(A)') 'in two input files.'
  write(6,'(A)') 'gjfname1='//TRIM(gjfname1)
  write(6,'(A)') 'gjfname2='//TRIM(gjfname2)
  stop
 end if

 allocate(coor1(3,natom), coor2(3,natom), coor(3,natom), elem(natom), nuc(natom))
 call read_elem_and_coor_from_gjf(gjfname1, natom, elem, nuc, coor1, charge, mult)
 call read_elem_and_coor_from_gjf(gjfname2, natom, elem, nuc, coor2, charge, mult)
 deallocate(nuc)

 ! in case that two geometries have different orientations, let's do DMSD first
 call rmsd(natom, coor1, coor2, rmsd_v, trans1, trans2, rotation)
 if(rmsd_v < 1d-2) then
  write(6,'(/,A)') 'Warning from subroutine geom_lin_intrplt: the initial geome&
                   &try is too similar'
  write(6,'(A,F8.4)') 'to the final geometry. RMSD=', rmsd_v
  write(6,'(A)') 'You should be aware of what you are doing. Anyway, this subro&
                 &utine will continue.'
 end if

 ! move coor1 into the frame of coor2
 forall(i = 1:natom) coor1(:,i) = coor1(:,i) - trans2

 do i = 1, n, 1
  write(gjfname,'(A,I0,A)') gjfname1(1:k-1)//'_',i,'.gjf'
  r = DBLE(i)/DBLE(n+1)
  coor = (1d0-r)*coor1 + r*coor2
  call write_gjf(gjfname, charge, mult, natom, elem, coor)
 end do ! for i

 deallocate(coor1, coor2, elem, coor)
end subroutine geom_lin_intrplt

