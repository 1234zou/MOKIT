! written by jxzou at 20230516: generate the wave function of the mirror of a
! molecule

module bas_rot
 implicit none
 real(kind=8), parameter :: root2 = DSQRT(2d0), root3 = DSQRT(3d0), &
  root5 = DSQRT(5d0), root7 = DSQRT(7d0), root10 = DSQRT(10d0), &
  root14 = DSQRT(14d0), root15 = DSQRT(15d0), root21 = DSQRT(21d0), &
  root30 = DSQRT(30d0), root35 = DSQRT(35d0), root70 = DSQRT(70d0), &
  root105 = DSQRT(105d0)

 ! Note: there is no need to declare arrays vec7f, vec9g, etc, we can simply use
 ! part of the array vec21h
 real(kind=8), parameter :: vec21h(3,21) = RESHAPE([1d0,0d0,0d0, 0d0,0d0,1d0, &
  1d0,root3,1d0, 1d0,0d0,1d0, 0d0,1d0,1d0, 1d0,1d0,0d0, root3,1d0,root5, &
  root5,1d0,root3, root5,root3,1d0, root3,-10d0,-root14, -root3,-root14,10d0, &
  -root10,7d0,-root15, -root10,-root15,7d0, 7d0,-root10,-root15, &
  7d0,-root15,-root10, -root14,root3,-10d0, -root14,-10d0,root3, & 
  -root21,5d0,root7, -root21,root7,5d0, root7,-root21,5d0, root7,5d0,-root21],&
  [3,21])

 real(kind=8), allocatable :: A_7f(:,:), A_10f(:,:), A_9g(:,:), A_15g(:,:), &
  A_11h(:,:), A_21h(:,:), invA_7f(:,:), invA_10f(:,:), invA_9g(:,:), &
  invA_15g(:,:), invA_11h(:,:), invA_21h(:,:)

contains

! Hard-coding of invA_7f, invA_9g, or higher is so tedious that typos will
!  probably occur, so here we solve the three inverse matrices numerically.
! The invA_5d is hard-coded in the subroutine get_rot_mat_5d.
subroutine get_invA(n)
 implicit none
 integer :: i
 integer, intent(in) :: n

 select case(n)
 case(7) ! 7F
  allocate(A_7f(7,7), invA_7f(7,7))
  forall(i = 1:7) A_7f(:,i) = get7f_vector(vec21h(:,i))
  call inverse(7, A_7f, invA_7f)
  deallocate(A_7f)
 case(10) ! 10F
  allocate(A_10f(10,10), invA_10f(10,10))
  forall(i = 1:10) A_10f(:,i) = get10f_vector(vec21h(:,i))
  call inverse(10, A_10f, invA_10f)
  deallocate(A_10f)
 case(9) ! 9G
  allocate(A_9g(9,9), invA_9g(9,9))
  forall(i = 1:9) A_9g(:,i) = get9g_vector(vec21h(:,i))
  call inverse(9, A_9g, invA_9g)
  deallocate(A_9g)
 case(15) ! 15G
  allocate(A_15g(15,15), invA_15g(15,15))
  forall(i = 1:15) A_15g(:,i) = get15g_vector(vec21h(:,i))
  call inverse(15, A_15g, invA_15g)
  deallocate(A_15g)
 case(11) ! 11H
  allocate(A_11h(11,11), invA_11h(11,11))
  forall(i = 1:11) A_11h(:,i) = get11h_vector(vec21h(:,i))
  call inverse(11, A_11h, invA_11h)
  deallocate(A_11h)
 case(21) ! 21H
  allocate(A_21h(21,21), invA_21h(21,21))
  forall(i = 1:21) A_21h(:,i) = get21h_vector(vec21h(:,i))
  call inverse(21, A_21h, invA_21h)
  deallocate(A_21h)
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
 integer :: i
 real(kind=8) :: vecs(3,6)
 real(kind=8), intent(in) :: rotation(3,3)
 real(kind=8), intent(out) :: rot6d(6,6)
 real(kind=8), parameter :: s = 4d0*DSQRT(DATAN(1d0)/15d0)
 real(kind=8), parameter :: t1 = 1d0/3d0, t2 = 1.2d0, t3 = -0.3d0
 real(kind=8), parameter :: vecs0(3,6) = &
  RESHAPE([1d0,0d0,0d0, 0d0,0d0,1d0, 1d0,1d0,0d0, 1d0,0d0,1d0, 0d0,1d0,1d0, &
           0d0,1d0,0d0], [3,6])
 ! The inverse of square matrix A is hard-coded here
 real(kind=8), parameter :: invA(6,6) = RESHAPE([root3,0d0,0d0,0d0,0d0,0d0, &
  0d0,0d0,0d0,0d0,0d0,root3, 0d0,root3,0d0,0d0,0d0,0d0, &
  -1d0,0d0,2d0,0d0,0d0,-1d0, -1d0,-1d0,0d0,2d0,0d0,0d0, &
  0d0,-1d0,0d0,0d0,2d0,-1d0], [6,6])
 real(kind=8), parameter :: ovlp(6,6) = RESHAPE([1d0,t1,t1,0d0,0d0,0d0, &
  t1,1d0,t1,0d0,0d0,0d0, t1,t1,1d0,0d0,0d0,0d0, 0d0,0d0,0d0,1d0,0d0,0d0, &
  0d0,0d0,0d0,0d0,1d0,0d0, 0d0,0d0,0d0,0d0,0d0,1d0], [6,6])
 real(kind=8), parameter :: invS(6,6) = RESHAPE([t2,t3,t3,0d0,0d0,0d0, &
  t3,t2,t3,0d0,0d0,0d0, t3,t3,t2,0d0,0d0,0d0, 0d0,0d0,0d0,1d0,0d0,0d0, &
  0d0,0d0,0d0,0d0,1d0,0d0, 0d0,0d0,0d0,0d0,0d0,1d0], [6,6])
 ! invS = (ovlp)^(-1)

 vecs = MATMUL(TRANSPOSE(rotation), vecs0)

 forall(i = 1:6) rot6d(:,i) = get6d_vector(vecs(:,i))

 rot6d = MATMUL(rot6d, s*invA)

 ! Note: the (x2,y2,z2) in 6D are not orthogonal, so we need the overlap matrix
 ! and its inverse
 rot6d = MATMUL(MATMUL(invS, rot6d), ovlp)
end subroutine get_rot_mat_6d

subroutine get_rot_mat_7f(rotation, rot7f)
 implicit none
 integer :: i
 real(kind=8) :: vecs(3,7)
 real(kind=8), intent(in) :: rotation(3,3)
 real(kind=8), intent(out) :: rot7f(7,7)

 vecs = MATMUL(TRANSPOSE(rotation), vec21h(:,1:7))
 forall(i = 1:7) rot7f(:,i) = get7f_vector(vecs(:,i))
 rot7f = MATMUL(rot7f, invA_7f)
 deallocate(invA_7f)
end subroutine get_rot_mat_7f

subroutine get_rot_mat_10f(rotation, rot10f)
 implicit none
 integer :: i
 real(kind=8) :: vecs(3,10)
 real(kind=8), intent(in) :: rotation(3,3)
 real(kind=8), intent(out) :: rot10f(10,10)
 real(kind=8), parameter :: t1 = 1d0/DSQRT(5d0), t2 = 1d0/3d0, t3 = 10d0/7d0, &
  t4 = -3d0*DSQRT(5d0)/14d0, t5 = 9d0/7d0, t6 = -3d0/14d0
 real(kind=8), parameter :: ovlp(10,10) = RESHAPE([&
  1d0,0d0,0d0,t1,0d0,0d0,t1,0d0,0d0,0d0, 0d0,1d0,0d0,0d0,t1,0d0,0d0,t1,0d0,0d0,&
  0d0,0d0,1d0,0d0,0d0,t1,0d0,0d0,t1,0d0, t1,0d0,0d0,1d0,0d0,0d0,t2,0d0,0d0,0d0,&
  0d0,t1,0d0,0d0,1d0,0d0,0d0,t2,0d0,0d0, 0d0,0d0,t1,0d0,0d0,1d0,0d0,0d0,t2,0d0,&
  t1,0d0,0d0,t2,0d0,0d0,1d0,0d0,0d0,0d0, 0d0,t1,0d0,0d0,t2,0d0,0d0,1d0,0d0,0d0,&
  0d0,0d0,t1,0d0,0d0,t2,0d0,0d0,1d0,0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0, &
  1d0], [10,10])
 real(kind=8), parameter :: invS(10,10) = RESHAPE([&
  t3,0d0,0d0,t4,0d0,0d0,t4,0d0,0d0,0d0, 0d0,t3,0d0,0d0,t4,0d0,0d0,t4,0d0,0d0, &
  0d0,0d0,t3,0d0,0d0,t4,0d0,0d0,t4,0d0, t4,0d0,0d0,t5,0d0,0d0,t6,0d0,0d0,0d0, &
  0d0,t4,0d0,0d0,t5,0d0,0d0,t6,0d0,0d0, 0d0,0d0,t4,0d0,0d0,t5,0d0,0d0,t6,0d0, &
  t4,0d0,0d0,t6,0d0,0d0,t5,0d0,0d0,0d0, 0d0,t4,0d0,0d0,t6,0d0,0d0,t5,0d0,0d0, &
  0d0,0d0,t4,0d0,0d0,t6,0d0,0d0,t5,0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0, &
  1d0], [10,10])

 vecs = MATMUL(TRANSPOSE(rotation), vec21h(:,1:10))
 forall(i = 1:10) rot10f(:,i) = get10f_vector(vecs(:,i))
 rot10f = MATMUL(rot10f, invA_10f)
 deallocate(invA_10f)

 ! Note: the 10F are not orthogonal, so we need the overlap matrix and its inverse
 rot10f = MATMUL(MATMUL(invS, rot10f), ovlp)
end subroutine get_rot_mat_10f

subroutine get_rot_mat_9g(rotation, rot9g)
 implicit none
 integer :: i
 real(kind=8) :: vecs(3,9)
 real(kind=8), intent(in) :: rotation(3,3)
 real(kind=8), intent(out) :: rot9g(9,9)

 vecs = MATMUL(TRANSPOSE(rotation), vec21h(:,1:9))
 forall(i = 1:9) rot9g(:,i) = get9g_vector(vecs(:,i))
 rot9g = MATMUL(rot9g, invA_9g)
 deallocate(invA_9g)
end subroutine get_rot_mat_9g

subroutine get_rot_mat_15g(rotation, rot15g)
 implicit none
 integer :: i
 real(kind=8) :: vecs(3,15)
 real(kind=8), intent(in) :: rotation(3,3)
 real(kind=8), intent(out) :: rot15g(15,15)
 real(kind=8), parameter :: t1 = DSQRT(105d0)/21d0, t2 = 3d0/35d0, &
  t3 = 1d0/DSQRT(105d0), t4 = 0.6d0, t5 = 1d0/DSQRT(5d0), t6 = 1d0/3d0, &
  t7 = 5d0/3d0, t8 = -DSQRT(105d0)/14d0, t9 = 5d0/24d0, t10=DSQRT(105d0)/84d0, &
  t11 = -5d0/6d0, t12 = -DSQRT(5d0)/6d0, t13 = 51d0/28d0, t14 = -5d0/28d0, &
  t15 = 4d0/3d0
 real(kind=8), parameter :: ovlp(15,15) = RESHAPE([&
  1d0,0d0,t1,0d0,t2,0d0,0d0,0d0,0d0,t1,0d0,t3,0d0,0d0,t2, &
  0d0,1d0,0d0,t4,0d0,0d0,0d0,0d0,0d0,0d0,t5,0d0,0d0,0d0,0d0, &
  t1,0d0,1d0,0d0,t1,0d0,0d0,0d0,0d0,t6,0d0,t6,0d0,0d0,t3, &
  0d0,t4,0d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0,t5,0d0,0d0,0d0,0d0, &
  t2,0d0,t1,0d0,1d0,0d0,0d0,0d0,0d0,t3,0d0,t1,0d0,0d0,t2, &
  0d0,0d0,0d0,0d0,0d0,1d0,0d0,t5,0d0,0d0,0d0,0d0,t4,0d0,0d0, &
  0d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0,t5,0d0,0d0,0d0,0d0,t5,0d0, &
  0d0,0d0,0d0,0d0,0d0,t5,0d0,1d0,0d0,0d0,0d0,0d0,t5,0d0,0d0, &
  0d0,0d0,0d0,0d0,0d0,0d0,t5,0d0,1d0,0d0,0d0,0d0,0d0,t4,0d0, &
  t1,0d0,t6,0d0,t3,0d0,0d0,0d0,0d0,1d0,0d0,t6,0d0,0d0,t1, &
  0d0,t5,0d0,t5,0d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,0d0, &
  t3,0d0,t6,0d0,t1,0d0,0d0,0d0,0d0,t6,0d0,1d0,0d0,0d0,t1, &
  0d0,0d0,0d0,0d0,0d0,t4,0d0,t5,0d0,0d0,0d0,0d0,1d0,0d0,0d0, &
  0d0,0d0,0d0,0d0,0d0,0d0,t5,0d0,t4,0d0,0d0,0d0,0d0,1d0,0d0, &
  t2,0d0,t3,0d0,t2,0d0,0d0,0d0,0d0,t1,0d0,t1,0d0,0d0,1d0], [15,15])
 real(kind=8), parameter :: invS(15,15) = RESHAPE([&
  t7,0d0,t8,0d0,t9,0d0,0d0,0d0,0d0,t8,0d0,t10,0d0,0d0,t9, &
  0d0,t7,0d0,t11,0d0,0d0,0d0,0d0,0d0,0d0,t12,0d0,0d0,0d0,0d0, &
  t8,0d0,t13,0d0,t8,0d0,0d0,0d0,0d0,t14,0d0,t14,0d0,0d0,t10, &
  0d0,t11,0d0,t7,0d0,0d0,0d0,0d0,0d0,0d0,t12,0d0,0d0,0d0,0d0, &
  t9,0d0,t8,0d0,t7,0d0,0d0,0d0,0d0,t10,0d0,t8,0d0,0d0,t9, &
  0d0,0d0,0d0,0d0,0d0,t7,0d0,t12,0d0,0d0,0d0,0d0,t11,0d0,0d0, &
  0d0,0d0,0d0,0d0,0d0,0d0,t15,0d0,t12,0d0,0d0,0d0,0d0,t12,0d0, &
  0d0,0d0,0d0,0d0,0d0,t12,0d0,t15,0d0,0d0,0d0,0d0,t12,0d0,0d0, &
  0d0,0d0,0d0,0d0,0d0,0d0,t12,0d0,t7,0d0,0d0,0d0,0d0,t11,0d0, &
  t8,0d0,t14,0d0,t10,0d0,0d0,0d0,0d0,t13,0d0,t14,0d0,0d0,t8, &
  0d0,t12,0d0,t12,0d0,0d0,0d0,0d0,0d0,0d0,t15,0d0,0d0,0d0,0d0, &
  t10,0d0,t14,0d0,t8,0d0,0d0,0d0,0d0,t14,0d0,t13,0d0,0d0,t8, &
  0d0,0d0,0d0,0d0,0d0,t11,0d0,t12,0d0,0d0,0d0,0d0,t7,0d0,0d0, &
  0d0,0d0,0d0,0d0,0d0,0d0,t12,0d0,t11,0d0,0d0,0d0,0d0,t7,0d0, &
  t9,0d0,t10,0d0,t9,0d0,0d0,0d0,0d0,t8,0d0,t8,0d0,0d0,t7], [15,15])

 vecs = MATMUL(TRANSPOSE(rotation), vec21h(:,1:15))
 forall(i = 1:15) rot15g(:,i) = get15g_vector(vecs(:,i))
 rot15g = MATMUL(rot15g, invA_15g)
 deallocate(invA_15g)

 ! Note: the 15G are not orthogonal, so we need the overlap matrix and its inverse
 rot15g = MATMUL(MATMUL(invS, rot15g), ovlp)
end subroutine get_rot_mat_15g

subroutine get_rot_mat_11h(rotation, rot11h)
 implicit none
 integer :: i
 real(kind=8) :: vecs(3,11)
 real(kind=8), intent(in) :: rotation(3,3)
 real(kind=8), intent(out) :: rot11h(11,11)

 vecs = MATMUL(TRANSPOSE(rotation), vec21h(:,1:11))
 forall(i = 1:11) rot11h(:,i) = get11h_vector(vecs(:,i))
 rot11h = MATMUL(rot11h, invA_11h)
 deallocate(invA_11h)
end subroutine get_rot_mat_11h

subroutine get_rot_mat_21h(rotation, rot21h)
 implicit none
 integer :: i
 real(kind=8) :: vecs(3,21)
 real(kind=8), intent(in) :: rotation(3,3)
 real(kind=8), intent(out) :: rot21h(21,21)
 real(kind=8), parameter :: t1 = DSQRT(21d0)/9d0, t2 = 1d0/7d0, &
  t3 = DSQRT(105d0)/63d0, t4 = DSQRT(21d0)/7d0, t5 = DSQRT(105d0)/21d0, &
  t6 = DSQRT(21d0)/35d0, t7 = 3d0/35d0, t8 = 1d0/3d0, t9 = 1d0/DSQRT(5d0), &
  t10 = 0.6d0, t11 = 21d0/11d0, t12 = -5d0*DSQRT(21d0)/22d0, t13 = 35d0/88d0, &
  t14 = DSQRT(105d0)/44d0, t15 = 70d0/33d0, t16 = -10d0*DSQRT(21d0)/33d0, &
  t17 = -2d0*DSQRT(105d0)/33d0, t18 = 5d0*DSQRT(21d0)/132d0, t19 = 35d0/264d0, &
  t20 = 115d0/44d0, t21 = -5d0/44d0, t22 = -7d0*DSQRT(5d0)/44d0, t23=20d0/11d0,&
  t24 = -15d0/22d0, t25 = 83d0/44d0
 real(kind=8), parameter :: ovlp(21,21) = RESHAPE([&
 1d0,0d0,t1,0d0,t2,0d0,0d0,0d0,0d0,0d0,0d0,t1,0d0,t3,0d0,0d0,0d0,0d0,t2,0d0,0d0,&
 0d0,1d0,0d0,t4,0d0,t2,0d0,0d0,0d0,0d0,0d0,0d0,t5,0d0,t6,0d0,0d0,0d0,0d0,t7,0d0,&
 t1,0d0,1d0,0d0,t4,0d0,0d0,0d0,0d0,0d0,0d0,t8,0d0,t9,0d0,0d0,0d0,0d0,t6,0d0,0d0,&
 0d0,t4,0d0,1d0,0d0,t1,0d0,0d0,0d0,0d0,0d0,0d0,t9,0d0,t8,0d0,0d0,0d0,0d0,t6,0d0,&
 t2,0d0,t4,0d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0,t6,0d0,t5,0d0,0d0,0d0,0d0,t7,0d0,0d0,&
 0d0,t2,0d0,t1,0d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0,t3,0d0,t1,0d0,0d0,0d0,0d0,t2,0d0,&
 0d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0,t5,0d0,t7,0d0,0d0,0d0,0d0,t4,0d0,t6,0d0,0d0,t2,&
 0d0,0d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0,t10,0d0,0d0,0d0,0d0,0d0,0d0,t10,0d0,0d0,0d0,0d0,&
 0d0,0d0,0d0,0d0,0d0,0d0,t5,0d0,1d0,0d0,t5,0d0,0d0,0d0,0d0,t9,0d0,t9,0d0,0d0,t3,&
 0d0,0d0,0d0,0d0,0d0,0d0,0d0,t10,0d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0,t10,0d0,0d0,0d0,0d0,&
 0d0,0d0,0d0,0d0,0d0,0d0,t7,0d0,t5,0d0,1d0,0d0,0d0,0d0,0d0,t6,0d0,t4,0d0,0d0,t2,&
 t1,0d0,t8,0d0,t6,0d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0,t9,0d0,0d0,0d0,0d0,t4,0d0,0d0,&
 0d0,t5,0d0,t9,0d0,t3,0d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0,t9,0d0,0d0,0d0,0d0,t5,0d0,&
 t3,0d0,t9,0d0,t5,0d0,0d0,0d0,0d0,0d0,0d0,t9,0d0,1d0,0d0,0d0,0d0,0d0,t5,0d0,0d0,&
 0d0,t6,0d0,t8,0d0,t1,0d0,0d0,0d0,0d0,0d0,0d0,t9,0d0,1d0,0d0,0d0,0d0,0d0,t4,0d0,&
 0d0,0d0,0d0,0d0,0d0,0d0,t4,0d0,t9,0d0,t6,0d0,0d0,0d0,0d0,1d0,0d0,t8,0d0,0d0,t1,&
 0d0,0d0,0d0,0d0,0d0,0d0,0d0,t10,0d0,t10,0d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,0d0,&
 0d0,0d0,0d0,0d0,0d0,0d0,t6,0d0,t9,0d0,t4,0d0,0d0,0d0,0d0,t8,0d0,1d0,0d0,0d0,t1,&
 t2,0d0,t6,0d0,t7,0d0,0d0,0d0,0d0,0d0,0d0,t4,0d0,t5,0d0,0d0,0d0,0d0,1d0,0d0,0d0,&
 0d0,t7,0d0,t6,0d0,t2,0d0,0d0,0d0,0d0,0d0,0d0,t5,0d0,t4,0d0,0d0,0d0,0d0,1d0,0d0,&
 0d0,0d0,0d0,0d0,0d0,0d0,t2,0d0,t3,0d0,t2,0d0,0d0,0d0,0d0,t1,0d0,t1,0d0,0d0,1d0],&
 [21,21])
 real(kind=8), parameter :: invS(21,21) = RESHAPE([&
 t11,0d0,t12,0d0,t13,0d0,0d0,0d0,0d0,0d0,0d0,t12,0d0,t14,0d0,0d0,0d0,0d0,t13,0d0,0d0,&
 0d0,t15,0d0,t16,0d0,t13,0d0,0d0,0d0,0d0,0d0,0d0,t17,0d0,t18,0d0,0d0,0d0,0d0,t19,0d0,&
 t12,0d0,t20,0d0,t16,0d0,0d0,0d0,0d0,0d0,0d0,t21,0d0,t22,0d0,0d0,0d0,0d0,t18,0d0,0d0,&
 0d0,t16,0d0,t20,0d0,t12,0d0,0d0,0d0,0d0,0d0,0d0,t22,0d0,t21,0d0,0d0,0d0,0d0,t18,0d0,&
 t13,0d0,t16,0d0,t15,0d0,0d0,0d0,0d0,0d0,0d0,t18,0d0,t17,0d0,0d0,0d0,0d0,t19,0d0,0d0,&
 0d0,t13,0d0,t12,0d0,t11,0d0,0d0,0d0,0d0,0d0,0d0,t14,0d0,t12,0d0,0d0,0d0,0d0,t13,0d0,&
 0d0,0d0,0d0,0d0,0d0,0d0,t15,0d0,t17,0d0,t19,0d0,0d0,0d0,0d0,t16,0d0,t18,0d0,0d0,t13,&
 0d0,0d0,0d0,0d0,0d0,0d0,0d0,t23,0d0,t24,0d0,0d0,0d0,0d0,0d0,0d0,t24,0d0,0d0,0d0,0d0,&
 0d0,0d0,0d0,0d0,0d0,0d0,t17,0d0,t25,0d0,t17,0d0,0d0,0d0,0d0,t22,0d0,t22,0d0,0d0,t14,&
 0d0,0d0,0d0,0d0,0d0,0d0,0d0,t24,0d0,t23,0d0,0d0,0d0,0d0,0d0,0d0,t24,0d0,0d0,0d0,0d0,&
 0d0,0d0,0d0,0d0,0d0,0d0,t19,0d0,t17,0d0,t15,0d0,0d0,0d0,0d0,t18,0d0,t16,0d0,0d0,t13,&
 t12,0d0,t21,0d0,t18,0d0,0d0,0d0,0d0,0d0,0d0,t20,0d0,t22,0d0,0d0,0d0,0d0,t16,0d0,0d0,&
 0d0,t17,0d0,t22,0d0,t14,0d0,0d0,0d0,0d0,0d0,0d0,t25,0d0,t22,0d0,0d0,0d0,0d0,t17,0d0,&
 t14,0d0,t22,0d0,t17,0d0,0d0,0d0,0d0,0d0,0d0,t22,0d0,t25,0d0,0d0,0d0,0d0,t17,0d0,0d0,&
 0d0,t18,0d0,t21,0d0,t12,0d0,0d0,0d0,0d0,0d0,0d0,t22,0d0,t20,0d0,0d0,0d0,0d0,t16,0d0,&
 0d0,0d0,0d0,0d0,0d0,0d0,t16,0d0,t22,0d0,t18,0d0,0d0,0d0,0d0,t20,0d0,t21,0d0,0d0,t12,&
 0d0,0d0,0d0,0d0,0d0,0d0,0d0,t24,0d0,t24,0d0,0d0,0d0,0d0,0d0,0d0,t23,0d0,0d0,0d0,0d0,&
 0d0,0d0,0d0,0d0,0d0,0d0,t18,0d0,t22,0d0,t16,0d0,0d0,0d0,0d0,t21,0d0,t20,0d0,0d0,t12,&
 t13,0d0,t18,0d0,t19,0d0,0d0,0d0,0d0,0d0,0d0,t16,0d0,t17,0d0,0d0,0d0,0d0,t15,0d0,0d0,&
 0d0,t19,0d0,t18,0d0,t13,0d0,0d0,0d0,0d0,0d0,0d0,t17,0d0,t16,0d0,0d0,0d0,0d0,t15,0d0,&
 0d0,0d0,0d0,0d0,0d0,0d0,t13,0d0,t14,0d0,t13,0d0,0d0,0d0,0d0,t12,0d0,t12,0d0,0d0,t11],&
 [21,21])

 vecs = MATMUL(TRANSPOSE(rotation), vec21h)
 forall(i = 1:21) rot21h(:,i) = get21h_vector(vecs(:,i))
 rot21h = MATMUL(rot21h, invA_21h)
 deallocate(invA_21h)

 ! Note: the 21H are not orthogonal, so we need the overlap matrix and its inverse
 rot21h = MATMUL(MATMUL(invS, rot21h), ovlp)
end subroutine get_rot_mat_21h

! input a 3D-vector (x,y,z), calculate the values of spherical harmonic 5D functions
pure function get5d_vector(v1) result(v2)
 implicit none
 real(kind=8) :: r2, v2(5)
 real(kind=8), intent(in) :: v1(3)

 r2 = DOT_PRODUCT(v1, v1)
 v2(1) = (3d0*v1(3)*v1(3) - r2)/root3
 v2(2) = 2d0*v1(1)*v1(3)
 v2(3) = 2d0*v1(2)*v1(3)
 v2(4) = v1(1)*v1(1) - v1(2)*v1(2)
 v2(5) = 2d0*v1(1)*v1(2)
 v2 = v2*0.25d0*DSQRT(3.75d0/DATAN(1d0))/r2
end function get5d_vector

! input a 3D-vector (x,y,z), calculate the values of Cartesian 6D functions
pure function get6d_vector(v1) result(v2)
 implicit none
 real(kind=8) :: x, y, z, x2, y2, z2, r2, v2(6)
 real(kind=8), intent(in) :: v1(3)

 x = v1(1); y = v1(2); z = v1(3)
 x2 = x*x; y2 = y*y; z2 = z*z
 r2 = x2 + y2 + z2
 v2 = [x2, y2, z2, root3*x*y, root3*x*z, root3*y*z]
 v2 = 0.5d0*v2*DSQRT(1.25d0/DATAN(1d0))/r2
end function get6d_vector

pure function get7f_vector(v1) result(v2)
 implicit none
 real(kind=8) :: x2, y2, z2, r2, r3, fivez2, v2(7)
 real(kind=8), intent(in) :: v1(3)

 x2 = v1(1)*v1(1); y2 = v1(2)*v1(2); z2 = v1(3)*v1(3)
 r2 = x2 + y2 + z2
 r3 = r2*DSQRT(r2)
 fivez2 = 5d0*z2
 v2(1) = root2*(fivez2 - 3d0*r2)*v1(3)
 v2(2) = root3*(fivez2 - r2)*v1(1)
 v2(3) = root3*(fivez2 - r2)*v1(2)
 v2(4) = root30*(x2 - y2)*v1(3)
 v2(5) = 2d0*root30*v1(1)*v1(2)*v1(3)
 v2(6) = root5*v1(1)*(x2 - 3d0*y2)
 v2(7) = root5*v1(2)*(3d0*x2 - y2)
 v2 = v2*0.125d0*DSQRT(3.5d0/DATAN(1d0))/r3
end function get7f_vector

pure function get10f_vector(v1) result(v2)
 implicit none
 real(kind=8) :: x, y, z, x2, y2, z2, r2, r3, v2(10)
 real(kind=8), intent(in) :: v1(3)

 x = v1(1); y = v1(2); z = v1(3)
 x2 = v1(1)*v1(1); y2 = v1(2)*v1(2); z2 = v1(3)*v1(3)
 r2 = x2 + y2 + z2
 r3 = r2*DSQRT(r2)
 v2 = [x2*x, y2*y, z2*z, x*y2, x2*y, x2*z, x*z2, y*z2, y2*z, x*y*z]
 v2(4:9) = root5*v2(4:9)
 v2(10) = root15*v2(10)
 v2 = 0.5d0*v2*DSQRT(1.75/DATAN(1d0))/r3
end function get10f_vector

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

pure function get15g_vector(v1) result(v2)
 implicit none
 integer :: i
 integer, parameter :: idx1(6) = [2,4,6,9,13,14]
 integer, parameter :: idx2(3) = [3,10,12]
 integer, parameter :: idx3(3) = [7,8,11]
 real(kind=8) :: x, y, z, x2, y2, z2, x3, y3, z3, r2, v2(15)
 real(kind=8), intent(in) :: v1(3)

 x = v1(1); y = v1(2); z = v1(3)
 x2 = x*x; y2 = y*y; z2 = z*z
 r2 = x2 + y2 + z2
 x3 = x2*x; y3 = y2*y; z3 = z2*z
 v2 = [z2*z2, y*z3, y2*z2, y3*z, y2*y2, x*z3, x*y*z2, x*y2*z, x*y3, x2*z2, &
       x2*y*z, x2*y2, x3*z, x3*y, x2*x2]
 forall(i = 1:6) v2(idx1(i)) = v2(idx1(i))*DSQRT(7d0)
 forall(i = 1:3) v2(idx2(i)) = v2(idx2(i))*DSQRT(105d0)/3d0
 forall(i = 1:3) v2(idx3(i)) = v2(idx3(i))*DSQRT(35d0)
 v2 = 0.75d0*v2/(r2*r2*DSQRT(DATAN(1d0)))
end function get15g_vector

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

pure function get21h_vector(v1) result(v2)
 implicit none
 integer :: i
 integer, parameter :: idx1(6) = [2,5,7,11,19,20]
 integer, parameter :: idx2(6) = [3,4,12,15,16,18]
 integer, parameter :: idx3(3) = [8,10,17]
 integer, parameter :: idx4(3) = [9,13,14]
 real(kind=8) :: x, y, z, x2, y2, z2, x3, y3, z3, x4, y4, z4, r2, r5, v2(21)
 real(kind=8), intent(in) :: v1(3)

 x = v1(1); y = v1(2); z = v1(3)
 x2 = x*x; y2 = y*y; z2 = z*z
 r2 = x2 + y2 + z2
 r5 = r2*r2*DSQRT(r2)
 x3 = x2*x; y3 = y2*y; z3 = z2*z
 x4 = x2*x2; y4 = y2*y2; z4 = z2*z2
 v2 = [z*z4, y*z4, y2*z3, y3*z2, y4*z, y4*y, x*z4, x*y*z3, x*y2*z2, x*y3*z, &
       x*y4, x2*z3, x2*y*z2, x2*y2*z, x2*y3, x3*z2, x3*y*z, x3*y2, x4*z, &
       x4*y, x4*x]
 forall(i = 1:6) v2(idx1(i)) = v2(idx1(i))*3d0
 forall(i = 1:6) v2(idx2(i)) = v2(idx2(i))*root21
 forall(i = 1:3) v2(idx3(i)) = v2(idx3(i))*root7*3d0
 forall(i = 1:3) v2(idx4(i)) = v2(idx4(i))*root105
 v2 = 0.25d0*v2*DSQRT(11d0/DATAN(1d0))/r5
end function get21h_vector

end module bas_rot

! Get the wave function of a rotated molecule. The original coordinates are
!  provided in fchname, while the new coordinates are provided in coor_file
!  (which is a .xyz or .gjf file).
! Note: the one-to-one correspondence of element symbols in two files will not
!  be checked, so it is assumed that you already ensure the correspondence.
! Assuming the input file is a.fch, the output file would be a_r.fch
subroutine rotate_atoms_wfn(fchname, coor_file)
 implicit none
 integer :: i, natom
 real(kind=8), allocatable :: coor(:,:)
 character(len=2), allocatable :: elem(:)
 character(len=240) :: new_fch
 character(len=240), intent(in) :: fchname, coor_file
!f2py intent(in) :: fchname, coor_file

 call find_specified_suffix(fchname, '.fch', i)
 ! the filename to store rotated molecule and wfn
 new_fch = fchname(1:i-1)//'_r.fch'

 call check_natom_eq_in_fch_and_gjf_or_xyz(fchname, coor_file, natom)
 allocate(elem(natom), coor(3,natom))
 call read_elem_and_coor_from_file(coor_file, natom, elem, coor)
 deallocate(elem)

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
  beta_coeff, tot_dm, spin_dm, check_uhf_in_fch, read_fch, coor0=>coor, &
  free_arrays_in_fch_content
 use bas_rot, only: get_invA, get_rot_mat_5d, get_rot_mat_6d, get_rot_mat_7f, &
  get_rot_mat_10f, get_rot_mat_9g, get_rot_mat_15g, get_rot_mat_11h, &
  get_rot_mat_21h
 implicit none
 integer :: i, j, nif1
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 integer, allocatable :: no_mark(:), mark(:,:)
 real(kind=8) :: rmsd_v, trans1(3), trans2(3), rotation(3,3), rot5d(5,5), &
  rot6d(6,6)
 real(kind=8), allocatable :: rot7f(:,:), rot10f(:,:), rot9g(:,:), rot15g(:,:),&
  rot11h(:,:), rot21h(:,:)
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

 if(rmsd_v > 0.03d0) then
  write(6,'(/,A)') 'Warning in subroutine rotate_atoms_wfn2: RMSD value > 0.03.'
  write(6,'(A)') 'Anyway, the program will continue, but you should be aware of&
                 & what you are'
  write(6,'(A,I0)') 'doing. fchname='//TRIM(fchname)//', natom=', natom
  write(6,'(A,F12.6)') 'RMSD=', rmsd_v
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

 if(ANY(IABS(shell_type)>5)) then
  write(6,'(/,A)') 'ERROR in subroutine rotate_atoms_wfn2: basis functions with&
                   & angular momentum'
  write(6,'(A)') 'higher than H are not supported currently.'
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
  allocate(rot7f(7,7))
  call get_rot_mat_7f(rotation, rot7f)
  allocate(rtmp(7,nif1), source=0d0)
  do i = 1, no_mark(4), 1
   j = mark(i,4)
   rtmp = coeff(j:j+6,:)
   call dgemm('N', 'N', 7, nif1, 7, 1d0,rot7f,7, rtmp,7, 0d0,coeff(j:j+6,:),7)
  end do ! for i
  deallocate(rot7f, rtmp)
 end if

 ! 10F
 if(no_mark(5) > 0) then
  call get_invA(10)
  allocate(rot10f(10,10))
  call get_rot_mat_10f(rotation, rot10f)
  allocate(rtmp(10,nif1), source=0d0)
  do i = 1, no_mark(5), 1
   j = mark(i,5)
   rtmp = coeff(j:j+9,:)
   call dgemm('N', 'N', 10, nif1, 10, 1d0,rot10f,10, rtmp,10, 0d0,coeff(j:j+9,:),10)
  end do ! for i
  deallocate(rot10f, rtmp)
 end if

 ! 9G
 if(no_mark(6) > 0) then
  call get_invA(9)
  allocate(rot9g(9,9))
  call get_rot_mat_9g(rotation, rot9g)
  allocate(rtmp(9,nif1), source=0d0)
  do i = 1, no_mark(6), 1
   j = mark(i,6)
   rtmp = coeff(j:j+8,:)
   call dgemm('N', 'N', 9, nif1, 9, 1d0,rot9g,9, rtmp,9, 0d0,coeff(j:j+8,:),9)
  end do ! for i
  deallocate(rot9g, rtmp)
 end if

 ! 15G
 if(no_mark(7) > 0) then
  call get_invA(15)
  allocate(rot15g(15,15))
  call get_rot_mat_15g(rotation, rot15g)
  allocate(rtmp(15,nif1), source=0d0)
  do i = 1, no_mark(7), 1
   j = mark(i,7)
   rtmp = coeff(j:j+14,:)
   call dgemm('N', 'N', 15, nif1, 15, 1d0,rot15g,15, rtmp,15, 0d0,coeff(j:j+14,:),15)
  end do ! for i
  deallocate(rot15g, rtmp)
 end if

 ! 11H
 if(no_mark(8) > 0) then
  call get_invA(11)
  allocate(rot11h(11,11))
  call get_rot_mat_11h(rotation, rot11h)
  allocate(rtmp(11,nif1), source=0d0)
  do i = 1, no_mark(8), 1
   j = mark(i,8)
   rtmp = coeff(j:j+10,:)
   call dgemm('N', 'N', 11, nif1, 11, 1d0,rot11h,11, rtmp,11, 0d0,coeff(j:j+10,:),11)
  end do ! for i
  deallocate(rot11h, rtmp)
 end if

 ! 21H
 if(no_mark(9) > 0) then
  call get_invA(21)
  allocate(rot21h(21,21))
  call get_rot_mat_21h(rotation, rot21h)
  allocate(rtmp(21,nif1), source=0d0)
  do i = 1, no_mark(9), 1
   j = mark(i,9)
   rtmp = coeff(j:j+20,:)
   call dgemm('N', 'N', 21, nif1, 21, 1d0,rot21h,21, rtmp,21, 0d0,coeff(j:j+20,:),21)
  end do ! for i
  deallocate(rot21h, rtmp)
 end if

 deallocate(no_mark, mark)
 allocate(alpha_coeff(nbf,nif), source=coeff(:,1:nif))
 if(is_uhf) allocate(beta_coeff(nbf,nif), source=coeff(:,nif+1:nif1))
 deallocate(coeff)

 ! currently the tot_dm and spin_dm are not rotated, so set to 0
 tot_dm = 0d0
 if(is_uhf) spin_dm = 0d0

 call write_fch(new_fch)
 call free_arrays_in_fch_content()
end subroutine rotate_atoms_wfn2

! Get the wavefunction (MO coefficients, actually) of a molecule after some
!  atoms are permuted (permuted coordinates are stored in coor_file).
subroutine permute_atoms_wfn(fchname, coor_file)
 implicit none
 integer :: i, natom
 integer, allocatable :: idx(:)
 real(kind=8), allocatable :: coor1(:,:), coor2(:,:)
 character(len=2), allocatable :: elem(:)
 character(len=240) :: new_fch
 character(len=240), intent(in) :: fchname, coor_file
!f2py intent(in) :: fchname, coor_file

 call find_specified_suffix(fchname, '.fch', i)
 ! the filename to store molecule and wfn after permutation
 new_fch = fchname(1:i-1)//'_p.fch'

 call check_natom_eq_in_fch_and_gjf_or_xyz(fchname, coor_file, natom)

 allocate(coor1(3,natom), coor2(3,natom), idx(natom))
 call read_coor_from_fch(fchname, natom, coor1)
 allocate(elem(natom))
 call read_elem_and_coor_from_file(coor_file, natom, elem, coor2)
 deallocate(elem)
 call get_atom_map_idx_from_coor(natom, coor1, coor2, idx)
 deallocate(coor1, coor2)

 call permute_atoms_wfn2(fchname, natom, idx, new_fch)
 deallocate(idx)
end subroutine permute_atoms_wfn

! Get the wavefunction (MO coefficients, actually) of a molecule after some
!  atoms are permuted (permuted indices are stored in the integer array idx).
! The updated wavefunction is stored in new_fch.
subroutine permute_atoms_wfn2(fchname, natm, idx, new_fch)
 use fch_content
 implicit none
 integer :: i, j, k, m
 integer, intent(in) :: natm
!f2py intent(in) :: natm
 integer, intent(in) :: idx(natm)
!f2py depend(natm) :: idx
!f2py intent(in) :: idx
 integer, allocatable :: i1(:), begin_idx(:), bas_idx(:), shl2atmp(:), &
  kfirst1(:,:), klast1(:,:)
 real(kind=8), allocatable :: prim_exp1(:), contr_coeff1(:), contr_coeff_sp1(:),&
  r2(:,:)
 character(len=240), intent(in) :: fchname, new_fch
!f2py intent(in) :: fchname, new_fch
 logical :: has_sp

 if(natm < 1) then
  write(6,'(/,A)') 'ERROR in subroutine permute_atoms_wfn2: input natm<1.'
  stop
 else if(natm==1 .and. idx(1)/=1) then
  write(6,'(/,A)') 'ERROR in subroutine permute_atoms_wfn2: natm=1 but idx/=1.'
  write(6,'(A)') 'Input parameters/arguments are wrong.'
  stop
 end if

 ! check whether the integer array idx = {1,2,3,...}
 allocate(i1(natm))
 forall(i = 1:natm) i1(i) = i
 if(ALL(i1 == idx)) then
  deallocate(i1)
  call copy_file(fchname, new_fch)
  return
 end if
 ! now idx is not {1,2,3,...}

 call check_uhf_in_fch(fchname, is_uhf)
 call read_fch(fchname, is_uhf)

 has_sp = .false.
 if(allocated(contr_coeff_sp)) has_sp = .true.

 if(natm /= natom) then
  write(6,'(/,A)') 'ERROR in subroutine permute_atoms_wfn2: natm/=natom.'
  write(6,'(A,2I5)') 'natm, natom=', natm, natom
  stop
 end if
 ! natom is used in the module fch_content, and natm is used for idx. They are
 ! supposed to be equal.

 ! permute arrays ielem and coor
 i1 = ielem
 allocate(r2(3,natom), source=coor)
 forall(i = 1:natom)
  ielem(i) = i1(idx(i))
  coor(:,i) = r2(:,idx(i))
 end forall
 deallocate(i1, r2)

 ! find beginning indices of primitive functions of each atom
 allocate(begin_idx(natom+1), source=0)
 j = 2; begin_idx(1) = 1
 if(ncontr > 1) begin_idx(2) = 1

 do i = 1, ncontr, 1
  begin_idx(j) = begin_idx(j) + prim_per_shell(i)
  if(i < ncontr) then
   if(j == shell2atom_map(i+1)) then
    j = j + 1
    begin_idx(j) = begin_idx(j-1)
   end if
  else ! i = ncontr
   exit
  end if
 end do ! for i

 if(begin_idx(natom+1) /= nprim+1) then
  write(6,'(/,A)') 'ERROR in subroutine permute_atoms_wfn2: begin_idx(natom+1) &
                   &/= nprim+1.'
  write(6,'(A,2I5)') 'begin_idx(natom+1), nprim=', begin_idx(natom+1), nprim
  stop
 end if

 ! permute arrays prim_exp, contr_coeff, and (possible) contr_coeff_sp1
 allocate(prim_exp1(nprim), source=prim_exp)
 allocate(contr_coeff1(nprim), source=contr_coeff)
 j = 1
 do i = 1, natom, 1
  k = begin_idx(idx(i))
  m = begin_idx(idx(i)+1) - k
  prim_exp(j:j+m-1) = prim_exp1(k:k+m-1)
  contr_coeff(j:j+m-1) = contr_coeff1(k:k+m-1)
  j = j + m
 end do ! for i
 deallocate(prim_exp1, contr_coeff1)

 if(has_sp) then
  allocate(contr_coeff_sp1(nprim), source=contr_coeff_sp)
  j = 1
  do i = 1, natom, 1
   k = begin_idx(idx(i))
   m = begin_idx(idx(i)+1) - k
   contr_coeff_sp(j:j+m-1) = contr_coeff_sp1(k:k+m-1)
   j = j + m
  end do ! for i
  deallocate(contr_coeff_sp1)
 end if

 ! find basis function permutation indices
 call get_bas_begin_idx_from_shltyp(ncontr, shell_type, shell2atom_map, natom, &
                                    begin_idx)
 allocate(bas_idx(nbf), source=0)
 allocate(i1(nbf))
 forall(i = 1:nbf) i1(i) = i
 j = 1
 do i = 1, natom, 1
  k = begin_idx(idx(i))
  m = begin_idx(idx(i)+1) - k
  bas_idx(j:j+m-1) = i1(k:k+m-1)
  j = j + m
 end do ! for i
 deallocate(i1)

 ! permute MO coefficients, total density matrix and (possible) spin density
 ! matrix
 allocate(r2(nbf,nbf), source=0d0)
 r2(:,1:nif) = alpha_coeff
 forall(i = 1:nbf) alpha_coeff(i,:) = r2(bas_idx(i),1:nif)
 r2 = tot_dm
 forall(i=1:nbf, j=1:nbf) tot_dm(i,j) = r2(bas_idx(i),bas_idx(j))

 if(is_uhf) then
  r2(:,1:nif) = beta_coeff
  forall(i = 1:nbf) beta_coeff(i,:) = r2(bas_idx(i),1:nif)
  if(allocated(spin_dm)) then
   r2 = spin_dm
   forall(i=1:nbf, j=1:nbf) spin_dm(i,j) = r2(bas_idx(i),bas_idx(j))
  end if
 end if

 deallocate(bas_idx, r2)

 ! permute arrays shell_type, prim_per_shell and shell2atom_map
 begin_idx = 0; begin_idx(1) = 1; begin_idx(natom+1) = ncontr+1
 j = 2
 do i = 2, ncontr, 1
  if(shell2atom_map(i-1)+1 == shell2atom_map(i)) then
   begin_idx(j) = i
   j = j + 1
  end if
 end do ! for i

 allocate(i1(ncontr), source=shell_type)
 allocate(bas_idx(ncontr), source=prim_per_shell)
 allocate(shl2atmp(ncontr), source=shell2atom_map)
 j = 1
 do i = 1, natom, 1
  k = begin_idx(idx(i))
  m = begin_idx(idx(i)+1) - k
  shell_type(j:j+m-1) = i1(k:k+m-1)
  prim_per_shell(j:j+m-1) = bas_idx(k:k+m-1)
  shell2atom_map(j:j+m-1) = shl2atmp(k:k+m-1)
  j = j + m
 end do ! for i

 deallocate(i1, bas_idx, begin_idx)
 allocate(i1(natom))
 forall(i = 1:natom) i1(idx(i)) = i
 shl2atmp = shell2atom_map
 forall(i = 1:ncontr) shell2atom_map(i) = i1(shl2atmp(i))
 deallocate(shl2atmp)

 if(LenNCZ > 0) then ! If ECP/PP is used, permute ECP data
  i1 = Lmax
  allocate(bas_idx(natom), source=LPSkip)
  allocate(kfirst1(natom,10), source=KFirst)
  allocate(klast1(natom,10), source=KLast)
  allocate(prim_exp1(natom), source=RNFroz)
  forall(i = 1:natom)
   Lmax(i) = i1(idx(i))
   LPSkip(i) = bas_idx(idx(i))
   KFirst(i,:) = kfirst1(idx(i),:)
   KLast(i,:) = klast1(idx(i),:)
   RNFroz(i) = prim_exp1(idx(i))
  end forall
  deallocate(i1, bas_idx, kfirst1, klast1, prim_exp1)
  ! There is no need to permute arrays NLP, CLP1, CLP2 and ZLP
 end if

 call write_fch(new_fch)
 call free_arrays_in_fch_content()
end subroutine permute_atoms_wfn2

! get the basis function beginning indices of each atom
subroutine get_bas_begin_idx_from_shltyp(ncontr, shell_type, shell2atom_map, &
                                         natom, begin_idx)
 implicit none
 integer :: i, j
 integer, parameter :: shltyp2nbas(-5:5) = [11,9,7,5,4,1,3,6,10,15,21]
 integer, intent(in) :: ncontr, natom
 integer, intent(in) :: shell_type(ncontr), shell2atom_map(ncontr)
 integer, intent(out) :: begin_idx(natom+1)

 begin_idx = 0
 j = 2; begin_idx(1) = 1
 if(ncontr > 1) begin_idx(2) = 1

 do i = 1, ncontr, 1
  begin_idx(j) = begin_idx(j) + shltyp2nbas(shell_type(i))
  if(i < ncontr) then
   if(j == shell2atom_map(i+1)) then
    j = j + 1
    begin_idx(j) = begin_idx(j-1)
   end if
  else ! i = ncontr
   exit
  end if
 end do ! for i
end subroutine get_bas_begin_idx_from_shltyp

! Calculate the RMSD value of two molecules (in .xyz/.gjf/.fch files).
! This subroutine assumes that the atomic labels are in one-to-one correspondence
! in two files.
subroutine rmsd_wrapper(fname1, fname2, reorder, rmsd_v)
 use periodic_table, only: write_xyz
 implicit none
 integer :: i, natom, natom1, natom2
 real(kind=8) :: trans1(3), trans2(3), rotation(3,3)
 real(kind=8), allocatable :: coor1(:,:), coor2(:,:)
 real(kind=8), intent(out) :: rmsd_v
!f2py intent(out) :: rmsd_v
 character(len=2), allocatable :: elem1(:), elem2(:)
 character(len=240) :: fname
 character(len=240), intent(in) :: fname1, fname2
!f2py intent(in) :: fname1, fname2
 logical :: relabel
 logical, intent(in), optional :: reorder
!f2py intent(in), optional :: reorder

 i = INDEX(fname2, '.', back=.true.)
 fname = fname2(1:i-1)//'_new.xyz'

 call read_natom_from_file(fname1, natom1)
 call read_natom_from_file(fname2, natom2)
 if(natom1 /= natom2) then
  write(6,'(/,A)') 'ERROR in subroutine rmsd_wrapper: the number of atoms are n&
                   &ot equal in two files.'
  write(6,'(A)') 'fname1='//TRIM(fname1)
  write(6,'(A)') 'fname2='//TRIM(fname2)
  stop
 end if

 natom = natom1
 allocate(elem1(natom), elem2(natom), coor1(3,natom), coor2(3,natom))
 call read_elem_and_coor_from_file(fname1, natom, elem1, coor1)
 call read_elem_and_coor_from_file(fname2, natom, elem2, coor2)

 relabel = .false.
 if(present(reorder)) relabel = reorder

 if(relabel) then ! reordering atomic labels requested
  call rmsd_reorder(natom, elem1, elem2, coor1, coor2)
 else             ! keep atomic labels unchanged

  if(.not. ALL(elem1 == elem2)) then
   write(6,'(/,A)') 'ERROR in subroutine rmsd_wrapper: elements in two files are&
                    & not identical.'
   write(6,'(A)') 'fname1='//TRIM(fname1)
   write(6,'(A)') 'fname2='//TRIM(fname2)
   deallocate(elem1, elem2, coor1, coor2)
   stop
  end if
  call rmsd(natom, coor2, coor1, rmsd_v, trans2, trans1, rotation)
 end if

 deallocate(elem1, coor1)
 call write_xyz(natom, elem2, coor2, fname)
 deallocate(elem2, coor2)
end subroutine rmsd_wrapper

! Calculate the RMSD value of two sets of Cartesian coordindates.
! This subroutine assumes that the atomic labels are in one-to-one correspondence
! in two files.
subroutine rmsd(natom, coor1, coor2, rmsd_v, trans1, trans2, rotation)
 implicit none
 integer :: i, lwork
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
 forall(i = 1:natom) coor1(:,i) = coor1(:,i) - trans1
end subroutine rmsd

! Calculate the RMSD value of two molecules, where permutations of atomic labels
! in the 2nd molecule will be performed.
subroutine rmsd_reorder(natom, elem1, elem2, coor1, coor2)
 use fch_content, only: elem2nuc, ram
 implicit none
 integer :: i
 integer, intent(in) :: natom
 integer, allocatable :: nuc1(:), nuc2(:), idx1(:), idx2(:), idx3(:)
 real(kind=8) :: trans1(3), trans2(3), rot1(3,3), rot2(3,3)
 real(kind=8), intent(in) :: coor1(3,natom)
 real(kind=8), intent(inout) :: coor2(3,natom)
 real(kind=8), allocatable :: aw(:), coor(:,:)
 character(len=2), intent(in) :: elem1(natom)
 character(len=2), intent(inout) :: elem2(natom)
 character(len=2), allocatable :: elem(:)
 character(len=240) :: gjfname1, gjfname2

 write(6,'(/,A)') 'ERROR in subroutine rmsd_reorder: not implemented yet.'
 stop
 allocate(nuc1(natom), nuc2(natom))
 forall(i = 1:natom) nuc1(i) = elem2nuc(elem1(i))
 forall(i = 1:natom) nuc2(i) = elem2nuc(elem2(i))

 allocate(idx1(natom), idx2(natom))
 call sort_int_array(natom, nuc1, .true., idx1)
 call sort_int_array(natom, nuc2, .true., idx2)
 deallocate(nuc2)
 ! nuc2 = nuc1 after calling sort_int_array, so it can be deallocated

 allocate(coor(3,natom), source=coor2)
 forall(i = 1:natom) coor2(:,i) = coor(:,idx2(i))
 forall(i = 1:natom) coor(:,i) = coor1(:,idx1(i))
 allocate(elem(natom), source=elem2)
 forall(i = 1:natom) elem2(i) = elem(idx2(i))
 deallocate(elem)

 call move_coor_to_com(natom, nuc1, coor, trans1)
 call move_coor_to_com(natom, nuc1, coor2, trans2)

 allocate(aw(natom))
 forall(i = 1:natom) aw(i) = ram(nuc1(i))
 call rot_to_principal_axis(natom, aw, coor, rot1)
 call rot_to_principal_axis(natom, aw, coor2, rot2)
 deallocate(aw)

! allocate(idx3(natom))
! call permute_atoms_of_same_elem(natom, nuc1, coor, coor2, idx3)

! gjfname1 = '1.gjf'
! gjfname2 = '2.gjf'
! call write_gjf(gjfname1, 0, 1, natom, elem2, coor)
! call write_gjf(gjfname2, 0, 1, natom, elem2, coor2)
end subroutine rmsd_reorder

! move the geometry center (centeroid) of a set of coordinates to the origin
! (0,0,0), return the corresponding translational vector and translated coordinates
subroutine move_coor_to_center(natom, coor, trans)
 implicit none
 integer :: i
 integer, intent(in) :: natom
 real(kind=8), intent(inout) :: coor(3,natom)
 real(kind=8), intent(out) :: trans(3)

 forall(i = 1:3) trans(i) = -SUM(coor(i,:))
 trans = trans/DBLE(natom)
 forall(i = 1:natom) coor(:,i) = coor(:,i) + trans
end subroutine move_coor_to_center

! move the center of mass (COM) of a set of coordinates to the origin (0,0,0),
! return the corresponding translational vector and translated coordinates
subroutine move_coor_to_com(natom, nuc, coor, trans)
 use fch_content, only: ram
 implicit none
 integer :: i
 integer, intent(in) :: natom
 integer, intent(in) :: nuc(natom)
 real(kind=8), intent(inout) :: coor(3,natom)
 real(kind=8), intent(out) :: trans(3)
 real(kind=8), allocatable :: aw(:)

 allocate(aw(natom))
 forall(i = 1:natom) aw(i) = ram(nuc(i))
 forall(i = 1:3) trans(i) = DOT_PRODUCT(aw, coor(i,:))
 trans = -trans/SUM(aw)
 deallocate(aw)
 forall(i = 1:natom) coor(:,i) = coor(:,i) + trans
end subroutine move_coor_to_com

! rotate the molecule into its principal axis
! Note: it is assumed that COM of this molecule has been moved to (0,0,0).
subroutine rot_to_principal_axis(natom, aw, coor, rot)
 implicit none
 integer :: i
 integer, intent(in) :: natom
 real(kind=8) :: w(3)
 real(kind=8), intent(in) :: aw(natom) ! relative atomic weight
 real(kind=8), intent(inout) :: coor(3,natom)
 real(kind=8), intent(out) :: rot(3,3)
 real(kind=8), allocatable:: coor2(:,:)

 call get_inertia_tensor(natom, aw, coor, rot)
 call diag_get_e_and_vec(3, rot, w)
 write(6,'(3F12.8)') w
 write(6,'(3F12.8)') rot
 allocate(coor2(3,natom), source=coor)
 coor = MATMUL(TRANSPOSE(rot), coor2)
 deallocate(coor2)
end subroutine rot_to_principal_axis

! Calculate the inertia tensor
subroutine get_inertia_tensor(natom, aw, coor, it)
 implicit none
 integer :: i, j
 integer, intent(in) :: natom
 real(kind=8), intent(in) :: aw(natom) ! relative atomic weight
 real(kind=8), intent(in) :: coor(3,natom)
 real(kind=8), intent(out) :: it(3,3)
 real(kind=8), allocatable :: coor2(:,:)

 allocate(coor2(3,natom))
 forall(i=1:3, j=1:natom) coor2(i,j) = coor(i,j)*coor(i,j)
 it(1,1) = DOT_PRODUCT(aw, coor2(2,:)+coor2(3,:))
 it(2,2) = DOT_PRODUCT(aw, coor2(1,:)+coor2(3,:))
 it(3,3) = DOT_PRODUCT(aw, coor2(1,:)+coor2(2,:))

 forall(i = 1:natom)
  coor2(1,i) = coor(1,i)*coor(2,i)
  coor2(2,i) = coor(1,i)*coor(3,i)
  coor2(3,i) = coor(2,i)*coor(3,i)
 end forall

 it(2,1) = DOT_PRODUCT(aw, coor2(1,:))
 it(3,1) = DOT_PRODUCT(aw, coor2(2,:))
 it(3,2) = DOT_PRODUCT(aw, coor2(3,:))
 deallocate(coor2)

 it(1,2) = it(2,1); it(1,3) = it(3,1); it(2,3) = it(3,2)
end subroutine get_inertia_tensor

! permute atoms of the same elements to achieve maximum overlap
! Note: the input nuc should have beend sorted in ascending order.
subroutine permute_atoms_of_same_elem(natom, nuc, coor1, coor2, idx)
 implicit none
 integer :: i, j, k, k1, k2, nelem
 integer, intent(in) :: natom
 integer, intent(in) :: nuc(natom)
 integer, intent(out) :: idx(natom)
 integer, allocatable :: natom_per_elem(:)
 real(kind=8), intent(in) :: coor1(3,natom)
 real(kind=8), intent(inout) :: coor2(3,natom)

 forall(i = 1:natom) idx(i) = i
 call find_nelem_in_nuc(natom, nuc, nelem)
 allocate(natom_per_elem(nelem))
 call find_natom_per_elem(natom, nuc, nelem, natom_per_elem)

 k1 = 1
 do i = 1, nelem, 1
  k = natom_per_elem(i)
  if(k == 1) cycle
  k2 = k1 + k - 1
  call find_max_match_for_points(k, coor1(:,k1:k2), coor2(:,k1:k2), idx(k1:k2))
  k1 = k1 + k
 end do ! for i
end subroutine permute_atoms_of_same_elem

! find maximum match/overlap for a set of points
subroutine find_max_match_for_points(n, coor1, coor2, idx)
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 integer, intent(inout) :: idx(n)
 integer, allocatable :: matched(:)
 real(kind=8) :: vtmp(3), sum_dis0, sum_dis
 real(kind=8), intent(in) :: coor1(3,n)
 real(kind=8), intent(inout) :: coor2(3,n)
 real(kind=8), external :: calc_sum_dis_for_points
 real(kind=8), allocatable :: dis0(:,:), dis(:,:)

 allocate(dis0(n,n))
 do i = 1, n, 1
  do j = 1, n, 1
   vtmp = coor1(:,i) - coor2(:,j)
   dis0(j,i) = DSQRT(DOT_PRODUCT(vtmp, vtmp))
  end do ! for j
 end do ! for i

 allocate(dis(n,n), matched(n))

 do i = 1, n, 1
  dis = dis0
  matched = 0

  call find_max_match_for_partial_points(i, n, n, dis, matched)
  call find_max_match_for_partial_points(1, i-1, n, dis, matched)
  sum_dis = calc_sum_dis_for_points(n, dis0, matched)
  write(6,'(10I3)') matched
  write(6,'(A,F18.8)') 'sum_dis = ', sum_dis

  if(i == 1) then
   sum_dis0 = sum_dis
  else
   if(sum_dis < sum_dis0) sum_dis0 = sum_dis
  end if
 end do ! for i

 stop
end subroutine find_max_match_for_points

subroutine find_max_match_for_partial_points(k1, k2, n, dis, matched)
 implicit none
 integer :: i, j(1), k(1)
 integer, intent(in) :: k1, k2, n
 integer, intent(inout) :: matched(n)
 real(kind=8), intent(inout) :: dis(n,n)

 do i = k1, k2, 1
  k = MAXLOC(dis(:,i))

  do while(.true.)
   j = MINLOC(dis(:,i))
   if(matched(i) > 0) then
    dis(j(1),i) = dis(k(1),i) + 1d0
    k = j
   else
    matched(i) = j(1)
    exit
   end if
  end do ! for while

 end do ! for i
end subroutine find_max_match_for_partial_points

function calc_sum_dis_for_points(n, dis, matched) result(sum_dis)
 implicit none
 integer :: i
 integer, intent(in) :: n
 integer, intent(in) :: matched(n)
 real(kind=8) :: sum_dis
 real(kind=8), intent(in) :: dis(n,n)

 sum_dis = 0d0
 do i = 1, n, 1
  sum_dis = sum_dis + dis(matched(i),i)
 end do ! for i
end function calc_sum_dis_for_points

! .chk -> .chk, a wrapper of {formchk, mirror_wfn, and unfchk}
subroutine mirror_c2c(chkname)
 use util_wrapper, only: formchk, unfchk
 implicit none
 integer :: i
 character(len=240) :: fchname1, fchname2
 character(len=240), intent(in) :: chkname
!f2py intent(in) :: chkname

 call formchk(chkname)
 i = INDEX(chkname, '.chk', back=.true.)
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
  is_uhf, coor, tot_dm, spin_dm, check_uhf_in_fch, read_fch, free_arrays_in_fch_content
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
 call find_specified_suffix(fchname, '.fch', i)
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
 call free_arrays_in_fch_content()
end subroutine mirror_wfn

! calculate weights of 1D Lagrange's interpolation
subroutine calc_1d_lagrange_w(n, x, weight)
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 real(kind=8) :: p, r
 real(kind=8), parameter :: diff = 1d-7
 real(kind=8), intent(in) :: x(n)
 real(kind=8), intent(out) :: weight(n-1)
 real(kind=8), allocatable :: numerator(:), denominator(:), x_xi(:)

 select case(n)
 case(1)
  write(6,'(/,A)') 'ERROR in subroutine calc_1d_lagrange_w: invalid n=1.'
  write(6,'(A)') 'n>1 is required.'
  stop
 case(2)
  weight(1) = 1d0 ! we have no choice but 1.0
  return
 case default
  weight = 0d0   ! initialization
 end select

 ! if some x(i) is extremely close to x(n), set its weight as 1.0 and return
! do i = 1, n-1, 1
!  r = x(n) - x(i)
!  if(DABS(r) < diff) then
!   weight(i) = 1d0
!   return
!  end if
! end do ! for i

 allocate(x_xi(n-1))
 forall(i = 1:n-1) x_xi(i) = x(n) - x(i)
 p = PRODUCT(x_xi)
 deallocate(x_xi)

 allocate(numerator(n-1))
 do i = 1, n-1, 1
  r = x(n) - x(i)

  if(DABS(r) > diff) then
   numerator(i) = p/r
  else
   numerator(i) = 1d0
   do j = 1, n-1, 1
    if(j == i) cycle
    numerator(i) = numerator(i)*(x(n) - x(j))
   end do ! for j
  end if
 end do ! for i

 allocate(denominator(n-1), source=1d0)
 do i = 1, n-1, 1
  do j = 1, n-1, 1
   if(j == i) cycle
   denominator(i) = denominator(i)*(x(i) - x(j))
  end do ! for j
 end do ! for i

 weight = numerator/denominator
 deallocate(numerator, denominator)
end subroutine calc_1d_lagrange_w

! Perform Grassmann interpolation to generate occupied MOs of a new geometry,
!  using known occupied MOs of known geometries.
! nmo: Usually the number of occupied MOs, i.e. na for R(O)HF, na+nb for UHF.
!      But can also be set as other numbers. E.g. na+npair for GVB.
! nfile: the number of .fch(k) files, i.e. the number of known geometries and
!        known MOs, plus a .gjf file with unknown MOs
! S: AO overlap matrices of each known geometry and the new geometry
!    S(:,:,nfile) is the overlap matrix of the the new geometry
! mo: occupied MOs of each known geometry
! new_mo: occupied MOs of the new geometry
subroutine mo_grassmann_intrplt(nbf, nmo, nfile, x, S, mo, new_mo)
 implicit none
 integer :: i, k, itmp(1)
 integer, intent(in) :: nbf, nmo, nfile
!f2py intent(in) :: nbf, nmo, nfile
 real(kind=8), intent(in) :: x(nfile), S(nbf,nbf,nfile), mo(nbf,nmo,nfile-1)
!f2py intent(in) :: x, S, mo
!f2py depend(nfile) :: x
!f2py depend(nbf,nfile) :: S
!f2py depend(nbf,nmo,nfile) :: mo
 real(kind=8), intent(out) :: new_mo(nbf,nmo)
!f2py intent(out) :: new_mo
!f2py depend(nbf,nmo) :: new_mo
 real(kind=8), parameter :: diff = 1d-7
 real(kind=8), allocatable :: sqrt_S(:,:), n_sqrt_S(:,:), mo_ref(:,:), mo_k(:,:),&
  u(:,:), vt(:,:), sv(:), sin_s(:,:), cos_s(:,:), u_sin_s(:,:), mo_ref_v(:,:), &
  weight(:)

! allocate(weight(nfile-1))
! call calc_1d_lagrange_w(nfile, x, weight)
 allocate(weight(2:nfile-1))
 call calc_1d_lagrange_w(nfile-1, x(2:), weight)
 write(6,'(/,A)') "Weights of Lagrange's interpolation:"
 write(6,'(5(1X,ES15.8))') weight

 ! In Ref 10.1063/5.0137775 and Ref 10.1063/5.0153440, the authors found that for
 !  RHF/UHF, it does not matter which geometry is set as the reference geometry.

 ! I've checked that is true. And choosing the geometry whose weight is closest
 !  to 1.0 will usually make the result slightly better, for RHF/UHF.
 ! For ROHF, I find the result is sensitive to the choice of reference geometry.
 !  So here we choose the geometry whose weight is closest to 1.0.
 !itmp = MINLOC(DABS(weight-1d0))
 !k = itmp(1)
 k = 1
 ! even if the maximum weight is 1.0 and other weights are 0, the following step
 ! is needed

 ! Calculate S^1/2, S^(-1/2) of the 1st geometry
 allocate(sqrt_S(nbf,nbf), n_sqrt_S(nbf,nbf))
 call mat_dsqrt(nbf, S(:,:,k), .false., sqrt_S, n_sqrt_S)

 ! C_ref' = (S^1/2)C_ref. MOs at orthogonal basis of the reference geometry
 allocate(mo_ref(nbf,nmo))
 call dsymm('L', 'L', nbf, nmo, 1d0, sqrt_S,nbf, mo(:,:,k),nbf, 0d0, mo_ref,nbf)
 allocate(mo_k(nbf,nmo), source=mo_ref)
 call grassmann_C2GAMMA(nbf, nmo, mo_ref, mo_k) ! GAMMA_1 stored in mo_k
 !new_mo = weight(k)*mo_k
 new_mo = 0d0

 ! GAMMA matrix of the new geometry would be stored in the array mo_k temporarily
 do i = 2, nfile-1, 1
  !if(i == k) cycle
  if(DABS(weight(i)) < diff) cycle

  call mat_dsqrt(nbf, S(:,:,i), .false., sqrt_S, n_sqrt_S)
  call dsymm('L', 'L', nbf, nmo, 1d0, sqrt_S, nbf, mo(:,:,i),nbf, 0d0, mo_k,nbf)
  call grassmann_C2GAMMA(nbf, nmo, mo_ref, mo_k) ! GAMMA_i stored in mo_k
  new_mo = new_mo + weight(i)*mo_k               ! linear combination of GAMMA_i
 end do ! for i
 deallocate(sqrt_S, n_sqrt_S, mo_k, weight)

 ! GAMMA = U(SIGMA)V^T
 allocate(u(nbf,nbf), vt(nmo,nmo), sv(nbf))
 call do_svd(nbf, nmo, new_mo, u, vt, sv)

 ! U*sin(SIGMA)
 allocate(sin_s(nbf,nmo), source=0d0)
 forall(i = 1:nmo) sin_s(i,i) = DSIN(sv(i))
 allocate(u_sin_s(nbf,nmo))
 call dgemm('N', 'N', nbf, nmo, nbf, 1d0, u, nbf, sin_s, nbf, 0d0, u_sin_s,nbf)
 deallocate(u, sin_s)

 ! (C_ref')V
 allocate(mo_ref_v(nbf,nmo))
 call dgemm('N', 'T', nbf, nmo, nmo, 1d0, mo_ref,nbf, vt,nmo, 0d0,mo_ref_v,nbf)
 deallocate(mo_ref)

 ! (C_ref')V*cos(SIGMA) + U*sin(SIGMA)
 allocate(cos_s(nmo,nmo), source=0d0)
 forall(i = 1:nmo) cos_s(i,i) = DCOS(sv(i))
 deallocate(sv)
 call dgemm('N','N', nbf,nmo,nmo, 1d0,mo_ref_v,nbf, cos_s,nmo, 1d0,u_sin_s,nbf)
 deallocate(mo_ref_v, cos_s)

 ! C_unk' = ((C_ref')V*cos(SIGMA) + U*sin(SIGMA))(V^T), unk: unknown
 call dgemm('N','N', nbf,nmo,nmo, 1d0, u_sin_s, nbf, vt, nmo, 0d0, new_mo, nbf)
 deallocate(u_sin_s, vt)

 ! C_unk = (S^(-1/2))(C_unk')
 allocate(sqrt_S(nbf,nbf), n_sqrt_S(nbf,nbf))
 call mat_dsqrt(nbf, S(:,:,nfile), .true., sqrt_S, n_sqrt_S)
 deallocate(sqrt_S)
 allocate(mo_ref(nbf,nmo))
 call dsymm('L', 'L', nbf, nmo, 1d0,n_sqrt_S,nbf, new_mo, nbf, 0d0, mo_ref,nbf)
 deallocate(n_sqrt_S)

 new_mo = mo_ref
 deallocate(mo_ref)
end subroutine mo_grassmann_intrplt

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

 k = INDEX(gjfname1, '.gjf', back=.true.)
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

! check whether natom in a .fch file and a .gjf/.xyz file are the same
subroutine check_natom_eq_in_fch_and_gjf_or_xyz(fchname, coor_file, natom)
 implicit none
 integer :: i, natom1
 integer, intent(out) :: natom
 character(len=240), intent(in) :: fchname, coor_file

 natom = 0 

 i = INDEX(coor_file, '.xyz', back=.true.)
 if(i > 0) then
  call read_natom_from_xyz(coor_file, natom)
 else ! not .xyz, assume it be .gjf
  call read_natom_from_gjf(coor_file, natom)
 end if

 call read_natom_from_fch(fchname, natom1)
 if(natom1 /= natom) then
  write(6,'(/,A)') 'ERROR in subroutine check_natom_eq_in_fch_and_gjf_or_xyz: t&
                   &he number of atoms'
  write(6,'(A)') 'is not equal to each other in two files.'
  write(6,'(A,I0)') 'fchname = '//TRIM(fchname)//', natom1 = ', natom1
  write(6,'(A,I0)') 'coor_file = '//TRIM(coor_file)//', natom = ', natom
  stop
 end if
end subroutine check_natom_eq_in_fch_and_gjf_or_xyz

! get the atomic correspondence indices of two geometries
subroutine get_atom_map_idx_from_coor(natom, coor1, coor2, idx)
 implicit none
 integer :: i, j
 integer, intent(in) :: natom
 integer, intent(out) :: idx(natom)
 real(kind=8) :: norm, rtmp1(3), rtmp2(3)
 real(kind=8), parameter :: thres = 1d-5
 real(kind=8), intent(in) :: coor1(3,natom), coor2(3,natom)
 logical, allocatable :: used1(:), used2(:)

 idx = 0 ! initialization
 allocate(used1(natom), used2(natom))
 used1 = .false.; used2 = .false.

 do i = 1, natom, 1
  rtmp1 = coor1(:,i)

  do j = 1, natom, 1
   if(used2(j)) cycle
   rtmp2 = coor2(:,j) - rtmp1
   norm = DOT_PRODUCT(rtmp2, rtmp2)
   if(norm < thres) then
    used1(i) = .true.; used2(j) = .true.
    idx(j) = i
    exit
   end if
  end do ! for j

 end do ! for i

 if(ANY(used1 .eqv. .false.) .or. ANY(used2 .eqv. .false.)) then
  write(6,'(/,A)') 'ERROR in subroutine get_atom_map_idx_from_coor: one-to-one &
                   &correspondence'
  write(6,'(A)') 'cannot be built for these two geometries.'
  write(6,'(A,I0)') 'natom=', natom
  stop
 end if

 deallocate(used1, used2)
end subroutine get_atom_map_idx_from_coor

