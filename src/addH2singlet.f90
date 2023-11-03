! written by jxzou at 20220214: add hydrogen atoms to make sure the whole
!  molecule is singlet.
!  The hydrogen atoms are added onto the principal axis of the molecule

program main
 implicit none
 integer :: i
 character(len=4) :: str
 character(len=240) :: fchname
 logical :: gvb_no

 i = iargc()
 if(.not. (i==1 .or. i==2)) then
  write(6,'(/,A)') ' ERROR in subroutine addH2singlet: wrong command line&
                   & arguments!'
  write(6,'(A)')   ' Example 1 (HF MOs) : addH2singlet a.fch'
  write(6,'(A,/)') ' Example 2 (GVB NOs): addH2singlet a.fch -gvb'
  stop
 end if

 call getarg(1, fchname)

 gvb_no = .false.
 if(i == 2) then
  call getarg(2, str)
  if(str /= '-gvb') then
   write(6,'(/,A)') 'ERROR in subroutine addH2singlet: wrong command line&
                    & arguments!'
   write(6,'(A)') "The 2nd argument can only be '-gvb'."
   stop
  else
   gvb_no = .true.
  end if
 end if

 call require_file_exist(fchname)
 call addH2singlet(fchname, gvb_no)
end program main

! add hydrogen atoms for a specified .fch file
subroutine addH2singlet(fchname, gvb_no)
 use fch_content
 implicit none
 integer :: i, j, nadd
 integer, allocatable :: ielem1(:), shell_type1(:), prim_per_shell1(:)
 integer, allocatable :: shell2atom_map1(:)
 real(kind=8), parameter :: dis = 50d0 ! Angstrom
 real(kind=8), parameter :: root2 = 0.5d0*DSQRT(2d0)
 real(kind=8) :: w(3), it(3,3) ! it: inertia tensor
 real(kind=8), allocatable :: prim_exp1(:), contr_coeff1(:), contr_coeff_sp1(:)
 real(kind=8), allocatable :: eigen_e_a1(:), coor1(:,:), coorH(:,:)
 real(kind=8), allocatable :: alpha_coeff1(:,:), beta_coeff1(:,:), coeff(:,:)
 character(len=240) :: new_fch
 character(len=240), intent(in) :: fchname
 logical :: uhf
 logical, intent(in) :: gvb_no

 call check_uhf_in_fch(fchname, uhf)
 call read_fch(fchname, uhf)
 if(uhf .and. gvb_no) then
  write(6,'(/,A)') 'ERROR in subroutine addH2singlet: both uhf and gvb_no are .True.'
  write(6,'(A)') "This is a UHF-type .fch file, you cannot specify '-gvb'."
  stop
 end if

 if(mult == 1) then
  write(6,'(/,A)') 'ERROR in subroutine addH2singlet: this is a singlet molecule.'
  write(6,'(A)') 'There is no need to add hydrogen atoms.'
  stop
 end if

 if(mult > 6) then
  write(6,'(/,A)') 'ERROR in subroutine addH2singlet: spin is too high.'
  write(6,'(A)') 'This subroutine does not support >=7 singly occupied &
                 &electrons.'
  stop
 end if

 nadd = na - nb
 allocate(ielem1(natom+nadd))
 ielem1(1:natom) = ielem
 ielem1(natom+1:natom+nadd) = 1
 ielem = ielem1
 deallocate(ielem1)

 allocate(shell_type1(ncontr+nadd))
 shell_type1(1:ncontr) = shell_type
 shell_type1(ncontr+1:ncontr+nadd) = 0
 shell_type = shell_type1
 deallocate(shell_type1)

 allocate(prim_per_shell1(ncontr+nadd))
 prim_per_shell1(1:ncontr) = prim_per_shell
 prim_per_shell1(ncontr+1:ncontr+nadd) = 2
 prim_per_shell = prim_per_shell1
 deallocate(prim_per_shell1)

 allocate(shell2atom_map1(ncontr+nadd))
 shell2atom_map1(1:ncontr) = shell2atom_map
 forall(i=1:nadd) shell2atom_map1(ncontr+i) = natom + i
 shell2atom_map = shell2atom_map1
 deallocate(shell2atom_map1)

 allocate(prim_exp1(nprim+2*nadd))
 prim_exp1(1:nprim) = prim_exp
 forall(i = 1:nadd) ! STO-2G basis set
  prim_exp1(nprim+2*i-1) = 1.309756377d0
  prim_exp1(nprim+2*i) = 0.2331359749d0
 end forall
 prim_exp = prim_exp1
 deallocate(prim_exp1)

 allocate(contr_coeff1(nprim+2*nadd))
 contr_coeff1(1:nprim) = contr_coeff
 forall(i = 1:nadd) ! STO-2G basis set
  contr_coeff1(nprim+2*i-1) = 0.4301284983d0
  contr_coeff1(nprim+2*i) = 0.6789135305d0
 end forall
 contr_coeff = contr_coeff1
 deallocate(contr_coeff1)

 if(allocated(contr_coeff_sp)) then
  allocate(contr_coeff_sp1(nprim+2*nadd))
  contr_coeff_sp1(1:nprim) = contr_coeff_sp
  contr_coeff_sp1(nprim+1:nprim+2*nadd) = 0d0
  contr_coeff_sp = contr_coeff_sp1
  deallocate(contr_coeff_sp1)
 end if

 allocate(eigen_e_a1(nif+nadd), source=0d0)
 if(gvb_no) then
  eigen_e_a1(1:na) = eigen_e_a(1:na)
  eigen_e_a1(na+1:na+nadd) = 1d0
  eigen_e_a1(na+nadd+1:) = eigen_e_a(na+1:)
  eigen_e_a = eigen_e_a1
 else
  eigen_e_a = eigen_e_a1
  if(uhf) eigen_e_b = eigen_e_a1
 end if
 deallocate(eigen_e_a1)

 ! compute the inertia tensor
 call get_inertia_tensor(natom, coor, it)
 call diag_get_e_and_vec(3, it, w)

 ! compute the geometry center
 forall(i = 1:3) w(i) = SUM(coor(i,:))/DBLE(natom)

 ! calculate the 6 positions of added hydrogen atoms
 allocate(coorH(3,6), source=0d0)
 coorH(:,1) = it(:,3)*dis
 coorH(:,3) = it(:,2)*dis
 coorH(:,5) = it(:,1)*dis
 forall(i = 1:3) coorH(:,2*i) = -coorH(:,2*i-1)
 forall(i = 1:6) coorH(:,i) = coorH(:,i) + w

 ! add hydrogen atoms onto the principal axis of the molecule
 allocate(coor1(3,natom+nadd))
 coor1(:,1:natom) = coor
 coor1(:,natom+1:natom+nadd) = coorH(:,1:nadd)
 coor = coor1
 deallocate(coor1, coorH)

 allocate(alpha_coeff1(nbf+nadd,nif+nadd), source=0d0)
 ! copy alpha occupied MOs
 alpha_coeff1(1:nbf,1:na) = alpha_coeff(:,1:na)
 ! add hydrogen orbitals into alpha MOs
 forall(i = 1:nadd) alpha_coeff1(nbf+i,na+i) = 1d0
 ! copy alpha virtual MOs
 alpha_coeff1(1:nbf,na+1+nadd:nif+nadd) = alpha_coeff(:,na+1:nif)

 if(uhf) then
  ! add hydrogen orbitals into beta MOs
  allocate(beta_coeff1(nbf+nadd,nif+nadd), source=0d0)
  beta_coeff1(1:nbf,1:nb) = beta_coeff(:,1:nb)
  forall(i = 1:nadd) beta_coeff1(nbf+i,nb+i) = 1d0
  beta_coeff1(1:nbf,nb+1+nadd:nif+nadd) = beta_coeff(:,nb+1:nif)
 else
  if(gvb_no) then   ! convert into GVB NOs for each pair
   allocate(coeff(nbf+nadd,2))
   do i = 1, nadd, 1
    j = na + nadd + 1 - i
    coeff(:,1) = alpha_coeff1(:,nb+i)
    coeff(:,2) = alpha_coeff1(:,j)
    alpha_coeff1(:,nb+i) = root2*(coeff(:,1) + coeff(:,2))
    alpha_coeff1(:,j)    = root2*(coeff(:,1) - coeff(:,2))
   end do ! for i
   deallocate(coeff)
  end if
 end if

 ! copy back to alpha_coeff, auto-reallocation of alpha_coeff
 alpha_coeff = alpha_coeff1
 deallocate(alpha_coeff1)
 if(uhf) then
  beta_coeff = beta_coeff1
  deallocate(beta_coeff1)
 end if

 ! update dimension parameters
 nbf = nbf + nadd
 nif = nif + nadd
 nb = na; nopen = 0
 ncontr = ncontr + nadd
 nprim = nprim + 2*nadd
 mult = 1
 natom = natom + nadd

 if(allocated(tot_dm)) deallocate(tot_dm)
 allocate(tot_dm(nbf,nbf), source=0d0)
 ! if this is UHF/UDFT wave funtion, update tot_dm
 if(uhf) then
  if(allocated(spin_dm)) deallocate(spin_dm)
  allocate(spin_dm(nbf,nbf), source=0d0)
  allocate(coeff(nbf,nbf)) ! use array coeff to store beta density matrix
  tot_dm = MATMUL(alpha_coeff(:,1:na), TRANSPOSE(alpha_coeff(:,1:na)))
  coeff = MATMUL(beta_coeff(:,1:nb), TRANSPOSE(beta_coeff(:,1:nb)))
  spin_dm = tot_dm - coeff
  tot_dm = tot_dm + coeff
  deallocate(coeff)
 end if

 if(allocated(mull_char)) deallocate(mull_char)
 allocate(mull_char(natom), source=0d0)

 ! create/generate a new .fch file
 i = INDEX(fchname, '.fch')
 new_fch = fchname(1:i-1)//'aH.fch'
 call write_fch(new_fch)
end subroutine addH2singlet

! compute the inertia tensor, mass not considered
subroutine get_inertia_tensor(natom, coor, it)
 implicit none
 integer :: i
 integer, intent(in) :: natom
 real(kind=8) :: rtmp(3)
 real(kind=8), intent(in) :: coor(3,natom)
 real(kind=8), intent(out) :: it(3,3) ! inertia tensor

 forall(i = 1:3) rtmp(i) = DOT_PRODUCT(coor(i,:), coor(i,:))
 it(1,1) = rtmp(2) + rtmp(3)
 it(2,2) = rtmp(3) + rtmp(1)
 it(3,3) = rtmp(1) + rtmp(2)
 it(2,1) = -DOT_PRODUCT(coor(1,:), coor(2,:))
 it(3,1) = -DOT_PRODUCT(coor(1,:), coor(3,:))
 it(3,2) = -DOT_PRODUCT(coor(2,:), coor(3,:))
 it(1,2) = it(2,1)
 it(1,3) = it(3,1)
 it(2,3) = it(3,2)
end subroutine get_inertia_tensor

