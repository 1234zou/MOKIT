! written by jxzou at 20220302: generate/calculate 1-e integrals

! solver AO-based overlap matrix (S) from condition (C^T)SC=I
! Note: this subroutine only applies to nbf=nif, i.e. no linear dependence
subroutine solve_ovlp_from_ctsc(nbf, C, S)
 implicit none
 integer :: i
 integer, intent(in) :: nbf
 real(kind=8), intent(in) :: C(nbf,nbf)
 real(kind=8), intent(out) :: S(nbf,nbf)
 real(kind=8), allocatable :: SC(:,:)

 forall(i = 1:nbf) S(i,i) = 1d0 ! unit matrix I
 allocate(SC(nbf,nbf), source=0d0)
 call solve_multi_lin_eqs(nbf, nbf, TRANSPOSE(C), nbf, S, SC)
 ! SC = X -> (C^T)S = X^T
 call solve_multi_lin_eqs(nbf, nbf, TRANSPOSE(C), nbf, TRANSPOSE(SC), S)
end subroutine solve_ovlp_from_ctsc

! solver AO-based Fock matrix (F) from condition (C^T)FC=E
subroutine solve_fock_from_ctfc(nbf, nif, C, E, F)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: C(nbf,nif), E(nif)
 real(kind=8), intent(out) :: F(nbf,nbf)
 real(kind=8), allocatable :: FC(:,:), E1(:,:)

 allocate(E1(nif,nif), source=0d0)
 forall(i = 1:nif) E1(i,i) = E(i) ! diagonal matrix
 allocate(FC(nbf,nif), source=0d0)
 call solve_multi_lin_eqs(nif, nbf, TRANSPOSE(C), nif, E1, FC)
 deallocate(E1)
 ! FC = X -> (C^T)(F^T) = X^T, (F^T) = F
 call solve_multi_lin_eqs(nif, nbf, TRANSPOSE(C), nif, TRANSPOSE(FC), F)
end subroutine solve_fock_from_ctfc

! explicitly calculate the AO-based overlap matrix (S) using basis set data
! in a given .fch(k) file
!TODO: implement this subroutine
subroutine calc_ovlp_using_fch(fchname, nbf1, S)
 use fch_content, only: nbf, ncontr, natom, shltyp2nbf, shell_type, &
  shell2atom_map, coor, read_fch
 implicit none
 integer :: i, j, k, m, i1, i2, i3, i4, Cnp2
 integer, intent(in) :: nbf1
 character(len=240), intent(in) :: fchname
 real(kind=8) :: dis0, coor1(3), coor2(3)
 real(kind=8), intent(out) :: S(nbf,nbf)
 real(kind=8), allocatable :: dis(:)

 S = 0d0
 call read_fch(fchname, .false.)
 if(nbf1 /= nbf) then
  write(6,'(A)') 'ERROR in subroutine calc_ovlp_using_fch: nbf1/=nbf.'
  write(6,'(A)') 'Input number of basis functions is not equal to that in &
                 &file '//TRIM(fchname)
  stop
 end if

 Cnp2 = natom*(natom-1)/2
 allocate(dis(Cnp2))
 ! calculate the square of distances between any two atoms
 do i = 1, natom, 1
  k = (2*natom-i)*(i-1)/2 - i
  coor1 = coor(:,i)
  do j = i+1, natom, 1
   coor2 = coor(:,j) - coor1
   dis(k+j) = DOT_PRODUCT(coor2,coor2)
  end do ! for j
 end do ! for i

 do i = 1, ncontr, 1   ! loop for each shell
  i1 = shell2atom_map(i); i2 = shltyp2nbf(shell_type(i))
  do j = i, ncontr, 1
   i3 = shell2atom_map(j); i4 = shltyp2nbf(shell_type(j))
   dis0 = 0d0
   if(i1 < i3) dis0 = dis((2*natom-i1)*(i1-1)/2 + i3-i1)
   do k = 1, i2, 1
    do m = 1, i4, 1

    end do ! for m
   end do ! for k
  end do ! for j
 end do ! for i

 call symmetrize_dmat(nbf, S)
 deallocate(dis)
end subroutine calc_ovlp_using_fch

! explicitly calculate the AO-based dipole moment matrix (D) using basis set
! data in a given .fch(k) file
!TODO: implement this subroutine
subroutine calc_dipole_using_fch(fchname, nbf, D)
 implicit none
 integer :: i
 integer, intent(in) :: nbf
 character(len=240), intent(in) :: fchname
 real(kind=8), intent(out) :: D(nbf,nbf,3) ! x,y,z 3 components

 D = 0d0
 forall(i = 1:3) call symmetrize_dmat(nbf, D(:,:,i))
end subroutine calc_dipole_using_fch

! symmetrize a double precision matrix
pure subroutine symmetrize_dmat(n, a)
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 real(kind=8), intent(inout) :: a(n)

 forall(i=1:n,j=1:n,j>i) a(i,j) = a(j,i)
end subroutine symmetrize_dmat

