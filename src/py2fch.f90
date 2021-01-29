! written by jxzou at 20171218: adjust the orders of d, f, g, h functions etc.
!  in the input file of PySCF, to the oeders of .fch(k) file in Gaussian

! updated by jxzou at 20180317: add a input parameter 'a' or 'b' to write the MO into the Alpha or Beta part in .fchk file
! updated by jxzou at 20180404: support the case that D functions preceding L functions
! updated by jxzou at 20180406: code optimization
! updated by jxzou at 20180520: support linear dependence
! updated by jxzou at 20190226: change data to parameter when declare constant arrays
! updated by jxzou at 20190411: pass the Sdiag in
! updated by jxzou at 20200326: renamed as py2fch
! updated by jxzou at 20200425: add input arrays ev (for eigenvalues or NOONs)
! updated by jxzou at 20210126: generate Total SCF Density using MOs and ONs

! This subroutine is designed to be imported as a module by Python.
!  For INTEL compiler, use
! ---------------------------------------------------------------------
!  f2py -m py2fch -c py2fch.f90 --fcompiler=intelem --compiler=intelem
! ---------------------------------------------------------------------
!  For GNU compiler, use
! ------------------------------
!  f2py -m py2fch -c py2fch.f90
! ------------------------------

!  to compile this file (a py2fch.so file will be generated). Then in
!  Python you can import the py2fch module.

! read the MOs in .fch(k) file and adjust its d,f,g etc. functions order
!  of PySCF to that of Gaussian
subroutine py2fch(fchname, nbf, nif, coeff2, Sdiag, ab, ev, gen_density)
 implicit none
 integer :: i, j, k, fid, fid1
 integer :: ncoeff, nbf, nif, nocc
!f2py intent(in) :: nbf, nif
 integer :: n6dmark,n10fmark,n15gmark,n21hmark
 integer :: n5dmark,n7fmark, n9gmark, n11hmark
 integer, allocatable :: shell_type(:), shell_to_atom_map(:)
 integer, parameter :: iout = 6
 ! mark the index where d, f, g, h functions begin
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 integer :: RENAME

 real(kind=8), parameter :: diff = 1.0d-6
 real(kind=8) :: coeff2(nbf,nif), Sdiag(nbf), ev(nif)
!f2py intent(in,copy) :: coeff2, Sdiag
!f2py intent(in) :: ev
!f2py depend(nbf,nif) :: coeff2
!f2py depend(nbf) :: Sdiag
!f2py depend(nif) :: ev
 real(kind=8), allocatable :: coeff(:), den(:,:)

 character(len=1) :: ab
!f2py intent(in) :: ab
 character(len=8) :: key0, key
 character(len=8), parameter :: key1 = 'Alpha MO'
 character(len=7), parameter :: key2 = 'Beta MO'
 character(len=8), parameter :: key3 = 'Alpha Or'
 character(len=7), parameter :: key4 = 'Beta Or'
 character(len=49) :: str
 character(len=240) :: fchname, buf, new_fchk
!f2py intent(in) :: fchname

 logical :: alive,  gen_density
!f2py intent(in) :: gen_density

! If gen_density = .True., generate total density using input MOs and eigenvalues
! In this case, the array ev should contain occupation numbers
 logical, allocatable :: eq1(:)

 inquire(file=TRIM(fchname),exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine py2fch: file does not exist!'
  write(iout,'(A)') 'Filename='//TRIM(fchname)
  stop
 end if

 ! normalize MO coefficients as Gaussian
 Sdiag = DSQRT(Sdiag)
 allocate(eq1(nbf), source=.false.)
 do i = 1, nbf, 1
  if( DABS(Sdiag(i)-1d0) < diff) eq1(i) = .true.
 end do

 do i = 1, nif, 1
  do k = 1, nbf, 1
   if(.not. eq1(k)) coeff2(k,i) = coeff2(k,i)*Sdiag(k)
  end do ! for k
 end do ! for i
 deallocate(eq1)
 ! normalize done

 buf = ' '
 new_fchk = ' '
 ncoeff = 0

 key0 = key3
 key = key1
 if(ab/='a' .and. ab/='A') then
  key = key2//' '
  key0 = key4//'b'
 end if

 i = index(fchname,'.fch')
 new_fchk = fchname(1:i-1)//'.p2f'
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 ! find and read Shell types
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(buf(1:11) == 'Shell types') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine py2fch: missing 'Shell types' section&
                   & in file "//TRIM(fchname)
  close(fid)
  return
 end if

 read(buf,'(A49,2X,I10)') str, k
 allocate(shell_type(2*k), source=0)
 read(fid,'(6(6X,I6))') (shell_type(i),i=1,k)
! read Shell types done

 ! find and read Shell to atom map
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(buf(1:13) == 'Shell to atom') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine py2fch: missing 'Shell to atom map'&
                   & section in file "//TRIM(fchname)
  close(fid)
  return
 end if

 allocate(shell_to_atom_map(2*k), source=0)
 read(fid,'(6(6X,I6))') (shell_to_atom_map(i),i=1,k)
 ! read Shell to atom map done

! first we adjust the basis functions in each MO according to the Shell to atom map
! this is to ensure that the order of basis functions in PySCF is converted to be that of Gaussian
 ! unsort the shell_type, shell_to_atom_map, MOs will be adjusted accordingly
 call unsort_shell_and_mo(k, shell_type, shell_to_atom_map, nbf, nif, coeff2)
 ! Note that k will be updated
 deallocate(shell_to_atom_map)
! adjust done

! then we adjust the basis functions in each MO according to the type of basis functions
 n5dmark  = 0 ; n6dmark  = 0
 n7fmark  = 0 ; n10fmark = 0
 n9gmark  = 0 ; n15gmark = 0
 n11hmark = 0 ; n21hmark = 0
 allocate(d_mark(k), source=0)
 allocate(f_mark(k), source=0)
 allocate(g_mark(k), source=0)
 allocate(h_mark(k), source=0)
 nbf = 0

 do i = 1, k, 1
  select case(shell_type(i))
  case( 0)   !'S'
   nbf = nbf + 1
  case( 1)   !'P'
   nbf = nbf + 3
  case(-1)   !'SP' or 'L'
   nbf = nbf + 4
  case(-2)   ! 5D
   n5dmark = n5dmark + 1
   d_mark(n5dmark) = nbf + 1
   nbf = nbf + 5
  case( 2)   ! 6D
   n6dmark = n6dmark + 1
   d_mark(n6dmark) = nbf + 1
   nbf = nbf + 6
  case(-3)   ! 7F
   n7fmark = n7fmark + 1
   f_mark(n7fmark) = nbf + 1
   nbf = nbf + 7
  case( 3)   ! 10F
   n10fmark = n10fmark + 1
   f_mark(n10fmark) = nbf + 1
   nbf = nbf + 10
  case(-4)   ! 9G
   n9gmark = n9gmark + 1
   g_mark(n9gmark) = nbf + 1
   nbf = nbf + 9
  case( 4)   ! 15G
   n15gmark = n15gmark + 1
   g_mark(n15gmark) = nbf + 1
   nbf = nbf + 15
  case(-5)   ! 11H
   n11hmark = n11hmark + 1
   h_mark(n11hmark) = nbf + 1
   nbf = nbf + 11
  case( 5)   ! 21H
   n21hmark = n21hmark + 1
   h_mark(n21hmark) = nbf + 1
   nbf = nbf + 21
  end select
 end do ! for i
 deallocate(shell_type)

 ! adjust the order of d, f, g etc. functions
 do i = 1,n5dmark,1
  call py2fch_permute_5d(nif,coeff2(d_mark(i):d_mark(i)+4,:))
 end do
 do i = 1,n6dmark,1
  call py2fch_permute_6d(nif,coeff2(d_mark(i):d_mark(i)+5,:))
 end do
 do i = 1,n7fmark,1
  call py2fch_permute_7f(nif,coeff2(f_mark(i):f_mark(i)+6,:))
 end do
 do i = 1,n10fmark,1
  call py2fch_permute_10f(nif,coeff2(f_mark(i):f_mark(i)+9,:))
 end do
 do i = 1,n9gmark,1
  call py2fch_permute_9g(nif,coeff2(g_mark(i):g_mark(i)+8,:))
 end do
 do i = 1,n15gmark,1
  call py2fch_permute_15g(nif,coeff2(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1,n11hmark,1
  call py2fch_permute_11h(nif,coeff2(h_mark(i):h_mark(i)+10,:))
 end do
 do i = 1,n21hmark,1
  call py2fch_permute_21h(nif,coeff2(h_mark(i):h_mark(i)+20,:))
 end do
 deallocate(d_mark, f_mark, g_mark, h_mark)
 ! adjustment finished

 ! write the MOs into the .fchk file
 rewind(fid)
 open(newunit=fid1,file=TRIM(new_fchk),status='replace')
 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(buf(1:8) == key0) exit ! Alpha/Beta Orbital Energies
 end do
 write(fid1,'(5(1X,ES15.8))') (ev(i),i=1,nif)

 do while(.true.) ! skip the Orbital Energies in old fch(k)
  read(fid,'(A)') buf
  if(buf(49:49) == '=') exit
 end do
 write(fid1,'(A)') TRIM(buf)

 if(buf(1:8) /= key) then
  do while(.true.)
   read(fid,'(A)') buf
   write(fid1,'(A)') TRIM(buf)
   if(buf(1:8) == key) exit  ! Alpha/Beta MO
  end do ! for while
 end if
 read(buf,'(A49,2X,I10)') str, ncoeff

 if(ncoeff /= nbf*nif) then
  write(iout,'(A)') 'ERROR in subroutine py2fch: inconsistent basis set in&
                   & .fch(k) file and in PySCF script. ncoeff/=nbf*nif!'
  write(iout,'(A)') 'fchname='//TRIM(fchname)
  write(iout,'(3(A,I0))') 'ncoeff=', ncoeff, ', nif=', nif, ', nbf=', nbf
  close(fid)
  close(fid1,status='delete')
  return
 end if

 allocate(coeff(ncoeff))
 coeff = RESHAPE(coeff2,(/ncoeff/))
 write(fid1,'(5(1X,ES15.8))') (coeff(i),i=1,ncoeff)
 deallocate(coeff)

 ! copy the rest of the .fchk file
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(49:49) == '=') exit
 end do
 ! here should be 'Orthonormal basis' or 'Total SCF Density'
 BACKSPACE(fid)

 if(gen_density) then
  if( ANY(ev<-0.1D0) ) then
   write(iout,'(A)') 'ERROR in subroutine py2fch: occupation numbers of some&
                    & orbitals < -0.1. Did you'
   write(iout,'(A)') 'mistake orbital energies for occupation numbers?'
   close(fid1)
   close(fid,status='delete')
   stop
  end if

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   write(fid1,'(A)') TRIM(buf)
   if(buf(1:11) == 'Total SCF D') exit
  end do ! for while

  if(i /= 0) then
   write(iout,'(A)') "ERROR in subroutine py2fch: no 'Total SCF D' found in&
                    & file "//TRIM(fchname)
   close(fid1)
   close(fid,status='delete')
   stop
  end if

  allocate(den(nbf,nbf), source=0d0)
  ! only den(j,i) j<=i will be assigned values
  do i = 1, nbf, 1
   do j = 1, i, 1
    do k = 1, nif, 1
     if(DABS(ev(k)) < 1D-7) cycle
     den(j,i) = den(j,i) + ev(k)*coeff2(j,k)*coeff2(i,k)
    end do ! for k
   end do ! for j
  end do ! for i

  write(fid1,'(5(1X,ES15.8))') ((den(k,i),k=1,i),i=1,nbf)
  deallocate(den)

  do while(.true.) ! skip density in the original .fch(k) file
   read(fid,'(A)') buf
   if(buf(49:49) == '=') then
    write(fid1,'(A)') TRIM(buf)
    exit
   end if
  end do ! for while
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while
 ! writing new fchk file done

 close(fid1)
 close(fid, status='delete')
 i = RENAME(TRIM(new_fchk),TRIM(fchname))
 return
end subroutine py2fch

! split the 'L' into 'S' and 'P'
subroutine split_L_func(k, shell_type, shell_to_atom_map, length)
 implicit none
 integer :: i, k0
 integer,intent(in) :: k
 integer,intent(inout) :: shell_type(2*k), shell_to_atom_map(2*k)
 integer,intent(out) :: length
 integer,allocatable :: temp1(:), temp2(:)

 k0 = 2*k
 length = k
 ! set initial values for arrays shell_type, assume 15 will not be used
 shell_type(k+1:k0) = 15
 i = 1
 do while(shell_type(i) /= 15)
  if(shell_type(i) /= -1) then
   i = i + 1
   cycle
  end if
  shell_type(i) = 0
  allocate( temp1(i+1 : k0-1), temp2(i+1 : k0-1) )
  temp1(i+1 : k0-1) = shell_type(i+1 : k0-1)
  shell_type(i+2 : k0) = temp1(i+1 : k0-1)
  temp2(i+1 : k0-1) = shell_to_atom_map(i+1 : k0-1)
  shell_to_atom_map(i+2 : k0) = temp2(i+1 : k0-1)
  deallocate(temp1, temp2)
  shell_type(i+1) = 1
  shell_to_atom_map(i+1) = shell_to_atom_map(i)
  i = i + 2
 end do

 length = i - 1
 shell_type(i : k0) = 0
 return
end subroutine split_L_func

! unsort the shell_type, shell_to_atom_map according to the order in Gaussian
! MOs will be adjusted accordingly
subroutine unsort_shell_and_mo(ilen, shell_type, shell_to_atom_map, nbf, nif, coeff2)
 implicit none
 integer :: i, j, k
 integer :: length, natom
 integer :: ibegin, iend, jbegin, jend
 integer, parameter :: ntype = 10
 integer, parameter :: num0(ntype) = (/0, 1, -2, 2, -3, 3, -4, 4, -5, 5/)
 integer, parameter :: num1(ntype) = (/1, 3, 5, 6, 7, 10, 9, 15, 11, 21/)
 !                                     S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer :: num(ntype)
 integer,intent(inout) :: ilen
 integer,intent(in) :: nbf, nif
 ! Note that this 'ilen' is not identical to that in fchk2py.f90.
 !  The -1 in array shell_type has not been split. So below we use 2*ilen.

 integer,intent(inout) :: shell_type(2*ilen), shell_to_atom_map(2*ilen)
 integer new_shell_type(2*ilen)
 integer new_shell_to_atom_map(2*ilen)
 integer,allocatable :: ith(:), new_ith(:), ith_bas(:), tmp_type(:)
 real(kind=8),intent(inout) :: coeff2(nbf,nif)

! split the 'L' into 'S' and 'P'
 new_shell_to_atom_map = shell_to_atom_map
 call split_L_func(ilen, shell_type, new_shell_to_atom_map, length)
 new_shell_type = shell_type

 ! find the number of atoms
 natom = shell_to_atom_map(ilen)

 allocate(ith(0:natom), new_ith(0:natom), ith_bas(0:natom))
 ith = 0
 new_ith = 0
 ith_bas = 0

 ! find the end position of each atom in array shell_to_atom_map
 do i = 1, natom, 1
  ith(i) = count(shell_to_atom_map==i) + ith(i-1)
  new_ith(i) = count(new_shell_to_atom_map==i) + new_ith(i-1)
 end do

 ! find the end position of basis functions of each atom
 do i = 1, natom, 1
  ibegin = new_ith(i-1) + 1
  iend = new_ith(i)
  allocate(tmp_type(ibegin:iend))
  tmp_type = 0
  tmp_type = new_shell_type(ibegin:iend)
  num = 0
  do j = 1, ntype, 1
   num(j) = count(tmp_type == num0(j))
  end do
  k = 0
  do j = 1, ntype, 1
   k = k + num(j)*num1(j)
  end do
  ith_bas(i) = ith_bas(i-1) + k
  deallocate(tmp_type)
 end do

 ! adjust the MOs in each atom
 do i = 1, natom, 1
  ibegin = new_ith(i-1) + 1
  iend = new_ith(i)
  jbegin = ith_bas(i-1) + 1
  jend = ith_bas(i)
  call sort_shell_type_in_each_atom2(iend-ibegin+1, new_shell_type(ibegin:iend))
  call unsort_mo_in_each_atom(iend-ibegin+1, shell_type(ibegin:iend), new_shell_type(ibegin:iend), &
       & jend-jbegin+1, nif, coeff2(jbegin:jend,1:nif))
 end do
 ! adjust the MOs in each atom done
 deallocate(ith, new_ith, ith_bas)

 ! update ilen
 ilen = length
 return
end subroutine unsort_shell_and_mo

! sort the shell_type within each atom
subroutine sort_shell_type_in_each_atom2(ilen, shell_type)
 implicit none
 integer :: i, tmp_type
 integer, intent(in) :: ilen
 integer, intent(inout) :: shell_type(ilen)
 logical :: sort_done

 sort_done = .false.
 do while(.not. sort_done)
  sort_done = .true.
  do i = 1, ilen-1, 1
    if(shell_type(i) == 0) cycle
    if(ABS(shell_type(i+1)) >= ABS(shell_type(i))) cycle
    sort_done = .false.
    tmp_type = shell_type(i+1)
    shell_type(i+1) = shell_type(i)
    shell_type(i) = tmp_type
  end do
 end do
 return
end subroutine sort_shell_type_in_each_atom2

! Unsort the MO within each atom. (Only for 'L' in Pople basis)
! Example:
!  PySCF:    1s, 2s, 3s, 2px, 2py, 2pz, 3px, 3py, 3pz
!  Gaussian: 1s, 2s, 2px, 2py, 2pz, 3s, 3px, 3py, 3pz
subroutine unsort_mo_in_each_atom(ilen1, shell_type, new_shell_type, ilen2, nif, coeff2)
 implicit none
 integer :: i, j, k, m   ! temporary variables
 integer :: ibegin
 integer, parameter :: ntype = 10
 integer :: ith_shell(ntype)
 integer, parameter :: num0(ntype) = (/0, 1, -2, 2, -3, 3, -4, 4, -5, 5/)
 integer, parameter :: num1(ntype) = (/1, 3, 5, 6, 7, 10, 9, 15, 11, 21/)
 !                                    S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer, parameter :: rnum(-5:5) = (/9, 7, 5, 3, 0, 1, 2, 4, 6, 8, 10/)

 integer, intent(in) :: ilen1, ilen2, nif
 integer, intent(in) :: shell_type(ilen1), new_shell_type(ilen1)
 integer, allocatable :: ith_bas(:)
 real(kind=8), intent(inout) :: coeff2(ilen2,nif)
 real(kind=8) :: new_coeff2(ilen2,nif)

 new_coeff2 = 0.0d0
 ! get the begin position of each type of basis functions within an atom
 ith_shell = 0
 k = rnum(new_shell_type(1))
 ith_shell(k) = 1
 do i = k+1, ntype, 1
  call get_1st_loc(num0(i), ith_shell(i), ilen1, new_shell_type)
 end do

 ! find the begin position of basis functions within an atom
 allocate(ith_bas(ilen1))
 ith_bas = 0
 ith_bas(1) = 1
 do i = 2, ilen1, 1
  ith_bas(i) = ith_bas(i-1) + num1(rnum(new_shell_type(i-1)))
 end do

 ! unsort the MO within an atom
 j = 0
 do i = 1, ilen1, 1
  k = rnum(shell_type(i))
  ibegin = ith_bas(ith_shell(k))
  m = num1(k) - 1
  j = j + 1
  new_coeff2(j:j+m,:) = coeff2(ibegin:ibegin+m,:)
  ith_shell(k) = ith_shell(k) + 1
  j = j + m
 end do
 deallocate(ith_bas)
 ! unsort done

 ! update array coeff2
 coeff2 = new_coeff2
 return
end subroutine unsort_mo_in_each_atom

subroutine get_1st_loc(inum, loc, ilen, a)
 implicit none
 integer :: i
 integer, intent(out) :: loc
 integer, intent(in) :: inum, ilen
 integer, intent(in) :: a(ilen)

 do i = 1, ilen, 1
  if(a(i) == inum) exit
 end do
 loc = i
 if(i == ilen+1) loc = 0
 return
end subroutine get_1st_loc

subroutine py2fch_permute_5d(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(5) = (/3, 4, 2, 5, 1/)
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(5,nif)
 real(kind=8) :: coeff2(5,nif)
! From: the order of spherical d functions in PySCF
! To: the order of spherical d functions in Gaussian
! 1    2    3    4    5
! d-2, d-1, d0 , d+1, d+2
! d0 , d+1, d-1, d+2, d-2

 coeff2 = 0d0
 forall(i = 1:5)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine py2fch_permute_5d

subroutine py2fch_permute_6d(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(6) = (/1, 4, 6, 2, 3, 5/)
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(6,nif)
 real(kind=8) :: coeff2(6,nif)
! From: the order of Cartesian d functions in PySCF
! To: the order of Cartesian d functions in Gaussian
! 1  2  3  4  5  6
! XX,XY,XZ,YY,YZ,ZZ
! XX,YY,ZZ,XY,XZ,YZ

 coeff2 = 0d0
 forall(i = 1:6)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine py2fch_permute_6d

subroutine py2fch_permute_7f(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(7) = (/4, 5, 3, 6, 2, 7, 1/)
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(7,nif)
 real(kind=8) :: coeff2(7,nif)
! From: the order of spherical f functions in PySCF
! To: the order of spherical f functions in Gaussian
! 1    2    3    4    5    6    7
! f-3, f-2, f-1, f0 , f+1, f+2, f+3
! f0 , f+1, f-1, f+2, f-2, f+3, f-3

 coeff2 = 0d0
 forall(i = 1:7)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine py2fch_permute_7f

subroutine py2fch_permute_10f(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(10) = (/1, 7, 10, 4, 2, 3, 6, 9, 8, 5/)
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(10,nif)
 real(kind=8) :: coeff2(10,nif)
! From: the order of Cartesian f functions in PySCF
! To: the order of Cartesian f functions in Gaussian
! 1   2   3   4   5   6   7   8   9   10
! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ

 coeff2 = 0d0
 forall(i = 1:10)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine py2fch_permute_10f

subroutine py2fch_permute_9g(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(9) = (/5, 6, 4, 7, 3, 8, 2, 9, 1/)
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(9,nif)
 real(kind=8) :: coeff2(9,nif)
! From: the order of spherical g functions in PySCF
! To: the order of spherical g functions in Gaussian
! 1    2    3    4    5    6    7    8    9
! g-4, g-3, g-2, g-1, g0 , g+1, g+2, g+3, g+4
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4

 coeff2 = 0d0
 forall(i = 1:9)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine py2fch_permute_9g

subroutine py2fch_permute_15g(nif,coeff)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(15,nif)
 real(kind=8) :: coeff2(15,nif)
! From: the order of Cartesian g functions in PySCF
! To: the order of Cartesian g functions in Gaussian
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX

 coeff2 = 0d0
 forall(i = 1:15)
  coeff2(i,:) = coeff(16-i,:)
 end forall
 coeff = coeff2
 return
end subroutine py2fch_permute_15g

subroutine py2fch_permute_11h(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(11) = (/6, 7, 5, 8, 4, 9, 3, 10, 2, 11, 1/)
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(11,nif)
 real(kind=8) :: coeff2(11,nif)
! From: the order of Cartesian h functions in PySCF
! To: the order of Cartesian h functions in Gaussian
! 1    2    3    4    5    6    7    8    9    10   11
! h-5, h-4, h-3, h-2, h-1, h0 , h+1, h+2, h+3, h+4, h+5
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5

 coeff2 = 0d0
 forall(i = 1:11)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine py2fch_permute_11h

subroutine py2fch_permute_21h(nif,coeff)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(21,nif)
 real(kind=8) :: coeff2(21,nif)
! From: the order of Cartesian h functions in PySCF
! To: the order of Cartesian h functions in Gaussian
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! xxxxx,xxxxy,xxxxz,xxxyy,xxxyz,xxxzz,xxyyy,xxyyz,xxyzz,xxzzz,xyyyy,xyyyz,xyyzz,xyzzz,xzzzz,yyyyy,yyyyz,yyyzz,yyzzz,yzzzz,zzzzz
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX

 coeff2 = 0d0
 forall(i = 1:21)
  coeff2(i,:) = coeff(22-i,:)
 end forall
 coeff = coeff2
 return
end subroutine py2fch_permute_21h

