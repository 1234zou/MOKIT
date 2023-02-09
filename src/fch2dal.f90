! written by jxzou at 20210613: adjust the orders of d,f,g, etc functions in
!  .fch(k) file, to the orders in Dalton

program main
 use util_wrapper, only: formchk, fch2inp_wrap
 implicit none
 integer :: i, system
 character(len=240) :: fchname, inpname
 logical :: sph

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in subroutine fch2dal: wrong command line argument!'
  write(6,'(A,/)') ' Example (R(O)HF, CAS): fch2dal a.fch'
  stop
 end if

 fchname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)

 i = index(fchname, '.fch', back=.true.)
 inpname = fchname(1:i-1)//'.inp'

 call fch2inp_wrap(fchname, .false., 0, 0) ! generate GAMESS .inp file

 call check_sph(fchname, sph)
 if(sph) then
  i = system('bas_gms2dal '//TRIM(inpname)//' -sph')
 else
  i = system('bas_gms2dal '//TRIM(inpname))
 end if

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine fch2dal: failed to call utility&
                  & bas_gms2dal. Two possible reasons:'
  write(6,'(A)')   '(1) The file '//TRIM(fchname)//' may be incomplete.'
  write(6,'(A)')   '(2) You forgot to compile the utility bas_gms2dal.'
  write(6,'(A,/)') '(3) This is a bug of the utility bas_gms2dal.'
  stop
 end if

 call delete_file(inpname)
 call fch2dal(fchname)
end program main

! read the MOs in .fch(k) file and adjust its d,f,g,h functions order of Gaussian
!  to that of Dalton
subroutine fch2dal(fchname)
 use fch_content, only: check_uhf_in_fch
 implicit none
 integer :: i, k, length, fid, fid1, nbf, nif, RENAME
 integer :: n6dmark,n10fmark,n15gmark,n21hmark
 integer :: n5dmark,n7fmark, n9gmark, n11hmark
 integer, allocatable :: shell_type(:), shell2atom_map(:)
 ! mark the index where d, f, g, h functions begin
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 real(kind=8), allocatable :: coeff(:,:)
 character(len=24) :: data_string
 character(len=240) :: buf, dalfile, dalfile1
 character(len=240), intent(in) :: fchname
 logical :: uhf, sph

 buf = ' '
 call check_uhf_in_fch(fchname, uhf)
 if(uhf) then
  write(6,'(/,A)') 'ERROR in subroutine fch2dal: UHF wave function not supported.'
  write(6,'(A)') 'Because there is no UHF method in Dalton. You can compute&
                 & ROHF instead.'
  stop
 end if

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(coeff(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

 if(ANY(shell_type<-1) .and. ANY(shell_type>1)) then
  write(6,'(A)') 'ERROR in subroutine fch2dal: mixed spherical harmonic/&
                 &Cartesian functions detected.'
  write(6,'(A)') 'You probably used a basis set like 6-31G(d) in Gaussian. Its&
                 & default setting is (6D,7F).'
  write(6,'(A)') "You need to add '5D 7F' or '6D 10F' keywords in Gaussian."
  stop
 else if( ANY(shell_type>1) ) then
  sph = .false.
 else
  sph = .true.
 end if

! first we adjust the basis functions in each MO according to the Shell to atom map
! this is to ensure that D comes after L functions
 ! split the 'L' into 'S' and 'P'
 call split_L_func(k, shell_type, shell2atom_map, length)

 ! sort the shell_type, shell_to_atom_map by ascending order
 ! MOs will be adjusted accordingly
 call sort_shell_and_mo(length, shell_type, shell2atom_map, nbf, nif, coeff)
 deallocate(shell2atom_map)
! adjust done

! then we adjust the basis functions in each MO according to the type of basis functions
 k = length  ! update k
 n6dmark = 0; n10fmark = 0; n15gmark = 0; n21hmark = 0
 n5dmark = 0; n7fmark = 0; n9gmark = 0; n11hmark = 0
 allocate(d_mark(k), source=0)
 allocate(f_mark(k), source=0)
 allocate(g_mark(k), source=0)
 allocate(h_mark(k), source=0)
 nbf = 0
 do i = 1, k, 1
  select case(shell_type(i))
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   nbf = nbf + 3
  case(-1)   ! SP or L
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
 end do
 deallocate(shell_type)

 ! adjust the order of d, f, etc. functions (share subroutines with fch2inporb)
 call fch2inporb_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, d_mark, &
                             f_mark, g_mark, h_mark, nbf, nif, coeff)
 call fch2inporb_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, d_mark, &
                              f_mark, g_mark, h_mark, nbf, nif, coeff)
! adjustment finished
 deallocate(d_mark, f_mark, g_mark, h_mark)

 i = index(fchname, '.fch', back=.true.)
 dalfile = fchname(1:i-1)//'.dal'
 dalfile1 = fchname(1:i-1)//'.dal1'
 open(newunit=fid,file=TRIM(dalfile),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(dalfile1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:9) == '.PUNCHOUT') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine fch2dal: no '' found in file "//TRIM(dalfile)
  close(fid)
  close(fid1,status='delete')
  stop
 end if
 close(fid,status='delete')

 write(fid1,'(A)') '.MOSTART'
 write(fid1,'(A)') 'FORM18'
 write(fid1,'(A)') '.PUNCHOUTPUTORBITALS'
 data_string = ' '
 call fdate(data_string)
 write(fid1,'(A)') '**MOLORB (punched by fch2dal of MOKIT '//TRIM(data_string)//')'
 do i = 1, nif, 1
  write(fid1,'(4F18.14)') (coeff(k,i),k=1,nbf)
 end do ! for i
 deallocate(coeff)

 write(fid1,'(A)') '**END OF INPUT'
 close(fid1)
 i = RENAME(TRIM(dalfile1), TRIM(dalfile))
end subroutine fch2dal

