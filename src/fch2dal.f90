! written by jxzou at 20210613: adjust the orders of d,f,g, etc functions in
!  .fch(k) file, to the orders in Dalton

program main
 use util_wrapper, only: formchk, fch2inp_wrap
 implicit none
 integer :: i, k, narg, icart
 character(len=4) :: str4
 character(len=26), parameter :: error_warn = 'ERROR in program fch2dal: '
 character(len=30) :: dftname
 character(len=240) :: fchname, inpname
 character(len=260) :: buf
 logical :: sph

 narg = iargc()
 if(.not. (narg==1 .or. narg==3)) then
  write(6,'(/,1X,A)') error_warn//'wrong command line arguments!'
  write(6,'(A)')   ' Example 1 (R(O)HF, CAS): fch2dal h2o.fch'
  write(6,'(A,/)') ' Example 2 (DFT):         fch2dal h2o.fch -dft "B3LYP"'
  stop
 end if

 str4 = ' '; dftname = ' '; fchname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 if(narg == 3) then
  call getarg(2, str4)
  if(str4 == '-dft') then
   call getarg(3, dftname)
  else
   write(6,'(/,A)') error_warn//'the 2nd argument can only be `-dft` currently.'
   write(6,'(A)') 'But got `'//str4//'`'
   stop
  end if
 end if

 ! if .chk file provided, convert into .fch file automatically
 k = LEN_TRIM(fchname)
 if(fchname(k-3:k) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:k-3)//'fch'
 end if

 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)

 call find_specified_suffix(fchname, '.fch', i)
 inpname = fchname(1:i-1)//'.inp'

 call fch2inp_wrap(fchname, .false., 0, 0, .true., .false., .false.)

 call find_icart_in_fch(fchname, .false., icart)
 if(icart == 2) then
  sph = .false.
  buf = 'bas_gms2dal '//TRIM(inpname)
  call run_command(TRIM(buf), .false., .false.)
 else
  sph = .true.
  buf = 'bas_gms2dal '//TRIM(inpname)//' -sph'
  call run_command(TRIM(buf), .false., .false.)
 end if

 call delete_file(inpname)
 call upper(dftname)
 select case(TRIM(dftname))
 case('B3LYP')
  dftname = 'B3LYPg'
  write(6,'(/,A)') 'Remark: B3LYP is specified by the user. B3LYPg will be writ&
                   &ten into .dal'
  write(6,'(A,/)') 'This functional is equal to B3LYP in Gaussian.'
 case('CAM-B3LYP')
  dftname = 'CAMB3LYP'
 case('PBEPBE')
  dftname = 'PBE'
 case('PBE1PBE')
  dftname = 'PBE0'
 end select

 call fch2dal(fchname, dftname, sph)
end program main

! read the MOs in .fch(k) file and adjust its d,f,g,h functions order of Gaussian
!  to that of Dalton
subroutine fch2dal(fchname, dftname, sph)
 implicit none
 integer :: i, j, k, length, fid, fid1, nbf, nif, RENAME
 integer :: n6dmark, n10fmark, n15gmark, n21hmark, n28imark
 integer :: n5dmark, n7fmark, n9gmark, n11hmark, n13imark
 integer, allocatable :: shell_type(:), shell2atom_map(:), idx(:)
 ! mark the index where d,f,g,h,i functions begin
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:), i_mark(:)
 real(kind=8), allocatable :: coeff0(:,:), coeff(:,:), norm(:)
 character(len=24) :: data_string
 character(len=30), intent(in) :: dftname
 character(len=240) :: buf, dalfile, dalfile1
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: sph
 logical :: uhf, dft

 buf = ' '
 call check_uhf_in_fch(fchname, uhf)
 if(uhf) then
  write(6,'(/,A)') 'ERROR in subroutine fch2dal: UHF wave function not supported.'
  write(6,'(A)') 'Because there is no UHF method in Dalton. You can compute ROH&
                 &F instead.'
  stop
 end if

 dft = .false.
 if(LEN_TRIM(dftname) > 0) dft = .true.

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(coeff(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

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
 allocate(d_mark(k), f_mark(k), g_mark(k), h_mark(k), i_mark(k))

 ! adjust the order of d, f, etc. functions (share subroutines with fch2inporb)
 allocate(idx(nbf))
 forall(i = 1:nbf) idx(i) = i ! initialization
 allocate(coeff0(nbf,nif), source=coeff)
 if(sph) then ! spherical harmonic
  call read_mark_from_shltyp_sph(k, shell_type, n5dmark, n7fmark, n9gmark, &
                 n11hmark, n13imark, d_mark, f_mark, g_mark, h_mark, i_mark)
  call fch2inporb_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, d_mark, &
                              f_mark, g_mark, h_mark, nbf, idx)
  forall(i=1:nif, j=1:nbf) coeff(j,i) = coeff0(idx(j),i)
 else         ! Cartesian-type basis
  allocate(norm(nbf), source=1d0)
  call read_mark_from_shltyp_cart(k, shell_type, n6dmark, n10fmark, n15gmark, &
                    n21hmark, n28imark, d_mark, f_mark, g_mark, h_mark, i_mark)
  call fch2inporb_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, d_mark, &
                               f_mark, g_mark, h_mark, nbf, idx, norm)
  forall(i=1:nif, j=1:nbf) coeff(j,i) = coeff0(idx(j),i)*norm(j)
  deallocate(norm)
 end if
 deallocate(shell_type, d_mark, f_mark, g_mark, h_mark, i_mark, coeff0, idx)
! adjustment finished

 i = INDEX(fchname, '.fch', back=.true.)
 dalfile = fchname(1:i-1)//'.dal'
 dalfile1 = fchname(1:i-1)//'.dt'

 open(newunit=fid,file=TRIM(dalfile),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(dalfile1),status='replace')

 if(dft) then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:3) == '.HF') exit
   write(fid1,'(A)') TRIM(buf)
  end do ! for while
  write(fid1,'(A,/,A)') '.DFT', TRIM(dftname)
  write(fid1,'(A,/,A)') '*DFT INPUT', '.ULTRAF'
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:9) == '.PUNCHOUT') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine fch2dal: no '.PUNCHOUT' found in file "&
                   //TRIM(dalfile)
  close(fid)
  close(fid1,status='delete')
  stop
 end if
 close(fid,status='delete')

 write(fid1,'(A)') '.PUNCHOUTPUTORBITALS'
 write(fid1,'(A)') '.MOSTART'
 write(fid1,'(A)') 'FORM18'

 ! For basis set linear dependency, Dalton require it to be printed as 0.0
 if(nif < nbf) then
  write(fid1,'(A)') '.DELETE'
  write(fid1,'(I0)') nbf-nif
  allocate(coeff0(nbf,nbf))
  coeff0(:,1:nif) = coeff
  coeff0(:,nif+1:nbf) = 0d0
 else
  allocate(coeff0(nbf,nif), source=coeff)
 end if
 deallocate(coeff)

 data_string = ' '
 call fdate(data_string)
 write(fid1,'(A)') '**MOLORB (punched by fch2dal of MOKIT '//TRIM(data_string)//')'
 do i = 1, nbf, 1
  write(fid1,'(4F18.14)') (coeff0(k,i),k=1,nbf)
 end do ! for i
 deallocate(coeff0)

 write(fid1,'(A)') '**END OF INPUT'
 close(fid1)
 i = RENAME(TRIM(dalfile1), TRIM(dalfile))
end subroutine fch2dal

! DFT ultrafine grid in Dalton is
! which may be used in the future

