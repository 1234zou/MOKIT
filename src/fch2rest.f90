! Transform MOs from Gaussian -> REST
! Note: the atoms of one element must share one basis set. REST does not support
!  H1 using STO-6G while H2 using cc-pVDZ.
! TODO: ghost atoms are not supported yet.

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i, k, disp_type ! 0/1/2/3 for none/D3/D3BJ/D4
 character(len=4) :: str4
 character(len=15) :: dftname
 character(len=27), parameter :: error_warn = 'ERROR in prorgam fch2rest: '
 character(len=240) :: fchname
 logical :: is_hf, rotype, untype

 i = iargc()
 if(.not. (i==1 .or. i==3)) then
  write(6,'(/,1X,A)') error_warn//'wrong command line arguments!'
  write(6,'(A)')  ' Example 1 (R(O)HF/UHF)   : fch2rest h2o.fch'
  write(6,'(A)')  " Example 2 (MP2)          : fch2rest h2o.fch -wft 'MP2'"
  write(6,'(A)')  "                          : fch2rest h2o.fch -dft 'MP2'"
  write(6,'(A)')  " Example 3 (DFT)          : fch2rest h2o.fch -dft 'B3LYP'"
  write(6,'(A)')  " Example 4 (DFT-D3)       : fch2rest h2o.fch -dft 'B3LYP D3'"
  write(6,'(A)')  "                            fch2rest h2o.fch -dft 'B3LYP D3BJ'"
  write(6,'(A)')  " Example 5 (DFT-D4)       : fch2rest h2o.fch -dft 'B3LYP D4'"
  write(6,'(A)')  " Example 6 (double hybrid): fch2rest h2o.fch -dft 'XYG3'"
  write(6,'(A)')  "                            fch2rest h2o.fch -dft 'XYGJOS'"
  write(6,'(A)')  ' Example 7 (find functional in .fch automatically):'
  write(6,'(A,/)')'                            fch2rest h2o.fch -dft auto'
  stop
 end if

 fchname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)
 str4 = ' '; dftname = ' '; disp_type = 0

 if(i == 3) then
  call getarg(2, str4)
  select case(str4)
  case('-dft')
  case('-wft')
   call getarg(3, dftname)
   if(TRIM(dftname) /= 'MP2') then
    write(6,'(/,A)') error_warn//'`MP2` can only be used after the `-wft` flag.'
    write(6,'(A)') 'But got `'//dftname//'`'
    stop
   end if
  case default
   write(6,'(/,A)') error_warn//"the 2nd argument can only be -dft or -wft."
   write(6,'(A)') "But got `"//str4//"`"
   stop
  end select

  call getarg(3, dftname)
  if(INDEX(dftname,'-')>0 .or. INDEX(dftname,',')>0) then
   write(6,'(/,A)') error_warn//"'-'/',' symbols are not allowed in the 3rd"
   write(6,'(A)') "argument. For example, XYGJ-OS should be written as 'XYGJOS'&
                  &; B3LYP-D3BJ"
   write(6,'(A)') "should be written as 'B3LYP D3BJ'."
   stop
  end if
  call lower(dftname)
  if(TRIM(dftname) == 'auto') then
   if(str4 /= '-dft') then
    write(6,'(/,A)') error_warn//'`auto` can only be used after the `-dft` flag.'
    write(6,'(A)') 'But got `'//str4//'`'
    stop
   end if
   call find_dftname_in_fch(fchname, dftname, is_hf, rotype, untype)
  else ! dftname is not 'auto'
   k = LEN_TRIM(dftname)
   i = INDEX(dftname(1:k), ' ')
   if(i > 0) then
    select case(dftname(i+1:k))
    case('d3')
     disp_type = 1
    case('d3bj')
     disp_type = 2
    case('d4')
     disp_type = 3
    case default
     write(6,'(/,A)') error_warn//'dispersion correction cannot be recognized.'
     write(6,'(A)') 'Currently only D3/D3BJ/D4 are allowed.'
     stop
    end select
    dftname(i+1:k) = ' '
   end if
  end if
 end if

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 call fch2rest(fchname, dftname, disp_type)
end program main

subroutine fch2rest(fchname, dftname, disp_type)
 use fch_content
 implicit none
 integer :: i, icart
 integer, intent(in) :: disp_type ! type of dispersion correction
 character(len=15), intent(in) :: dftname
 character(len=240), intent(in) :: fchname
 character(len=240) :: inpname, dirname
 logical :: uhf, ghf, sph
 logical, allocatable :: ghost(:) ! size natom

 call find_specified_suffix(fchname, '.fch', i)
 inpname = fchname(1:i-1)//'.in'
 dirname = fchname(1:i-1)//'-basis'

 call check_ghf_in_fch(fchname, ghf) ! determine whether GHF
 if(ghf) then
  write(6,'(/,A)') 'ERROR in subroutine fch2rest: GHF not supported currently.'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 call read_fch(fchname, uhf)         ! read content in .fch(k) file

 allocate(ghost(natom))
 do i = 1, natom, 1
  if(iatom_type(i) == 1000) then
   ghost(i) = .true.
  else
   ghost(i) = .false.
  end if
 end do ! for i

 sph = .true.
 call find_icart_from_shell_type(.false., ncontr, shell_type, icart)
 if(icart == 2) sph = .false.

 call write_rest_in_and_basis(inpname, dftname, disp_type, charge, mult, natom,&
                              elem, coor, sph, uhf, ghost)
 deallocate(ghost)

 call gen_rest_bas_dir(dirname)
 call free_arrays_in_fch_content()

 call rest_fch2pchk(fchname, uhf)
end subroutine fch2rest

! Auto-detect the basis set directory of the REST program, and use it later for
! auxiliary basis set. Since recent versions of REST does not require the absolute
! path of the auxiliary basis set, this subroutine is seldom used.
subroutine find_rest_basis_set_pool(path)
 implicit none
 integer :: i, k, fid, SYSTEM
 character(len=30) :: tmpfile
 character(len=240) :: buf, home, rest_home
 character(len=480), intent(out) :: path

 path = ' ' ! initialization

 ! If $REST_HOME is defined by the user, use $REST_HOME/rest/basis-set-pool
 rest_home = ' '
 call getenv('REST_HOME', rest_home)
 k = LEN_TRIM(rest_home)
 if(k > 0) then
  path = rest_home(1:k)//'/rest/basis-set-pool'
  return
 end if

 ! Otherwise check the path from `which rest`
 call get_a_random_int(k)
 write(tmpfile,'(A,I0)') 'rest_bas.', k
 i = SYSTEM('which rest >'//TRIM(tmpfile)//' 2>&1')

 open(newunit=fid,file=TRIM(tmpfile),status='old',position='rewind')
 read(fid,'(A)') buf
 close(fid, status='delete')

 ! if `rest` is not found, or not installed, return empty path
 if(buf(1:26) == '/usr/bin/which: no rest in') return

 k = LEN_TRIM(buf)
 path(1:k) = buf(1:k)
 if(buf(k-7:k) == 'bin/rest') path = buf(1:k-8)//'share/rest/basis-set-pool'

 if(path(1:1) == '~') then
  call getenv('HOME', home)
  path = TRIM(home)//TRIM(path(2:))
 end if
end subroutine find_rest_basis_set_pool

subroutine write_rest_in_and_basis(inpname, dftname, disp_type, charge, mult, &
                                   natom, elem, coor, sph, uhf, ghost)
 implicit none
 integer :: i, fid
 integer, intent(in) :: disp_type, charge, mult, natom
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=2), intent(in) :: elem(natom)
 character(len=15), intent(in) :: dftname
 character(len=240), intent(in) :: inpname
 character(len=240) :: basename
 logical, intent(in) :: sph, uhf, ghost(natom)
 
 call find_specified_suffix(inpname, '.in', i)
 basename = inpname(1:i-1)

 open(newunit=fid,file=TRIM(inpname),status='replace')
 write(fid,'(A)') '[ctrl]'
 write(fid,'(2X,A)') 'print_level = 2'
 write(fid,'(2X,A)') 'num_threads = 4'
 write(fid,'(2X,A)') 'basis_path = "./'//TRIM(basename)//'-basis"'
 write(fid,'(2X,A)') 'auxbas_path = "def2-SV(P)-JKFIT"'
 write(fid,'(2X,A)') '#auxbas_path = "def2-universal-JKFIT"'
 write(fid,'(2X,A)') 'chkfile = "'//TRIM(basename)//'.pchk"'
 write(fid,'(2X,A,I0)') 'charge = ', charge
 write(fid,'(2X,A,I0)') 'spin = ', mult

 if(LEN_TRIM(dftname) == 0) then
  write(fid,'(2X,A)') 'xc = "hf"'
 else
  write(fid,'(2X,A)') 'xc = "'//TRIM(dftname)//'"'
 end if

 if(uhf) then
  write(fid,'(2X,A)') 'spin_polarization = true'
 else
  if(mult > 1) write(fid,'(2X,A)') 'spin_polarization = false'
 end if

 selectcase(disp_type)
 case(0) ! no dispersion correction
 case(1) ! D3
  write(fid,'(2X,A)') 'empirical_dispersion = "d3"'
 case(2) ! D3BJ
  write(fid,'(2X,A)') 'empirical_dispersion = "d3bj"'
 case(3) ! D4
  write(fid,'(2X,A)') 'empirical_dispersion = "d4"'
 case default
  write(6,'(/,A)') 'ERROR in subroutine write_rest_in_and_basis: disp_type out &
                   &of range!'
  write(6,'(A,I0)') 'Currently only disp_type=1,2,3 are allowed. But got ',disp_type
  close(fid)
  stop
 end select

 if(.not. sph) write(fid,'(2X,A)') 'basis_type = "Cartesian"'
 write(fid,'(2X,A)') 'max_scf_cycle = 32'
 write(fid,'(2X,A)') 'outputs = ["fchk"]'
 write(fid,'(2X,A)') '# Try options below if not converged but close to converg&
                     &ence at the beginning'
 write(fid,'(2X,A)') '#start_diis_cycle = 10'
 write(fid,'(2X,A)') '#scf_acc_eev = 1e-4'
 write(fid,'(A)') '[geom]'
 write(fid,'(2X,A)') 'name = "'//TRIM(basename)//'"'
 write(fid,'(2X,A)') 'unit = "angstrom"'
 write(fid,'(2X,A)') "position = '''"

 do i = 1, natom, 1
  if(ghost(i)) then
   write(fid,'(A4,2X,3(1X,F17.8))') 'X-'//elem(i), coor(1:3,i)
  else
   write(fid,'(A2,2X,3(1X,F17.8))') elem(i), coor(1:3,i)
  end if
 end do ! for i

 write(fid,'(A)') "'''"
 close(fid)
end subroutine write_rest_in_and_basis

subroutine rest_fch2pchk(fchname, uhf)
 implicit none
 integer :: i, k, fid
 character(len=240) :: pchk, pyname, outname
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: uhf

 call find_specified_suffix(fchname, '.fch', i)
 write(pchk,'(A,I0,A)') fchname(1:i-1)//'.pchk'
 call get_a_random_int(k)
 write(pyname,'(A,I0,A)') fchname(1:i-1)//'_', k, '.py'
 write(outname,'(A,I0,A)') fchname(1:i-1)//'_', k, '.out'

 open(newunit=fid,file=TRIM(pyname),status='replace')
 write(fid,'(A)') 'from pyscf.scf.chkfile import dump_scf'
 write(fid,'(A)') 'from mokit.lib.gaussian import load_mol_from_fch'
 write(fid,'(A)') 'from mokit.lib.rwwfn import ('
 write(fid,'(4X,A)') 'read_nbf_and_nif_from_fch,'
 write(fid,'(4X,A)') 'read_na_and_nb_from_fch,'
 write(fid,'(4X,A)') 'read_eigenvalues_from_fch,'
 if(uhf) then
  write(fid,'(4X,A)') 'get_occ_from_na_nb2'
 else
  write(fid,'(4X,A)') 'get_occ_from_na_nb'
 end if
 write(fid,'(A)') ')'
 write(fid,'(A)') 'from mokit.lib.fch2py import fch2py'
 write(fid,'(A)') 'import numpy as np'

 write(fid,'(/,A)') "fchname = '"//TRIM(fchname)//"'"
 write(fid,'(A)') "chkfile = '"//TRIM(pchk)//"'"
 write(fid,'(A)') 'mol = load_mol_from_fch(fchname)'
 write(fid,'(A)') 'nbf, nif = read_nbf_and_nif_from_fch(fchname)'
 write(fid,'(A)') 'na, nb = read_na_and_nb_from_fch(fchname)'
 write(fid,'(A)') 'e_tot = 0e0'
 if(uhf) then
  write(fid,'(A)') "ene_a = read_eigenvalues_from_fch(fchname, nif, 'a')"
  write(fid,'(A)') "ene_b = read_eigenvalues_from_fch(fchname, nif, 'b')"
  write(fid,'(A)') 'mo_ene = np.array((ene_a, ene_b))'
  write(fid,'(A)') "coeff_a = fch2py(fchname, nbf, nif, 'a')"
  write(fid,'(A)') "coeff_b = fch2py(fchname, nbf, nif, 'b')"
  write(fid,'(A)') 'mo_coeff = np.array((coeff_a, coeff_b))'
  write(fid,'(A)') 'mo_occ = get_occ_from_na_nb2(nif, na, nb)'
 else
  write(fid,'(A)') "mo_ene = read_eigenvalues_from_fch(fchname, nif, 'a')"
  write(fid,'(A)') "mo_coeff = fch2py(fchname, nbf, nif, 'a')"
  write(fid,'(A)') 'mo_occ = get_occ_from_na_nb(nif, na, nb)'
 end if
 write(fid,'(A)') 'dump_scf(mol,chkfile,e_tot,mo_ene,mo_coeff,mo_occ,False)'
 close(fid)

 call submit_pyscf_job(pyname, .false.)
 call delete_files(2, [pyname, outname])
 call remove_dir('__pycache__')
end subroutine rest_fch2pchk

