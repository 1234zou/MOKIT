! Transform MOs from Gaussian -> REST
! Note: the atoms of one element must share one basis set. REST does not support
!  H1 using STO-6G while H2 using cc-pVDZ.
! TODO: ghost atoms are not supported yet.

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i, k, disp_type ! 0/1/2/3 for none/D3/D3BJ/D4
 character(len=4) :: str4
 character(len=16) :: verstr
 character(len=27), parameter :: error_warn = 'ERROR in prorgam fch2rest: '
 character(len=30) :: dftname, dftname1
 character(len=240) :: fchname
 logical :: new_format ! .true. means REST >= 2026.1.1

 i = iargc()
 if(i < 1) then
  write(6,'(/,1X,A)') error_warn//'wrong command line arguments!'
  write(6,'(A)')  ' Example 1 (R(O)HF/UHF)   : fch2rest h2o.fch'
  write(6,'(A)')  ' Example 2 (MP2)          : fch2rest h2o.fch -wft "MP2"'
  write(6,'(A)')  '                          : fch2rest h2o.fch -dft "MP2"'
  write(6,'(A)')  ' Example 3 (DFT)          : fch2rest h2o.fch -dft "B3LYP"'
  write(6,'(A)')  ' Example 4 (DFT-D3)       : fch2rest h2o.fch -dft "B3LYP D3"'
  write(6,'(A)')  '                            fch2rest h2o.fch -dft "B3LYP D3BJ"'
  write(6,'(A)')  ' Example 5 (DFT-D4)       : fch2rest h2o.fch -dft "B3LYP D4"'
  write(6,'(A)')  ' Example 6 (double hybrid): fch2rest h2o.fch -dft "XYG3"'
  write(6,'(A)')  '                            fch2rest h2o.fch -dft "XYGJOS"'
  write(6,'(A)')  ' Example 7 (find functional in .fch automatically):'
  write(6,'(A)')  '                            fch2rest h2o.fch -dft auto'
  write(6,'(A)')  ' Example 8 (REST < 2026.1.1 reads the legacy chkfile format):'
  write(6,'(A,/)')'                            fch2rest h2o.fch -ver 2025.02'
  stop
 end if

 fchname = ' '; dftname = ' '; disp_type = 0; new_format = .true.
 call getarg(1, fchname)
 call require_file_exist(fchname)

 k = 2
 do while(k <= i)
  str4 = ' '
  call getarg(k, str4)
  select case(TRIM(str4))
  case('-dft','-wft')
   if(k == i) then
    write(6,'(/,A)') error_warn//'the flag `'//TRIM(str4)//'` needs a value.'
    stop
   end if
   call getarg(k+1, dftname)
   if(TRIM(str4) == '-wft') then
    if(TRIM(dftname) /= 'MP2') then
     write(6,'(/,A)') error_warn//'unrecognized dftname='//TRIM(dftname)
     stop
    end if
   end if
   call check_dftname_in_fch2rest(fchname, str4, dftname, disp_type)
  case('-ver')
   if(k == i) then
    write(6,'(/,A)') error_warn//'the flag `'//TRIM(str4)//'` needs a value.'
    stop
   end if
   verstr = ' '
   call getarg(k+1, verstr)
   call is_new_rest_chkfile_format(verstr, new_format)
  case default
   write(6,'(/,A)') error_warn//"the flags can only be -dft, -wft or -ver."
   write(6,'(A)') "But got `"//TRIM(str4)//"`"
   stop
  end select
  k = k + 2
 end do ! for while

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 dftname1 = dftname
 call upper(dftname1)
 select case(TRIM(dftname1))
 case('M052X')
  dftname = 'M05-2X'
 case('M062X')
  dftname = 'M06-2X'
 case('M06L')
  dftname = 'M06-L'
 case('MN15L')
  dftname = 'MN15-L'
 case('XDHPBE0')
  dftname = 'xDH-PBE0'
 case('RXDH7')
  dftname = 'R-xDH7'
 case('XYGJ-OS')
  dftname = 'XYGJOS'
 end select

 call fch2rest(fchname, dftname, disp_type, new_format)
end program main

! Check the method name given by the user after the `-dft` or `-wft` flag.
! `str4` is the flag name, which is required to be '-dft' when dftname='auto'.
! Dispersion corrections D3/D3BJ/D4 are recognized and stripped from dftname.
subroutine check_dftname_in_fch2rest(fchname, str4, dftname, disp_type)
 implicit none
 integer :: i, k
 integer, intent(inout) :: disp_type ! type of dispersion correction
 character(len=4), intent(in) :: str4
 character(len=30), intent(inout) :: dftname
 character(len=240), intent(in) :: fchname
 character(len=27), parameter :: error_warn = 'ERROR in subroutine check_df&
                                              &tname_in_fch2rest: '
 logical :: is_hf, rotype, untype

 if(INDEX(dftname,',') > 0) then
  write(6,'(/,A)') error_warn//'the symbol "," is not allowed in the method name.'
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
  if(k >= 3) then
   if(dftname(k-2:k)=='-d3' .or. dftname(k-2:k)=='-d4') dftname(k-2:k-2) = ' '
  end if
  if(k >= 5) then
   if(dftname(k-4:k)=='-d3bj') dftname(k-4:k-4) = ' '
  end if
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
end subroutine check_dftname_in_fch2rest

! Determine whether the REST program of the given version reads the new chkfile
! format (with the 'molecule/basis4elem' dataset). Only REST >= 2026.1.1 reads
! the new format; older versions need the legacy chkfile format. Leading 'v'
! of a version tag is tolerated, and the patch level is optional.
subroutine is_new_rest_chkfile_format(version, new_format)
 implicit none
 integer :: i, k, maj, mino, pat, ios
 character(len=16), intent(in) :: version
 character(len=16) :: buf
 logical, intent(out) :: new_format

 buf = version
 k = LEN_TRIM(buf)
 if(k == 0) then
  write(6,'(/,A)') 'ERROR in subroutine is_new_rest_chkfile_format: empty v&
                   &ersion string.'
  stop
 end if
 if(buf(1:1) == 'v' .or. buf(1:1) == 'V') then
  buf = buf(2:k)//' '
  k = k - 1
 end if

 do i = 1, k, 1
  if(buf(i:i) == '.') buf(i:i) = ' '
 end do ! for i

 maj = 0; mino = 0; pat = 0
 read(buf(1:k),*,iostat=ios) maj, mino, pat
 if(ios /= 0) then
  maj = 0; mino = 0; pat = 0
  read(buf(1:k),*,iostat=ios) maj, mino
  if(ios /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine is_new_rest_chkfile_format: unr&
                    &ecognized REST version `'//TRIM(version)//'`.'
   stop
  end if
 end if

 if(maj > 2026) then
  new_format = .true.
 else if(maj == 2026) then
  if(mino > 1) then
   new_format = .true.
  else if(mino == 1) then
   new_format = (pat >= 1)
  else
   new_format = .false.
  end if
 else
  new_format = .false.
 end if
end subroutine is_new_rest_chkfile_format

subroutine fch2rest(fchname, dftname, disp_type, new_format)
 use fch_content
 implicit none
 integer :: i, icart
 integer, intent(in) :: disp_type ! type of dispersion correction
 character(len=30), intent(in) :: dftname
 character(len=30), parameter :: error_warn = 'ERROR in subroutine fch2rest: '
 character(len=240), intent(in) :: fchname
 character(len=240) :: inpname, dirname, basjson, geomjson
 logical, intent(in) :: new_format
 logical :: uhf, ghf, sph, sfx2c
 logical, allocatable :: ghost(:) ! size natom

 call find_specified_suffix(fchname, '.fch', i)
 inpname = fchname(1:i-1)//'.in'
 dirname = fchname(1:i-1)//'-basis'

 call check_ghf_in_fch(fchname, ghf) ! determine whether GHF
 if(ghf) then
  write(6,'(/,A)') error_warn//'GHF is unsupported currently.'
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

 if(ANY(ghost)) then
  write(6,'(/,A)') error_warn//'ghost atoms are unsupported currently.'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 sfx2c = .false.
 call find_irel_in_fch(fchname, irel)
 select case(irel)
 case(-1) ! non-relativistic calculation
 case(-3) ! sfX2C, sfX2C1e
  sfx2c = .true.
 case default
  write(6,'(/,A)') error_warn//'relativistic type cannot be recognized.'
  write(6,'(A,I0)')'Currently only NONE/sfX2C are supported. But got irel=',irel
  stop
 end select

 sph = .true.
 call find_icart_from_shell_type(.false., ncontr, shell_type, icart)
 if(icart == 2) sph = .false.
 if(.not. sph) then
  write(6,'(/,A)') 'ERROR in subroutine fch2rest: Cartesian-type basis function&
                   &s (6D 10F) are'
  write(6,'(A)') 'not fully tested. It can be used in the near future. Currentl&
                 &y you can write'
  write(6,'(A)') '`5D 7F` in .gjf file, in order to use spherical harmonic type&
                 & basis functions.'
  !write(6,'(/,A)') 'Warning from subroutine fch2rest: Cartesian-type basis func&
  !                 &tions (6D 10F)'
  !write(6,'(A)') 'detected. You must REST >= , otherwise the REST result would &
  !               &be incorrect.'
  stop
 end if

 call write_rest_in_and_basis(inpname, dftname, disp_type, charge, mult, natom,&
                              elem, coor, sph, uhf, sfx2c, ghost, new_format)
 deallocate(ghost)

 if(new_format) then
  call find_specified_suffix(fchname, '.fch', i) ! i is altered by the loops above
  basjson = fchname(1:i-1)//'.basis4elem.json'
  geomjson = fchname(1:i-1)//'.geom.json'
  call write_rest_basis4elem_json(fchname, sph, basjson)
  call write_rest_geom_json(fchname, geomjson)
 else ! legacy REST reads the basis set from local JSON files
  call gen_rest_bas_dir(dirname)
 end if

 call free_arrays_in_fch_content()
 call rest_fch2pchk(fchname, mult, uhf, charge, new_format)
end subroutine fch2rest

! Write the basis set data of each atom into a JSON file, in the format of
! Vec<Basis4Elem> required by the 'molecule/basis4elem' dataset of a REST
! chkfile (see fileop/chkfile.rs of REST). The contraction coefficients are
! normalized in the libcint convention (BasCell::basis_normalization of REST):
!   coefficients = c*N(l,a)/sqrt(S), native_coefficients = c
! where N(l,a) = sqrt(2^(l+2)*(2a)^(l+1.5)/((2l+1)!!*sqrt(pi))) is the radial
! normalization factor of the primitive GTO r^l*exp(-a*r^2), and S is the
! self-overlap of the primitive-normalized contraction.
subroutine write_rest_basis4elem_json(fchname, sph, jsonname)
 use fch_content
 use basis_data, only: ncol, nline, bas4atom, init_bas4atom_for_an_atom, clear_bas4atom
 implicit none
 integer :: i, j, k, m, iatom, i1, i2, highest, fid, nao, ioff, nl, nc, df
 integer :: n1, n2, iecp
 character(len=16) :: str16
 character(len=240), intent(in) :: fchname, jsonname
 real(kind=8) :: fac, s_ovlp
 real(kind=8), allocatable :: gto_norm(:) ! radial norm of each primitive
 real(kind=8), allocatable :: coeff1(:,:) ! primitive-normalized contractions
 real(kind=8), allocatable :: coeff2(:,:) ! libcint-normalized contractions
 logical, intent(in) :: sph
 logical, allocatable :: ecp(:) ! size natom

 allocate(ecp(natom), source=.false.)
 if(LenNCZ > 0) then
  where(LPSkip == 0) ecp = .true.
 end if

 open(newunit=fid,file=TRIM(jsonname),status='replace')
 write(fid,'(A)') '['

 iatom = 1; i1 = 1; i2 = 1; ioff = 0 ! initialization
 do while(.true.)
  call init_bas4atom_for_an_atom(iatom, ncontr, nprim, shell_type, prim_per_shell,&
   shell2atom_map, prim_exp, contr_coeff, contr_coeff_sp, i1, i2, highest)

  nao = 0 ! number of AO basis functions of this atom
  do i = 0, highest, 1
   if(sph) then
    nao = nao + (2*i+1)*ncol(i)
   else
    nao = nao + (i+1)*(i+2)/2*ncol(i)
   end if
  end do ! for i

  write(fid,'(A)') ' {'
  write(fid,'(A)') '  "electron_shells": ['

  do i = 0, highest, 1
   nl = nline(i); nc = ncol(i)
   df = 1 ! (2*i+1)!!
   do j = 3, 2*i+1, 2
    df = df*j
   end do ! for j

   allocate(gto_norm(nl), coeff1(nc,nl), coeff2(nc,nl))
   fac = 2d0**(i+2)/(dble(df)*dsqrt(4d0*datan(1d0)))
   do j = 1, nl, 1
    gto_norm(j) = dsqrt(fac*(2d0*bas4atom(i)%prim_exp(j))**(dble(i)+1.5d0))
   end do ! for j
   do k = 1, nc, 1
    coeff1(k,1:nl) = bas4atom(i)%coeff(k,1:nl)*gto_norm(1:nl)
   end do ! for k
   do k = 1, nc, 1
    s_ovlp = 0d0
    do j = 1, nl, 1
     do m = 1, nl, 1
      s_ovlp = s_ovlp + coeff1(k,j)*coeff1(k,m)/(fac*(bas4atom(i)%prim_exp(j)&
               + bas4atom(i)%prim_exp(m))**(dble(i)+1.5d0))
     end do ! for m
    end do ! for j
    coeff2(k,1:nl) = coeff1(k,1:nl)/dsqrt(s_ovlp)
   end do ! for k

   write(fid,'(16X,A)') '{'
   write(fid,'(20X,A)') '"function_type": "gto",'
   write(fid,'(20X,A)') '"region": null,'
   write(fid,'(20X,A,/,21X,I4,/,20X,A)') '"angular_momentum": [',i,'],'

   write(fid,'(20X,A)') '"exponents": ['
   do j = 1, nl-1, 1
    call dp2str16(bas4atom(i)%prim_exp(j), str16, m)
    write(fid,'(24X,A)') str16(1:m)//','
   end do ! for j
   call dp2str16(bas4atom(i)%prim_exp(nl), str16, m)
   write(fid,'(24X,A)') str16(1:m)
   write(fid,'(20X,A)') '],' ! exponents

   call write_coeff_json(fid, 'coefficients', coeff2, nc, nl, .false.)
   call write_coeff_json(fid, 'native_coefficients', bas4atom(i)%coeff, nc, nl, .true.)
   deallocate(gto_norm, coeff1, coeff2)

   if(i < highest) then
    write(fid,'(16X,A)') '},'
   else
    write(fid,'(16X,A)') '}'
   end if
  end do ! for i

  write(fid,'(A)') '  ],' ! electron_shells
  write(fid,'(A)') '  "references": null,'

  if(ecp(iatom)) then
   write(fid,'(A,I0,A)') '  "ecp_electrons": ', IDNINT(RNFroz(iatom)), ','
   write(fid,'(A)') '  "ecp_potentials": ['
   k = Lmax(iatom); nl = COUNT(KFirst(iatom,:) > 0)

   do iecp = 1, nl, 1
    n1 = KFirst(iatom,iecp); n2 = KLast(iatom,iecp)
    write(fid,'(16X,A)') '{'
    if(iecp == 1) then
     write(fid,'(20X,A,I0,A)') '"angular_momentum": [', k, '],'
    else
     write(fid,'(20X,A,I0,A)') '"angular_momentum": [', iecp-2, '],'
    end if
    write(fid,'(20X,A)') '"ecp_type": "scalar_ecp",'
    write(fid,'(20X,A)',advance='no') '"r_exponents": ['
    do j = n1, n2-1, 1
     write(fid,'(I0,A)',advance='no') NLP(j), ', '
    end do ! for j
    write(fid,'(I0,A)') NLP(n2), '],'
    write(fid,'(20X,A)',advance='no') '"gaussian_exponents": ['
    do j = n1, n2-1, 1
     call dp2str16(ZLP(j), str16, m)
     write(fid,'(A)',advance='no') str16(1:m)//', '
    end do ! for j
    call dp2str16(ZLP(n2), str16, m)
    write(fid,'(A)') str16(1:m)//'],'
    write(fid,'(20X,A)',advance='no') '"coefficients": [['
    do j = n1, n2-1, 1
     call dp2str16(CLP(j), str16, m)
     write(fid,'(A)',advance='no') str16(1:m)//', '
    end do ! for j
    call dp2str16(CLP(n2), str16, m)
    write(fid,'(A)') str16(1:m)//']]'
    if(iecp < nl) then
     write(fid,'(16X,A)') '},'
    else
     write(fid,'(16X,A)') '}'
    end if
   end do ! for iecp

   write(fid,'(A)') '  ],' ! ecp_potentials
  else
   write(fid,'(A)') '  "ecp_electrons": null,'
   write(fid,'(A)') '  "ecp_potentials": null,'
  end if

  write(fid,'(A,I0,A,I0,A)') '  "global_index": [', ioff, ',', nao, ']'
  if(iatom < natom) then
   write(fid,'(A)') ' },'
  else
   write(fid,'(A)') ' }'
  end if

  ioff = ioff + nao
  call clear_bas4atom()
  if(iatom == natom) exit
  iatom = iatom + 1
 end do ! for while

 write(fid,'(A)') ']'
 close(fid)
 deallocate(ecp)
end subroutine write_rest_basis4elem_json

! Print one array of contraction coefficients into a JSON array of arrays.
subroutine write_coeff_json(fid, keyname, coeff, nc, nl, last)
 implicit none
 integer, intent(in) :: fid, nc, nl
 character(len=*), intent(in) :: keyname
 real(kind=8), intent(in) :: coeff(nc,nl)
 logical, intent(in) :: last ! whether this is the last key of the shell object
 integer :: j, k, m
 character(len=16) :: str16

 write(fid,'(20X,A)') '"'//TRIM(keyname)//'": ['
 do k = 1, nc, 1
  write(fid,'(24X,A)',advance='no') '['
  do j = 1, nl-1, 1
   call dp2str16(coeff(k,j), str16, m)
   write(fid,'(A)',advance='no') str16(1:m)//', '
  end do ! for j
  call dp2str16(coeff(k,nl), str16, m)
  if(k < nc) then
   write(fid,'(A)') str16(1:m)//'],'
  else
   write(fid,'(A)') str16(1:m)//']'
  end if
 end do ! for k
 if(last) then
  write(fid,'(20X,A)') ']'
 else
  write(fid,'(20X,A)') '],'
 end if
end subroutine write_coeff_json

! Write the geometry (elements and Cartesian coordinates) into a JSON file,
! in the format required by the 'molecule/geom' dataset of a REST chkfile
! (see GeomCell::geom_to_json of REST). Note that the positions in a REST
! chkfile are always stored in Bohr, regardless of the value of 'unit'.
subroutine write_rest_geom_json(fchname, jsonname)
 use fch_content
 use phys_cons, only: Bohr_const
 implicit none
 integer :: i, j, k, fid, m, natm3
 character(len=16) :: str16
 character(len=240) :: basename
 character(len=240), intent(in) :: fchname, jsonname

 call find_specified_suffix(fchname, '.fch', k)
 basename = fchname(1:k-1)

 open(newunit=fid,file=TRIM(jsonname),status='replace')
 write(fid,'(A)') '{'
 write(fid,'(A)') ' "name": "'//TRIM(basename)//'",'
 write(fid,'(A)',advance='no') ' "elem": ['
 do i = 1, natom-1, 1
  write(fid,'(A)',advance='no') '"'//TRIM(elem(i))//'", '
 end do ! for i
 write(fid,'(A)') '"'//TRIM(elem(natom))//'"],'
 write(fid,'(A)') ' "unit": "angstrom",'
 write(fid,'(A)',advance='no') ' "position": ['
 natm3 = 3*natom
 do i = 1, natm3-1, 1
  j = MOD(i-1,3) + 1 ! x/y/z component
  k = (i-1)/3 + 1    ! atom index
  call dp2str16(coor(j,k)/Bohr_const, str16, m)
  write(fid,'(A)',advance='no') str16(1:m)//', '
 end do ! for i
 call dp2str16(coor(3,natom)/Bohr_const, str16, m)
 write(fid,'(A)') str16(1:m)//'],'
 write(fid,'(A)') ' "ghost_bs_elem": [],'
 write(fid,'(A)') ' "ghost_bs_pos": []'
 write(fid,'(A)') '}'
 close(fid)
end subroutine write_rest_geom_json

! Auto-detect the basis set directory of the REST program, and use it later for
! auxiliary basis set. Since recent versions of REST does not require the absolute
! path of the auxiliary basis set, this subroutine is seldom used.
subroutine find_rest_basis_set_pool(path)
 implicit none
 integer :: k, fid
 character(len=30) :: tmpfile
 character(len=240) :: home, rest_home
 character(len=260) :: buf
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
 buf = 'which rest >'//TRIM(tmpfile)//' 2>&1'
 call run_command(TRIM(buf), .false., .false.)

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
                                   natom, elem, coor, sph, uhf, sfx2c, ghost, &
                                   new_format)
 implicit none
 integer :: i, fid
 integer, intent(in) :: disp_type, charge, mult, natom
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=2), intent(in) :: elem(natom)
 character(len=30) :: dftname1
 character(len=30), intent(in) :: dftname
 character(len=240), intent(in) :: inpname
 character(len=240) :: basename
 logical :: mp_or_dh
 logical, intent(in) :: sph, uhf, sfx2c, ghost(natom), new_format

 mp_or_dh = .false. ! not MP2 or Double-hybrid functional
 dftname1 = dftname
 call upper(dftname1)
 select case(TRIM(dftname1))
 case('MP2','XYGJOS','XYG3','XYG7','XDH-PBE0','SBGE2','ZRPS','SCSRPA','R-XDH7',&
      'RPA@PBE','RPA@B3LYP')
  mp_or_dh = .true.
 end select

 call find_specified_suffix(inpname, '.in', i)
 basename = inpname(1:i-1)

 open(newunit=fid,file=TRIM(inpname),status='replace')
 write(fid,'(A)') '[ctrl]'
 write(fid,'(2X,A)') 'print_level = 2'
 write(fid,'(2X,A)') 'num_threads = 4'
 if(new_format) then
  ! the basis set is taken from 'molecule/basis4elem' of the guessfile
  write(fid,'(2X,A)') 'basis_path = "chkfile"'
 else
  write(fid,'(2X,A)') 'basis_path = "./'//TRIM(basename)//'-basis"'
 end if
 write(fid,'(2X,A)') 'auxbas_path = "def2-universal-JKFIT"'
 write(fid,'(2X,A)') '#auxbas_path = "def2-SV(P)-JKFIT"'
 write(fid,'(2X,A)') 'guessfile = "'//TRIM(basename)//'.pchk"'
 write(fid,'(2X,A,I0)') 'charge = ', charge
 write(fid,'(2X,A,I0)') 'spin = ', mult

 if(LEN_TRIM(dftname) == 0) then
  write(fid,'(2X,A)') 'xc = "hf"'
 else
  write(fid,'(2X,A)') 'xc = "'//TRIM(dftname)//'"'
 end if

 ! Note:
 ! 1) By default, REST does not freeze any core orbitals in MP2/DH calculations.
 ! 2) It seems that currently REST does not adopt the most commonly used frozen-
 !  core option (i.e. directly sets the number of frozen core orbitals or has a
 !  set of built-in frozen core settings). So here we set `21` which resembles
 !  the number of frozen core orbitals in other quantum chemistry programs.
 if(mp_or_dh) write(fid,'(2X,A)') 'frozen_core_postscf = 21'

 if(uhf) then
  write(fid,'(2X,A)') 'spin_polarization = true'
 else
  if(mult > 1) write(fid,'(2X,A)') 'spin_polarization = false'
 end if
 if(sfx2c) write(fid,'(2X,A)') 'rel = "sfx2c"'

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

subroutine rest_fch2pchk(fchname, mult, uhf, charge, new_format)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: mult, charge
 character(len=240) :: pchk, pyname, outname, basjson, geomjson
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: uhf, new_format

 call find_specified_suffix(fchname, '.fch', i)
 write(pchk,'(A,I0,A)') fchname(1:i-1)//'.pchk'
 call get_a_random_int(k)
 write(pyname,'(A,I0,A)') fchname(1:i-1)//'_', k, '.py'
 write(outname,'(A,I0,A)') fchname(1:i-1)//'_', k, '.out'
 basjson = fchname(1:i-1)//'.basis4elem.json'
 geomjson = fchname(1:i-1)//'.geom.json'

 open(newunit=fid,file=TRIM(pyname),status='replace')
 if(new_format) then
  write(fid,'(A)') 'from mokit.lib.dump_chk import dump_scf_for_rest'
 else
  write(fid,'(A)') 'from mokit.lib.dump_chk import dump_scf_no_mol'
 end if
 write(fid,'(A)') 'from mokit.lib.rwwfn import ('
 write(fid,'(4X,A)') 'read_nbf_and_nif_from_fch,'
 write(fid,'(4X,A)') 'read_na_and_nb_from_fch,'
 write(fid,'(4X,A)') 'read_eigenvalues_from_fch,'
 if((.not.uhf) .and. mult==1) then
  write(fid,'(4X,A)') 'get_occ_from_na_nb'
 else
  write(fid,'(4X,A)') 'get_occ_from_na_nb2'
 end if
 write(fid,'(A)') ')'
 write(fid,'(A)') 'from mokit.lib.fch2py import fch2py'
 write(fid,'(A)') 'import numpy as np'

 write(fid,'(/,A)') "fchname = '"//TRIM(fchname)//"'"
 write(fid,'(A)') "chkfile = '"//TRIM(pchk)//"'"
 if(new_format) then
  write(fid,'(A)') "basis4elem_json = open('"//TRIM(basjson)//"').read()"
  write(fid,'(A)') "geom_json = open('"//TRIM(geomjson)//"').read()"
 end if
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
  if(mult == 1) then
   write(fid,'(A)') 'mo_occ = get_occ_from_na_nb(nif, na, nb)'
  else
   write(fid,'(A)') 'mo_occ = get_occ_from_na_nb2(nif, na, nb)'
  end if
 end if

 if(new_format) then
  if(uhf) then
   write(fid,'(A)') 'spin_channel = 2'
  else
   write(fid,'(A)') 'spin_channel = 1'
  end if
  write(fid,'(A)') 'num_elec = float(np.sum(mo_occ))'
  write(fid,'(A)') 'dump_scf_for_rest(chkfile,e_tot,mo_ene,mo_coeff,mo_occ,'
  write(fid,'(A)') '  basis4elem_json,geom_json,"spheric",nbf,nif,spin_channel,'
  write(fid,'(A,I0,A,I0,A)') '  spin=',mult,', charge=',charge,', num_elec=num_elec)'
 else
  write(fid,'(A)') 'dump_scf_no_mol(chkfile,e_tot,mo_ene,mo_coeff,mo_occ)'
 end if
 close(fid)

 call submit_pyscf_job(pyname, .false.)
 if(new_format) then
  call delete_files(4, [pyname, outname, basjson, geomjson])
 else
  call delete_files(2, [pyname, outname])
 end if
 call remove_dir('__pycache__')
end subroutine rest_fch2pchk

