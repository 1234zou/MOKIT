! written by jxzou at 20171120: generate .inp file (GAMESS) from .fch(k) file (Gaussian)
! updated by jxzou at 20171224: use '-gvb npair' in command line to permute the active orbtials
! updated by jxzou at 20180211: support UHF type MOs
! updated by jxzou at 20180620: support open shell orbitals (doublet, triplet, ...)
! updated by jxzou at 20180824: fix the bug: forgot to move open shell orbitals when npair=1
! updated by jxzou at 20190226: change data to parameter when declare constant arrays
! updated by jxzou at 20200322: renamed as fch2inp; simplify code; support ECP/PP
! updated by jxzou at 20200811: add coupling coefficients for high spin GVB (MULT>3)
! updated by jxzou at 20201113: support for spherical harmonic functions
!                               (expand 5D,7F,9G,11H to 6D,10F,15G,21H)
! updated by jxzou at 20201118: detect DKH/RESC keywords in .fch(k) file
! updated by jxzou at 20210407: remove '-uhf', add automatic determination

! If '-gvb [npair]' is specified, the orbitals in active space will be permuted to the
!  order of Gamess. In this case, you must specify the argument [npair] in command line.

! After '-gvb [npair]' is specified, you can also append '-open [nopen0]' if your system
!  is truly open shell, i.e., doublet, triplet and so on.

! When ISPHER=1 used in GAMESS, the resulting electronic energy is the same as
! '5D 7F' in Gaussian, but the MO coefficients in .dat file are expanded on Car-
! tesian functions. So, if a '5D 7F' .fch(k) file is provided, this subroutine
! will expand the MO coefficients from spherical harmonic functions to Cartesian
! functions.

! Warning: GAMESS does not support negative contraction coefficients for a CGTO
! which contains only one primitive Gaussian function. Even if -1.0 is provided
! in the .inp file, it will be automatically set to 1.0. For example, ORCA built
! -in basis set DKH-def2-TZVP for N has a negative contraction coefficient on the
! 5-th S basis function.

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: narg, k, npair, nopen0, itype
 character(len=4) :: string
 character(len=8) :: arg2, arg4, arg6
 character(len=26), parameter :: error_warn = 'ERROR in program fch2inp: '
 character(len=240) :: fchname
 logical :: fc ! whether to print the $MOFRZ section for frozen core MOs
 logical :: no_vec ! whether to print the $VEC section

 narg = iargc()
 if(narg<1 .or. narg>6) then
  write(6,'(/,1X,A)') error_warn//'wrong command line arguments!'
  write(6,'(A)')  ' Example 1 (R(O)HF/UHF/CAS): fch2inp h2o.fch'
  write(6,'(A)')  ' Example 2 (SF-CIS)        : fch2inp high_spin.fch -sfcis'
  write(6,'(A)')  ' Example 3 (SF-TDDFT)      : fch2inp high_spin.fch -sf'
  write(6,'(A)')  ' Example 4 (MRSF-CIS)      : fch2inp triplet.fch -mrsfcis'
  write(6,'(A)')  ' Example 5 (MRSF-TDDFT)    : fch2inp triplet.fch -mrsf'
  write(6,'(A)')  ' Example 6 (GVB)           : fch2inp h2o.fch -gvb [Npair]'
  write(6,'(A)')  ' Example 7 (frozen core)   : fch2inp h2o.fch -gvb [Npair] -fc'
  write(6,'(A)')  ' Example 8 (ROGVB)         : fch2inp h2o.fch -gvb [Npair] -open [Nopen]'
  write(6,'(A)')  ' Example 9 (frozen core)   : fch2inp h2o.fch -gvb [Npair] -open [Nopen] -fc'
  write(6,'(A,/)')' Example10 (no $VEC)       : fch2inp h2o.fch -novec'
  stop
 end if

 fchname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 k = LEN_TRIM(fchname)
 if(fchname(k-3:k) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:k-3)//'fch'
 end if

 npair = 0; nopen0 = 0; itype = 0; fc = .false.; no_vec = .false.
 string = ' '; arg2 = ' '; arg4 = ' '; arg6 = ' '

 if(narg > 1) then
  call getarg(2, arg2)
  select case(TRIM(arg2))
  case('-sfcis')
   itype = 1
  case('-sf')
   itype = 2
  case('-mrsfcis')
   itype = 3
  case('-mrsf')
   itype = 4
  case('-gvb')
   itype = 5
  case('-novec')
   no_vec = .true.
   if(narg > 2) then
    write(6,'(/,A)') error_warn//'-novec is incompatible with other arguments!'
    stop
   end if
  case default
   write(6,'(/,A)') error_warn//'the 2nd argument is wrong!'
   write(6,'(A)') 'It can only be one of -gvb/-sf/-mrsf/-novec'
   stop
  end select

  if(narg > 2) then
   if(TRIM(arg2) /= '-gvb') then
    write(6,'(/,A)') error_warn//'when there are more than 2 arguments, the'
    write(6,'(A)') '2nd argument can only be "-gvb". But got "'//arg2//'"'
    stop
   end if
   call getarg(3, string)
   read(string,*) npair

   if(narg == 4) then ! fch2inp h2o.fch -gvb [Npair] -fc
    call getarg(4, arg4)
    if(TRIM(arg4) == '-fc') then
     fc = .true.
    else
     write(6,'(/,A)') error_warn//'when there are only 4 arguments specified,'
     write(6,'(A)') 'the 4th argument can only be "-fc". But got "'//arg4//'"'
     stop
    end if
   else if(narg > 4) then
    call getarg(4, arg4)
    if(TRIM(arg4) /= '-open') then
     write(6,'(/,A)') error_warn//'when there are more than 4 arguments specified,'
     write(6,'(A)') 'the 4th argument can only be "-open". But got "'//arg4//'"'
     stop
    end if
    call getarg(5, string)
    read(string,*) nopen0
    if(narg == 6) then
     call getarg(6, arg6)
     if(TRIM(arg6) == '-fc') then
      fc = .true.
     else
      write(6,'(/,A)') error_warn//'when there are 6 arguments specified, the 6-th'
      write(6,'(A)') 'argument can only be "-fc". But got "'//arg6//'"'
      stop
     end if
    end if
   end if
  end if
 end if

 call fch2inp(fchname, no_vec, fc, itype, npair, nopen0)
end program main

! Generate GAMESS .inp file from Gaussian .fch(k) file.
subroutine fch2inp(fchname, no_vec, fc, itype, npair, nopen0)
 use fch_content
 implicit none
 integer :: i, j, k, m, n, n1, n2, nd, nf, ng, nh, ni, icart, fid
 integer :: ncore   ! the number of core MOs
 integer :: nif1    ! new nif, where nif is number of MOs
 integer :: nbf1    ! new nbf, where nbf is number of basis functions
 integer :: itype1  ! -3/-2/-1/0/1/2/3/4/5 for GHF/ROHF/RHF/UHF/SF-CIS/SF-TDDFT/
 !                                             MRSF-CIS/MRSF-TDDFT/GVB
 integer, intent(in) :: itype, npair, nopen0
 ! here nopen0 used since nopen already used in module fch_content
 integer, allocatable :: order(:), d_mark(:), f_mark(:), g_mark(:), h_mark(:), &
  i_mark(:)
 character(len=1) :: str = ' '
 character(len=1), parameter :: am_type(-1:6) = ['L','S','P','D','F','G','H','I']
 character(len=1), parameter :: am_type1(0:6) = ['s','p','d','f','g','h','i']
 real(kind=8), allocatable :: temp_coeff(:,:), open_coeff(:,:)
 character(len=29), parameter :: error_warn = 'ERROR in subroutine fch2inp: '
 character(len=240), intent(in) :: fchname
 character(len=240) :: inpname = ' '
 logical :: uhf, ghf, ecp, so_ecp, sph, X2C, DIIS
 logical, intent(in) :: no_vec, fc

 itype1 = itype; ecp = .false.; so_ecp = .false.; X2C = .false.; DIIS = .false.

 call find_specified_suffix(fchname, '.fch', i)
 inpname = fchname(1:i-1)//'.inp'
 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)

 sph = .true.
 call find_icart_in_fch(fchname, .false., icart)
 if(icart == 2) sph = .false.

 call find_irel_in_fch(fchname, irel)
 select case(irel)
 case(-3) ! X2C, i.e. sf-x2c1e
  X2C = .true.
  write(6,'(/,A)') 'Warning in subroutine fch2inp: X2C detected. But GAMESS doe&
                   &s not support'
  write(6,'(A)') 'X2C. DKH2 keywords will be printed into GAMESS .inp file.'
 case(-2,-1,2) ! RESC/none/DKH2
 case(0)  ! DKH0
  write(6,'(/,A)') 'Warning in subroutine fch2inp: DKH0 detected.'
  write(6,'(A)') 'But GAMESS does not support this DKH 0-th order correction.'
  write(6,'(A)') 'DKH2 keywords will be printed into GAMESS .inp file.'
 case(4) ! DKH4
  write(6,'(/,A)') 'Warning in subroutine fch2inp: DKHSO detected.'
  write(6,'(A)') 'But GAMESS does not support this DKH 4-th order correction.'
  write(6,'(A)') 'DKH2 keywords will be printed into GAMESS .inp file.'
 case default
  write(6,'(/,A)') error_warn//'irel out of range.'
  write(6,'(A,I0)') 'irel=', irel
  stop
 end select

 call check_ghf_in_fch(fchname, ghf) ! determine whether GHF
 if(ghf) then
  itype1 = -3
  write(6,'(/,A)') 'Warning in subroutine fch2inp: GHF detected in file '//&
                    TRIM(fchname)
  write(6,'(A)') 'GAMESS does not support GHF currently. But fch2inp will be&
                 & continued.'
  write(6,'(A)') 'Orbitals in the generated .inp file are meanningless.'
 end if

 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 call read_fch(fchname, uhf)         ! read content in .fch(k) file

 if((itype==1 .or. itype==2) .and. mult<3) then
  write(6,'(/,A)') error_warn//'SF-CIS/SF-TDDFT in GAMESS is supposed to be based'
  write(6,'(A)') 'on a high-spin ROHF/UHF reference, with spin multiplicity >=3&
                 &. But the spin'
  write(6,'(A,I0)') 'multiplicity in your .fch file is: ', mult
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 if((itype==3 .or. itype==4) .and. mult/=3) then
  write(6,'(/,A)') error_warn//'MRSF-CIS/MRSF-TDDFT in GAMESS can only be based'
  write(6,'(A,I0)') 'on a triplet ROHF reference! The spin multiplicity in your&
                    & .fch(k) file is ', mult
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 if(uhf) then
  nif1 = 2*nif
  if(itype==3 .or. itype==4) then
   write(6,'(/,A)') error_warn//'UHF and MRSF both activated.'
   write(6,'(A)') 'Please use an ROHF .fch(k) file.'
   stop
  else if(itype == 5) then
   write(6,'(/,A)') error_warn//'UHF and GVB both activated.'
   write(6,'(A)') 'Did you provide a wrong .fch(k) file, or specify wrong argum&
                  &ents?'
   write(6,'(A)') 'fchname='//TRIM(fchname)
   stop
  end if
 else
  nif1 = nif
  if(itype==0 .and. (.not.ghf)) then
   if(mult == 1) then
    itype1 = -1 ! RHF
   else
    itype1 = -2 ! ROHF
   end if
  end if
 end if

 if(ghf) then
  if(ANY(DABS(CLP2) > 1d-4)) so_ecp = .true. ! SOHF/SODFT
 end if

 if(itype1==5 .and. nopen0/=nopen) then
  k = nopen0 - nopen
  if(MOD(IABS(k),2) /= 0) then
   write(6,'(/,A)') error_warn//'odd number of singly occupied orbitals is changed!'
   write(6,'(A)') 'When you change the number of singly occupied orbitals, only&
                  & an even number is'
   write(6,'(A)') 'allowed, e.g. from 0 to 2/4/... (from singlet to triplet/qui&
                  &ntuplet/...). If you'
   write(6,'(A)') 'want to change an odd number of singly occupied orbitals (e.&
                  &g. from singlet to'
   write(6,'(A)') 'doublet), you need to firstly use the Python API to modify t&
                  &he charge and spin'
   write(6,'(A)') 'multiplicity in .fch, for example'
   write(6,'(A)') '```'
   write(6,'(A)') 'from mokit.lib.rwwfn import modify_charge_and_mult_in_fch'
   write(6,'(A)') "modify_charge_and_mult_in_fch(fchname='h2o.fch',charge=1,mult=2)"
   write(6,'(A)') '```'
   write(6,'(A)') 'then use the `fch2inp` utility.'
   stop
  end if
  na = na + k/2
  nb = nb - k/2
  mult = mult + k
 end if

 ncore = na - npair - nopen0
 ! arrays eigen_e_a and eigen_e_b are useless for GAMESS inp file
 deallocate(eigen_e_a)
 if(allocated(eigen_e_b)) deallocate(eigen_e_b)

 if(sph) then
  nbf1 = nbf + COUNT(shell_type==-2) + 3*COUNT(shell_type==-3) + &
           & 6*COUNT(shell_type==-4) + 10*COUNT(shell_type==-5)
  ! [6D,10F,15G,21H] - [5D,7F,9G,11H] = [1,3,6,10]
 else
  nbf1 = nbf
 end if

 ! in GAMESS inp format file, the contr_coeff_sp is zero when there is no 'L'/'SP'
 if(.not. allocated(contr_coeff_sp)) allocate(contr_coeff_sp(nprim), source=0d0)

 if(LenNCZ > 0) ecp = .true.
 if(natom>2 .and. ALL(elem=='H ')) DIIS = .true.
 ! By default, SOSCF is good for most systems; but for a system containing only
 ! hydrogen atoms, DIIS is better than SOSCF

 ! create the GAMESS .inp file and print the keywords information
 call creat_gamess_inp_head(inpname, charge, mult, ncore, npair, nopen0, nif, &
                       nbf, na, nb, itype1, irel, fc, uhf, ecp, sph, X2C, DIIS)
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')

 ! for ghost atoms (0 charge, has basis function), make ielem(i) negative,
 ! which can be recognized by GAMESS
 do i = 1, natom, 1
  if(iatom_type(i) == 1000) ielem(i) = -ielem(i)
 end do ! for i

 ! print basis sets into the .inp file
 write(fid,'(A2,2X,I3,A1,3(1X,F18.8))') elem(1), ielem(1), '.', coor(1:3,1)
 k = 0
 do i = 1, ncontr, 1
  m = shell2atom_map(i)
  if(m > 1) then
   if(shell2atom_map(i-1) == m-1) then
    write(fid,'(A2,2X,I3,A1,3(1X,F18.8))') elem(m), ielem(m), '.', coor(1:3,m)
   end if
  end if

  m = shell_type(i); n = prim_per_shell(i)
  if(m < -1) m = -m
  write(fid,'(3X,A1,1X,I2)') am_type(m), n
  do j = k+1, k+n, 1
   write(fid,'(2X,I2,3(2X,ES15.8))') j-k, prim_exp(j), contr_coeff(j), &
                                     contr_coeff_sp(j)
  end do ! for j

  if(i == ncontr) then
   write(fid,'(/)',advance='no')
  else ! i < ncontr
   if(shell2atom_map(i+1) == shell2atom_map(i)+1) write(fid,'(/)',advance='no')
  end if

  k = k + n
 end do ! for i

 k = shell2atom_map(ncontr)
 if(k < natom) then
  write(6,'(/,A)') error_warn//'shell2atom_map(ncontr)<natom.'
  write(6,'(A)') 'Internal inconsistency. Stop and check.'
  write(6,'(2(A,I0))') 'k=', k, ', natom=', natom
  close(fid)
  stop
 end if
 deallocate(ielem, elem, coor, shell2atom_map, prim_per_shell, prim_exp, &
            contr_coeff, contr_coeff_sp)
 write(fid,'(A)') ' $END'

 ! print ECP/PP (if any) into the .inp file
 if(ecp) then
  write(fid,'(A)') ' $ECP'

  do i = 1, natom, 1
   if(LPSkip(i) /= 0) then
    write(fid,'(I0,A9)') i, ' ECP-NONE'
    cycle
   else
    write(fid,'(I0,A8,2X,I3,2X,I2)') i,'-ECP GEN', NINT(RNFroz(i)), LMax(i)
    str = am_type1(LMax(i))
    do j = 1, 10, 1
     n1 = KFirst(i,j); n2 = KLast(i,j)
     if(n1 == 0) exit
     if(j == 1) then
      write(fid,'(I0,5X,A)') n2-n1+1,'----- '//str//'-ul potential -----'
     else
      write(fid,'(I0,5X,A)') n2-n1+1,'----- '//am_type1(j-2)//'-'//str//' potential -----'
     end if
     if(so_ecp) then ! GHF (SOHF, SODFT)
      ! Note: SO-ECP of GHF is not supported in GAMESS. I add SO-ECP as the 3rd
      ! column. If in the future GAMESS supports SO-ECP of GHF, then this may be
      ! updated
      do n = n1, n2, 1
       if(j==1 .or. j==2) then
        write(fid,'(3X,ES15.8,3X,I2,2X,ES15.8)') CLP(n),NLP(n),ZLP(n)
       else
        write(fid,'(3X,ES15.8,3X,I2,2(2X,ES15.8))') CLP(n),NLP(n),ZLP(n),CLP2(n)
       end if
      end do ! for n
     else            ! R(O)HF, UHF
      do n = n1, n2, 1
       write(fid,'(3X,ES15.8,3X,I2,2X,ES15.8)') CLP(n), NLP(n), ZLP(n)
      end do ! for n
     end if
    end do ! for j
   end if
  end do ! for i

  write(fid,'(A)') ' $END'
  deallocate(KFirst, KLast, Lmax, LPSkip, NLP, RNFroz, CLP, CLP2, ZLP)
 end if

 if(no_vec) then ! do not print the $VEC section
  close(fid)
  deallocate(alpha_coeff)
  if(allocated(beta_coeff)) deallocate(beta_coeff)
  return
 end if
 allocate(temp_coeff(nbf1,nif1))

 if(sph) then ! spherical harmonic functions, expanding
  if(uhf) then
   allocate(open_coeff(nbf,nif1))
   open_coeff(:,1:nif) = alpha_coeff
   open_coeff(:,nif+1:2*nif) = beta_coeff
  else
   allocate(open_coeff(nbf,nif))
   open_coeff = alpha_coeff
  end if
  call mo_sph2cart(ncontr, shell_type, nbf, nbf1, nif1, open_coeff, temp_coeff)
  deallocate(open_coeff)
 else         ! Cartesian functions
  temp_coeff(:,1:nif) = alpha_coeff
  if(uhf) temp_coeff(:,nif+1:nif1) = beta_coeff
 end if
 deallocate(alpha_coeff)
 if(allocated(beta_coeff)) deallocate(beta_coeff)

 ! record the indices of Cartesian f, g and h functions
 allocate(d_mark(ncontr), f_mark(ncontr), g_mark(ncontr), h_mark(ncontr), &
          i_mark(ncontr))
 call read_mark_from_shltyp_cart(ncontr, shell_type, nd, nf, ng, nh, ni, &
                                 d_mark, f_mark, g_mark, h_mark, i_mark)
 deallocate(d_mark, i_mark, shell_type)
 ! done recording

 ! adjust the order of Cartesian f, g, h functions
 do i = 1, nf, 1
  call fch2inp_permute_10f(nif1, temp_coeff(f_mark(i)+3:f_mark(i)+8,:))
 end do
 do i = 1, ng, 1
  call fch2inp_permute_15g(nif1, temp_coeff(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1, nh, 1
  call fch2inp_permute_21h(nif1, temp_coeff(h_mark(i):h_mark(i)+20,:))
 end do
 deallocate(f_mark, g_mark, h_mark)
 ! adjustment finished

 ! if nbf1 > nbf, the following two lines are auto-reallocation
 allocate(alpha_coeff(nbf1,nif), source=temp_coeff(:,1:nif))
 if(uhf) allocate(beta_coeff(nbf1,nif), source=temp_coeff(:,nif+1:2*nif))
 deallocate(temp_coeff)
 nbf = nbf1 ! update nbf

 ! if active orbitals in GAMESS order are required, permute them
 if(itype == 5) then ! GVB
  if(npair > 1) then
   allocate(order(2*npair), source=0)
   allocate(temp_coeff(nbf,2*npair), source=0d0)

   if(nopen0 > 0) open_coeff = alpha_coeff(:,na-nopen0+1:na)
   forall(i = 1:npair)
    order(2*i-1) = i
    order(2*i) = 2*npair + 1 - i + nopen0
   end forall
   forall(i = 1:2*npair) temp_coeff(:,i) = alpha_coeff(:,order(i)+ncore)
   forall(i = 1:2*npair) alpha_coeff(:,i+ncore+nopen0) = temp_coeff(:,i)
   deallocate(order, temp_coeff)
   if(nopen0 > 0) then
    alpha_coeff(:,ncore+1:ncore+nopen0) = open_coeff(:,1:nopen0)
    deallocate(open_coeff)
   end if

  else if(npair==1 .and. nopen0>0) then
   open_coeff = alpha_coeff(:,na-nopen0+1:na)
   alpha_coeff(:,na) = alpha_coeff(:,ncore+1)
   alpha_coeff(:,ncore+1:na-1) = open_coeff
   deallocate(open_coeff)
  end if
 end if

 write(fid,'(1X,A)') '$VEC' ! print MOs

 ! print Alpha MOs
 call prt_gms_mo(fid, nbf, nif, alpha_coeff)
 deallocate(alpha_coeff)

 ! print Beta MOs (if any)
 if(uhf) then
  call prt_gms_mo(fid, nbf, nif, beta_coeff)
  deallocate(beta_coeff)
 end if

 write(fid,'(1X,A)') '$END'
 close(fid)
end subroutine fch2inp

! create the GAMESS .inp file and print the keywords information
subroutine creat_gamess_inp_head(inpname, charge, mult, ncore, npair, nopen, &
               nif, nbf, na, nb, itype, irel, fcgvb, uhf, ecp, sph, X2C, DIIS)
 implicit none
 integer :: i, j, k, m, n, fid
 integer, intent(in) :: charge, mult, ncore, npair, nopen, nif, nbf, na, nb, &
  itype, irel
 character(len=2), allocatable :: ideg(:)
 character(len=43), parameter :: error_warn='ERROR in subroutine creat_gamess_i&
                                            &np_head: '
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: fcgvb, uhf, ecp, sph, X2C, DIIS

 open(newunit=fid,file=TRIM(inpname),status='replace')
 write(fid,'(A)',advance='no') ' $CONTRL SCFTYP='

 select case(itype)
 case(-3) ! GHF
  write(fid,'(A)',advance='no') 'GHF' ! unsupported in GAMESS, print anyway
 case(-2) ! ROHF
  write(fid,'(A)',advance='no') 'ROHF'
 case(-1) ! RHF
  write(fid,'(A)',advance='no') 'RHF'
 case(0)  ! UHF
  write(fid,'(A)',advance='no') 'UHF'
 case(1,2,3,4) ! SF/MRSF
  if(uhf) then
   write(fid,'(A)',advance='no') 'UHF'
  else
   write(fid,'(A)',advance='no') 'ROHF'
  end if
 case(5) ! GVB
  write(fid,'(A)',advance='no') 'GVB'
 case default
  write(6,'(/,A)') error_warn//'itype out of range.'
  write(6,'(A,I0)') 'itype=', itype
  close(fid)
  stop
 end select

 write(fid,'(2(A,I0))',advance='no') ' RUNTYP=ENERGY ICHARG=',charge,' MULT=',mult
 write(fid,'(A)',advance='no') ' NOSYM=1 ICUT=11'
 if(ecp) write(fid,'(A)',advance='no') ' PP=READ'

 if(itype == 5) then ! GVB
  write(fid,'(/,A)',advance='no') '  MAXIT=500'
 else
  write(fid,'(/,A)',advance='no') '  MAXIT=200'
 end if

 if(irel == -2) then
  write(fid,'(A)',advance='no') ' RELWFN=RESC'
 else if(irel /= -1) then
  write(fid,'(A)',advance='no') ' RELWFN=DK'
 end if

 if(sph) write(fid,'(A)',advance='no') ' ISPHER=1'

 select case(itype)
 case(1) ! SF-CIS
  write(fid,'(A)',advance='no') ' CITYP=SFCIS'
 case(2) ! SF-TDDFT
  write(fid,'(A)',advance='no') ' DFTTYP=BHHLYP TDDFT=SPNFLP'
 case(3) ! MRSF-CIS
  write(fid,'(A)',advance='no') ' DFTTYP=NONE TDDFT=MRSF'
 case(4) ! MRSF-TDDFT
  write(fid,'(A)',advance='no') ' DFTTYP=BHHLYP TDDFT=MRSF'
 end select

 if(X2C) then
  write(fid,'(A)') ' $END X2C'
 else
  write(fid,'(A)') ' $END'
 end if

 write(fid,'(A)') ' $SYSTEM MWORDS=300 $END'

 if(irel == 0) write(fid,'(A)') ' $RELWFN NORDER=1 $END' ! DKH 0-th

 if(itype == 5) then
  write(fid,'(2(A,I0))',advance='no') ' $SCF NCO=',ncore,' NPAIR=',npair
  if(nopen > 0) then
   write(fid,'(A,I0,A)',advance='no') ' NSETO=',nopen,' NO(1)=1'
   if(nopen > 1) then
    allocate(ideg(nopen-1))
    ideg = ',1'
    write(fid,'(15A2)',advance='no') ideg
    deallocate(ideg)
    if(nopen > 16) then
     write(fid,'(/,A)') error_warn//'nopen>16 cannot be dealt'
     write(fid,'(A)') 'with currently. Please contact MOKIT developers to updat&
                      &e code.'
     close(fid)
     stop
    end if
   end if
  end if
  if(mult < 4) then
   write(fid,'(A)',advance='no') ' DIRSCF=.T.'
   if(DIIS) write(fid,'(A)',advance='no') ' DIIS=.T. SOSCF=.F.'
  else ! mult >=4, i.e. >=3 singly occpuied orb
   write(fid,'(A)') ' DIRSCF=.T. COUPLE=.T.'
   if(DIIS) write(fid,'(2X,A)') 'DIIS=.T. SOSCF=.F.'
   call prt_gvb_couple_coeff(fid, ncore, nopen)
  end if
  write(fid,'(A)') ' $END'
 else
  write(fid,'(A)',advance='no') ' $SCF DIRSCF=.T.'
  if(irel==0 .or. irel==2 .or. irel==4) then
   write(fid,'(A)',advance='no') ' DIIS=.T. SOSCF=.F.'
  end if
  ! HDOMO: the highest doubly occupied MO
  ! LSOMO: the lowest singly occupied MO
  ! The ROHF/ROKS Fock operator is not unique. When transforming converged MOs
  ! from other programs to GAMESS, E(LSOMO) < E(HDOMO) happens somtimes. In such
  ! case GAMESS will exchange these two special MOs, and thus make SCF oscillated.
  ! To avoid this case, `MOM=.T.` is added. If there is no orbital energy
  ! inversion occurs, `MOM=.T.` will not take effect.
  if(itype>0 .and. itype<5) write(fid,'(A)',advance='no') ' MOM=.T.'
  write(fid,'(A)') ' FDIFF=.F. $END'
 end if

 ! freeze GVB doubly occupied MOs if the user requests that
 if(itype==5 .and. ncore>0 .and. fcgvb) then
  if(ncore <= 20) then        ! 0 < ncore <= 20
   write(fid,'(A)',advance='no') ' $MOFRZ FRZ=.T. IFRZ(1)=1'
   do i = 2, ncore, 1
    write(fid,'(A,I0)',advance='no') ',', i
   end do ! for i
   write(fid,'(A)') ' $END'
  else if(ncore <= 1000) then ! 20 < ncore <= 1000
   write(fid,'(A)') ' $MOFRZ FRZ=.T.'
   k = ncore/17
   do i = 1, k, 1
    m = (i-1)*17 + 1
    write(fid,'(2X,2(A,I0))',advance='no') 'IFRZ(',m,')=',m
    write(fid,'(16(A,I0))') (',',m+j,j=1,16)
   end do ! for i
   n = ncore - 17*k
   if(n > 0) then
    m = k*17 + 1
    write(fid,'(2X,2(A,I0))',advance='no') 'IFRZ(',m,')=',m
    write(fid,'(16(A,I0))') (',',m+j-1,j=2,n)
   end if
   write(fid,'(A)') ' $END'
  else                        ! 1000 < ncore
   write(6,'(/,A)') error_warn//'ncore>1000 is not supported by fch2inp currently.'
   write(6,'(A,I0)') 'ncore=', ncore
   close(fid)
   stop
  end if
 end if

 write(fid,'(A,I0,A)') ' $GUESS GUESS=MOREAD NORB=', nif, ' $END'
 select case(itype)
 case(1)
  write(fid,'(A)') ' $CIS NSTATE=5 $END'
 case(2,4)
  write(fid,'(A)') ' $DFT NRAD0=99 NLEB0=590 NRAD=99 NLEB=590 $END'
  write(fid,'(A)') ' $TDDFT NSTATE=5 NRAD=99 NLEB=590 $END'
 case(3)
  write(fid,'(A)') ' $TDDFT NSTATE=5 NRAD=99 NLEB=590 $END'
 end select

 write(fid,'(A)') ' $DATA'
 write(fid,'(A)',advance='no') 'GAMESS inp file produced by MOKIT'
 if(itype /= 5) write(fid,'(2(A,I0))',advance='no') ',na=',na,',nb=',nb
 write(fid,'(2(A,I0))') ',nif=',nif,',nbf=',nbf
 write(fid,'(A)') 'C1   1'
 close(fid)
end subroutine creat_gamess_inp_head

! print GAMESS format of MOs into a given file ID
subroutine prt_gms_mo(fid, nbf, nif, coeff)
 implicit none
 integer :: i, j, k, nline, nleft
 integer, intent(in) :: fid, nbf, nif
 real(kind=8), intent(in) :: coeff(nbf,nif)

 nline = nbf/5
 nleft = nbf - 5*nline

 do i = 1, nif, 1
  k = MOD(i,100)

  do j = 1, nline, 1
   write(fid,'(I2,I3,5ES15.8)') k, MOD(j,1000), coeff(5*j-4:5*j,i)
  end do ! for j

  if(nleft > 0) then
   write(fid,'(I2,I3,5ES15.8)') k, MOD(j,1000), coeff(5*j-4:nbf,i)
  end if
 end do ! for i
end subroutine prt_gms_mo

subroutine fch2inp_permute_10f(nif,coeff)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 integer, parameter :: order(6) = [2, 3, 1, 6, 4, 5]
 real(kind=8), intent(inout) :: coeff(6,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in Gamess
!                     1    2    3    4    5    6
! From: XXX,YYY,ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ,XYZ
! To:   XXX,YYY,ZZZ, XXY, XXZ, XYY, YYZ, XZZ, YZZ,XYZ

 allocate(coeff2(6,nif))
 forall(i = 1:6) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2inp_permute_10f

subroutine fch2inp_permute_15g(nif,coeff)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 integer, parameter :: order(15) = [15, 5, 1, 14, 13, 9, 4, 6, 2, 12, 10, 3, &
                                    11, 8, 7]
 real(kind=8), intent(inout) :: coeff(15,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian g functions in Gaussian
! To: the order of Cartesian g functions in Gamess
!       1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! From: ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
! To:   XXXX,YYYY,ZZZZ,XXXY,XXXZ,XYYY,YYYZ,XZZZ,YZZZ,XXYY,XXZZ,YYZZ,XXYZ,XYYZ,XYZZ

 allocate(coeff2(15,nif))
 forall(i = 1:15) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2inp_permute_15g

subroutine fch2inp_permute_21h(nif,coeff)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 integer, parameter :: order(21) = [21, 6, 1, 20, 19, 11, 5, 7, 2, 18, 16, 15, &
                                    4, 12, 3, 17, 10, 8, 14, 13, 9]
 real(kind=8),intent(inout) :: coeff(21,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian h functions in Gaussian
! To: the order of Cartesian h functions in Gamess
!       1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! From: ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX
! To:   XXXXX,YYYYY,ZZZZZ,XXXXY,XXXXZ,XYYYY,YYYYZ,XZZZZ,YZZZZ,XXXYY,XXXZZ,XXYYY,YYYZZ,XXZZZ,YYZZZ,XXXYZ,XYYYZ,XYZZZ,XXYYZ,XXYZZ,XYYZZ

 allocate(coeff2(21,nif))
 forall(i = 1:21) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2inp_permute_21h

