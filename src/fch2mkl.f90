! written by jxzou at 20200113: convert .fch(k) file (Gaussian) to .mkl file (Molekel, ORCA)
! updated by jxzou at 20200202: ECP/PP supported
! updated by jxzou at 20200304: Pople-type basis sets supported
! updated by jxzou at 20200322: move read_fch to read_fch.f90
! updated by jxzou at 20200622: fix the bug in F, G, H (some multiply by -1); add 1 more digit for MOs
! updated by jxzou at 20200802: add $CHARGES section to .mkl file
! updated by jxzou at 20201118: detect DKH/RESC keywords in .fch(k) file
! updated by jxzou at 20201221: add 'DelGTO' for elements Rb~Rn
! updated by jxzou at 20210407: remove '-uhf', use automatic determination
! updated by jxzou at 20210730: add 'sthresh 1e-6' in %scf, in case of linear dependence

! The 'Shell types' array in Gaussian .fch file:
!
!   Spherical     |     Cartesian
! -6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6
!  I  H  G  F  D  L  S  P  D  F  G  H  I
!
! 'L' is 'SP' in Pople-type basis sets

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i, k
 integer :: itype ! 0/1/2 for default/SF-TDDFT/DFT
 character(len=4) :: str4
 character(len=30) :: str30
 character(len=240) :: fchname

 itype = 0; str4 = ' '; str30 = ' '; fchname = ' '

 i = iargc()
 if(i<1 .or. i>3) then
  write(6,'(/,A)')' ERROR in program fch2mkl: wrong command line arguments!'
  write(6,'(A)')  ' Example 1 (R(O)HF/UHF/CAS): fch2mkl h2o.fch'
  write(6,'(A)')  ' Example 2 (SF-TDDFT)      : fch2mkl O2_T.fch -sf'
  write(6,'(A)')  " Example 3 (R(O)KS/UKS)    : fch2mkl h2o.fch -dft 'B3LYP D3BJ'"
  write(6,'(A)')  "                             fch2mkl h2o.fch -dft 'HSE06 D3zero'"
  write(6,'(A)')  "                             fch2mkl h2o.fch -dft 'wB97M-V'"
  write(6,'(A)')  " Example 4 (DLPNO-DH)      : fch2mkl h2o.fch -dft 'DLPNO-B2PLYP D3BJ'"
  write(6,'(A,/)')"                             fch2mkl h2o.fch -dft 'DLPNO-wB97X-2 D3'"
  stop
 end if

 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 k = LEN_TRIM(fchname)
 if(fchname(k-3:k) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:k-4)//'.fch'
 end if

 if(i > 1) then
  call getarg(2, str4)
  select case(TRIM(str4))
  case('-sf')
   itype = 1
  case('-dft')
   itype = 2
   if(i == 2) then
    write(6,'(/,A)') 'ERROR in program fch2mkl: you should specify the DFT name.'
    write(6,'(A)') "Example: fch2mkl h2o.fch -dft 'B3LYP D3BJ'"
    stop
   end if
   call getarg(3, str30)
  case default
   write(6,'(/,A)') 'ERROR in program fch2mkl: the 2nd command line argument is&
                    & wrong!'
   stop
  end select
 end if

 call fch2mkl(fchname, itype, str30)
end program main

! convert .fch(k) file (Gaussian) to .mkl file (Molekel, ORCA)
subroutine fch2mkl(fchname, itype, dftname)
 use fch_content
 implicit none
 integer :: i, j, k, m, n, n1, n2, am, fid1, fid2
 integer :: ndmark, nfmark, ngmark, nhmark, nimark
 integer, intent(in) :: itype
 integer, parameter :: ndh = 32
 integer, parameter :: list(10) = [2,3,4,5,6,7,8,9,10,1]
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:), i_mark(:)
 real(kind=8), allocatable :: coeff(:,:)
 character(len=1) :: str = ' '
 character(len=1), parameter :: am_type(0:6) = ['S','P','D','F','G','H','I']
 character(len=1), parameter :: am_type1(0:6) = ['s','p','d','f','g','h','i']
 character(len=16), parameter :: dhname(ndh) = ['B2PLYP          ', &
  'MPW2PLYP        ','B2GP-PLYP       ','B2K-PLYP        ','B2T-PLYP        ',&
  'PWPB95          ','PBE-QIDH        ','PBE0-DH         ','DSD-BLYP        ',&
  'DSD-PBEP86      ','DSD-PBEB95      ','WB2PLYP         ','WB2GP-PLYP      ',&
  'RSX-QIDH        ','RSX-0DH         ','WB88PP86        ','WPBEPP86        ',&
  'WB97X-2         ','SCS/SOS-B2PLYP21','SCS-PBE-QIDH    ','SOS-PBE-QIDH    ',&
  'SCS-B2GP-PLYP21 ','SOS-B2GP-PLYP21 ','SCS/SOS-WB2PLYP ','SCS-WB2GP-PLYP  ',&
  'SOS-WB2GP-PLYP  ','SCS-RSX-QIDH    ','SOS-RSX-QIDH    ','SCS-WB88PP86    ',&
  'SOS-WB88PP86    ','SCS-WPBEPP86    ','SOS-WPBEPP86    ']
 character(len=30) :: dftname1
 character(len=30), intent(in) :: dftname
 character(len=240) :: mklname, inpname
 character(len=240), intent(in) :: fchname
 logical :: uhf, ecp, composite, hse06, m052x, b972, d3zero, d3bj

 ecp = .false.; composite = .false.; hse06 = .false.; m052x = .false.
 b972 = .false.; d3zero = .false.; d3bj = .false.

 if(itype == 2) then
  dftname1 = ADJUSTL(dftname)
  call upper(dftname1)
  i = LEN_TRIM(dftname1)
  if(i > 3) then ! B97-3c, r2SCAN-3c, HF-3c, etc.
   if(dftname1(i-2:i) == '-3C') composite = .true.
   if(dftname1(1:4) == 'B972') b972 = .true.
  end if
  if(i > 4) then
   if(dftname1(i-3:i) == 'D3BJ') d3bj = .true.
   if(dftname1(1:5) == 'HSE06') hse06 = .true.
   if(dftname1(1:5) == 'M052X') m052x = .true.
  end if
  if(i > 6) then
   if(dftname1(i-5:i) == 'D3ZERO') d3zero = .true.
  end if
  if(dftname1(1:5)=='B3LYP' .and. dftname1(6:7)/='/G') then
   write(6,'(/,A)') 'Remark: it seems that B3LYP is used. Note that B3LYP in Ga&
                    &ussian is equi-'
   write(6,'(A)') 'valent to B3LYP/G in ORCA. If you want a close/closer energy&
                  & result between'
   write(6,'(A)') 'Gaussian/ORCA B3LYP, you need to specify B3LYP/G when using &
                  &fch2mkl.'
  end if
 end if

 call find_specified_suffix(fchname, '.fch', i)
 mklname = fchname(1:i-1)//'_o.mkl'
 inpname = fchname(1:i-1)//'_o.inp'

 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)
 call find_irel_in_fch(fchname, irel)
 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 call read_fch(fchname, uhf) ! read content in .fch(k) file
 ecp = .false.
 if(LenNCZ > 0) ecp = .true.

 ! check if any Cartesian functions
 if( ANY(shell_type > 1) ) then
  write(6,'(/,A)') 'ERROR in subroutine fch2mkl: Cartesian functions detected i&
                   &n file '//TRIM(fchname)//'.'
  write(6,'(A)') "ORCA supports only spherical harmonic functions. You need to &
                 &add '5D 7F' keywords"
  write(6,'(A)') 'in Gaussian.'
  stop
 else if( ANY(shell_type < -6) ) then
  write(6,'(/,A)') 'ERROR in subroutine fch2mkl: angular momentum too high! not&
                   & supported for fch2mkl.'
  stop
 end if

 if(itype==1) then
  if(mult /= 3) then
   write(6,'(/,A)') 'ERROR in subroutine fch2mkl: SF-TDDFT in ORCA can only be &
                    &based on a triplet'
   write(6,'(A,I0)') 'UHF/UKS reference. The spin multiplicity in your provided&
                    & .fch(k) file is ', mult
   stop
  end if
  if(.not. uhf) then
   write(6,'(/,A)') 'ERROR in subroutine fch2mkl: SF-TDDFT in ORCA can only be &
                    &based on a triplet'
   write(6,'(A,I0)') 'UHF/UKS reference. It seems that you provide an ROHF/ROKS&
                     & .fch(k) file.'
   stop
  end if
 end if

 ! update MO coefficients
 if(uhf) then ! UHF
  k = 2*nif
  allocate(coeff(nbf,k))
  coeff(:,1:nif) = alpha_coeff
  coeff(:,nif+1:) = beta_coeff
 else         ! R(O)HF
  k = nif
  allocate(coeff(nbf,k), source=alpha_coeff)
 end if

 ! find F+3, G+3 and H+3 functions, multiply them by -1
 allocate(d_mark(ncontr), f_mark(ncontr), g_mark(ncontr), h_mark(ncontr), &
          i_mark(ncontr))
 call read_mark_from_shltyp_sph(ncontr, shell_type, ndmark, nfmark, ngmark, &
                      nhmark, nimark, d_mark, f_mark, g_mark, h_mark, i_mark)
 deallocate(d_mark)
 call update_mo_using_mark_orca(nbf, k, nfmark, ngmark, nhmark, nimark, ncontr,&
                                f_mark, g_mark, h_mark, i_mark, coeff)
 deallocate(f_mark, g_mark, h_mark, i_mark)

 if(uhf) then ! UHF
  alpha_coeff = coeff(:,1:nif)
  beta_coeff = coeff(:,nif+1:)
 else         ! R(O)HF
  alpha_coeff = coeff
 end if
 deallocate(coeff)
 ! update MO coefficients done

 ! print elements and coordinates into .mkl file
 open(newunit=fid1,file=TRIM(mklname),status='replace')
 write(fid1,'(A)') '$MKL'
 write(fid1,'(A)') '#'
 write(fid1,'(A)') '# MKL format file produced by MOKIT'
 write(fid1,'(A)') '#'
 write(fid1,'(A)') '$CHAR_MULT'
 write(fid1,'(I0,1X,I0)') charge, mult
 write(fid1,'(A,/)') '$END'

 write(fid1,'(A)') '$COORD'
 do i = 1, natom, 1
  write(fid1,'(I3,1X,3F17.9)') ielem(i), (coor(j,i), j=1,3)
 end do ! for i
 write(fid1,'(A,/)') '$END'

 write(fid1,'(A)') '$CHARGES'
 do i = 1, natom, 1
  write(fid1,'(1X,A3)') '0.0'
 end do ! for i
 write(fid1,'(A,/)') '$END'

 ! print basis sets into .mkl file (Note: mkl file contains no ECP/PP data)
 write(fid1,'(A)') '$BASIS'
 k = 0
 do i = 1, ncontr, 1
  m = shell2atom_map(i)
  if(m > 1) then
   if(shell2atom_map(i-1) == m-1) write(fid1,'(A2)') '$$'
  end if

  m = shell_type(i); n = prim_per_shell(i)

  if(m /= -1) then ! m<-1 or m=0,1

   m = IABS(m)
   write(fid1,'(I2,1X,A)') 2*m+1, am_type(m)//' 1.0'
   do j = k+1, k+n, 1
    write(fid1,'(2(2X,ES15.8))') prim_exp(j), contr_coeff(j)
   end do ! for j

  else ! m = -1, 'L' or 'SP' in Pople-type basis sets

   write(fid1,'(I2,1X,A)') 1, am_type(0)//' 1.0'
   do j = k+1, k+n, 1
    write(fid1,'(2(2X,ES15.8))') prim_exp(j), contr_coeff(j)
   end do ! for j

   write(fid1,'(I2,1X,A)') 3, am_type(1)//' 1.0'
   do j = k+1, k+n, 1
    write(fid1,'(2(2X,ES15.8))') prim_exp(j), contr_coeff_sp(j)
   end do ! for j

  end if

  k = k + n
 end do ! for i

 write(fid1,'(/,A,/)') '$END'

 ! print coordinates, basis sets and ECP/PP data into .inp file
 open(newunit=fid2,file=TRIM(inpname),status='replace')
 write(fid2,'(A)') '%pal nprocs 4 end'
 write(fid2,'(A)') '%maxcore 1000'

 if(itype == 0) then ! HF
  if(uhf) then
   write(fid2,'(A)',advance='no') '! UHF'
  else
   if(nopen == 0) then
    write(fid2,'(A)',advance='no') '! RHF'
   else ! nopen > 0
    write(fid2,'(A)',advance='no') '! ROHF'
   end if
  end if
  write(fid2,'(A)') ' VeryTightSCF noTRAH'
 else                ! KS-DFT
  if(uhf) then
   write(fid2,'(A)',advance='no') '! UKS'
  else
   if(nopen == 0) then
    write(fid2,'(A)',advance='no') '! RKS'
   else ! nopen > 0
    write(fid2,'(A)',advance='no') '! ROKS'
   end if
  end if
  select case(itype)
  case(1)
   write(fid2,'(A)',advance='no') ' BHANDHLYP'
  case(2)
   if(hse06 .or. b972) then
    if(d3bj) then
     write(fid2,'(A)',advance='no') ' D3BJ'
    else if(d3zero) then
     write(fid2,'(A)',advance='no') ' D3ZERO'
    end if
   else if(m052x) then
    if(d3bj) then
     write(6,'(/,A)') 'ERROR in subroutine fch2mkl: M052X cannot be combined wi&
                      &th D3BJ.'
     close(fid1)
     close(fid2)
     stop
    else if(d3zero) then
     write(fid2,'(A)',advance='no') ' D3ZERO'
    end if
   else
    write(fid2,'(A)',advance='no') ' '//TRIM(dftname)
   end if
   do i = 1, ndh, 1
    if(INDEX(TRIM(dftname1), TRIM(dhname(i))) > 0) then
     write(fid2,'(A)',advance='no') ' def2-TZVPP/C'
     exit
    end if
   end do ! for i
  end select

  if(INDEX(TRIM(dftname1), 'DLPNO') > 0) then
   write(fid2,'(A)',advance='no') ' TightPNO'
  end if
  if(composite) then
   write(fid2,'(A)') ' defgrid3 TightSCF noTRAH'
  else
   write(fid2,'(A)') ' def2/J RIJCOSX defgrid3 TightSCF noTRAH'
  end if

  if(hse06 .or. m052x .or. b972) then
   write(fid2,'(A)') '%method'
   write(fid2,'(A)') ' method dft'
   if(hse06) then
    write(fid2,'(A)') ' functional hyb_gga_xc_hse06'
    if(d3bj) then
     write(fid2,'(A)') ' D3S6 1.0'
     write(fid2,'(A)') ' D3A1 0.383'
     write(fid2,'(A)') ' D3S8 2.31'
     write(fid2,'(A)') ' D3A2 5.685'
    else if(d3zero) then
     write(fid2,'(A)') ' D3S6 1.0'
     write(fid2,'(A)') ' D3RS6 1.129'
     write(fid2,'(A)') ' D3S8 0.109'
     write(fid2,'(A)') ' D3alpha6 14'
    end if
   end if
   if(m052x) then
    write(fid2,'(A)') ' exchange hyb_mgga_x_m05_2x'
    write(fid2,'(A)') ' correlation mgga_c_m05_2x'
    if(d3zero) then
     write(fid2,'(A)') ' D3S6 1.0'
     write(fid2,'(A)') ' D3RS6 1.417'
     write(fid2,'(A)') ' D3S8 0.0'
     write(fid2,'(A)') ' D3alpha6 14'
    end if
   end if
   if(b972) then
    write(fid2,'(A)') ' functional hyb_gga_xc_b97_2'
    if(d3bj) then
     write(fid2,'(A)') ' D3S6 1.0'
     write(fid2,'(A)') ' D3A1 0.0'
     write(fid2,'(A)') ' D3S8 0.9448'
     write(fid2,'(A)') ' D3A2 5.4603'
    else if(d3zero) then
     write(fid2,'(A)') ' D3S6 1.0'
     write(fid2,'(A)') ' D3RS6 1.7066'
     write(fid2,'(A)') ' D3S8 2.4661'
     write(fid2,'(A)') ' D3alpha6 14'
    end if
   end if
   write(fid2,'(A)') 'end'
  end if
 end if

 select case(irel)
 case(-3) ! X2C, i.e. sf-x2c1e
  write(6,'(/,A)') 'Remark from subroutine fch2mkl: X2C detected. X2C is suppor&
                   &ted since ORCA 6,'
  write(6,'(A)') 'please make sure that your ORCA version is appropriate.'
 case(-2) ! RESC
  write(6,'(/,A)') 'ERROR in subroutine fch2mkl: RESC keyword detected in file &
                   &'//TRIM(fchname)//'.'
  write(6,'(A)') 'But ORCA does not support RESC.'
  stop
 case(-1,2) ! none/DKH2
 case(0) ! DKH0
  write(6,'(/,A)') 'Warning in subroutine fch2mkl: DKH0 detected in file '//&
                    TRIM(fchname)
  write(6,'(A)') 'But ORCA does not support DKH0. DKH2 keyword will be printed &
                 &into ORCA input file.'
 case(4) ! DKHSO, DKH4 with SO
  write(6,'(/,A)') 'Warning in subroutine fch2mkl: DKHSO detected in file '//&
                    TRIM(fchname)//'.'
  write(6,'(A)') 'But ORCA does not support DKHSO. DKH2 keyword will be printed&
                 & into ORCA input file.'
 case default
  write(6,'(/,A)') 'ERROR in subroutine fch2mkl: irel out of range!'
  write(6,'(A)') 'This type of Hamiltonian is not supported.'
  write(6,'(A,I0)') 'irel=', irel
  stop
 end select

 if(irel > 0) then
  write(fid2,'(A,/,A)') '%rel',' method DKH'
  write(fid2,'(A,I0,/,A)') ' order ', irel, 'end'
 else if(irel == -3) then
  write(fid2,'(A,/,A,/,A)') '%rel',' method X2C','end'
 end if

 if(itype == 1) then
  write(fid2,'(A)') '%tddft'
  write(fid2,'(A)') ' sf true'
  write(fid2,'(A)') ' nroots 5'
  write(fid2,'(A)') 'end'
 end if

 write(fid2,'(A)') '%scf'
 write(fid2,'(A)') ' Thresh 1e-12'
 write(fid2,'(A)') ' Tcut 1e-14'
 write(fid2,'(A)') ' CNVDamp False'
 ! ORCA default is `sthresh 1e-7`
 if(nif < nbf) write(fid2,'(A)') ' sthresh 1e-6'
 write(fid2,'(A)') 'end'

 do i = 1, natom, 1
  if(elem(i) == 'Bq') elem(i) = 'X '
 end do ! for i

 if(composite) then ! composite method
  write(fid2,'(A,I0,1X,I0)') '* xyz ', charge, mult
  do i = 1, natom, 1
   if(iatom_type(i) == 1000) str = ':'
   write(fid2,'(A,3(1X,F17.9))') TRIM(elem(i))//str, coor(:,i)
  end do ! for i
  write(fid2,'(A)') '*'
  if(ecp) deallocate(KFirst, KLast, Lmax, LPSkip, NLP, RNFroz, CLP, ZLP)
 else               ! not composite method
  write(fid2,'(A)') '%coords'
  write(fid2,'(A)') ' Units = angs'
  write(fid2,'(A,I0)') ' Charge = ', charge
  write(fid2,'(A,I0)') ' Mult = ', mult
  write(fid2,'(A)') ' Coords'
  k = 0

  do i = 1, ncontr, 1
   m = shell2atom_map(i)
 
   if(m == 1) then
    if(i == 1) then
     if(iatom_type(1) == 1000) str = ':'
     write(fid2,'(1X,A,3(1X,F17.9))') TRIM(elem(1))//str, coor(:,1)
     if(ielem(1)>36 .and. ielem(1)<87) write(fid2,'(2X,A)') 'DelECP'
     write(fid2,'(2X,A)') 'NewGTO'
    end if
 
   else ! m > 1
    if(shell2atom_map(i-1) == m-1) then
     write(fid2,'(2X,A)') 'end'   ! print GTO end of last atom
 
     if(ecp) then   ! print ECP/PP data of last atom
      if(LPSkip(m-1) == 0) then
       write(fid2,'(2X,A)') 'NewECP'
       write(fid2,'(3X,A,1X,I3)') 'N_core', NINT(RNFroz(m-1))
       write(fid2,'(3X,A)') 'lmax '//am_type1(LMax(m-1))
       am = 0
       do j = 1, 10, 1
        n1 = KFirst(m-1,list(j)); n2 = KLast(m-1,list(j))
        if(n1 == 0) cycle
        am = am + 1
        write(fid2,'(3X,A1,1X,I1)') am_type1(am-1), n2-n1+1
        do n = n1, n2, 1
         write(fid2,'(3X,I2,2(1X,ES15.8),1X,I1)') n-n1+1, ZLP(n), CLP(n), NLP(n)
        end do ! for n
       end do ! for j
 
       write(fid2,'(2X,A)') 'end'  ! in accord with 'NewECP'
      end if
     end if         ! print ECP/PP data done
 
     ! print coordinates of the current atom
     str = ' '
     if(iatom_type(m) == 1000) str = ':'
     write(fid2,'(1X,A,3(1X,F17.9))') TRIM(elem(m))//str, coor(:,m)
     if(ielem(m)>36 .and. ielem(m)<87) write(fid2,'(2X,A)') 'DelECP'
     write(fid2,'(2X,A)') 'NewGTO'
    end if
   end if
 
   m = shell_type(i); n = prim_per_shell(i)
 
   if(m /= -1) then
    m = IABS(m)
    write(fid2,'(4X,A1,1X,I3)') am_type(m), n
    do j = k+1, k+n, 1
     write(fid2,'(2X,I3,2(2X,ES15.8))') j-k,prim_exp(j), contr_coeff(j)
    end do ! for j
 
   else ! m = -1, 'L' or 'SP'
    write(fid2,'(4X,A1,1X,I3)') 'S', n
    do j = k+1, k+n, 1
     write(fid2,'(2X,I3,3(2X,ES15.8))') j-k, prim_exp(j), contr_coeff(j)
    end do ! for j
    write(fid2,'(4X,A1,1X,I3)') 'P', n
    do j = k+1, k+n, 1
     write(fid2,'(2X,I3,3(2X,ES15.8))') j-k, prim_exp(j), contr_coeff_sp(j)
    end do ! for j
   end if
 
   k = k + n
  end do ! for i
 
  write(fid2,'(2X,A)') 'end'  ! in accord with 'NewGTO'
  k = shell2atom_map(ncontr)
  if(k < natom) then
   do i = k+1, natom, 1
    write(fid2,'(1X,A2,3(1X,F17.9))') 'H:', coor(:,i)
    write(fid2,'(2X,A)') 'NewGTO S 1 1 1e6 1 end'
   end do ! for i
  end if
 
  if(ecp) then   ! print ECP/PP data of the last atom
   if(LPSkip(natom) == 0) then
    write(fid2,'(2X,A)') 'NewECP'
    write(fid2,'(3X,A,1X,I3)') 'N_core', NINT(RNFroz(natom))
    write(fid2,'(3X,A)') 'lmax '//am_type1(LMax(natom))
    am = 0
    do j = 1, 10, 1
     n1 = KFirst(natom,list(j)); n2 = KLast(natom,list(j))
     if(n1 == 0) cycle
     am = am + 1
     write(fid2,'(3X,A1,1X,I1)') am_type1(am-1), n2-n1+1
     do n = n1, n2, 1
      write(fid2,'(3X,I2,2(1X,ES15.8),1X,I1)') n-n1+1, ZLP(n), CLP(n), NLP(n)
     end do ! for n
    end do ! for j
 
    write(fid2,'(2X,A)') 'end'  ! in accord with 'NewECP'
   end if
 
   deallocate(KFirst, KLast, Lmax, LPSkip, NLP, RNFroz, CLP, ZLP)
  end if         ! print ECP/PP data done

  write(fid2,'(1X,A)') 'end'  ! in accord with ' Coords'
  write(fid2,'(A,/)') 'end'   ! in accord with '%coords'
 end if

 close(fid2)
 deallocate(ielem, elem, coor, shell_type, prim_per_shell, shell2atom_map, &
            prim_exp, contr_coeff)
 if(allocated(contr_coeff_sp)) deallocate(contr_coeff_sp)

 ! print Alpha MO and corresponding energies into .mkl file
 write(fid1,'(A)') '$COEFF_ALPHA'
 call prt_mo_and_e_in_mkl(fid1, nbf, nif, alpha_coeff, eigen_e_a)
 write(fid1,'(A,/)') '$END'
 deallocate(alpha_coeff, eigen_e_a)

 ! print Alpha orbital occupation numbers into .mkl file
 write(fid1,'(A)') '$OCC_ALPHA'
 allocate(eigen_e_a(nif), source=0d0)
 if(uhf) then
  forall(i = 1:na) eigen_e_a(i) = 1d0
 else ! .not. uhf
  if(nopen == 0) then
   forall(i = 1:na) eigen_e_a(i) = 2d0
  else ! nopen > 0
   forall(i = 1:nb) eigen_e_a(i) = 2d0
   forall(i = nb+1:na) eigen_e_a(i) = 1d0
  end if
 end if
 write(fid1,'(5(F12.7,1X))') (eigen_e_a(i), i=1,nif)
 deallocate(eigen_e_a)
 write(fid1,'(A,/)') '$END'

 if(uhf) then ! print Beta MOs and corresponding energies
  write(fid1,'(A)') '$COEFF_BETA'
  call prt_mo_and_e_in_mkl(fid1, nbf, nif, beta_coeff, eigen_e_b)
  write(fid1,'(A,/)') '$END'
  deallocate(beta_coeff, eigen_e_b)

  ! print Beta orbital occupation numbers (if any) into .mkl file
  write(fid1,'(A)') '$OCC_BETA'
  allocate(eigen_e_b(nif), source=0d0)
  forall(i = 1:nb) eigen_e_b(i) = 1d0
  write(fid1,'(5(F12.7,1X))') (eigen_e_b(i), i=1,nif)
  deallocate(eigen_e_b)
  write(fid1,'(A,/)') '$END'
 end if

 close(fid1)
end subroutine fch2mkl

