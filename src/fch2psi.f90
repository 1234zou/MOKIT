! written by jxzou at 20210210: transfer MOs from Gaussian to Psi4
! updated by jxzou at 20210215: supported ECP/PP

! This utility/subroutine will generate PSI4 input file from Gaussian .fch(k)
! file (with coordinates, basis sets and MOs written in)

module root_param_fch2psi
 implicit none
 real(kind=8), parameter :: root3 = DSQRT(3d0)
 real(kind=8), parameter :: root5 = DSQRT(5d0)
 real(kind=8), parameter :: root7 = DSQRT(7d0)
 real(kind=8), parameter :: root15 = DSQRT(15d0)
 real(kind=8), parameter :: root21 = DSQRT(21d0)
 real(kind=8), parameter :: root35 = DSQRT(35d0)
 real(kind=8), parameter :: root35_3 = DSQRT(35d0/3d0)
 real(kind=8), parameter :: root63 = 3d0*DSQRT(7d0)
 real(kind=8), parameter :: root105 = DSQRT(105d0)
end module root_param_fch2psi

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i
 character(len=4) :: str4
 character(len=15) :: dftname
 character(len=240) :: fchname

 i = iargc()
 if(.not. (i==1 .or. i==3)) then
  write(6,'(/,A)') ' ERROR in subroutine fch2psi: wrong command line arguments!'
  write(6,'(A)')   ' Example (R(O)HF/UHF/CAS): fch2psi a.fch'
  write(6,'(A,/)') " Example            (DFT): fch2psi a.fch -dft 'wB97M-D3BJ'"
  stop
 end if

 str4 = ' '; dftname = ' '; fchname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 if(i == 3) then
  call getarg(2, str4)
  if(str4 /= '-dft') then
   write(6,'(/,A)') "ERROR in subroutine fch2psi: the 2nd argument can only be &
                    &'-dft'."
   stop
  end if
  call getarg(3, dftname)
  if(TRIM(dftname) == 'scf') then
   write(6,'(/,A)') "ERROR in subroutine fch2psi: the 3rd argument cannot be 's&
                    &cf'. Please"
   write(6,'(A)') 'specify a proper density functional name.'
   stop
  end if
 end if

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 call fch2psi(fchname, dftname)
end program main

! transfer MOs from Gaussian to Psi4
subroutine fch2psi(fchname, dftname)
 use util_wrapper, only: fch2inp_wrap
 implicit none
 integer :: i, nbf0, nbf, nif, ncontr, icart
 integer :: n3pmark, n6dmark, n10fmark, n15gmark, n21hmark
 integer, allocatable :: shell_type(:), shl2atm(:)
 integer, allocatable :: p_mark(:), d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 real(kind=8), allocatable :: coeff(:,:)
 character(len=15), intent(in) :: dftname
 character(len=240) :: gms_inp, inpname, fileA, fileB
 character(len=240), intent(in) :: fchname
 logical :: sph, uhf

 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)
 call check_uhf_in_fch(fchname, uhf)
 call fch2inp_wrap(fchname, .false., 0, 0, .false., .false.)
 sph = .true.
 call find_icart_in_fch(fchname, .false., icart)
 if(icart == 2) sph = .false.

 call find_specified_suffix(fchname, '.fch', i)
 gms_inp = fchname(1:i-1)//'.inp'
 inpname = fchname(1:i-1)//'_psi.inp'
 fileA = fchname(1:i-1)//'.A'
 fileB = fchname(1:i-1)//'.B'

 call bas_gms2psi(gms_inp, dftname, sph)
 call delete_file(gms_inp)

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf

 ! read MO Coefficients
 if(uhf) then
  allocate(coeff(nbf,2*nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_fch(fchname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 else
  allocate(coeff(nbf,nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)
 end if

 call read_ncontr_from_fch(fchname, ncontr)
 allocate(shell_type(ncontr), source=0)
 allocate(shl2atm(ncontr), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, ncontr, shell_type, shl2atm)
 deallocate(shl2atm)

 n3pmark = 0; n6dmark = 0; n10fmark = 0; n15gmark = 0; n21hmark = 0
 allocate(p_mark(ncontr), source=0)
 allocate(d_mark(ncontr), source=0)
 allocate(f_mark(ncontr), source=0)
 allocate(g_mark(ncontr), source=0)
 allocate(h_mark(ncontr), source=0)
 nbf = 0
 do i = 1, ncontr, 1
  select case(shell_type(i))
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   n3pmark = n3pmark + 1
   p_mark(n3pmark) = nbf + 1
   nbf = nbf + 3
  case(-1)   ! SP or L
   n3pmark = n3pmark + 1
   p_mark(n3pmark) = nbf + 2
   nbf = nbf + 4
  case(-2)   ! 5D
   nbf = nbf + 5
  case( 2)   ! 6D
   n6dmark = n6dmark + 1
   d_mark(n6dmark) = nbf + 1
   nbf = nbf + 6
  case(-3)   ! 7F
   nbf = nbf + 7
  case( 3)   ! 10F
   n10fmark = n10fmark + 1
   f_mark(n10fmark) = nbf + 1
   nbf = nbf + 10
  case(-4)   ! 9G
   nbf = nbf + 9
  case( 4)   ! 15G
   n15gmark = n15gmark + 1
   g_mark(n15gmark) = nbf + 1
   nbf = nbf + 15
  case(-5)   ! 11H
   nbf = nbf + 11
  case( 5)   ! 21H
   n21hmark = n21hmark + 1
   h_mark(n21hmark) = nbf + 1
   nbf = nbf + 21
  case default
   write(6,'(/,A)') 'ERROR in subroutine fch2psi: angular momentum too high!'
   write(6,'(3(A,I0))') 'ncontr=', ncontr, ', i=', i, ', nbf=', nbf
   stop
  end select
 end do ! for i
 deallocate(shell_type)

 if(nbf0 /= nbf) then
  write(6,'(/,A)') 'ERROR in subroutine fch2psi: inconsistent nbf.'
  write(6,'(2(A,I0))') 'nbf0=', nbf0, ', nbf=', nbf
  stop
 end if

 ! adjust MO coefficients according to the order of AOs (p,)
 if(sph) then
  do i = 1, n3pmark, 1
   call fch2psi_permute_3p(nif, coeff(p_mark(i):p_mark(i)+2,:))
  end do ! for i
 end if
 do i = 1, n6dmark, 1
  call fch2psi_permute_6d(nif, coeff(d_mark(i):d_mark(i)+5,:))
 end do ! for i
 do i = 1, n10fmark, 1
  call fch2psi_permute_10f(nif, coeff(f_mark(i):f_mark(i)+9,:))
 end do ! for i
 do i = 1, n15gmark, 1
  call fch2psi_permute_15g(nif, coeff(g_mark(i):g_mark(i)+14,:))
 end do ! for i
 do i = 1, n21hmark, 1
  call fch2psi_permute_21h(nif, coeff(h_mark(i):h_mark(i)+20,:))
 end do ! for i
 deallocate(p_mark, d_mark, f_mark, g_mark, h_mark)

 if(uhf) nif = nif/2
 call write_mo_into_psi_mat(fileA, nbf, nif, coeff(:,1:nif))
 if(uhf) call write_mo_into_psi_mat(fileB, nbf, nif, coeff(:,nif+1:2*nif))
 deallocate(coeff)

 if(nbf > 500) then
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Warning: more than 500 basis functions. If you find PSI4 comp&
                 &utation is slow,'
  write(6,'(A)') 'you can modify'
  write(6,'(A)') ' scf_type pk => scf_type df'
  write(6,'(A)') '                df_basis_scf def2-universal-jkfit'
  write(6,'(A)') 'in the input file.'
  write(6,'(A)') REPEAT('-',79)
 end if
end subroutine fch2psi

! this subroutine is used only when spherical harmonic functions are used
subroutine fch2psi_permute_3p(nif, coeff)
 implicit none
 integer :: i
 integer, parameter :: order(3) = [3,1,2]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(3,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical p functions in Gaussian
! To: the order of spherical p functions in PSI4
! 1    2    3
! +1, -1,  0
! 0 , +1, -1

 allocate(coeff2(3,nif), source=0d0)
 forall(i = 1:3) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2psi_permute_3p

! The basis order of 6D, 10F, 15G, 21H of PSI4 is the same as that of (OpenMolcas).
! But their normalization factors are not the same.
subroutine fch2psi_permute_6d(nif, coeff)
 use root_param_fch2psi, only: root3
 implicit none
 integer :: i, j
 integer, parameter :: order(6) = [1, 4, 5, 2, 6, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(6,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian d functions in Gaussian
! To: the order of Cartesian d functions in PSI4
! 1  2  3  4  5  6
! XX,YY,ZZ,XY,XZ,YZ
! XX,XY,XZ,YY,YZ,ZZ

 forall(j=4:6, i=1:nif) coeff(j,i) = coeff(j,i)*root3

 allocate(coeff2(6,nif), source=0d0)
 forall(i = 1:6) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2psi_permute_6d

subroutine fch2psi_permute_10f(nif, coeff)
 use root_param_fch2psi, only: root5, root15
 implicit none
 integer :: i, j
 integer, parameter :: order(10) = [1, 5, 6, 4, 10, 7, 2, 9, 8, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(10,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in PSI4
! 1   2   3   4   5   6   7   8   9   10
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ

 forall(j=4:9, i=1:nif) coeff(j,i) = coeff(j,i)*root5
 coeff(10,:) = coeff(10,:)*root15

 allocate(coeff2(10,nif), source=0d0)
 forall(i = 1:10) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2psi_permute_10f

subroutine fch2psi_permute_15g(nif, coeff)
 use root_param_fch2psi, only: root7, root35, root35_3
 implicit none
 integer :: i, j
 integer, intent(in) :: nif
 real(kind=8), parameter :: ratio(15) = [1d0,root7,root35_3,root7,1d0,root7,&
  root35,root35,root7,root35_3,root35,root35_3,root7,root7,1d0]
 real(kind=8), intent(inout) :: coeff(15,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian g functions in Gaussian
! To: the order of Cartesian g functions in PSI4
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
! xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz

 forall(j=1:15, i=1:nif) coeff(j,i) = coeff(j,i)*ratio(j)

 allocate(coeff2(15,nif), source=0d0)
 forall(i = 1:15) coeff2(i,:) = coeff(16-i,:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2psi_permute_15g

subroutine fch2psi_permute_21h(nif, coeff)
 use root_param_fch2psi, only: root21, root63, root105
 implicit none
 integer :: i, j
 integer, intent(in) :: nif
 real(kind=8), parameter :: ratio(21) = [1d0,3d0,root21,root21,3d0,1d0,3d0,&
  root63,root105,root63,3d0,root21,root105,root105,root21,root21,root63,root21,&
  3d0,3d0,1d0]
 real(kind=8), intent(inout) :: coeff(21,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian h functions in Gaussian
! To: the order of Cartesian h functions in PSI4
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX
! xxxxx,xxxxy,xxxxz,xxxyy,xxxyz,xxxzz,xxyyy,xxyyz,xxyzz,xxzzz,xyyyy,xyyyz,xyyzz,xyzzz,xzzzz,yyyyy,yyyyz,yyyzz,yyzzz,yzzzz,zzzzz

 forall(j=1:21, i=1:nif) coeff(j,i) = coeff(j,i)*ratio(j)

 allocate(coeff2(21,nif), source=0d0)
 forall(i = 1:21) coeff2(i,:) = coeff(22-i,:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2psi_permute_21h

! convert a GAMESS .inp file into a PSI4 input file
subroutine bas_gms2psi(inpname, dftname, sph)
 use pg, only: natom, nuc, elem, coor, ntimes, all_ecp, ecp_exist
 implicit none
 integer :: i, j, k, m, n, nline, rel, charge, mult, isph, fid1, fid2
 real(kind=8) :: rtmp(3)
 character(len=1), parameter :: am(0:6) = ['s','p','d','f','g','h','i']
 character(len=2) :: stype
 character(len=15) :: dftname1
 character(len=15), intent(in) :: dftname
 character(len=240) :: buf, inpname1, fileA, fileB
 character(len=240), intent(in) :: inpname
 logical :: uhf, ghf, X2C
 logical, intent(in) :: sph
 logical, allocatable :: ghost(:)

 call find_specified_suffix(inpname, '.inp', i)
 inpname1 = inpname(1:i-1)//'_psi.inp'
 fileA = inpname(1:i-1)//'.A'
 fileB = inpname(1:i-1)//'.B'

 call check_X2C_in_gms_inp(inpname, X2C)
 call read_charge_mult_isph_from_gms_inp(inpname, charge, mult, isph, uhf, ghf,&
                                         ecp_exist)
 call read_natom_from_gms_inp(inpname, natom)
 allocate(elem(natom), coor(3,natom), ntimes(natom), nuc(natom), ghost(natom))
 call read_elem_nuc_coor_from_gms_inp(inpname, natom, elem, nuc, coor, ghost)
 deallocate(nuc, ghost)
 call calc_ntimes(natom, elem, ntimes)

 open(newunit=fid2,file=TRIM(inpname1),status='replace')
 write(fid2,'(A)') '# generated by fch2psi of MOKIT'
 write(fid2,'(A)') 'memory 4 GB'
 write(fid2,'(/,A)') 'molecule mymol {'
 write(fid2,'(A)') 'symmetry C1'
 write(fid2,'(A)') 'no_reorient'
 write(fid2,'(I0,1X,I0)') charge, mult

 do i = 1, natom, 1
  write(fid2,'(A,I0,3(2X,F15.8))') TRIM(elem(i)), ntimes(i), coor(1:3,i)
 end do ! for i

 deallocate(coor)
 write(fid2,'(A)') '}'
 write(fid2,'(/,A)') 'basis mybas {'

 do i = 1, natom, 1
  write(fid2,'(2(A,I0))') ' assign '//TRIM(elem(i)),ntimes(i),' gen',i
 end do ! for i

 call read_all_ecp_from_gms_inp(inpname)
 call goto_data_section_in_gms_inp(inpname, fid1)
 read(fid1,'(A)') buf
 read(fid1,'(A)') buf

 do i = 1, natom, 1
  read(fid1,'(A)') buf ! buf contains elem(i), nuc(i) and coor(1:3,i)
  write(fid2,'(A,I0,A)') '[gen', i, ']'
  if(sph) then
   write(fid2,'(A)') 'spherical'
  else
   write(fid2,'(A)') 'cartesian'
  end if
  write(fid2,'(A)') '****'
  write(fid2,'(A)') TRIM(elem(i))//' 0'

  do while(.true.)
   read(fid1,'(A)') buf
   if(LEN_TRIM(buf) == 0) then
    write(fid2,'(A)') '****'
    exit
   end if

   read(buf,*) stype, nline
   if(stype == 'L') then
    stype = 'SP'; k = 3
   else
    k = 2
   end if
   write(fid2,'(A,1X,I3,A)') TRIM(stype), nline, '  1.00'

   rtmp = 0d0
   do j = 1, nline, 1
    read(fid1,*) m, rtmp(1:k)
    write(fid2,'(3(2X,ES16.9))') rtmp(1:k)
   end do ! for j
  end do ! for while

  if(all_ecp(i)%ecp) then
   write(fid2,'(A)') TRIM(elem(i))//'   0'
   m = all_ecp(i)%highest
   write(fid2,'(A,2(1X,I3))') TRIM(elem(i))//'-ECP', m, all_ecp(i)%core_e

   do j = 0, m, 1
    if(j == 0) then
     write(fid2,'(A)') am(m)//' potential'
    else
     write(fid2,'(A)') am(j-1)//'-'//am(m)//' potential'
    end if
    n = all_ecp(i)%potential(j+1)%n
    write(fid2,'(I3)') n
    do k = 1, n, 1
     write(fid2,'(I1,2(1X,ES16.9))') all_ecp(i)%potential(j+1)%col2(k), &
      all_ecp(i)%potential(j+1)%col3(k), all_ecp(i)%potential(j+1)%col1(k)
    end do ! for k
   end do ! for j
   write(fid2,'(A)') '****'
  end if
 end do ! for i

 deallocate(all_ecp)
 write(fid2,'(A)') '}'

 if(X2C) then
  write(fid2,'(/,A)') 'basis mybas1 {'
  do i = 1, natom, 1
   write(fid2,'(2(A,I0))') ' assign '//TRIM(elem(i)),ntimes(i),' gen',i
  end do ! for i

  call goto_data_section_in_gms_inp(inpname, fid1)
  read(fid1,'(A)') buf
  read(fid1,'(A)') buf

  do i = 1, natom, 1
   read(fid1,'(A)') buf ! buf contains elem(i), nuc(i) and coor(1:3,i)
   write(fid2,'(A,I0,A)') '[gen', i, ']'
   write(fid2,'(A)') 'spherical'
   write(fid2,'(A)') '****'
   write(fid2,'(A)') TRIM(elem(i))//' 0'
 
   do while(.true.)
    read(fid1,'(A)') buf
    if(LEN_TRIM(buf) == 0) then
     write(fid2,'(A)') '****'
     exit
    end if
 
    read(buf,*) stype, nline
    if(stype == 'L') then
     stype = 'SP'; k = 3
    else
     k = 2
    end if
 
    rtmp = 0d0
    do j = 1, nline, 1
     write(fid2,'(A)') TRIM(stype)//'   1  1.00'
     read(fid1,*) m, rtmp(1:k)
     write(fid2,'(3(2X,ES16.9))') rtmp(1:k)
    end do ! for j
   end do ! for while
  end do ! for i

  write(fid2,'(A)') '}'
  write(fid2,'(/,A)') 'set relativistic x2c'
  write(fid2,'(A)') 'set basis mybas'
  write(fid2,'(A)') 'set basis_relativistic mybas1'
 end if

 close(fid1)
 deallocate(ntimes, elem)

 call check_DKH_in_gms_inp(inpname, rel)
 select case(rel)
 case(-2) ! nothing
 case(-1) ! RESC
  write(6,'(/,A)') 'ERROR in subroutine bas_gms2psi: RESC keywords detected.'
  write(6,'(A)') 'But RESC is not supported in PSI4.'
  stop
 case(0,1,2,4)  ! DKH0/1/2/4
  if(.not. X2C) then
   write(fid2,'(A)') 'set relativistic dkh'
   write(fid2,'(A)') 'set basis_relativistic mybas'
   if(rel /= 2) write(fid2,'(A,I0)') 'set DKH_order ', rel
  end if
 case default
  write(6,'(/,A)') 'ERROR in subroutine bas_gms2psi: rel out of range!'
  write(6,'(A,I0)') 'rel=', rel
  close(fid2,status='delete')
  stop
 end select

 write(fid2,'(/,A)') 'set {'
 write(fid2,'(A)') ' scf_type pk'
 write(fid2,'(A)') ' s_tolerance 1e-6'
 write(fid2,'(A)') ' e_convergence 1e5'
 write(fid2,'(A)') ' d_convergence 1e5'

 if(uhf) then
  write(fid2,'(A)') ' reference uhf'
 else
  if(mult == 1) then
   write(fid2,'(A)') ' reference rhf'
  else
   write(fid2,'(A)') ' reference rohf'
  end if
 end if
 write(fid2,'(A)') '}'

 write(fid2,'(/,A)') "scfenergy, scf_wfn = energy('scf', return_wfn=True)"
 write(fid2,'(A)') '# This scf makes every array allocated. Its energy is usele&
                   &ss.'
 write(fid2,'(/,A)') "scf_wfn.Ca().load('"//TRIM(fileA)//"')"
 if(uhf) write(fid2,'(A)') "scf_wfn.Cb().load('"//TRIM(fileB)//"')"
 write(fid2,'(A)') 'scf_wfn.to_file(scf_wfn.get_scratch_filename(180))'

 write(fid2,'(/,A)') 'set {'
 write(fid2,'(A)') ' guess read'
 write(fid2,'(A)') ' e_convergence 1e-8'
 write(fid2,'(A)') ' d_convergence 1e-6'
 if(LEN_TRIM(dftname) == 0) then
  write(fid2,'(A)') '}'
  write(fid2,'(/,A)') "scfenergy = energy('scf')"
 else
  write(fid2,'(A)') ' dft_radial_points 99'
  write(fid2,'(A)') ' dft_spherical_points 590'
  write(fid2,'(A)') '}'
  dftname1 = dftname
  call upper(dftname1)
  select case(TRIM(dftname1))
  case('PBEPBE')
   write(fid2,'(/,A)') "scfenergy = energy('PBE')"
  case('PBE1PBE')
   write(fid2,'(/,A)') "scfenergy = energy('PBE0')"
  case('TPSSTPSS')
   write(fid2,'(/,A)') "scfenergy = energy('TPSS')"
  case default
   write(fid2,'(/,A)') "scfenergy = energy('"//TRIM(dftname)//"')"
  end select
 end if
 close(fid2)
end subroutine bas_gms2psi

