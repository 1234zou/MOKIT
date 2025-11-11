! written by jxzou at 20171108: transfer MOs from GAMESS .dat/.inp file to Gaussian .fch(k) file
! updated by jxzou at 20180222: support UHF type MOs
! updated by jxzou at 20190105: support ROGVB type MOs
! updated by jxzou at 20190226: change data to parameter when declare constant arrays
! updated by jxzou at 20200326: renamed to dat2fch; simplify code
! updated by jxzou at 20200620: support CASCI/CASSCF NOs, simplify code
! updated by jxzou at 20201115: support for spherical harmonic functions
!                               (expand 5D,7F,9G,11H to 6D,10F,15G,21H)
! updated by jxzou at 20201212: if -no, read $OCCNO for CASCI, write to .fch
! updated by jxzou at 20210413: remove '-uhf', add automatic determination

! The orders of Cartesian f, g and h functions in .dat file will be permuted.
! Note: an initial .fch(k) file must be provided, and MOs in it will be replaced.

! The order of Cartesian functions can be known by adding keyword 'pop=reg' in .gjf file.
! In Gamess output files(.out, .gms), the order of Cartesian functions is printed without
!  extra keyword needed.

program main
 implicit none
 integer :: i, j, npair, nopen, idx2, itype
 ! itype=0/1/2/3 for R(O)HF/UHF/GVB/Some kind of NOs
 character(len=4) :: str4
 character(len=5) :: str5
 character(len=240) :: datname, fchname
 logical :: alive, uhf

 i = iargc()
 select case(i)
 case(1,2,4,6)
 case default
  write(6,'(/,A)')' ERROR in program dat2fch: wrong command line arguments!'
  write(6,'(A)')  ' Example 1 (R(O)HF/UHF/CAS): dat2fch a.dat'
  write(6,'(A)')  ' Example 2 (R(O)HF/UHF/CAS): dat2fch a.dat a.fch'
  write(6,'(A)')  ' Example 3 (GVB)           : dat2fch a.dat a.fch -gvb 4'
  write(6,'(A)')  ' Example 4 (ROGVB)         : dat2fch a.dat a.fch -gvb 4 -open 2'
  write(6,'(A,/)')' Example 5 (CAS NOs)       : dat2fch a.dat a.fch -no 10'
  stop
 end select

 npair = 0; nopen = 0; idx2 = 0; itype = 0; alive = .false.
 str4 = ' '; str5 = ' '; datname = ' '; fchname = ' '

 call getarg(1, datname)
 call require_file_exist(datname)

 if(i == 1) then
  call find_specified_suffix(datname, '.dat', j)
  fchname = datname(1:j-1)//'.fch'
 else
  call getarg(2, fchname)
 end if

 inquire(file=TRIM(fchname),exist=alive)
 if(alive) then
  write(6,'(/,A)') 'Remark from program dat2fch: file '//TRIM(fchname)//' exists.'
  write(6,'(A)') 'It will be used directly.'
 else
  write(6,'(/,A)') 'Warning from program dat2fch: file '//TRIM(fchname)//' does&
                   & not exist,'
  write(6,'(A)') 'the program is trying to generate one from scratch...'
  call gen_fch_from_gms_inp(datname, fchname)
 end if

 call check_uhf_in_fch(fchname, uhf)
 if(uhf) itype = 1

 if(i > 2) then
  call getarg(3, str4)
  str4 = ADJUSTL(str4)

  select case(TRIM(str4))
  case('-gvb')
   itype = 2
   call getarg(4, str5)
   read(str5,*) npair
   if(i == 6) then ! '-open' is only valid in GVB, means ROGVB
    call getarg(6, str5)
    read(str5,*) nopen
   end if
  case('-no')
   itype = 3
   call getarg(4, str5)
   read(str5,*) idx2
   if(idx2 < 2) then
    write(6,'(/,A,I0)') 'ERROR in subroutine dat2fch: invalid idx2=',idx2
    stop
   end if
  case default
   write(6,'(/,A)') 'ERROR in subroutine dat2fch: the 3rd argument is wrong!'
   write(6,'(A)') "It must be '-gvb' or '-no'."
   stop
  end select
 end if

 if(uhf .and. itype>1) then
  write(6,'(/,A)') "ERROR in subroutine dat2fch: '-gvb'/'-no' argument is imc&
                   &ompatible with"
  write(6,'(A)') 'a UHF-type .fch file. An R(O)HF-type .fch file is required &
                 &in such case.'
  stop
 end if

 call dat2fch(datname, fchname, itype, npair, nopen, idx2)
 if(.not. alive) write(6,'(/,A,/)') 'Done generation.'
end program main

! transform MOs in .dat file into .fchk file
subroutine dat2fch(datname, fchname, itype, npair, nopen, idx2)
 use parameter_sph2cart, only: rd, rf, rg, rh
 implicit none
 integer :: i, j, k, datid, nline, nleft
 integer :: nbf, nif, na, nb
 ! nbf: the number of basis functions
 ! nif: the number of independent functions, i.e. the number of MOs
 ! na/nb: the number of alpha/beta electrons
 integer :: nbf1, nif1, ncontr, nd, nf, ng, nh, ni
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:), i_mark(:)
 ! mark the index where f,g,h,i functions begin
 integer, intent(in) :: itype, npair, nopen, idx2
 integer, allocatable :: order(:), shltyp(:), shl2atm(:)
 real(kind=8), allocatable :: alpha_coeff(:,:), beta_coeff(:,:)
 real(kind=8), allocatable :: all_coeff(:,:), temp_coeff(:,:)
 real(kind=8), allocatable :: noon(:) ! Natural Orbital Occupation Number
 character(len=5) :: str1
 character(len=30) :: str2
 character(len=240), intent(in) :: datname, fchname
 character(len=240) :: fchname1, buf
 logical :: sph, noon_exist

 str1 = ' '; str2 = ' '; buf = ' '; nline = 0; nleft = 0
 call find_specified_suffix(fchname, '.fch', i)
 fchname1 = fchname(1:i-1)//'.d2f'
 call read_na_and_nb_from_fch(fchname, na, nb)

 ! check nopen ?= na-nb
 if(itype==2 .and. nopen/=na-nb) then
  write(6,'(/,A)') 'Warning in subroutine dat2fch: nopen /= na-nb detected.'
  write(6,'(A)') 'You should check if anything is wrong in your .fch(k) file. O&
                 &r maybe your'
  write(6,'(A)') 'command line arguments are incorrect. Anyway, continue...'
  write(6,'(3(A,I0))') 'nopen=', nopen, ', na=', na, ', nb=', nb
 end if
 ! check done

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)

 ! read array shell_type from .fch file
 call read_ncontr_from_fch(fchname, ncontr)
 allocate(shltyp(ncontr), shl2atm(ncontr))
 call read_shltyp_and_shl2atm_from_fch(fchname, ncontr, shltyp, shl2atm)

 ! check if any spherical functions
 if(ANY(shltyp<-1) .and. ANY(shltyp>1)) then
  write(6,'(/,A)') 'ERROR in subroutine dat2fch: mixed spherical harmonic/Carte&
                   &sian functions detected.'
  write(6,'(A)') 'You probably used a basis set like 6-31G(d) in Gaussian. Its &
                 &default setting is (6D,7F).'
  write(6,'(A)') "You need to add '5D 7F' or '6D 10F' keywords in Gaussian."
  stop
 else if( ANY(shltyp>1) ) then
  sph = .false.
 else
  sph = .true.
 end if

 if(sph) then ! spherical harmonic
  nbf1 = nbf + COUNT(shltyp==-2) + 3*COUNT(shltyp==-3) + &
       & 6*COUNT(shltyp==-4) + 10*COUNT(shltyp==-5) + 15*COUNT(shltyp==-6)
  ! [6D,10F,15G,21H,28I] - [5D,7F,9G,11H,13I] = [1,3,6,10,15]
 else         ! Cartesian basis functions
  nbf1 = nbf
 end if

 call read_cart_nbf_nif_from_dat(datname, .true., i, j)
 if(i /= nbf1) then
  write(6,'(/,A)') 'ERROR in subroutine dat2fch: inconsistent nbf between .fch &
                   &and .dat file.'
  write(6,'(A)') 'datname='//TRIM(datname)
  write(6,'(A)') 'fchname='//TRIM(fchname)
  write(6,'(2(A,I0))') 'i=', i, ', nbf1=', nbf1
  stop
 end if

 if(sph) then   ! use the array shl2atm to hold Cartesian shltyp
  do i = 1, ncontr, 1
   select case(shltyp(i))
   case(-1:5)
    shl2atm(i) = shltyp(i)
   case(-5:-2)
    shl2atm(i) = -shltyp(i)
   case default
    write(6,'(/,A)') 'ERROR in subroutine dat2fch: shltyp(i) out of range.'
    write(6,'(2(A,I0))') 'i=', i, ', shltyp(i)=', shltyp(i)
    stop
   end select
  end do ! for i
 else
  shl2atm = shltyp
 end if

 ! If there is more than one '$VEC' section in .dat file, by default only
 ! the last one will be read.
 ! E.g. in the .dat file of a MCSCF(CASSCF) job, the 1st '$VEC' holds natural
 ! orbitals (only doubly occupied orbitals and active orbitals) and the 2nd
 ! '$VEC' holds CASSCF orbitals (not natural orbitals)
 ! If '-no' is specified, the NOs in the first '$VEC' section will also be read

 if(itype == 3) then
  ! only the doubly occupied MOs and active space NOs are held in .dat file,
  ! so we change nif to the number of NOs
  nif1 = idx2
  allocate(noon(nif))
  call read_on_from_gms_dat(datname, nif1, noon(1:nif1), noon_exist)
  if(.not. noon_exist) deallocate(noon)

  open(newunit=datid,file=TRIM(datname),status='old',position='rewind')
  do while(.true.)
   read(datid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:5)=='$VEC' .or. buf(2:5)=='$vec' .or. buf(2:5)=='$Vec') exit
  end do ! for while
 else ! itype /= 3
  if(itype == 1) then
   nif1 = 2*nif
  else
   nif1 = nif
  end if
  open(newunit=datid,file=TRIM(datname),status='old',position='append')
  do while(.true.)
   BACKSPACE(datid)
   BACKSPACE(datid)
   read(datid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:5)=='$VEC' .or. buf(2:5)=='$vec' .or. buf(2:5)=='$Vec') exit
  end do ! for while
 end if

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine dat2fch: no '$VEC' found in file "//&
                    TRIM(datname)
  close(datid)
  stop
 end if

 allocate(all_coeff(nbf1,nif1), source=0d0)

 ! read Alpha MO
 nline = nbf1/5
 nleft = nbf1 - nline*5
 do i = 1, min(nif1,nif), 1
  k = 1
  do j = 1, nline, 1
   read(datid,'(A)') buf
   buf = buf(6:)
   read(buf,'(5ES15.8)') all_coeff(k:k+4,i)
   k = k + 5
  end do ! for j
  if(nleft > 0) then
   read(datid,'(A)') buf
   buf = buf(6:)
   str1 = ' '
   write(str1,'(I5)') nleft
   str1 = ADJUSTL(str1)
   str2 = '('//TRIM(str1)//'ES15.8)'
   read(buf,TRIM(str2)) all_coeff(k:nbf1,i)
  end if
 end do ! for i

 ! if '-uhf' is specified, read Beta MO
 if(itype == 1) then
  do i = 1, nif, 1
   k = 1
   do j = 1, nline, 1
    read(datid,'(A)') buf
    buf = buf(6:)
    read(buf,'(5ES15.8)') all_coeff(k:k+4,nif+i)
    k = k + 5
   end do ! for j
   if(nleft > 0) then
    read(datid,'(A)') buf
    buf = buf(6:)
    str1 = ' '
    write(str1,'(I5)') nleft
    str1 = ADJUSTL(str1)
    str2 = '('//TRIM(str1)//'ES15.8)'
    read(buf,TRIM(str2)) all_coeff(k:nbf1,nif+i)
   end if
  end do ! for i
 end if

 allocate(d_mark(ncontr), f_mark(ncontr), g_mark(ncontr), h_mark(ncontr), &
         i_mark(ncontr))
 call read_mark_from_shltyp_cart(ncontr, shl2atm, nd, nf, ng, nh, ni, d_mark, &
                                 f_mark, g_mark, h_mark, i_mark)
 deallocate(d_mark, i_mark, shl2atm)

 ! adjust the order of Cartesian f, g, h functions
 do i = 1, nf, 1
  call dat2fch_permute_10f(nif1, all_coeff(f_mark(i)+3:f_mark(i)+8,:))
 end do
 do i = 1, ng, 1
  call dat2fch_permute_15g(nif1, all_coeff(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1, nh, 1
  call dat2fch_permute_21h(nif1, all_coeff(h_mark(i):h_mark(i)+20,:))
 end do
 deallocate(f_mark, g_mark, h_mark)
 ! adjust Cartesian functions finished

 if(sph) then ! spherical harmonic functions, linear comb.
  allocate(temp_coeff(nbf,nif1), source=0d0)
  nbf = 0; j = 0
  do i = 1, ncontr, 1
   select case(shltyp(i))
   case( 0) ! S
    temp_coeff(nbf+1,:) = all_coeff(j+1,:)
    nbf = nbf + 1; j= j + 1
   case( 1) ! P
    temp_coeff(nbf+1:nbf+3,:) = all_coeff(j+1:j+3,:)
    nbf = nbf + 3; j = j + 3
   case(-1) ! L
    temp_coeff(nbf+1:nbf+4,:) = all_coeff(j+1:j+4,:)
    nbf = nbf + 4; j = j + 4
   case(-2) ! 5D
    call solve_multi_lin_eqs(6,5,rd,nif1,all_coeff(j+1:j+6,:),temp_coeff(nbf+1:nbf+5,:))
    nbf = nbf + 5; j= j + 6
   case(-3) ! 7F
    call solve_multi_lin_eqs(10,7,rf,nif1,all_coeff(j+1:j+10,:),temp_coeff(nbf+1:nbf+7,:))
    nbf = nbf + 7; j = j + 10
   case(-4) ! 9G
    call solve_multi_lin_eqs(15,9,rg,nif1,all_coeff(j+1:j+15,:),temp_coeff(nbf+1:nbf+9,:))
    nbf = nbf + 9; j = j + 15
   case(-5) ! 11H
    call solve_multi_lin_eqs(21,11,rh,nif1,all_coeff(j+1:j+21,:),temp_coeff(nbf+1:nbf+11,:))
    nbf = nbf + 11; j = j + 21
   end select
  end do ! for i

  deallocate(shltyp, all_coeff)
  allocate(alpha_coeff(nbf,MIN(nif,nif1)))
  alpha_coeff = temp_coeff(:,1:MIN(nif,nif1))
  if(itype == 1) then
   allocate(beta_coeff(nbf,nif))
   beta_coeff = temp_coeff(:,nif+1:nif1)
  end if
  deallocate(temp_coeff)

 else ! Cartesian functions

  deallocate(shltyp)
  allocate(alpha_coeff(nbf,MIN(nif,nif1)))
  alpha_coeff = all_coeff(:,1:MIN(nif,nif1))
  if(itype == 1) then
   allocate(beta_coeff(nbf,nif))
   beta_coeff = all_coeff(:,nif+1:nif1)
  end if
  deallocate(all_coeff)
 end if

 ! if '-gvb' is specified, permute the active orbitals
 ! if nopen > 0, singly occupied orbitals are considered
 if(itype==2 .and. npair>0) then
  k = 2*npair + nopen
  allocate(order(k), source=0)
  allocate(temp_coeff(nbf,k), source=0d0)

  forall(i = 1:npair) order(i) = 2*i - 1 + nopen
  if(nopen > 0) then
   forall(i = npair+1:npair+nopen)
    order(i) = i - npair
   end forall
  end if
  forall(i = npair+nopen+1:k) order(i) = 2*(k+1-i) + nopen

  j = na - nopen - npair
  order = order + j
  forall(i = 1:k) temp_coeff(:,i) = alpha_coeff(:,order(i))
  forall(i = 1:k) alpha_coeff(:,i+j) = temp_coeff(:,i)
  deallocate(order, temp_coeff)
 end if
 ! permute done

 if(itype == 3) then
  ! 1:idx2 are CASCI/CASSCF doubly occupied MOs and NOs
  ! copy virtual MOs from input .fch(k) file
  allocate(beta_coeff(nbf,nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', beta_coeff)
  beta_coeff(:,1:idx2) = alpha_coeff
  deallocate(alpha_coeff)
  allocate(alpha_coeff(nbf,nif), source=beta_coeff)
  deallocate(beta_coeff)
  if(noon_exist) call write_eigenvalues_to_fch(fchname,nif,'a',noon,.true.)
 end if

 ! output the MOs to .fch(k) file
 call write_mo_into_fch(fchname, nbf, nif, 'a', alpha_coeff)
 deallocate(alpha_coeff)

 ! write Beta MOs, if any
 if(itype == 1) then
  call write_mo_into_fch(fchname, nbf, nif, 'b', beta_coeff)
  deallocate(beta_coeff)
 end if
end subroutine dat2fch

subroutine dat2fch_permute_10f(nif,coeff)
 implicit none
 integer :: i, nif
 integer, parameter :: order(6) = [3, 1, 2, 5, 6, 4]
 real(kind=8), intent(inout) :: coeff(6,nif)
 real(kind=8) :: coeff2(6,nif)
! From: the order of Cartesian f functions in Gamess
! To: the order of Cartesian f functions in Gaussian
!                     1    2    3    4    5    6
! From: XXX,YYY,ZZZ, XXY, XXZ, XYY, YYZ, XZZ, YZZ, XYZ
! To:   XXX,YYY,ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ

 coeff2 = 0d0
 forall(i = 1:6) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
end subroutine dat2fch_permute_10f

subroutine dat2fch_permute_15g(nif,coeff)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 integer, parameter :: order(15) = [3,9,12,7,2,8,15,14,6,11,13,10,5,4,1]
 real(kind=8), intent(inout) :: coeff(15,nif)
 real(kind=8) :: coeff2(15,nif)
! From: the order of Cartesian g functions in Gamess
! To: the order of Cartesian g functions in Gaussian
!       1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! From: XXXX,YYYY,ZZZZ,XXXY,XXXZ,XYYY,YYYZ,XZZZ,YZZZ,XXYY,XXZZ,YYZZ,XXYZ,XYYZ,XYZZ
! To:   ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX

 coeff2 = 0d0
 forall(i = 1:15) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
end subroutine dat2fch_permute_15g

subroutine dat2fch_permute_21h(nif,coeff)
 implicit none
 integer :: i, nif
 integer, parameter :: order(21) = [3,9,15,13,7,2,8,18,21,17,6,14,20,19,12,11,16,10,5,4,1]
 real(kind=8), intent(inout) :: coeff(21,nif)
 real(kind=8) :: coeff2(21,nif)
! From: the order of Cartesian h functions in Gamess
! To: the order of Cartesian h functions in Gaussian
!       1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! From: XXXXX,YYYYY,ZZZZZ,XXXXY,XXXXZ,XYYYY,YYYYZ,XZZZZ,YZZZZ,XXXYY,XXXZZ,XXYYY,YYYZZ,XXZZZ,YYZZZ,XXXYZ,XYYYZ,XYZZZ,XXYYZ,XXYZZ,XYYZZ
! To:   ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX

 coeff2 = 0d0
 forall(i = 1:21) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
end subroutine dat2fch_permute_21h

! Generate a .fch file from a GAMESS .dat/.inp file. This requires the basis
! data and ECP/PP are well written in the .dat file.
subroutine gen_fch_from_gms_inp(datname, fchname)
 use fch_content
 use pg, only: natom0=>natom, ecp_exist, all_ecp
 use mkl_content, only: natom1=>natom, ncontr1=>ncontr,shell_type1=>shell_type,&
  shl2atm1=>shl2atm, all_pg, read_all_pg_from_gms_inp, un_normalized_all_pg, &
  merge_s_and_p_into_sp
 implicit none
 integer :: i, ne, isph, nbf1
 character(len=240), intent(in) :: datname, fchname
 logical :: alive, ghf, has_sp
 logical, allocatable :: ghost(:)

 call read_irel_from_gms_inp(datname, irel)
 call read_natom_from_gms_inp(datname, natom)
 allocate(mull_char(natom))
 call read_mul_char_from_dat(datname, natom, mull_char, alive)

 allocate(elem(natom), ielem(natom), coor(3,natom), ghost(natom))
 call read_elem_nuc_coor_from_gms_inp(datname, natom, elem, ielem, coor, ghost)

 allocate(iatom_type(natom), source=0)
 forall(i=1:natom, ghost(i)) iatom_type(i) = 1000
 deallocate(ghost)

 natom0 = natom
 call read_charge_mult_isph_from_gms_inp(datname, charge, mult, isph, is_uhf, &
                                         ghf, ecp_exist)
 if(ghf) then
  write(6,'(/,A)') 'ERROR in subroutine gen_fch_from_gms_inp: GHF-type is unsup&
                   &ported.'
  stop
 end if
 if(isph /= 1) then
  write(6,'(/,A)') 'ERROR in subroutine gen_fch_from_gms_inp: ISPHER/=1 in file&
                   & '//TRIM(datname)
  stop
 end if

 ne = SUM(ielem) - charge
 if(ecp_exist) then
  ! Do not change ielem here. It will be considered in subroutine write_fch.
  call read_all_ecp_from_gms_inp(datname)
  ne = ne - SUM(all_ecp(:)%core_e)
  call find_LenNCZ_in_all_ecp(natom, all_ecp, LenNCZ)
  allocate(KFirst(natom,10), KLast(natom,10), Lmax(natom), LPSkip(natom), &
   NLP(LenNCZ), RNFroz(natom), CLP(LenNCZ), ZLP(LenNCZ))
  call all_ecp2ecp_arrays(all_ecp, natom, LenNCZ, KFirst, KLast, Lmax, LPSkip, &
                          NLP, RNFroz, CLP, ZLP)
  deallocate(all_ecp)
 end if
 nopen = mult - 1
 na = (ne + nopen)/2
 nb = ne - na

 call read_ncontr_from_gms_inp(datname, ncontr)
 allocate(shell_type(ncontr), shell2atom_map(ncontr))
 call read_shltyp_and_shl2atm_from_gms_inp(datname, ncontr, shell_type, &
                                           shell2atom_map)
 ! the shell_type read from .inp/.dat does not have the difference of sph/Cart,
 ! now check this
 if(isph == 1) then
  forall(i=1:ncontr, (shell_type(i)>1)) shell_type(i) = -shell_type(i)
 end if

 natom1 = natom
 ncontr1 = ncontr
 allocate(shell_type1(ncontr), source=shell_type)
 allocate(shl2atm1(ncontr), source=shell2atom_map)
 call read_all_pg_from_gms_inp(datname)
 call un_normalized_all_pg()
 call merge_s_and_p_into_sp()
 if(ncontr1 < ncontr) then
  deallocate(shell_type, shell2atom_map)
  allocate(shell_type(ncontr1), source=shell_type1)
  allocate(shell2atom_map(ncontr1), source=shl2atm1)
  ncontr = ncontr1
 end if

 allocate(prim_per_shell(ncontr))
 call find_nprim_from_all_pg(ncontr, prim_per_shell, nprim, has_sp)
 call all_pg2prim_exp_and_contr_coeff(has_sp)
 deallocate(all_pg, shell_type1, shl2atm1)

 ! read the number of basis functions under Cartesian-type basis functions
 call read_cart_nbf_nif_from_dat(datname, .false., nbf1, nif)
 if(is_uhf) then
  if(MOD(nif,2) /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine gen_fch_from_gms_inp: nif is not even!'
   write(6,'(A)') 'Strange case. File='//TRIM(datname)
   stop
  end if
  nif = nif/2
 end if
 allocate(eigen_e_a(nif), source=0d0)

 ! convert to the number of basis functions using spherical harmonic functions
 nbf = nbf1 - COUNT(shell_type==-2) - 3*COUNT(shell_type==-3) - &
     & 6*COUNT(shell_type==-4) - 10*COUNT(shell_type==-5) - 15*COUNT(shell_type==-6)

 allocate(alpha_coeff(nbf,nif), source=0d0)
 allocate(tot_dm(nbf,nbf), source=0d0)
 if(is_uhf) then
  allocate(eigen_e_b(nif), source=0d0)
  allocate(beta_coeff(nbf,nif), source=0d0)
  allocate(spin_dm(nbf,nbf), source=0d0)
 end if

 call write_fch(fchname)
 call free_arrays_in_fch_content()
end subroutine gen_fch_from_gms_inp

