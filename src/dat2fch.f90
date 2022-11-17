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
 use fch_content, only: check_uhf_in_fch
 implicit none
 integer :: i, npair, nopen, idx1, idx2
 character(len=4) :: gvb_or_uhf_or_cas, string
 character(len=240) :: datname, fchname
 logical :: uhf

 npair = 0; nopen = 0; idx1 = 0; idx2 = 0
 gvb_or_uhf_or_cas = ' '
 string = ' '
 datname = ' '
 fchname = ' '

 i = iargc()
 if(.not. (i==2 .or. i==4 .or. i==5 .or. i==6) ) then
  write(6,'(/,A)') ' ERROR in subroutine dat2fch: wrong command line arguments!'
  write(6,'(A)')   ' Example 1 (for R(O)HF, UHF, CAS): dat2fch a.dat a.fch'
  write(6,'(A)')   ' Example 2 (for GVB)             : dat2fch a.dat a.fch -gvb 4'
  write(6,'(A)')   ' Example 3 (for ROGVB)           : dat2fch a.dat a.fch -gvb 4 -open 2'
  write(6,'(A,/)') ' Example 4 (for CAS NOs)         : dat2fch a.dat a.fch -no 5 10'
  stop
 end if

 call getarg(1,datname)
 call require_file_exist(datname)

 call getarg(2,fchname)
 call require_file_exist(fchname)
 call check_uhf_in_fch(fchname, uhf)

 if(i > 2) then
  call getarg(3,gvb_or_uhf_or_cas)
  gvb_or_uhf_or_cas = ADJUSTL(gvb_or_uhf_or_cas)

  select case(TRIM(gvb_or_uhf_or_cas))
  case('-gvb')
   call getarg(4,string)
   read(string,*) npair
   ! '-open' is only valid in GVB, means ROGVB
   if(i == 6) then
    call getarg(6,string)
    read(string,*) nopen
   end if
  case('-no')
   call getarg(4,string)
   read(string,*) idx1
   call getarg(5,string)
   read(string,*) idx2
   if(idx1<1 .or. idx2<2) then
    write(6,'(/,A)') 'ERROR in subroutine dat2fch: invalid idx1 and/or idx2.'
    write(6,'(2(A,I0))') 'idx1=', idx1, ', idx2=', idx2
    stop
   end if
  case default
   write(6,'(/,A)') 'ERROR in subroutine dat2fch: the 3rd argument is wrong!'
   write(6,'(A)') "It must be '-gvb' or '-no'."
   stop
  end select
 else
  if(uhf) gvb_or_uhf_or_cas = '-uhf'
 end if

 call dat2fch(datname, fchname, gvb_or_uhf_or_cas, npair, nopen, idx1, idx2)
 stop
end program main

! transform MOs in .dat file into .fchk file
subroutine dat2fch(datname, fchname, gvb_or_uhf_or_cas, npair, nopen, idx1, idx2)
 use r_5D_2_6D, only: rd, rf, rg, rh
 use fch_content, only: read_mark_from_shltyp_cart
 implicit none
 integer :: i, j, k, datid, nline, nleft
 integer :: nbf, nif, na, nb
 ! nbf: the number of basis functions
 ! nif: the number of independent functions, i.e. the number of MOs
 ! na/nb: the number of alpha/beta electrons
 integer :: nbf1, nif1, ncontr, nd, nf, ng, nh
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 ! mark the index where f,g,h functions begin
 integer, intent(in) :: npair, nopen, idx1, idx2
 integer, allocatable :: order(:), shltyp(:), shl2atm(:)
 real(kind=8), allocatable :: alpha_coeff(:,:), beta_coeff(:,:)
 real(kind=8), allocatable :: all_coeff(:,:), temp_coeff(:,:)
 real(kind=8), allocatable :: noon(:) ! Natural Orbital Occupation Number
 character(len=4), intent(in) :: gvb_or_uhf_or_cas
 character(len=5) :: str1
 character(len=30) :: str2
 character(len=240), intent(in) :: datname, fchname
 character(len=240) :: fchname1, buf
 logical :: sph, noon_exist

 str1 = ' '   ! initialization
 str2 = ' '
 buf = ' '
 nline = 0
 nleft = 0

 fchname1 = TRIM(fchname)//'_d2f.tmp'
 call read_na_and_nb_from_fch(fchname, na, nb)

 ! check nopen ?= na - nb
 if(gvb_or_uhf_or_cas=='-gvb' .and. nopen/=na-nb) then
  write(6,'(/,A)') 'Warning in subroutine dat2fch: nopen /= na-nb detected.'
  write(6,'(A)') 'You should check if anything is wrong in .fch(k) file&
                   & or your command line arguments.'
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
  write(6,'(A)') 'ERROR in subroutine dat2fch: mixed spherical harmonic/&
                   &Cartesian functions detected.'
  write(6,'(A)') 'You probably used a basis set like 6-31G(d) in Gaussian. Its&
                   & default setting is (6D,7F).'
  write(6,'(A)') "You need to add '5D 7F' or '6D 10F' keywords in Gaussian."
  stop
 else if(ANY(shltyp<-1)) then
  sph = .true.
 else
  sph = .false.
 end if

 if(sph) then ! spherical harmonic
  nbf1 = nbf + COUNT(shltyp==-2) + 3*COUNT(shltyp==-3) + &
           & 6*COUNT(shltyp==-4) + 10*COUNT(shltyp==-5)
  ! [6D,10F,15G,21H] - [5D,7F,9G,11H] = [1,3,6,10]
 else         ! Cartesian
  nbf1 = nbf
 end if

 call read_nbf_from_dat(datname, i)
 if(i /= nbf1) then
  write(6,'(/,A)') 'ERROR in subroutine dat2fch: inconsistent nbf between&
                     & .fch and .dat file.'
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
    write(6,'(A)') 'ERROR in subroutine dat2fch: shltyp(i) out of range.'
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

 if(gvb_or_uhf_or_cas(1:3) == '-no') then
  ! only the doubly occupied MOs and active space NOs are held in .dat file,
  ! so we change nif to the number of NOs
  nif1 = idx2
  allocate(noon(nif), source=0d0)
  call read_on_from_dat(datname, nif1, noon(1:nif1), noon_exist)
  if(.not. noon_exist) deallocate(noon)

  open(newunit=datid,file=TRIM(datname),status='old',position='rewind')
  do while(.true.)
   read(datid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:5)=='$VEC' .or. buf(2:5)=='$vec' .or. buf(2:5)=='$Vec') exit
  end do ! for while

 else
  if(gvb_or_uhf_or_cas == '-uhf') then
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
  write(6,'(A)') "ERROR in subroutine dat2fch: No '$VEC' section in file "//TRIM(datname)
  close(datid)
  stop
 end if

 allocate(all_coeff(nbf1,nif1), source=0.0d0)

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
 if(gvb_or_uhf_or_cas == '-uhf') then
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

 allocate(d_mark(ncontr), f_mark(ncontr), g_mark(ncontr), h_mark(ncontr))
 call read_mark_from_shltyp_cart(ncontr, shl2atm, nd, nf, ng, nh, d_mark, &
                                 f_mark, g_mark, h_mark)
 deallocate(d_mark, shl2atm)

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
  allocate(alpha_coeff(nbf,min(nif,nif1)))
  alpha_coeff = temp_coeff(:,1:min(nif,nif1))
  if(gvb_or_uhf_or_cas == '-uhf') then
   allocate(beta_coeff(nbf,nif))
   beta_coeff = temp_coeff(:,nif+1:nif1)
  end if
  deallocate(temp_coeff)

 else ! Cartesian functions

  deallocate(shltyp)
  allocate(alpha_coeff(nbf,min(nif,nif1)))
  alpha_coeff = all_coeff(:,1:min(nif,nif1))
  if(gvb_or_uhf_or_cas == '-uhf') then
   allocate(beta_coeff(nbf,nif))
   beta_coeff = all_coeff(:,nif+1:nif1)
  end if
  deallocate(all_coeff)
 end if

 ! if '-gvb' is specified, permute the active orbitals
 ! if nopen > 0, singly occupied orbitals are considered
 if(gvb_or_uhf_or_cas=='-gvb' .and. npair>0) then
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

 if(gvb_or_uhf_or_cas(1:3) == '-no') then
  ! 1:idx2 are CASCI/CASSCF doubly occupied MOs and NOs
  ! copy virtual MOs from input .fch(k) file
  allocate(beta_coeff(nbf,nif), source=0d0)
  call read_mo_from_fch(fchname, nbf, nif, 'a', beta_coeff)
  beta_coeff(:,1:idx2) = alpha_coeff
  alpha_coeff = beta_coeff ! auto-reallocation
  deallocate(beta_coeff)
  if(noon_exist) call write_eigenvalues_to_fch(fchname,nif,'a',noon,.true.)
 end if

 ! output the MOs to .fch(k) file
 call write_mo_into_fch(fchname, nbf, nif, 'a', alpha_coeff)
 deallocate(alpha_coeff)

 ! write Beta MOs, if any
 if(gvb_or_uhf_or_cas == '-uhf') then
  call write_mo_into_fch(fchname, nbf, nif, 'b', beta_coeff)
  deallocate(beta_coeff)
 end if
end subroutine dat2fch

subroutine dat2fch_permute_10f(nif,coeff)
 implicit none
 integer i, nif
 integer, parameter :: order(6) = (/3, 1, 2, 5, 6, 4/)
 real(kind=8), intent(inout) :: coeff(6,nif)
 real(kind=8) coeff2(6,nif)
! From: the order of Cartesian f functions in Gamess
! To: the order of Cartesian f functions in Gaussian
!                     1    2    3    4    5    6
! From: XXX,YYY,ZZZ, XXY, XXZ, XYY, YYZ, XZZ, YZZ, XYZ
! To:   XXX,YYY,ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ

 coeff2 = 0.0d0
 forall(i = 1:6)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine dat2fch_permute_10f

subroutine dat2fch_permute_15g(nif,coeff)
 implicit none
 integer i
 integer, intent(in) :: nif
 integer, parameter :: order(15) = (/3, 9, 12, 7, 2, 8, 15, 14, 6, 11, 13, 10, 5, 4, 1/)
 real(kind=8), intent(inout) :: coeff(15,nif)
 real(kind=8) coeff2(15,nif)
! From: the order of Cartesian g functions in Gamess
! To: the order of Cartesian g functions in Gaussian
!       1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! From: XXXX,YYYY,ZZZZ,XXXY,XXXZ,XYYY,YYYZ,XZZZ,YZZZ,XXYY,XXZZ,YYZZ,XXYZ,XYYZ,XYZZ
! To:   ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX

 coeff2 = 0.0d0
 forall(i = 1:15)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine dat2fch_permute_15g

subroutine dat2fch_permute_21h(nif,coeff)
 implicit none
 integer i, nif
 integer, parameter :: order(21) = (/3, 9, 15, 13, 7, 2, 8, 18, 21, 17, 6, 14, &
                                     20, 19, 12, 11, 16, 10, 5, 4, 1/)
 real(kind=8), intent(inout) :: coeff(21,nif)
 real(kind=8) coeff2(21,nif)
! From: the order of Cartesian h functions in Gamess
! To: the order of Cartesian h functions in Gaussian
!       1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! From: XXXXX,YYYYY,ZZZZZ,XXXXY,XXXXZ,XYYYY,YYYYZ,XZZZZ,YZZZZ,XXXYY,XXXZZ,XXYYY,YYYZZ,XXZZZ,YYZZZ,XXXYZ,XYYYZ,XYZZZ,XXYYZ,XXYZZ,XYYZZ
! To:   ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX

 coeff2 = 0.0d0
 forall(i = 1:21)
  coeff2(i,:) = coeff(order(i),:)
 end forall
 coeff = coeff2
 return
end subroutine dat2fch_permute_21h

! get the integer value after the keyword flag
subroutine get_int_after_flag(buf, flag, k)
 implicit none
 integer :: i
 integer, intent(out) :: k
 character(len=1), intent(in) :: flag
 character(len=240), intent(in) :: buf

 i = index(buf, flag)
 if(i == 0) then
  write(*,'(A)') "ERROR in subroutine get_int_after_flag: no keyword '"//flag &
               & //"' found in the string '"//TRIM(buf)//"'."
  stop
 end if

 k = 0
 read(buf(i+1:),*) k
 return
end subroutine get_int_after_flag

