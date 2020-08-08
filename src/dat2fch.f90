! written by jxzou at 20171108
! updated by jxzou at 20180222: support UHF type MOs
! updated by jxzou at 20190105: support ROGVB type MOs
! updated by jxzou at 20190226: change data to parameter when declare constant arrays
! updated by jxzou at 20200326: renamed to dat2fch; simplify code
! updated by jxzou at 20200620: support CASCI/CASSCF NOs, simplify code

! Copy the MOs in .dat or .inp file (Gamess) to .fch(k) file (Gaussian).
! The orders of Cartesian f, g and h functions in .dat file will be permuted.
! Note: an initial .fch(k) file must be provided, and MOs in it will be replaced.

! The order of Cartesian functions can be known by adding keyword 'pop=reg' in .gjf file.
! In Gamess output files(.out, .gms), the order of Cartesian functions is printed without
!  extra keyword needed.
! Limitation: only deal with Cartesian basis functions

program main
 implicit none
 integer :: i, npair, nopen, idx1, idx2
 integer, parameter :: iout = 6
 character(len=4) :: gvb_or_uhf_or_cas, string
 character(len=240) :: datname, fchname

 npair = 0
 nopen = 0
 gvb_or_uhf_or_cas = ' '
 string = ' '
 datname = ' '
 fchname = ' '

 i = iargc()
 if(i<1 .or. i>6) then
  write(iout,'(/,A)') ' ERROR in subroutine dat2fch: wrong command line parameter!'
  write(iout,'(/,A)') ' Example 1 (for RHF, ROHF and CASSCF): ./dat2fch a.dat a.fchk'
  write(iout,'(/,A)') ' Example 2 (for GVB): ./dat2fch a.dat a.fchk -gvb 4'
  write(iout,'(/,A)') ' Example 3 (for ROGVB): ./dat2fch a.dat a.fchk -gvb 4 -open 2'
  write(iout,'(/,A)') ' Example 4 (for UHF): ./dat2fch a.dat a.fchk -uhf'
  write(iout,'(/,A,/)') ' Example 5 (for CAS NOs): ./dat2fch a.dat a.fchk -no 5 10'
  stop
 end if

 call getarg(1,datname)
 call getarg(2,fchname)

 if(i > 2) then
  call getarg(3,gvb_or_uhf_or_cas)
  gvb_or_uhf_or_cas = ADJUSTL(gvb_or_uhf_or_cas)
  if(.not. (gvb_or_uhf_or_cas=='-gvb' .or. gvb_or_uhf_or_cas=='-uhf' &
       .or. gvb_or_uhf_or_cas=='-no') ) then
   write(iout,'(A)') 'ERROR in subroutine dat2fch: the 3rd argument is wrong!'
   write(iout,'(A)') "It must be '-gvb', '-uhf' or '-no'."
   stop
  end if

  if(gvb_or_uhf_or_cas == '-gvb') then
   call getarg(4,string)
   read(string,*) npair
   ! '-open' is only valid in GVB, means ROGVB
   if(i == 6) then
    call getarg(6,string)
    read(string,*) nopen
   end if
  else if(gvb_or_uhf_or_cas(1:3) == '-no') then
   call getarg(4,string)
   read(string,*) idx1
   call getarg(5,string)
   read(string,*) idx2
   if(idx1<1 .or. idx2<2) then
    write(iout,'(A)') 'ERROR in subroutine dat2fch: invalid idx1 and/or idx2.'
    write(iout,'(2(A,I0))') 'idx1=', idx1, ', idx2=', idx2
    stop
   end if
  end if
 end if

 call dat2fch(datname, fchname, gvb_or_uhf_or_cas, npair, nopen, idx1, idx2)
 stop
end program main

! transform MOs in .dat file into .fchk file
subroutine dat2fch(datname, fchname, gvb_or_uhf_or_cas, npair, nopen, idx1, idx2)
 implicit none
 integer :: i, j, k, datid, nline, nleft, RENAME
 integer :: nbf, nif, nif0, nalpha, nbeta
 ! nbf: the number of basis functions
 ! nif: the number of independent functions, i.e. the number of MOs
 ! nalpha: the number of alpha electrons
 ! nbeta: the number of beta electrons
 integer :: n10fmark, n15gmark, n21hmark
 integer, allocatable :: f_mark(:), g_mark(:), h_mark(:) ! mark the index where f,g,h functions begin
 integer, parameter :: iout = 6
 integer, intent(in) :: npair, nopen, idx1, idx2
 integer, allocatable :: order(:)
 real(kind=8), allocatable :: alpha_coeff(:,:), beta_coeff(:,:)
 real(kind=8), allocatable :: temp_coeff(:,:)
 character(len=4), intent(in) :: gvb_or_uhf_or_cas
 character(len=5) :: str1
 character(len=30) :: str2
 character(len=240), intent(in) :: datname, fchname
 character(len=240) :: fchname1, buf

 str1 = ' '   ! initialization
 str2 = ' '
 buf = ' '
 nline = 0
 nleft = 0

 fchname1 = TRIM(fchname)//'_d2f.tmp'
 call read_na_and_nb_from_fch(fchname, nalpha, nbeta)

 ! check nopen ?= nalpha - nbeta
 if(gvb_or_uhf_or_cas=='-gvb' .and. nopen/=nalpha-nbeta) then
  write(iout,'(A)') 'Warning in subroutine dat2fch: nopen /= nalpha-nbeta detected.'
  write(iout,'(A)') 'You should check if anything is wrong in .fch(k) file&
                   & or your command line arguments.'
 end if
 ! check done

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)

 ! Note: if there is more than one '$VEC' section in .dat file, by default only the last one will be read
 ! e.g. in the .dat file of a MCSCF(CASSCF) job, the 1st '$VEC' holds natural orbitals (only doubly
 ! occupied orbitals and active orbitals) and the 2nd '$VEC' holds CASSCF orbitals (not natural orbitals)
 ! If '-no' is specified, the NOs in the first '$VEC' section will also be read

 if(gvb_or_uhf_or_cas(1:3) == '-no') then
  ! only the doubly occupied MOs and active space NOs are held in .dat file,
  ! so we change nif to the number of NOs
  nif0 = nif
  nif = idx2 - idx1 + 1

  open(newunit=datid,file=TRIM(datname),status='old',position='rewind')
  do while(.true.)
   read(datid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:5)=='$VEC' .or. buf(2:5)=='$vec' .or. buf(2:5)=='$Vec') exit
  end do ! for while

 else
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
  write(iout,'(A)') "ERROR in subroutine dat2fch: No '$VEC' section in file "//TRIM(datname)//'.'
  close(datid)
  stop
 end if

 ! read Alpha MO
 allocate(alpha_coeff(nbf,nif), source=0.0d0)
 nline = nbf/5
 nleft = nbf - nline*5
 do i = 1, nif, 1
  k = 1
  do j = 1, nline, 1
   read(datid,'(A)') buf
   buf = buf(6:)
   read(buf,'(5ES15.8)') alpha_coeff(k:k+4,i)
   k = k + 5
  end do ! for j
  if(nleft > 0) then
   read(datid,'(A)') buf
   buf = buf(6:)
   str1 = ' '
   write(str1,'(I5)') nleft
   str1 = ADJUSTL(str1)
   str2 = '('//TRIM(str1)//'ES15.8)'
   read(buf,TRIM(str2)) alpha_coeff(k:nbf,i)
  end if
 end do ! for i

 ! if '-uhf' is specified, read Beta MO
 if(gvb_or_uhf_or_cas == '-uhf') then
  allocate(beta_coeff(nbf,nif), source=0.0d0)
  do i = 1,nif,1
   k = 1
   do j = 1,nline,1
    read(datid,'(A)') buf
    buf = buf(6:)
    read(buf,'(5ES15.8)') beta_coeff(k:k+4,i)
    k = k + 5
   end do ! for j
   if(nleft > 0) then
    read(datid,'(A)') buf
    buf = buf(6:)
    str1 = ' '
    write(str1,'(I5)') nleft
    str1 = ADJUSTL(str1)
    str2 = '('//TRIM(str1)//'ES15.8)'
    read(buf,TRIM(str2)) beta_coeff(k:nbf,i)
   end if
  end do ! for i
 end if

 ! find the $DATA section and record the indices of Cartesian f, g and h functions
 rewind(datid)
 do while(.true.)
  read(datid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:6)=='$DATA' .or. buf(2:6)=='$data' .or. buf(2:6)=='$Data') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine dat2fch: no '$DATA' section in the input file:"
  write(iout,'(A)') TRIM(datname)
  close(datid)
  stop
 end if

 ! skip the point group and Title Card lines
 read(datid,'(A)') buf
 read(datid,'(A)') buf

 n10fmark = 0
 n15gmark = 0
 n21hmark = 0
 allocate(f_mark(nbf), source=0)
 allocate(g_mark(nbf), source=0)
 allocate(h_mark(nbf), source=0)
 nbf = 0
 do while(.true.)
  read(datid,'(A)') buf
  buf = ADJUSTL(buf)
  if(buf(1:1)=='$' .or. buf(1:1)==' ') exit ! meet '$END' or blank line
  do while(.true.)
   read(datid,'(A)') buf
   buf = ADJUSTL(buf)
   if(buf(1:1)=='$' .or. buf(1:1)==' ') exit ! meet '$END' or blank line
   select case(buf(1:1))
   case('S')
    nbf = nbf + 1
   case('P')
    nbf = nbf + 3
   case('L') ! 'L' is 'SP' in Gaussian format
    nbf = nbf + 4
   case('D')
    nbf = nbf + 6
   case('F')
    n10fmark = n10fmark + 1
    f_mark(n10fmark) = nbf + 1
    nbf = nbf + 10
   case('G')
    n15gmark = n15gmark + 1
    g_mark(n15gmark) = nbf + 1
    nbf = nbf + 15
   case('H')
    n21hmark = n21hmark + 1
    h_mark(n21hmark) = nbf + 1
    nbf = nbf + 21
   end select
   call get_int_after_flag(buf,' ',nline)
   do i = 1, nline, 1  ! skip n lines
    read(datid,'(A)') buf
   end do
  end do
 end do
 close(datid)
 ! done recording

 ! adjust the order of Cartesian f, g, h functions
 do i = 1, n10fmark, 1
  call dat2fch_permute_10f(nif, alpha_coeff(f_mark(i)+3:f_mark(i)+8,:))
 end do
 do i = 1, n15gmark, 1
  call dat2fch_permute_15g(nif, alpha_coeff(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1, n21hmark, 1
  call dat2fch_permute_21h(nif, alpha_coeff(h_mark(i):h_mark(i)+20,:))
 end do

 if(gvb_or_uhf_or_cas == '-uhf') then
  do i = 1, n10fmark, 1
   call dat2fch_permute_10f(nif, beta_coeff(f_mark(i)+3:f_mark(i)+8,:))
  end do
  do i = 1, n15gmark, 1
   call dat2fch_permute_15g(nif, beta_coeff(g_mark(i):g_mark(i)+14,:))
  end do
  do i = 1, n21hmark, 1
   call dat2fch_permute_21h(nif, beta_coeff(h_mark(i):h_mark(i)+20,:))
  end do
 end if
 deallocate(f_mark, g_mark, h_mark)
 ! adjust Cartesian functions finished

 ! if '-gvb' is specified, permute the active orbitals
 ! if nopen > 0, singly occupied orbitals are considered
 if(gvb_or_uhf_or_cas=='-gvb' .and. npair>0) then

  k = 2*npair + nopen
  allocate(order(k), temp_coeff(nbf,k))
  order = 0
  temp_coeff = 0.0d0

  forall(i = 1:npair)
   order(i) = 2*i - 1 + nopen
  end forall
  if(nopen > 0) then
   forall(i = npair+1:npair+nopen)
    order(i) = i - npair
   end forall
  end if
  forall(i = npair+nopen+1:k)
   order(i) = 2*(k+1-i) + nopen
  end forall

  j = nalpha - nopen - npair
  forall(i = 1:k)
   order(i) = order(i) + j
  end forall
  forall(i = 1:k)
   temp_coeff(1:nbf,i) = alpha_coeff(1:nbf,order(i))
  end forall
  forall(i = 1:k)
   alpha_coeff(1:nbf,i+j) = temp_coeff(1:nbf,i)
  end forall
  deallocate(order,temp_coeff)

 end if
 ! permute done

 if(gvb_or_uhf_or_cas(1:3) == '-no') then
  ! 1:idx2 are CASCI/CASSCF doubly occupied MOs and NOs
  ! copy virtual MOs from input .fch(k) file
  allocate(beta_coeff(nbf,nif0), source=0.0d0)
  call read_mo_from_fch(fchname, nbf, nif0, 'a', beta_coeff)
  beta_coeff(:,1:idx2) = alpha_coeff
  alpha_coeff = beta_coeff ! auto-reallocation
  deallocate(beta_coeff)
  nif = nif0 ! update nif
 end if

 ! output the MOs to .fch(k) file
 call write_mo_into_fch(fchname, nbf, nif, 'a', alpha_coeff)
 deallocate(alpha_coeff)

 ! write Beta MOs, if any
 if(gvb_or_uhf_or_cas == '-uhf') then
  call write_mo_into_fch(fchname, nbf, nif, 'b', beta_coeff)
  deallocate(beta_coeff)
 end if

 return
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

