! written by jxzou at 20171108
! updated by jxzou at 20180222: support UHF type MOs
! updated by jxzou at 20190105: support ROGVB type MOs
! updated by jxzou at 20190226: change data to parameter when declare constant arrays
! updated by jxzou at 20200326: renamed to dat2fch; simplify code

! Copy the MOs in .dat or .inp file (Gamess) to .fch(k) file (Gaussian).
! The orders of Cartesian f, g and h functions in .dat file will be permuted.
! Note: an initial .fch(k) file must be provided, and MOs in it will be replaced.

! The order of Cartesian functions can be known by adding keyword 'pop=reg' in .gjf file.
! In Gamess output files(.out, .gms), the order of Cartesian functions is printed without
!  extra keyword needed.
! Limitation: only deal with Cartesian basis functions

program main
 implicit none
 integer i, npair, nopen
 integer, parameter :: iout = 6
 character(len=4) gvb_or_uhf, string
 character(len=240) datname, fchname

 npair = 0
 nopen = 0
 gvb_or_uhf = ' '
 string = ' '
 datname = ' '
 fchname = ' '

 i = iargc()
 if(.not. (i==2 .or. i==3 .or. i==4 .or. i==6) ) then
  write(iout,'(/,A)') ' ERROR in subroutine dat2fch: wrong command line parameter!'
  write(iout,'(/,A)') ' Example 1 (for RHF, ROHF and CASSCF): ./dat2fch a.dat a.fchk'
  write(iout,'(/,A)') ' Example 2 (for GVB): ./dat2fch a.dat a.fchk -gvb 4'
  write(iout,'(/,A)') ' Example 3 (for ROGVB): ./dat2fch a.dat a.fchk -gvb 4 -open 2'
  write(iout,'(/,A,/)') ' Example 4 (for UHF): ./dat2fch a.dat a.fchk -uhf'
  stop
 end if

 call getarg(1,datname)
 call getarg(2,fchname)

 if(i > 2) then

  call getarg(3,gvb_or_uhf)
  if(.not. (gvb_or_uhf=='-gvb' .or. gvb_or_uhf=='-uhf') ) then
   write(iout,'(A)') 'ERROR in subroutine dat2fch: the 3rd argument in command line is wrong!'
   write(iout,'(A)') "It must be '-gvb' or '-uhf'."
   stop
  end if

  if(gvb_or_uhf=='-gvb') then
   call getarg(4,string)
   read(string,*) npair
   ! '-open' is only valid in GVB, means ROGVB
   if(i == 6) then
    call getarg(6,string)
    read(string,*) nopen
   end if
  end if

 end if

 call dat2fch(datname, fchname, gvb_or_uhf, npair, nopen)
 stop
end program main

! transform MOs in .dat file into .fchk file
subroutine dat2fch(datname, fchname, gvb_or_uhf, npair, nopen)
 implicit none
 integer i, j, k
 integer RENAME
 integer ncoeff, nbf, nif, nalpha, nbeta
 ! nbf: the number of basis functions
 ! nif: the number of independent functions, i.e. the number of MOs
 ! nalpha: the number of alpha electrons
 ! nbeta: the number of beta electrons
 integer nline, nleft
 integer n10fmark, n15gmark, n21hmark
 integer, allocatable :: f_mark(:), g_mark(:), h_mark(:) ! mark the index where f,g,h functions begin
 integer, parameter :: datid = 11, fchid = 12, fchid1 = 13
 integer, parameter :: iout = 6
 integer, intent(in) :: npair, nopen
 integer, allocatable :: order(:)
 real(kind=8), allocatable :: alpha_coeff(:), alpha_coeff2(:,:)
 real(kind=8), allocatable :: beta_coeff(:), beta_coeff2(:,:)
 real(kind=8), allocatable :: temp_coeff(:,:)
 character(len=4), intent(in) :: gvb_or_uhf
 character(len=5) str1
 character(len=30) str2
 character(len=240), intent(in) :: datname, fchname
 character(len=240) fchname1, buffer

 str1 = ' '   ! initialization
 str2 = ' '
 buffer = ' '
 ncoeff = 0
 nline = 0
 nleft = 0

 fchname1 = TRIM(fchname)//'_d2f.tmp'
 open(unit=fchid,file=TRIM(fchname),status='old',position='rewind')

 ! find nalpha and nbeta
 do while(.true.)
  read(fchid,'(A)',iostat=i) buffer
  if(buffer(1:17) == 'Number of alpha e') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine dat2fch:'
  write(iout,'(A)') "'Number of alpha e' not found in file: "//TRIM(fchname)//'.'
  close(fchid)
  return
 end if
 BACKSPACE(fchid)
 read(fchid,'(A49,2X,I10)') buffer, nalpha
 read(fchid,'(A49,2X,I10)') buffer, nbeta
 ! nalpha and nbeta found

 ! check nopen ?= nalpha - nbeta
 if(gvb_or_uhf=='-gvb' .and. nopen/=nalpha-nbeta) then
  write(iout,'(A)') 'Warning in subroutine dat2fch:'
  write(iout,'(A)') 'nopen /= nalpha-nbeta detected. You should check if sth is wrong in .fchk file&
                   & or your command line parameters.'
 end if
 ! check done

 ! find nbf and nif
 do while(.true.)
  read(fchid,'(A)',iostat=i) buffer
  if(buffer(1:17) == 'Number of basis f') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine dat2fch:'
  write(iout,'(A)') "'Number of basis f' not found in the input file: "//TRIM(fchname)//'.'
  close(fchid)
  return
 end if
 BACKSPACE(fchid)
 read(fchid,'(A49,2X,I10)') buffer, nbf
 read(fchid,'(A49,2X,I10)') buffer, nif
 ! nbf and nif found

 ! find strings 'Alpha MO'
 do while(.true.)
  read(fchid,'(A)') buffer
  if(buffer(1:8) == 'Alpha MO') exit
 end do

 ! find ncoeff
 BACKSPACE(fchid)
 read(fchid,'(A49,2X,I10)') buffer, ncoeff
 ! ncoeff found

 if(nbf*nif /= ncoeff) then
  write(iout,'(A)') 'ERROR in subroutine dat2fch:'
  write(iout,'(A)') 'nbf*nif /= ncoeff detected in file '//TRIM(fchname)//'. Cannot deal with that.'
  close(fchid)
  return
 end if

 open(unit=datid,file=TRIM(datname),status='old',position='append')
 ! let's start from the end of the .dat file, and read upward

 ! Note: if there is more than one '$VEC' section in .dat file, only the last one will be read

 ! e.g. in the .dat file of a MCSCF(CASSCF) job, the 1st '$VEC' holds natural orbitals (only in active
 !  space) and the 2nd '$VEC' holds CASSCF orbitals (all MOs, maybe not natural orbitals)

 ! If you want to transform the natural orbitals into .fchk file, you need to use 'dat_no2fchk'

 do while(.true.)
  BACKSPACE(datid)
  BACKSPACE(datid)
  read(datid,'(A)',iostat=i) buffer
  if(i /= 0) exit
  if(buffer(2:5)=='$VEC' .or. buffer(2:5)=='$vec' .or. buffer(2:5)=='$Vec') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine dat2fch. No '$VEC' section in file "//TRIM(datname)//'.'
  close(datid)
  close(fchid)
  return
 end if

 ! read Alpha MO
 allocate(alpha_coeff2(nbf,nif))
 alpha_coeff2 = 0.0d0
 nline = nbf/5
 nleft = nbf - nline*5
 do i = 1, nif, 1
  k = 1
  do j = 1, nline, 1
   read(datid,'(A)') buffer
   buffer = buffer(6:)
   read(buffer,'(5ES15.8)') alpha_coeff2(k:k+4,i)
   k = k + 5
  end do
  if(nleft > 0) then
   read(datid,'(A)') buffer
   buffer = buffer(6:)
   str1 = ' '
   write(str1,'(I5)') nleft
   str1 = ADJUSTL(str1)
   str2 = '('//TRIM(str1)//'ES15.8)'
   read(buffer,TRIM(str2)) alpha_coeff2(k:nbf,i)
  end if
 end do

 ! if '-uhf' is specified, read Beta MO
 if(gvb_or_uhf == '-uhf') then
  allocate(beta_coeff2(nbf,nif))
  beta_coeff2 = 0.0d0
  do i = 1,nif,1
   k = 1
   do j = 1,nline,1
    read(datid,'(A)') buffer
    buffer = buffer(6:)
    read(buffer,'(5ES15.8)') beta_coeff2(k:k+4,i)
    k = k + 5
   end do
   if(nleft > 0) then
    read(datid,'(A)') buffer
    buffer = buffer(6:)
    str1 = ' '
    write(str1,'(I5)') nleft
    str1 = ADJUSTL(str1)
    str2 = '('//TRIM(str1)//'ES15.8)'
    read(buffer,TRIM(str2)) beta_coeff2(k:nbf,i)
   end if
  end do
 end if

 ! find the $DATA section and record the indices of Cartesian f, g and h functions
 rewind(datid)
 do while(.true.)
  read(datid,'(A)',iostat=i) buffer
  if(i /= 0) exit
  if(buffer(2:6)=='$DATA' .or. buffer(2:6)=='$data' .or. buffer(2:6)=='$Data') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine dat2fch: no '$DATA' section in the input file:"
  write(iout,'(A)') TRIM(datname)
  close(datid)
  close(fchid)
  return
 end if

 ! skip the point group and Title Card lines
 do i = 1,2
  read(datid,'(A)') buffer
 end do

 n10fmark = 0
 n15gmark = 0
 n21hmark = 0
 allocate(f_mark(nbf), g_mark(nbf), h_mark(nbf))
 f_mark = 0
 g_mark = 0
 h_mark = 0
 nbf = 0
 do while(.true.)
  read(datid,'(A)') buffer
  buffer = ADJUSTL(buffer)
  if(buffer(1:1)=='$' .or. buffer(1:1)==' ') exit ! meet '$END' or blank line
  do while(.true.)
   read(datid,'(A)') buffer
   buffer = ADJUSTL(buffer)
   if(buffer(1:1)=='$' .or. buffer(1:1)==' ') exit ! meet '$END' or blank line
   select case(buffer(1:1))
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
   call get_int_after_flag(buffer,' ',nline)
   do i = 1, nline, 1  ! skip n lines
    read(datid,'(A)') buffer
   end do
  end do
 end do
 close(datid)
 ! done recording

 ! adjust the order of Cartesian f, g, h functions
 do i = 1, n10fmark, 1
  call dat2fch_permute_10f(nif, alpha_coeff2(f_mark(i)+3:f_mark(i)+8,:))
 end do
 do i = 1, n15gmark, 1
  call dat2fch_permute_15g(nif, alpha_coeff2(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1, n21hmark, 1
  call dat2fch_permute_21h(nif, alpha_coeff2(h_mark(i):h_mark(i)+20,:))
 end do

 if(gvb_or_uhf == '-uhf') then
  do i = 1, n10fmark, 1
   call dat2fch_permute_10f(nif, beta_coeff2(f_mark(i)+3:f_mark(i)+8,:))
  end do
  do i = 1, n15gmark, 1
   call dat2fch_permute_15g(nif, beta_coeff2(g_mark(i):g_mark(i)+14,:))
  end do
  do i = 1, n21hmark, 1
   call dat2fch_permute_21h(nif, beta_coeff2(h_mark(i):h_mark(i)+20,:))
  end do
 end if
 deallocate(f_mark, g_mark, h_mark)
 ! adjust Cartesian functions finished

 ! if '-gvb' is specified, permute the active orbitals
 ! if nopen > 0, singly occupied orbitals are considered
 if(gvb_or_uhf=='-gvb' .and. npair>0) then

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
   temp_coeff(1:nbf,i) = alpha_coeff2(1:nbf,order(i))
  end forall
  forall(i = 1:k)
   alpha_coeff2(1:nbf,i+j) = temp_coeff(1:nbf,i)
  end forall
  deallocate(order,temp_coeff)

 end if
 ! permute done

 ! output the MOs to .fch(k) file
 open(unit=fchid1,file=TRIM(fchname1),status='replace')
 rewind(fchid)
 do while(.true.)
  read(fchid,'(A)') buffer
  write(fchid1,'(A)') TRIM(buffer)
  if(buffer(1:8) == 'Alpha MO') exit
 end do

 allocate(alpha_coeff(ncoeff))
 alpha_coeff = 0.0d0
 alpha_coeff = RESHAPE(alpha_coeff2,(/ncoeff/))
 write(fchid1,'(5(1X,ES15.8))') (alpha_coeff(i),i=1,ncoeff)
 deallocate(alpha_coeff, alpha_coeff2)

 ! write Beta MOs, if any
 if(gvb_or_uhf == '-uhf') then
  do while(.true.) ! skip the original Beta MO in .fchk file
   read(fchid,'(A)') buffer
   if(buffer(1:7) == 'Beta MO') exit
  end do
  write(fchid1,'(A)') TRIM(buffer)
  allocate(beta_coeff(ncoeff))
  beta_coeff = 0.0d0
  beta_coeff = RESHAPE(beta_coeff2,(/ncoeff/))
  write(fchid1,'(5(1X,ES15.8))') (beta_coeff(i),i=1,ncoeff)
  deallocate(beta_coeff, beta_coeff2)
 end if
 ! output MOs finished

 ! copy the rest of the .fch(k) file
 do while(.true.)
  read(fchid,'(A)') buffer
  if(index(buffer,'=') /= 0) exit
 end do
 BACKSPACE(fchid)

 do while(.true.)
  read(fchid,'(A)',iostat=i) buffer
  if(i /= 0) exit
  write(fchid1,'(A)') TRIM(buffer)
 end do
 ! copy finished

 close(fchid,status='delete')
 close(fchid1)
 i = RENAME(TRIM(fchname1),TRIM(fchname))
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
subroutine get_int_after_flag(buffer, flag, k)
 implicit none
 integer, intent(out) :: k
 character(len=1), intent(in) :: flag
 character(len=240), intent(in) :: buffer
 character(len=240) buf

 buf = ' '
 buf = buffer

 k = index(buf, flag)
 if(k == 0) then
  write(*,'(A)') "ERROR in subroutine get_int_after_flag: no keyword '"//flag &
               & //"' found in '"//TRIM(buf)//"'."
  stop
 end if
 buf(1:k) = ' '
 buf = ADJUSTL(buf)
 read(buf,*) k
 return
end subroutine get_int_after_flag

