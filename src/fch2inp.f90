! written by jxzou at 20171120
! updated by jxzou at 20171224: use '-gvb npair' in command line to permute the active orbtials
! updated by jxzou at 20180211: support UHF type MOs
! updated by jxzou at 20180620: support open shell orbitals (doublet, triplet, ...)
! updated by jxzou at 20180824: fix the bug: forgot to move open shell orbitals when npair=1
! updated by jxzou at 20190226: change data to parameter when declare constant arrays
! updated by jxzou at 20200322: renamed as fch2inp; simplify code; support ECP/PP

! generate .inp file (GAMESS) from .fch(k) file (Gaussian)

! If '-gvb [npair]' is specified, the orbitals in active space will be permuted to the
!  order of Gamess. In this case, you must specify the argument [npair] in command line.

! After '-gvb [npair]' is specified, you can also append '-open [nopen0]' if your system
!  is truely open shell, i.e., doublet, triplet and so on.

! If '-uhf' is specified, the .fch file must include the Beta MO part

! The order of Cartesian functions in Gaussian can be acquired by adding keyword
! 'pop=reg' in .gjf file. In GAMESS output files(.out, .gms), the order of Cartesian
! functions is printed without extra keyword needed.

program main
 implicit none
 integer :: i, npair, nopen0
 integer, parameter :: iout = 6
 character(len=4) :: string, gvb_or_uhf
 character(len=240) :: fchname

 npair = 0; nopen0 = 0
 string = ' '; gvb_or_uhf = ' '
 i = iargc()

 select case(i)
 case(1,2,3,5)
 case default
  write(iout,'(/,1X,A)') 'ERROR in subroutine fch2inp: wrong command line arguments!'
  write(iout,'(1X,A)') 'Example 1 (R(O)HF, CAS): fch2inp a.fch'
  write(iout,'(1X,A)') 'Example 2 (UHF)        : fch2inp a.fch -uhf'
  write(iout,'(1X,A)') 'Example 3 (GVB)        : fch2inp a.fch -gvb [npair]'
  write(iout,'(1X,A,/)') 'Example 4 (ROGVB)      : fch2inp a.fch -gvb [npair] -open [nopen]'
  stop
 end select

 call getarg(1,fchname)

 if(i > 1) then
  call getarg(2,gvb_or_uhf)
  if(.not. (gvb_or_uhf=='-gvb' .or. gvb_or_uhf== '-uhf')) then
   write(iout,'(A)') 'ERROR in subroutine fch2inp: the 2nd argument in command line is wrong!'
   write(iout,'(A)') "It must be '-gvb' or '-uhf'."
   stop
  end if

  if(gvb_or_uhf == '-gvb') then
   if(i == 3) then
    call getarg(3,string)
    read(string,*) npair
   else
    call getarg(3,string)
    read(string,*) npair
    call getarg(5,string)
    read(string,*) nopen0
   end if
  end if
 end if

 call fch2inp(fchname, gvb_or_uhf, npair, nopen0)
 stop
end program main

! generate .inp file (GAMESS) from .fch(k) file (Gaussian)
subroutine fch2inp(fchname, gvb_or_uhf, npair, nopen0)
 use fch_content
 implicit none
 integer :: i, j, k, m, n, n1, n2, nline, nleft, fid
 integer :: ncore   ! the number of core MOs
 integer :: n10fmark, n15gmark, n21hmark
 integer, intent(in) :: npair, nopen0
 ! here nopen0 used since nopen already used in module fch_content
 integer, allocatable :: f_mark(:), g_mark(:), h_mark(:)
 ! mark the index where f, g, h functions begin
 integer, allocatable :: order(:)
 ! six types of angular momentum
 character(len=1) :: str = ' '
 character(len=1), parameter :: am_type(-1:5) = ['L','S','P','D','F','G','H']
 character(len=1), parameter :: am_type1(0:5) = ['s','p','d','f','g','h']
 real(kind=8), allocatable :: temp_coeff(:,:), open_coeff(:,:)
 character(len=4), intent(in) :: gvb_or_uhf
 character(len=240), intent(in) :: fchname
 character(len=240) :: inpname = ' '
 logical :: uhf, ecp, sph

 i = INDEX(fchname,'.fch',back=.true.)
 if(i == 0) then
  write(iout,'(A)') "ERROR in subroutine fch2inp: input filename does not&
                     & contain '.fch' suffix!"
  write(iout,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if
 inpname = fchname(1:i-1)//'.inp'

 uhf = .false.
 if(gvb_or_uhf == '-uhf') uhf = .true.
 call read_fch(fchname, uhf) ! read content in .fch(k) file

 ncore = na - npair - nopen0
 ! arrays eigen_e_a and eigen_e_b is useless for GAMESS inp file
 deallocate(eigen_e_a)
 if(allocated(eigen_e_b)) deallocate(eigen_e_b)

 ! check if any spherical functions
 sph = .false.
 if( ANY(shell_type < -1) ) then
  write(iout,'(A)') 'Warning in subroutine fch2inp: spherical functions detected&
                   & in file '//TRIM(fchname)//'.'
  write(iout,'(A)') "GAMESS supports only Cartesian functions. You need to add&
                  & '6D 10F' keywords in Gaussian."
  sph = .true.
 end if

 ! in GAMESS inp format file, the contr_coeff_sp is zero when there is no 'L'/'SP'
 if(.not. allocated(contr_coeff_sp)) allocate(contr_coeff_sp(nprim), source=0.0d0)

 ! create the GAMESS .inp file and print the keywords information
 ecp = .false.
 if(LenNCZ > 0) ecp = .true.
 call creat_gamess_inp_head(inpname, charge, mult, ncore, npair, nopen0, nif, gvb_or_uhf, ecp)
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')

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
   write(fid,'(2X,I2,3(2X,ES15.8))') j-k, prim_exp(j), contr_coeff(j), contr_coeff_sp(j)
  end do ! for j

  if(i == ncontr) then
   write(fid,'(/)',advance='no')
  else ! i < ncontr
   if(shell2atom_map(i+1) == shell2atom_map(i)+1) write(fid,'(/)',advance='no')
  end if

  k = k + n
 end do ! for i

 deallocate(ielem, elem, coor)
 deallocate(shell2atom_map, prim_per_shell, prim_exp, contr_coeff, contr_coeff_sp)
 write(fid,'(A)') ' $END'

 ! print ECP/PP (if any) into the .inp file
 if(ecp) then
  write(fid,'(A)') ' $ECP'

  do i = 1, natom, 1
   if(LPSkip(i) /= 0) then
    write(fid,'(I0,A9)') i, ' ECP-NONE'
    cycle
   else
    write(fid,'(I0,A8,2X,I3,2X,I2)') i,'-ECP GEN', INT(RNFroz(i)), LMax(i)
    str = am_type1(LMax(i))
    do j = 1, 10, 1
     n1 = KFirst(i,j); n2 = KLast(i,j)
     if(n1 == 0) exit
     if(j == 1) then
      write(fid,'(I0,5X,A)') n2-n1+1, '----- '//str//'-ul potential -----'
     else
      write(fid,'(I0,5X,A)') n2-n1+1, '----- '//am_type1(j-2)//'-'//str//' potential -----'
     end if
     do n = n1, n2, 1
      write(fid,'(3X,ES15.8,3X,I2,2X,ES15.8)') CLP(n), NLP(n), ZLP(n)
     end do ! for n
    end do ! for j
   end if
  end do ! for i

  write(fid,'(A)') ' $END'
  deallocate(KFirst, KLast, Lmax, LPSkip, NLP, RNFroz, CLP, ZLP)
 end if

 if(sph) then
  write(iout,'(A)') 'Only the basis set and ECP(if any) are printed.'
  deallocate(shell_type, alpha_coeff)
  stop
 end if

 ! record the indices of Cartesian f, g and h functions
 n10fmark = 0
 n15gmark = 0
 n21hmark = 0
 allocate(f_mark(ncontr), source=0)
 allocate(g_mark(ncontr), source=0)
 allocate(h_mark(ncontr), source=0)

 nbf = 0
 do i = 1, ncontr, 1
  if(shell_type(i) == 0) then       !'S'
   nbf = nbf + 1
  else if(shell_type(i) == 1) then  !'P'
   nbf = nbf + 3
  else if(shell_type(i) == -1) then !'SP' or 'L'
   nbf = nbf + 4
  else if(shell_type(i) == 2) then  ! 6D
   nbf = nbf + 6
  else if(shell_type(i) == 3) then  ! 10F
   n10fmark = n10fmark + 1
   f_mark(n10fmark) = nbf + 1
   nbf = nbf + 10
  else if(shell_type(i) == 4) then  ! 15G
   n15gmark = n15gmark + 1
   g_mark(n15gmark) = nbf + 1
   nbf = nbf + 15
  else if(shell_type(i) == 5) then  ! 21H
   n21hmark = n21hmark + 1
   h_mark(n21hmark) = nbf + 1
   nbf = nbf + 21
  end if
 end do
 deallocate(shell_type)
 ! done recording

 ! adjust the order of Cartesian f, g, h functions
 do i = 1,n10fmark,1
  call fch2inp_permute_10f(nif,alpha_coeff(f_mark(i)+3:f_mark(i)+8,:))
 end do
 do i = 1,n15gmark,1
  call fch2inp_permute_15g(nif,alpha_coeff(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1,n21hmark,1
  call fch2inp_permute_21h(nif,alpha_coeff(h_mark(i):h_mark(i)+20,:))
 end do
 if(gvb_or_uhf == '-uhf') then
  do i = 1,n10fmark,1
   call fch2inp_permute_10f(nif,beta_coeff(f_mark(i)+3:f_mark(i)+8,:))
  end do
  do i = 1,n15gmark,1
   call fch2inp_permute_15g(nif,beta_coeff(g_mark(i):g_mark(i)+14,:))
  end do
  do i = 1,n21hmark,1
   call fch2inp_permute_21h(nif,beta_coeff(h_mark(i):h_mark(i)+20,:))
  end do
 end if
 deallocate(f_mark, g_mark, h_mark)
 ! adjustment finished

 ! if active orbitals in GAMESS order are required, permute them
 if(gvb_or_uhf == '-gvb') then
  if(npair > 1) then
   allocate(order(2*npair), source=0)
   allocate(temp_coeff(nbf,2*npair), source=0.0d0)

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

 ! output MOs to the .inp file
 write(fid,'(1X,A4)') '$VEC'

 ! print Alpha MOs
 nline = nbf/5
 nleft = nbf - 5*nline
 do i = 1,nif, 1
  k = MOD(i,100)
  do j = 1, nline, 1
   write(fid,'(I2,I3,5ES15.8)') k, MOD(j,1000), alpha_coeff(5*j-4:5*j,i)
  end do
  if(nleft > 0) then
   write(fid,'(I2,I3,5ES15.8)') k, MOD(j,1000), alpha_coeff(5*j-4:nbf,i)
  end if
 end do
 deallocate(alpha_coeff)

 ! print Beta MOs (if any)
 if(gvb_or_uhf == '-uhf') then
  do i = 1,nif, 1
   k = MOD(i,100)
   do j = 1, nline, 1
    write(fid,'(I2,I3,5ES15.8)') k, MOD(j,1000), beta_coeff(5*j-4:5*j,i)
   end do ! for j
   if(nleft > 0) then
    write(fid,'(I2,I3,5ES15.8)') k, MOD(j,1000), beta_coeff(5*j-4:nbf,i)
   end if
  end do ! for i
  deallocate(beta_coeff)
 end if

 write(fid,'(1X,A4)') '$END'
 close(fid)
 return
end subroutine fch2inp

! create the GAMESS .inp file and print the keywords information
subroutine creat_gamess_inp_head(inpname, charge, mult, ncore, npair, nopen, nif, gvb_or_uhf, ecp)
 implicit none
 integer :: fid
 integer, intent(in) :: charge, mult, ncore, npair, nopen, nif
 character(len=2), allocatable :: ideg(:)
 character(len=4), intent(in) :: gvb_or_uhf
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: ecp

 open(newunit=fid,file=TRIM(inpname),status='replace')
 write(fid,'(A)',advance='no') ' $CONTRL SCFTYP='

 select case(gvb_or_uhf)
 case('-uhf')
  write(fid,'(A)',advance='no') 'UHF'
 case('-gvb')
  write(fid,'(A)',advance='no') 'GVB'
 case default
  write(fid,'(A)',advance='no') 'RHF'
 end select

 write(fid,'(2(A,I0))',advance='no') ' RUNTYP=ENERGY ICHARG=',charge,' MULT=',mult
 write(fid,'(A)',advance='no') ' NOSYM=1 ICUT=10'
 if(ecp) write(fid,'(A)',advance='no') ' PP=READ'
 select case(gvb_or_uhf)
 case('-gvb')
  write(fid,'(/,A)') '  MAXIT=500 $END'
 case default
  write(fid,'(/,A)') '  MAXIT=200 $END'
 end select
 write(fid,'(A)') ' $SYSTEM MWORDS=2000 $END'

 select case(gvb_or_uhf)
 case('-gvb')
  write(fid,'(2(A,I0))',advance='no') ' $SCF NCO=',ncore,' NPAIR=',npair
  if(nopen > 0) then
   write(fid,'(A,I0,A)',advance='no') ' NSETO=',nopen,' NO(1)=1'
   if(nopen > 1) then
    allocate(ideg(nopen-1)); ideg = ',1'
    write(fid,'(10A2)',advance='no') ideg
    deallocate(ideg)
   end if
  end if
  write(fid,'(A)') ' DIRSCF=.TRUE. $END'
 case default
  write(fid,'(A)') ' $SCF DIRSCF=.TRUE. $END'
 end select

 write(fid,'(A,I0,A)') ' $GUESS GUESS=MOREAD NORB=', nif, ' $END'
 write(fid,'(A)') ' $DATA'
 write(fid,'(A)') 'GAMESS inp format file produced by MOKIT'
 write(fid,'(A)') 'C1   1'
 close(fid)
 return
end subroutine creat_gamess_inp_head

subroutine fch2inp_permute_10f(nif,coeff)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 integer, parameter :: order(6) = (/2, 3, 1, 6, 4, 5/)
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
 return
end subroutine fch2inp_permute_10f

subroutine fch2inp_permute_15g(nif,coeff)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 integer, parameter :: order(15) = (/15, 5, 1, 14, 13, 9, 4, 6, 2, 12, 10, 3, 11, 8, 7/)
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
 return
end subroutine fch2inp_permute_15g

subroutine fch2inp_permute_21h(nif,coeff)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 integer, parameter :: order(21) = (/21, 6, 1, 20, 19, 11, 5, 7, 2, 18, 16, &
                                     15, 4, 12, 3, 17, 10, 8, 14, 13, 9/)
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
 return
end subroutine fch2inp_permute_21h

