! written by jxzou at 20240125: generate Gaussian .fch file (from scratch) using
! a given .molden file
! TODO: extended to more formats of .molden files; support 6D 10F.

! The expansion matrices of spherical harmonic functions to Cartesian functions.
! The permutation matrices are included in these expansion matrices, so there
! is no need to permute MO coefficients after multiplying these expansion matrices.
module molden_sph2cart
 implicit none
 real(kind=8), parameter :: r1=1d0/DSQRT(12d0), r2=1d0/DSQRT(3d0), &
  r3=1d0/DSQRT(15d0), r4=DSQRT(3d0)/10d0, r5=1d0/DSQRT(40d0), r6=1d0/DSQRT(200d0),&
  r7=DSQRT(2d0)/5d0, r8=1d0/DSQRT(20d0), r9=1d0/DSQRT(24d0), r10=DSQRT(0.075d0),&
  r11=DSQRT(3d0/2240d0), r12=1d0/DSQRT(105d0), r13=1d0/DSQRT(2176d0), &
  r14=1d0/DSQRT(136d0), r15=DSQRT(3d0/392d0), r16=DSQRT(2d0/147d0), &
  r17=DSQRT(3d0/1960d0), r18=1d0/DSQRT(588d0), r19=DSQRT(3d0/245d0), &
  r20=1d0/DSQRT(336d0), r21=3d0/DSQRT(980d0), r22=1d0/DSQRT(168d0), &
  r23=DSQRT(3d0/280d0), r24=1d0/DSQRT(84d0), r25=1d0/DSQRT(192d0), &
  r26=3d0/DSQRT(560d0)
 real(kind=8), parameter :: rd(6,5) = RESHAPE([-r1, -r1, r2, 0d0, 0d0, 0d0,&
  0d0, 0d0, 0d0, 0d0,  r2, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,  r2, &
  0.5d0, -0.5d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,  r2, 0d0, 0d0], [6,5])
 real(kind=8), parameter :: &
  rf(10,7) = RESHAPE([0d0, 0d0,  r3, 0d0, 0d0, -r4, 0d0, 0d0, -r4, 0d0,&
                      -r5, 0d0, 0d0, -r6, 0d0, 0d0,  r7, 0d0, 0d0, 0d0,&
                      0d0, -r5, 0d0, 0d0, -r6, 0d0, 0d0,  r7, 0d0, 0d0,&
                      0d0, 0d0, 0d0, 0d0, 0d0,  r8, 0d0, 0d0, -r8, 0d0,&
                      0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,  r3,&
                       r9, 0d0, 0d0,-r10, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,&
                      0d0, -r9, 0d0, 0d0, r10, 0d0, 0d0, 0d0, 0d0, 0d0],[10,7])
 real(kind=8), parameter :: &
  rg(15,9) = RESHAPE([r11, r11, r12, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r13,-r14,-r14, 0d0, 0d0,0d0,&
                      0d0, 0d0, 0d0, 0d0,-r15, 0d0, 0d0, r16, 0d0, 0d0, 0d0, 0d0, 0d0,-r17,0d0,&
                      0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r15, 0d0, r16, 0d0, 0d0, 0d0,-r17, 0d0,0d0,&
                     -r20, r20, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r21,-r21, 0d0, 0d0,0d0,&
                      0d0, 0d0, 0d0,-r18, 0d0,-r18, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,r19,&
                      0d0, 0d0, 0d0, 0d0, r22, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r23,0d0,&
                      0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r22, 0d0, 0d0, 0d0, 0d0, 0d0, r23, 0d0,0d0,&
                      r25, r25, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r26, 0d0, 0d0, 0d0, 0d0,0d0,&
                      0d0, 0d0, 0d0, r24, 0d0,-r24, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,0d0],[15,9])
end module molden_sph2cart

program main
 implicit none
 integer :: i, iprog
 character(len=10) :: sprog
 character(len=240) :: molden
 logical :: alive

 i = iargc()
 if(i /= 2) then
  write(6,'(/,A)') ' ERROR in subroutine molden2fch: wrong command line arguments!'
  write(6,'(A)')   ' Example1 (R(O)HF, UHF, CAS): molden2fch a.molden -orca'
  write(6,'(A,/)') ' Example2 (R(O)HF, UHF, CAS): molden2fch a.molden -tm'
  stop
 end if

 iprog = 0; molden = ' '; sprog = ' '
 call getarg(1, molden)
 call getarg(2, sprog)

 inquire(file=TRIM(molden),exist=alive)
 if(.not. alive) then
  write(6,'(/,A)') 'ERROR in program molden2fch: input molden file does not exi&
                   &st!'
  write(6,'(A)') 'filename='//TRIM(molden)
  stop
 end if

 select case(TRIM(sprog))
 case('-orca')
  iprog = 1
 case('-tm')
  iprog = 2
 case default
  write(6,'(/,A)') 'ERROR in program fch2tm: this molden format cannot be recog&
                   &nized.'
  write(6,'(A)') 'molden format strongly depends on quantum chemistry programs.&
                 & So this utility'
  write(6,'(A)') 'requires you to specify a program.'
  stop
 end select

 call molden2fch(molden, iprog)
end program main

subroutine molden2fch(molden, iprog)
 use fch_content
 use mkl_content, only: natom1=>natom, ncontr1=>ncontr, shell_type1=>shell_type,&
  shl2atm1=>shl2atm, all_pg, un_normalized_all_pg, merge_s_and_p_into_sp
 use molden_sph2cart, only: rd, rf, rg
 implicit none
 integer :: i, j, ne0, ne, nbf1, nif1
 integer :: ndmark, nfmark, ngmark, nhmark
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 integer, intent(in) :: iprog ! 1/2 for ORCA/Turbomole
 real(kind=8), allocatable :: coeff(:,:), tmp_coeff(:,:), occ_a(:), occ_b(:)
 character(len=3), parameter :: fname = 'mos'
 character(len=240) :: fchname
 character(len=240), intent(in) :: molden
 logical :: has_sp

 i = INDEX(molden, '.molden', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in program molden2fch: illegal filename. No '.molden'&
                   & string found."
  stop
 end if

 fchname = molden(1:i-1)//'.fch'
 call read_natom_from_molden(molden, natom)
 allocate(ielem(natom), coor(3,natom))
 call read_nuc_and_coor_from_molden(molden, natom, ielem, coor)

 call read_ncontr_from_molden(molden, ncontr)
 allocate(shell_type(ncontr), shell2atom_map(ncontr))
 call read_shltyp_and_shl2atm_from_molden(molden, ncontr, shell_type, shell2atom_map)
 natom1 = natom
 ncontr1 = ncontr
 allocate(shell_type1(ncontr), source=shell_type)
 allocate(shl2atm1(ncontr), source=shell2atom_map)
 call read_all_pg_from_molden(molden)
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

 call read_nbf_and_nif_from_molden(molden, nbf, nif)
 call check_uhf_in_molden(molden, is_uhf)
 nif1 = nif
 if(is_uhf) nif1 = 2*nif
 allocate(coeff(nbf,nif1), eigen_e_a(nif), occ_a(nif))
 call read_mo_from_molden(molden, nbf, nif, 'a', coeff(:,1:nif), eigen_e_a, &
                          occ_a)
 if(is_uhf) then
  allocate(eigen_e_b(nif), occ_b(nif))
  call read_mo_from_molden(molden, nbf, nif, 'b', coeff(:,nif+1:nif1), eigen_e_b, &
                           occ_b)
 end if

 allocate(f_mark(ncontr), g_mark(ncontr), h_mark(ncontr))
 select case(iprog)
 case(1) ! ORCA
  ! find F+3, G+3 and H+3 functions, multiply them by -1
  call read_bas_mark_from_shltyp(ncontr, shell_type, nfmark, ngmark, nhmark, &
                                 f_mark, g_mark, h_mark)
  call update_mo_using_bas_mark(nbf, nif1, nfmark, ngmark, nhmark, f_mark, &
                                g_mark, h_mark, coeff)
 case(2) ! Turbomole
  ! It seems that Cartesian-type functions (6D,10F) are used in any .molden file
  ! generated by Turbomole, no matter that (5D,7F) or (6D,10F) is used in the
  ! file `control`. This is very similar to GAMESS. When spherical harmonic
  ! functions are used, MO coefficients in .molden are expanded from spherical
  ! harmonic functions, and they can be transformed back.
  ! [6D,10F,15G] - [5D,7F,9G] = [1,3,6]
  nbf1 = nbf - COUNT(shell_type==-2) - 3*COUNT(shell_type==-3) - &
         6*COUNT(shell_type==-4)
  allocate(d_mark(ncontr))
  allocate(shell_type1(ncontr), source=shell_type)
  forall(i = 1:ncontr, shell_type1(i)<-1) shell_type1(i) = -shell_type1(i)
  call read_mark_from_shltyp_cart(ncontr, shell_type1, ndmark, nfmark, ngmark, &
                                  nhmark, d_mark, f_mark, g_mark, h_mark)
  allocate(tmp_coeff(nbf1,nif1), source=0d0)
  nbf = 0; j = 0
  do i = 1, ncontr, 1
   select case(shell_type1(i))
   case( 0) ! S
    tmp_coeff(nbf+1,:) = coeff(j+1,:)
    nbf = nbf + 1; j= j + 1
   case( 1) ! P
    tmp_coeff(nbf+1:nbf+3,:) = coeff(j+1:j+3,:)
    nbf = nbf + 3; j = j + 3
   case(-1) ! L
    tmp_coeff(nbf+1:nbf+4,:) = coeff(j+1:j+4,:)
    nbf = nbf + 4; j = j + 4
   case(2) ! 6D -> 5D
    call solve_multi_lin_eqs(6,5,rd,nif1,coeff(j+1:j+6,:),tmp_coeff(nbf+1:nbf+5,:))
    nbf = nbf + 5; j= j + 6
   case(3) ! 10F -> 7F
    call solve_multi_lin_eqs(10,7,rf,nif1,coeff(j+1:j+10,:),tmp_coeff(nbf+1:nbf+7,:))
    nbf = nbf + 7; j = j + 10
   case(4) ! 15G -> 9G
    call solve_multi_lin_eqs(15,9,rg,nif1,coeff(j+1:j+15,:),tmp_coeff(nbf+1:nbf+9,:))
    nbf = nbf + 9; j = j + 15
   case default
    write(6,'(/,A)') 'ERROR in subroutine fch2tm: tm2molden does not support an&
                     &gular momentum'
    write(6,'(A)') 'higher than g functions. This molden file is suspicious.'
    stop
   end select
  end do ! for i
  deallocate(d_mark, shell_type1, coeff)
  allocate(coeff(nbf1,nif1), source=tmp_coeff)
  deallocate(tmp_coeff)
 case default
  write(6,'(/,A)') 'ERROR in subroutine fch2tm: iprog out of range!'
  write(6,'(A,I0)') 'iprog=', iprog
  stop
 end select

 deallocate(f_mark, g_mark, h_mark)
 allocate(alpha_coeff(nbf,nif), source=coeff(:,1:nif))
 if(is_uhf) then
  allocate(beta_coeff(nbf,nif), source=coeff(:,nif+1:nif1))
  na = NINT(SUM(occ_a))
  nb = NINT(SUM(occ_b))
  nopen = na - nb
  deallocate(occ_b)
 else
  call find_nopen_from_occ(nif, occ_a, nopen)
  na = NINT((SUM(occ_a) + DBLE(nopen))*0.5d0)
  nb = NINT((SUM(occ_a) - DBLE(nopen))*0.5d0)
 end if
 deallocate(coeff, occ_a)
 ne = na + nb
 ne0 = SUM(ielem)
 charge = ne0 - ne
 mult = nopen + 1
 write(6,'(/,A)') 'Warning from subroutine molden2fch: molden file does not inc&
                  &lude spin multiplicity.'
 write(6,'(A)') 'It will be guessed according to the occupation numbers. If the&
                & guessed spin is wrong,'
 write(6,'(A)') 'you need to modify it in file '//TRIM(fchname)

 call write_fch(fchname)
 call free_arrays_in_fch_content()
end subroutine molden2fch

! read atomic numbers and Cartesian coordinates from a given .molden file
subroutine read_nuc_and_coor_from_molden(molden, natom, nuc, coor)
 use phys_cons, only: Bohr_const
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: natom
 integer, intent(out) :: nuc(natom)
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=2) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: molden

 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:7) == '[Atoms]') exit
 end do ! for while

 do i = 1, natom, 1
  read(fid,*) str, j, nuc(i), coor(:,i) ! coor in Bohr
 end do ! for i

 close(fid)
 coor = coor*Bohr_const
end subroutine read_nuc_and_coor_from_molden

! find the array size of shell_type from a given .molden file
subroutine read_ncontr_from_molden(molden, ncontr)
 implicit none
 integer :: i, fid
 integer, intent(out) :: ncontr
 integer, external :: detect_ncol_in_buf
 character(len=240) :: buf
 character(len=240), intent(in) :: molden

 ncontr = 0
 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:5) == '[GTO]') exit
 end do ! for while

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:1) == '[') exit
  i = detect_ncol_in_buf(buf)
  if(i == 3) ncontr = ncontr + 1
 end do ! for while

 close(fid)
end subroutine read_ncontr_from_molden

! read array shell_type from a given .molden file
subroutine read_shltyp_and_shl2atm_from_molden(molden, ncontr, shltyp, shl2atm)
 implicit none
 integer :: i, j, k, fid, iatom
 integer, intent(in) :: ncontr
 integer, intent(out) :: shltyp(ncontr), shl2atm(ncontr)
 integer, external :: detect_ncol_in_buf
 character(len=1) :: ang
 character(len=240) :: buf
 character(len=240), intent(in) :: molden

 ang = ' '
 shltyp = 0; shl2atm = 0
 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:5) == '[GTO]') exit
 end do ! for while

 iatom = 1
 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, ncontr, 1
  shl2atm(i) = iatom
  read(buf,*,iostat=j) ang, k
  if(j /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine read_shltyp_and_shl2atm_from_molden: w&
                    &rong format in'
   write(6,'(A)') 'file '//TRIM(molden)
   write(6,'(A)') 'buf='//TRIM(buf)
   close(fid)
   stop
  end if

  select case(ang)
  case('s')
   shltyp(i) = 0
  case('p')
   shltyp(i) = 1
  case('d')
   shltyp(i) = -2
  case('f')
   shltyp(i) = -3
  case('g')
   shltyp(i) = -4
  case('h')
   shltyp(i) = -5
  case default
   write(6,'(/,A)') 'ERROR in subroutine read_shltyp_and_shl2atm_from_molden: u&
                    &nsupported angular'
   write(6,'(A)') 'momentum='//ang
   close(fid)
   stop
  end select

  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) then
    iatom = iatom + 1
    read(fid,'(A)') buf
    if(buf(1:1) == '[') exit
    read(fid,'(A)') buf
    exit
   else
    k = detect_ncol_in_buf(buf)
    if(k == 3) exit
   end if
  end do ! for while
 end do ! for i

 close(fid)
end subroutine read_shltyp_and_shl2atm_from_molden

! read nbf and nif from a given .molden file
subroutine read_nbf_and_nif_from_molden(molden, nbf, nif)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: molden
 logical :: uhf

 call check_uhf_in_molden(molden, uhf)
 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'Occup=') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nbf_and_nif_from_molden: key 'Occu&
                   &p=' not found in"
  write(6,'(A)') 'file '//TRIM(molden)
  close(fid)
  stop
 end if

 nbf = 0
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'Sym=')>0 .or. INDEX(buf,'Ene=')>0) exit
  if(LEN_TRIM(buf) == 0) exit
  nbf = nbf + 1
 end do ! for while

 nif = 1

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit

  do i = 1, 3
   read(fid,'(A)') buf
   if(INDEX(buf,'Occup=') > 0) exit
  end do ! for i

  do i = 1, nbf, 1
   read(fid,'(A)') buf
  end do ! for i
  nif = nif + 1
 end do ! for while

 close(fid)
 call check_uhf_in_molden(molden, uhf)
 if(uhf) nif = nif/2
end subroutine read_nbf_and_nif_from_molden

! find the number of tags/lines after [MO] and before MO coefficients
subroutine find_ntag_before_mo_in_molden(molden, ntag)
 implicit none
 integer :: fid
 integer, intent(out) :: ntag
 character(len=240) :: buf
 character(len=240), intent(in) :: molden

 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == '[MO]') exit
 end do ! for while

 do ntag = 1, 4
  read(fid,'(A)') buf
  if(INDEX(buf,'Occup=') > 0) exit
 end do

 close(fid)
end subroutine find_ntag_before_mo_in_molden

! check whether this is a UHF-type .molden file, i.e. there are two sets of MOs
subroutine check_uhf_in_molden(molden, uhf)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: molden
 logical, intent(out) :: uhf

 uhf = .false.
 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == '[MO]') exit
 end do ! for while

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  if(buf(1:1) == '[') exit
  if(INDEX(buf,'Spin=')>0 .and. INDEX(buf,'Beta')>0) then
   uhf = .true.
   exit
  end if
 end do ! for while

 close(fid)
end subroutine check_uhf_in_molden

! read MO coefficients from a given .molden file
subroutine read_mo_from_molden(molden, nbf, nif, ab, coeff, ev, occ)
 implicit none
 integer :: i, j, k, m, ntag, fid
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(out) :: coeff(nbf,nif), ev(nif), occ(nif)
 character(len=1), intent(in) :: ab
 character(len=5) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: molden
 logical :: uhf, oppo_spin

 coeff = 0d0; ev = 0d0; occ = 0d0
 call find_ntag_before_mo_in_molden(molden, ntag)
 call check_uhf_in_molden(molden, uhf)
 select case(ab)
 case('a')
  str = 'Alpha'
 case('b')
  str = 'Beta '
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_mo_from_molden: ab out of range.'
  write(6,'(A)') "Only 'a' or 'b' is allowed. But got '"//ab//"'."
  stop
 end select

 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == '[MO]') exit
 end do ! for while

 ! The beta MO coefficients are after all alpha MO coefficients in a molden
 ! generated by ORCA.
 ! The alpha/beta MO coefficients may occur alternatively in a molden generated
 ! by Turbomole.

 if(uhf) then ! UHF
  m = 0
  do i = 1, 2*nif, 1
   oppo_spin = .false.
   do j = 1, ntag, 1
    read(fid,'(A)') buf
    if(INDEX(buf,'Spin=') > 0) then
     if(INDEX(buf,TRIM(str)) == 0) oppo_spin = .true.
    else
     k = INDEX(buf, 'Ene=')
     if(k > 0) read(buf(k+4:),*) ev(m+1)
    end if
   end do ! for j

   if(oppo_spin) then
    do j = 1, nbf, 1
     read(fid,'(A)') buf
    end do ! for j
   else
    m = m + 1
    k = INDEX(buf, 'Occup=')
    if(k > 0) read(buf(k+6:),*) occ(m)
    do j = 1, nbf, 1
     read(fid,*) k, coeff(j,m)
    end do ! for j
   end if

   if(m == nif) exit
  end do ! for i

 else         ! R(O)HF
  do i = 1, nif, 1
   do j = 1, ntag, 1
    read(fid,'(A)') buf
    k = INDEX(buf, 'Ene=')
    if(k > 0) read(buf(k+4:),*) ev(i)
   end do ! for j
   k = INDEX(buf, 'Occup=')
   if(k > 0) read(buf(k+6:),*) occ(i)
   do j = 1, nbf, 1
    read(fid,*) k, coeff(j,i)
   end do ! for j
  end do ! for i
 end if

 close(fid)
end subroutine read_mo_from_molden

! read the basis set data of all atoms from .molden file
subroutine read_all_pg_from_molden(molden)
 use mkl_content, only: natom, shl2atm, all_pg
 implicit none
 integer :: i, j, k, nc, nline, fid
 integer, external :: detect_ncol_in_buf
 character(len=1) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: molden

 if(natom == 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_all_pg_from_molden: natom = 0.'
  stop
 end if

 if(.not. allocated(shl2atm)) then
  write(6,'(/,A)') 'ERROR in subroutine read_all_pg_from_mkl: array shl2atm is &
                   &not allocated.'
  stop
 end if
 allocate(all_pg(natom))

 do i = 1, natom, 1
  nc = COUNT(shl2atm == i)
  all_pg(i)%nc = nc
  allocate(all_pg(i)%prim_gau(nc))
  all_pg(i)%prim_gau(:)%ncol = 2
  all_pg(i)%prim_gau(:)%nline = 0
 end do ! for i

 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:5) == '[GTO]') exit
 end do ! for while

 do i = 1, natom, 1
  nc = all_pg(i)%nc
  read(fid,'(A)') buf
  read(fid,'(A)') buf

  do j = 1, nc, 1
   read(buf,*) str, k
   ! turn lower case to upper case
   all_pg(i)%prim_gau(j)%stype(1:1) = ACHAR(IACHAR(str)-32)

   nline = 0
   do while(.true.)
    read(fid,'(A)') buf
    if(LEN_TRIM(buf) == 0) exit
    if(detect_ncol_in_buf(buf) == 3) exit
    nline = nline + 1
   end do ! for while

   all_pg(i)%prim_gau(j)%nline = nline
   allocate(all_pg(i)%prim_gau(j)%coeff(nline,2), source=0d0)
   if(LEN_TRIM(buf) == 0) exit
  end do ! for j
 end do ! for i

 rewind(fid)
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:5) == '[GTO]') exit
 end do ! for while

 do i = 1, natom, 1
  nc = all_pg(i)%nc
  read(fid,'(A)') buf

  do j = 1, nc, 1
   read(fid,*) str, k
   nline = all_pg(i)%prim_gau(j)%nline

   do k = 1, nline, 1
    read(fid,*) all_pg(i)%prim_gau(j)%coeff(k,:)
   end do ! for k
  end do ! for j
  read(fid,'(A)') buf
 end do ! for i

 close(fid)
end subroutine read_all_pg_from_molden

subroutine find_nopen_from_occ(nif, occ, nopen)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 integer, intent(out) :: nopen
 real(kind=8), parameter :: thres = 1d-5
 real(kind=8), intent(in) :: occ(nif)

 nopen = 0
 do i = 1, nif, 1
  if(DABS(occ(i) - 1d0) < thres) nopen = nopen + 1 
 end do ! for i
end subroutine find_nopen_from_occ

