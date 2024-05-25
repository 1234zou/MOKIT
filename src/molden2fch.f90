! written by jxzou at 20240125: generate Gaussian .fch file (from scratch) using
! a given .molden file
! TODO: extended to more formats of .molden files; support 6D 10F.

! The expansion matrices of spherical harmonic functions to Cartesian functions.
! The permutation matrices are included in these expansion matrices, so there
! is no need to permute MO coefficients after multiplying these expansion matrices.
module molden_sph2cart
 implicit none
 real(kind=8), parameter :: r1=1d0/DSQRT(12d0), r2=1d0/DSQRT(3d0), &
  r3=1d0/DSQRT(15d0), r4=DSQRT(0.03d0), r5=1d0/DSQRT(40d0), r6=DSQRT(0.005d0),&
  r7=DSQRT(0.08d0), r8=DSQRT(0.05d0), r9=1d0/DSQRT(24d0), r10=DSQRT(0.075d0),&
  r11=DSQRT(3d0/2240d0), r12=1d0/DSQRT(105d0), r13=1d0/DSQRT(2176d0), &
  r14=1d0/DSQRT(136d0), r15=DSQRT(3d0/392d0), r16=DSQRT(2d0/147d0), &
  r17=DSQRT(3d0/1960d0), r18=1d0/DSQRT(588d0), r19=DSQRT(3d0/245d0), &
  r20=1d0/DSQRT(336d0), r21=3d0/DSQRT(980d0), r22=1d0/DSQRT(168d0), &
  r23=DSQRT(3d0/280d0), r24=1d0/DSQRT(84d0), r25=1d0/DSQRT(192d0), &
  r26=3d0/DSQRT(560d0), r27=DSQRT(0.15d0), r28=DSQRT(0.4d0), r29=DSQRT(0.375d0),&
  r30=DSQRT(3d0/560d0), r31=DSQRT(3d0/35d0), r32=DSQRT(3d0/56d0), &
  r33=DSQRT(2d0/21d0), r34=DSQRT(3d0/28d0), r35=DSQRT(3d0/7d0), r36=DSQRT(0.1875d0),&
  r37=DSQRT(5d0/1344d0), r38=DSQRT(5d0/336d0), r39=DSQRT(5d0/189d0), &
  r40=1d0/DSQRT(945d0), r41=1d0/DSQRT(4032d0), r42=1d0/DSQRT(1008d0), &
  r43=1d0/DSQRT(28d0), r44=1d0/DSQRT(63d0), r45=1d0/12d0, r46=1d0/6d0, &
  r47=1d0/3d0, r48=1d0/DSQRT(3456d0), r49=1d0/DSQRT(864d0), r50=1d0/DSQRT(54d0),&
  r51=1d0/DSQRT(384d0), r52=1d0/DSQRT(6d0), r53=1d0/DSQRT(1920d0), &
  r54=DSQRT(5d0/96d0), r55=DSQRT(5d0/384d0), r56=DSQRT(0.75d0), &
  r57=DSQRT(0.625d0), r58=DSQRT(1.125d0), r59=DSQRT(0.45d0), r60=DSQRT(1.2d0),&
  r61=DSQRT(0.140625d0), r62=DSQRT(0.3125d0), r63=DSQRT(0.546875d0), &
  r64=DSQRT(5d0/28d0), r65=DSQRT(1.25d0), r66=DSQRT(45d0/56d0), r67=DSQRT(10d0/7d0),&
  r68=DSQRT(27d0/560d0), r69=DSQRT(1.6875d0), r70=DSQRT(27d0/35d0), &
  r71=DSQRT(27d0/28d0), r72=DSQRT(9d0/56d0), r73=DSQRT(9d0/7d0)
 real(kind=8), parameter :: rd(6,5) = RESHAPE([-r1, -r1, r2, 0d0, 0d0, 0d0,&
  0d0, 0d0, 0d0, 0d0,  r2, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,  r2, &
  0.5d0, -0.5d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,  r2, 0d0, 0d0], [6,5])
 real(kind=8), parameter :: rd1(6,5) = RESHAPE([-r1, -r1, r2, 0d0, 0d0, 0d0,&
  0d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0, &
  0.5d0, -0.5d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0], [6,5])
 real(kind=8), parameter :: rd2(6,5) = RESHAPE([-0.5d0,-0.5d0,1d0,0d0,0d0,0d0,&
  0d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0, &
  r56,-r56, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0], [6,5])
 real(kind=8), parameter :: &
  rf(10,7) = RESHAPE([0d0, 0d0,  r3, 0d0, 0d0, -r4, 0d0, 0d0, -r4, 0d0,&
                      -r5, 0d0, 0d0, -r6, 0d0, 0d0,  r7, 0d0, 0d0, 0d0,&
                      0d0, -r5, 0d0, 0d0, -r6, 0d0, 0d0,  r7, 0d0, 0d0,&
                      0d0, 0d0, 0d0, 0d0, 0d0,  r8, 0d0, 0d0, -r8, 0d0,&
                      0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,  r3,&
                       r9, 0d0, 0d0,-r10, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,&
                      0d0, -r9, 0d0, 0d0, r10, 0d0, 0d0, 0d0, 0d0, 0d0],[10,7])
 real(kind=8), parameter :: &
  rf1(10,7) = RESHAPE([0d0, 0d0,  r3, 0d0, 0d0,-r27, 0d0, 0d0,-r27, 0d0,&
                       -r5, 0d0, 0d0, -r5, 0d0, 0d0, r28, 0d0, 0d0, 0d0,&
                       0d0, -r5, 0d0, 0d0, -r5, 0d0, 0d0, r28, 0d0, 0d0,&
                       0d0, 0d0, 0d0, 0d0, 0d0,0.5d0, 0d0, 0d0,-0.5d0, 0d0,&
                       0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0,&
                        r9, 0d0, 0d0,-r29, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,&
                       0d0, -r9, 0d0, 0d0, r29, 0d0, 0d0, 0d0, 0d0, 0d0],[10,7])
 real(kind=8), parameter :: &
  rf2(10,7) = RESHAPE([0d0, 0d0, 1d0, 0d0, 0d0,-r59, 0d0, 0d0,-r59, 0d0,&
                      -r29, 0d0, 0d0,-r10, 0d0, 0d0, r60, 0d0, 0d0, 0d0,&
                       0d0,-r29, 0d0, 0d0,-r10, 0d0, 0d0, r60, 0d0, 0d0,&
                       0d0, 0d0, 0d0, 0d0, 0d0, r56, 0d0, 0d0,-r56, 0d0,&
                       0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0,&
                       r57, 0d0, 0d0,-r58, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,&
                       0d0,-r57, 0d0, 0d0, r58, 0d0, 0d0, 0d0, 0d0, 0d0],[10,7])
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
 real(kind=8), parameter :: &
  rg1(15,9) = RESHAPE([r11, r11, r12, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r30,-r31,-r31, 0d0, 0d0,0d0,&
                       0d0, 0d0, 0d0, 0d0,-r32, 0d0, 0d0, r33, 0d0, 0d0, 0d0, 0d0, 0d0,-r32,0d0,&
                       0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r32, 0d0, r33, 0d0, 0d0, 0d0,-r32, 0d0,0d0,&
                      -r20, r20, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r34,-r34, 0d0, 0d0,0d0,&
                       0d0, 0d0, 0d0,-r24, 0d0,-r24, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,r35,&
                       0d0, 0d0, 0d0, 0d0,  r9, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r29,0d0,&
                       0d0, 0d0, 0d0, 0d0, 0d0, 0d0, -r9, 0d0, 0d0, 0d0, 0d0, 0d0, r29, 0d0,0d0,&
                       r25, r25, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r36, 0d0, 0d0, 0d0, 0d0,0d0,&
                       0d0, 0d0, 0d0,  r1, 0d0, -r1, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,0d0],[15,9])
 real(kind=8), parameter :: &
  rg2(15,9) = RESHAPE([r61, r61, 1d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r68,-r70,-r70, 0d0, 0d0,0d0,&
                       0d0, 0d0, 0d0, 0d0,-r66, 0d0, 0d0, r67, 0d0, 0d0, 0d0, 0d0, 0d0,-r72,0d0,&
                       0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r66, 0d0, r67, 0d0, 0d0, 0d0,-r72, 0d0,0d0,&
                      -r62, r62, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r71,-r71, 0d0, 0d0,0d0,&
                       0d0, 0d0, 0d0,-r64, 0d0,-r64, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,r73,&
                       0d0, 0d0, 0d0, 0d0, r57, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r58,0d0,&
                       0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r57, 0d0, 0d0, 0d0, 0d0, 0d0, r58, 0d0,0d0,&
                       r63, r63, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r69, 0d0, 0d0, 0d0, 0d0,0d0,&
                       0d0, 0d0, 0d0, r65, 0d0,-r65, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,0d0],[15,9])
 !real(kind=8), parameter :: &
 ! rh1(21,11) = RESHAPE([0d0, 0d0, r37, 0d0, 0d0, 0d0, 0d0, r38, 0d0,-r39,0d0, 0d0, 0d0,0d0,0d0,0d0,r37, 0d0,-r39,0d0,r40,&
 !                       r41, 0d0, 0d0, r42, 0d0,-r43, 0d0, 0d0, 0d0, 0d0,r41, 0d0,-r43,0d0,r44,0d0,0d0, 0d0, 0d0,0d0,0d0,&
 !                       0d0, r41, 0d0, 0d0, 0d0, 0d0, r42, 0d0,-r43, 0d0,0d0, 0d0, 0d0,0d0,0d0,r41,0d0,-r43, 0d0,r44,0d0,&
 !                       0d0, 0d0,-r45, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r46,0d0, 0d0, 0d0,0d0,0d0,0d0,r45, 0d0,-r46,0d0,0d0,&
 !                       0d0, 0d0, 0d0, 0d0,-r46, 0d0, 0d0, 0d0, 0d0, 0d0,0d0,-r46, 0d0,r47,0d0,0d0,0d0, 0d0, 0d0,0d0,0d0,&
 !                      -r48, 0d0, 0d0, r49, 0d0, r50, 0d0, 0d0, 0d0, 0d0,r51, 0d0,-r52,0d0,0d0,0d0,0d0, 0d0, 0d0,0d0,0d0,&
 !                       0d0,-r51, 0d0, 0d0, 0d0, 0d0,-r49, 0d0, r52, 0d0,0d0, 0d0, 0d0,0d0,0d0,r48,0d0,-r50, 0d0,0d0,0d0,&
 !                       0d0, 0d0, r25, 0d0, 0d0, 0d0, 0d0,-r36, 0d0, 0d0,0d0, 0d0, 0d0,0d0,0d0,0d0,r25, 0d0, 0d0,0d0,0d0,&
 !                       0d0, 0d0, 0d0, 0d0,  r1, 0d0, 0d0, 0d0, 0d0, 0d0,0d0, -r1, 0d0,0d0,0d0,0d0,0d0, 0d0, 0d0,0d0,0d0,&
 !                       r53, 0d0, 0d0,-r54, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,r55, 0d0, 0d0,0d0,0d0,0d0,0d0, 0d0, 0d0,0d0,0d0,&
 !                       0d0, r55, 0d0, 0d0, 0d0, 0d0,-r54, 0d0, 0d0, 0d0,0d0, 0d0, 0d0,0d0,0d0,r53,0d0, 0d0, 0d0,0d0,0d0],[21,11])
end module molden_sph2cart

program main
 implicit none
 integer :: i, j, iprog
 integer, parameter :: n = 13
 character(len=3) :: str3
 character(len=7) :: str7
 character(len=7), parameter :: sprog(n) = ['-bagel ','-cfour ','-cp2k  ', &
  '-dalton','-et    ','-molcas','-molpro','-mrcc  ','-nwchem','-orca  ', &
  '-psi4  ','-pyscf ','-tm    ']
 character(len=240) :: molden
 logical :: alive, natorb

 i = iargc()
 if(i<2 .or. i>3) then
  write(6,'(/,A)') ' ERROR in subroutine molden2fch: wrong command line arguments!'
  do j = 1, n, 1
   write(6,'(A,I2,A)') ' Example',j,': molden2fch a.molden '//TRIM(sprog(j))
  end do ! for j
  write(6,'(/,A)') " If you are transferring NOs, you can append a '-no' argume&
                   &nt, e.g."
  write(6,'(A,/)') '            molden2fch a.molden -orca -no'
  stop
 end if

 str3 = ' '; str7 = ' '; molden = ' '
 call getarg(1, molden)
 inquire(file=TRIM(molden),exist=alive)
 if(.not. alive) then
  write(6,'(/,A)') 'ERROR in program molden2fch: input molden file does not exi&
                   &st!'
  write(6,'(A)') 'Filename='//TRIM(molden)
  stop
 end if

 call getarg(2, str7)
 str7 = ADJUSTL(str7)
 do iprog = 1, n, 1
  if(TRIM(sprog(iprog)) == TRIM(str7)) exit
 end do ! for i

 if(iprog == n+1) then
  write(6,'(/,A)') 'ERROR in program molden2fch: this molden format cannot be r&
                   &ecognized.'
  write(6,'(A)') 'molden format strongly depends on quantum chemistry programs.&
                 & So this utility'
  write(6,'(A)') 'requires you to specify a program name.'
  stop
 end if

 natorb = .false.
 if(i == 3) then
  call getarg(3, str3)
  if(str3 /= '-no') then
   write(6,'(/,A)') "ERROR in program molden2fch: the 3rd argument can only be &
                    &'-no'."
   stop
  end if
  natorb = .true.
 end if

 call molden2fch(molden, iprog, natorb)
end program main

subroutine molden2fch(molden, iprog, natorb)
 use fch_content
 use mkl_content, only: natom1=>natom, ncontr1=>ncontr, shell_type1=>shell_type,&
  shl2atm1=>shl2atm, all_pg, un_normalized_all_pg, merge_s_and_p_into_sp
 use molden_sph2cart, only: rd, rd1, rd2, rf, rf1, rf2, rg, rg1, rg2
 use phys_cons, only: au2cm_1
 implicit none
 integer :: i, j, k, ne0, ne, nbf1, nif1
 integer :: ndmark, nfmark, ngmark, nhmark, nimark
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:), i_mark(:)
 integer, intent(in) :: iprog ! label of quantum chemistry packages
 real(kind=8), allocatable :: coeff(:,:), tmp_coeff(:,:), occ_a(:), occ_b(:), &
  frcoor(:,:), freq(:), all_mode(:,:)
 character(len=3), parameter :: fname = 'mos'
 character(len=240) :: fchname
 character(len=240), intent(in) :: molden
 logical :: sph, has_sp, has_freq
 logical, intent(in) :: natorb

 i = INDEX(molden, '.molden', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in program molden2fch: illegal filename. No '.molden'&
                   & string found."
  write(6,'(A)') 'Please rename it to xxx.molden.'
  stop
 end if

 call check_sph_in_molden(molden, sph)

 fchname = molden(1:i-1)//'.fch'
 call read_natom_from_molden(molden, natom)
 allocate(ielem(natom), coor(3,natom))
 call read_nuc_and_coor_from_molden(molden, natom, ielem, coor)
 ! It seems that the number of core electrons are not taken into consideration
 ! in the 3rd column of the [ATOMS] section.
 call read_ncontr_from_molden(molden, ncontr)
 allocate(shell_type(ncontr), shell2atom_map(ncontr))
 call read_shltyp_and_shl2atm_from_molden(molden, ncontr, shell_type, &
                                          shell2atom_map)
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
 if(natorb) eigen_e_a = occ_a

 ! occupation numbers in CFOUR RHF molden need to be *2
 if(iprog==2 .and. (.not.is_uhf) .and. (.not.natorb)) occ_a = 2d0*occ_a

 if(is_uhf) then
  allocate(eigen_e_b(nif), occ_b(nif))
  call read_mo_from_molden(molden, nbf, nif, 'b', coeff(:,nif+1:nif1), eigen_e_b, &
                           occ_b)
 end if

 allocate(f_mark(ncontr), g_mark(ncontr), h_mark(ncontr), i_mark(ncontr))
 select case(iprog)
 case(1,4,5,8,9,11,12)
  write(6,'(/,A)') 'Not implemented yet.'
  stop
 case(2) ! CFOUR
  nbf1 = nbf - COUNT(shell_type==-2) - 3*COUNT(shell_type==-3) - &
         6*COUNT(shell_type==-4) - 10*COUNT(shell_type==-5)
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
    nbf = nbf + 1; j = j + 1
   case( 1) ! P
    tmp_coeff(nbf+1:nbf+3,:) = coeff(j+1:j+3,:)
    nbf = nbf + 3; j = j + 3
   case(-1) ! L
    tmp_coeff(nbf+1:nbf+4,:) = coeff(j+1:j+4,:)
    nbf = nbf + 4; j = j + 4
   case(2) ! 6D -> 5D
    call solve_multi_lin_eqs(6,5,rd1,nif1,coeff(j+1:j+6,:),tmp_coeff(nbf+1:nbf+5,:))
    nbf = nbf + 5; j= j + 6
   case(3) ! 10F -> 7F
    call solve_multi_lin_eqs(10,7,rf1,nif1,coeff(j+1:j+10,:),tmp_coeff(nbf+1:nbf+7,:))
    nbf = nbf + 7; j = j + 10
   case(4) ! 15G -> 9G
    call solve_multi_lin_eqs(15,9,rg1,nif1,coeff(j+1:j+15,:),tmp_coeff(nbf+1:nbf+9,:))
    nbf = nbf + 9; j = j + 15
   !case(5) ! 21G -> 11G
   ! call solve_multi_lin_eqs(21,11,rh1,nif1,coeff(j+1:j+21,:),tmp_coeff(nbf+1:nbf+11,:))
   ! nbf = nbf + 11; j = j + 21
   case default
    write(6,'(/,A)') 'ERROR in subroutine molden2fch: the CFOUR-type molden fil&
                     &e is problematic'
    write(6,'(A)') 'when there are basis functions >=h. This is probably a bug &
                   &of CFOUR.'
    stop
   end select
  end do ! for i
  deallocate(d_mark, shell_type1, coeff)
  allocate(coeff(nbf1,nif1), source=tmp_coeff)
  deallocate(tmp_coeff)
 case(3) ! CP2K, nothing to do
 case(6) ! (Open)Molcas, nothing to do
  if(ANY(shell_type < -4)) then
    write(6,'(/,A)') 'ERROR in subroutine molden2fch: OpenMolcas cannot generat&
                     &e .molden file'
    write(6,'(A)') 'for angular momentum >=h. This molden file is suspicious.'
    stop
  end if
 case(7) ! Molpro
  nbf1 = nbf - COUNT(shell_type==-2) - 3*COUNT(shell_type==-3) - &
         6*COUNT(shell_type==-4) - 10*COUNT(shell_type==-5)
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
    nbf = nbf + 1; j = j + 1
   case( 1) ! P
    tmp_coeff(nbf+1:nbf+3,:) = coeff(j+1:j+3,:)
    nbf = nbf + 3; j = j + 3
   case(-1) ! L
    tmp_coeff(nbf+1:nbf+4,:) = coeff(j+1:j+4,:)
    nbf = nbf + 4; j = j + 4
   case(2) ! 6D -> 5D
    call solve_multi_lin_eqs(6,5,rd2,nif1,coeff(j+1:j+6,:),tmp_coeff(nbf+1:nbf+5,:))
    nbf = nbf + 5; j= j + 6
   case(3) ! 10F -> 7F
    call solve_multi_lin_eqs(10,7,rf2,nif1,coeff(j+1:j+10,:),tmp_coeff(nbf+1:nbf+7,:))
    nbf = nbf + 7; j = j + 10
   case(4) ! 15G -> 9G
    call solve_multi_lin_eqs(15,9,rg2,nif1,coeff(j+1:j+15,:),tmp_coeff(nbf+1:nbf+9,:))
    nbf = nbf + 9; j = j + 15
   case default
    write(6,'(/,A)') 'ERROR in subroutine molden2fch: Molpro cannot generate .m&
                     &olden file for'
    write(6,'(A)') 'angular momentum >=h. This molden file is suspicious.'
    stop
   end select
  end do ! for i
  deallocate(d_mark, shell_type1, coeff)
  allocate(coeff(nbf1,nif1), source=tmp_coeff)
  deallocate(tmp_coeff)
 case(10) ! ORCA
  ! find F+3, G+3, H+3, I+3 functions, multiply them by -1
  call read_bas_mark_from_shltyp(ncontr, shell_type, nfmark, ngmark, nhmark, &
                                 nimark, f_mark, g_mark, h_mark, i_mark)
  call update_mo_using_bas_mark(nbf, nif1, nfmark, ngmark, nhmark, nimark, &
                                ncontr, f_mark, g_mark, h_mark, i_mark, coeff)
 case(13) ! Turbomole
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
    nbf = nbf + 1; j = j + 1
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
    write(6,'(/,A)') 'ERROR in subroutine molden2fch: tm2molden does not suppor&
                     &t angular momentum'
    write(6,'(A)') 'higher than g functions. This molden file is suspicious.'
    stop
   end select
  end do ! for i
  deallocate(d_mark, shell_type1, coeff)
  allocate(coeff(nbf1,nif1), source=tmp_coeff)
  deallocate(tmp_coeff)
 case default
  write(6,'(/,A)') 'ERROR in subroutine molden2fch: iprog out of range!'
  write(6,'(A,I0)') 'iprog=', iprog
  stop
 end select

 deallocate(f_mark, g_mark, h_mark, i_mark)
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

 if(iprog == 13) then ! Turbomole
  call check_freq_in_molden(molden, has_freq)
  if(has_freq) then
   k = 3*natom
   allocate(frcoor(3,natom), freq(k), all_mode(k,k))
   call read_freq_from_molden(molden, natom, nmode, frcoor, freq, all_mode)
   frcoor = DABS(frcoor - coor)
   if(SUM(frcoor) > 1d-3) then
    write(6,'(/,A)') 'ERROR in subroutine molden2fch: Cartesian coordinates are &
                     &inconsistent in'
    write(6,'(A)') 'file '//TRIM(molden)
    stop
   end if
   deallocate(frcoor)
   ! It seems that the Hessian matrix cannot be regenerated here
   !do i = 1, k, 1
   ! j = i/3
   ! if(i - 3*j > 0) j = j + 1
   ! all_mode(i,:) = all_mode(i,:)*DSQRT(ram(ielem(j)))
   !end do ! for i
   !do i = 1, k, 1
   ! all_mode(:,i) = all_mode(:,i)/DSQRT(DOT_product(all_mode(:,i),all_mode(:,i)))
   !end do ! for i
   allocate(norm_mode(k,nmode), source=all_mode(:,k-nmode+1:k))
   deallocate(all_mode)
   allocate(vibe2(14*nmode), source=0d0)
   vibe2(1:nmode) = freq(k-nmode+1:k)
   deallocate(freq)
  end if
 end if

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
 logical :: coor_au

 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf(1:8),'[Atoms]')>0 .or. buf(1:7)=='[ATOMS]') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_nuc_and_coor_from_molden: no [ATOM&
                   &S] found in'
  write(6,'(A)') 'file '//TRIM(molden)
  close(fid)
  stop
 end if

 coor_au = .true.
 call upper(buf)
 if(INDEX(buf,'ANGS') > 0) coor_au = .false.

 do i = 1, natom, 1
  read(fid,*) str, j, nuc(i), coor(:,i) ! coor in Bohr
 end do ! for i

 close(fid)
 if(coor_au) coor = coor*Bohr_const
end subroutine read_nuc_and_coor_from_molden

! find the array size of shell_type from a given .molden file
subroutine read_ncontr_from_molden(molden, ncontr)
 implicit none
 integer :: i, j, fid
 integer, intent(out) :: ncontr
 integer, external :: detect_ncol_in_buf
 character(len=240) :: buf
 character(len=240), intent(in) :: molden

 ncontr = 0
 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(INDEX(buf(1:6),'[GTO]') > 0) exit
 end do ! for while

 do while(.true.)
  read(fid,'(A)') buf
  if(INDEX(buf(1:2),'[') > 0) exit
  i = detect_ncol_in_buf(buf)
  if(i == 2) then
   buf = ADJUSTL(buf)
   j = IACHAR(buf(1:1))
   if((j>64 .and. j<91) .or. (j>96 .and. j<123)) ncontr = ncontr + 1
  else if(i == 3) then
   ncontr = ncontr + 1
  end if
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

 ang = ' '; shltyp = 0; shl2atm = 0
 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(INDEX(buf(1:6),'[GTO]') > 0) exit
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
  case('i')
   shltyp(i) = -6
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
    if(INDEX(buf(1:2),'[') > 0) exit
    read(fid,'(A)') buf
    exit
   else
    k = detect_ncol_in_buf(buf)
    if(k == 3) exit
    if(k == 2) then
     buf = ADJUSTL(buf)
     j = IACHAR(buf(1:1))
     if((j>64 .and. j<91) .or. (j>96 .and. j<123)) exit
    end if
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
  if(INDEX(buf(1:10),'Occup=') > 0) exit
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
  if(buf(1:1) == '[') exit

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
  if(INDEX(buf(1:5),'[MO]') > 0) exit
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
  if(INDEX(buf(1:5),'[MO]') > 0) exit
 end do ! for while

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  if(INDEX(buf(1:2),'[') > 0) exit
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
  if(INDEX(buf(1:5),'[MO]') > 0) exit
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
 use mkl_content, only: natom, shl2atm, all_pg, del_zero_coeff_in_prim_gau
 implicit none
 integer :: i, j, k, m, nc, nline, fid
 integer, external :: detect_ncol_in_buf
 character(len=1) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: molden

 if(natom == 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_all_pg_from_molden: natom = 0.'
  stop
 end if

 if(.not. allocated(shl2atm)) then
  write(6,'(/,A)') 'ERROR in subroutine read_all_pg_from_molden: array shl2atm &
                   &is not allocated.'
  write(6,'(A)') 'Internal inconsistency.'
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
  if(INDEX(buf(1:6),'[GTO]') > 0) exit
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
    k = detect_ncol_in_buf(buf)
    if(k == 3) exit
    if(k == 2) then
     buf = ADJUSTL(buf)
     m = IACHAR(buf(1:1))
     if((m>64 .and. m<91) .or. (m>96 .and. m<123)) exit
    end if
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
  if(INDEX(buf(1:6),'[GTO]') > 0) exit
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

   call del_zero_coeff_in_prim_gau(all_pg(i)%prim_gau(j))
  end do ! for j
  read(fid,'(A)') buf
 end do ! for i

 close(fid)
 nline = all_pg(1)%prim_gau(1)%nline
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

subroutine align_mo(nbf, nif1, coeff, nbf1, tmp_coeff)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nif1, nbf1
 real(kind=8) :: r1(6), r2(6)
 real(kind=8), intent(in) :: coeff(nbf, nif1)
 real(kind=8), intent(inout) :: tmp_coeff(nbf1,nif1)

 do i = 1, nif1, 1
  r1 = coeff(1:6,i) - tmp_coeff(1:6,i)
  r2 = coeff(1:6,i) + tmp_coeff(1:6,i)
  if(SUM(DABS(r1)) > SUM(DABS(r2))) tmp_coeff(:,i) = -tmp_coeff(:,i)
 end do ! for i
end subroutine align_mo

! check spherical harmonic or Cartesian functions are used in a .molden file
subroutine check_sph_in_molden(molden, sph)
 implicit none
 integer :: fid
 character(len=5) :: buf ! no need to be a long string
 character(len=240), intent(in) :: molden
 logical, intent(out) :: sph

 sph = .true.
 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(INDEX(buf,'[5D') > 0) exit
  if(INDEX(buf,'[6D') > 0) then
   sph = .false.
   exit
  end if
  if(INDEX(buf,'[MO]') > 0) exit
 end do ! for while

 close(fid)
end subroutine check_sph_in_molden

! check whether there is frequency data in the molden file
subroutine check_freq_in_molden(molden, has_freq)
 implicit none
 integer :: i, fid
 character(len=6) :: str6
 character(len=240), intent(in) :: molden
 logical, intent(out) :: has_freq

 has_freq = .false.
 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) str6
  if(i /= 0) exit
  if(str6(1:6) == '[FREQ]') then
   has_freq = .true.
   exit
  end if
 end do ! for while

 close(fid)
end subroutine check_freq_in_molden

! read frequency data from the molden file
subroutine read_freq_from_molden(molden, natom, nmode, coor, e, ev)
 use phys_cons, only: Bohr_const
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 integer, intent(out) :: nmode
 real(kind=8), parameter :: thres = 1d-3
 real(kind=8), intent(out) :: coor(3,natom), e(3*natom), ev(3*natom,3*natom)
 character(len=2) :: elem
 character(len=240) :: buf
 character(len=240), intent(in) :: molden

 coor = 0d0; e = 0d0; ev = 0d0; k = 3*natom
 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:6) == '[FREQ]') exit
 end do ! for while
 read(fid,*) e

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:10) == '[FR-COORD]') exit
 end do ! for while

 do i = 1, natom, 1
  read(fid,*) elem, coor(:,i)
 end do ! for i

 read(fid,'(A)') buf ! [FR-NORM-COORD]
 do i = 1, k, 1
  read(fid,'(A)') buf ! vibration    i
  read(fid,*) ev(:,i)
 end do ! for i

 close(fid)
 coor = coor*Bohr_const

 do i = k, 1, -1
  if(DABS(e(i)) < thres) exit
 end do ! for i
 nmode = k - i
end subroutine read_freq_from_molden

