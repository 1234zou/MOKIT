! written by jxzou at 20180817: adjust the orders of d,f,g, etc functions in
!  .fch(k) file, to the orders in Molcas and OpenMolcas

! updated by jxzou at 20180826: support Cartesian H functions
! updated by jxzou at 20190226: change data to parameter when declare constant arrays
! updated by jxzou at 20190624: able to transform UHF orbitals
! updated by jxzou at 20200504: combined with rwwfn.f90
! updated by jxzou at 20200809: combined with util_wrapper.f90
! updated by jxzou at 20210407: remove '-uhf', add automatic determination

program main
 use util_wrapper, only: formchk, fch2inp_wrap
 implicit none
 integer :: k, narg
 character(len=4) :: str4
 character(len=29), parameter :: error_warn = 'ERROR in program fch2inporb: '
 character(len=30) :: dftname
 character(len=240) :: fchname, gms_inp, molcas_inp
 character(len=260) :: buf
 logical :: sph, prt_no

 narg = iargc()
 if(narg<1 .or. narg>3) then
  write(6,'(/,1X,A)') error_warn//'wrong command line arguments!'
  write(6,'(A)')  ' Example 1 (R(O)HF/UHF/CAS): fch2inporb a.fch'
  write(6,'(A)')  ' Example 2 (CAS NO)        : fch2inporb a.fch -no'
  write(6,'(A,/)')' Example 3 (DFT)           : fch2inporb a.fch -dft "B3LYP"'
  stop
 end if

 prt_no = .false.; dftname = ' '; fchname = ' '; str4 = ' '
 gms_inp = ' '; molcas_inp = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 k = LEN_TRIM(fchname)
 if(fchname(k-3:k) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:k-3)//'fch'
 end if

 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)

 call find_specified_suffix(fchname, '.fch', k)
 gms_inp = fchname(1:k-1)//'.inp'
 molcas_inp = fchname(1:k-1)//'.input'

 select case(narg)
 case(1) ! do nothing
 case(2)
  call getarg(2, str4)
  if(TRIM(str4) == '-no') then
   prt_no = .true.
  else
   write(6,'(/,A)') error_warn//'the 2nd argument is wrong!'
   write(6,'(A)') '`-no` is expected.'
   stop
  end if
 case(3)
  call getarg(2, str4)
  if(TRIM(str4) == '-dft') then
   call getarg(3, dftname)
  else
   write(6,'(/,A)') error_warn//'the 2nd argument is wrong!'
   write(6,'(A)') '`-dft` is expected.'
   stop
  end if
 case default
  write(6,'(/,A)') error_warn//'wrong command line arguments!'
  stop
 end select

 ! ->prt_no, <-sph
 call fch2inporb(fchname, prt_no, sph)
 call fch2inp_wrap(fchname, .false., 0, 0, .true., .false., .false.)

 buf = 'bas_gms2molcas '//TRIM(gms_inp)
 if(sph) buf = TRIM(buf)//' -sph'
 call run_command(TRIM(buf), .false., .false.)
 call delete_file(gms_inp)
 if(LEN_TRIM(dftname) > 0) call add_dftname2molcas_inp(molcas_inp, dftname)
end program main

! nbf: the number of basis functions
! nif: the number of independent functions, i.e., the number of MOs

! read the MOs in .fch(k) file and adjust its d,f,g, etc. functions order
!  of Gaussian to that of Molcas
subroutine fch2inporb(fchname, prt_no, sph)
 implicit none
 integer :: i, j, k, m, length, icart, orbid
 integer :: nalpha, nbeta, nbf, nif, nbf0, nbf1
 integer :: n6dmark, n10fmark, n15gmark, n21hmark, n28imark
 integer :: n5dmark, n7fmark, n9gmark, n11hmark, n13imark
 integer(kind=4) :: getpid, hostnm
 integer, allocatable :: shell_type(:), shell2atom_map(:)
 integer, allocatable :: idx(:), d_mark(:), f_mark(:), g_mark(:), h_mark(:), &
  i_mark(:)
 character(len=8) :: hostname
 character(len=24) :: data_string
 character(len=240), intent(in) :: fchname
 character(len=240) :: orbfile ! .INPORB file of (Open)Molcas
 real(kind=8), allocatable :: coeff(:,:), coeff0(:,:), norm(:), occ_num(:)
 logical :: uhf
 logical, intent(in) :: prt_no
 logical, intent(out) :: sph

 sph = .true.
 call find_icart_in_fch(fchname, .false., icart)
 if(icart == 2) sph = .false.

 call read_na_and_nb_from_fch(fchname, nalpha, nbeta)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf

 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF

 ! read MO Coefficients
 if(uhf) then
  allocate(coeff(nbf,2*nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_fch(fchname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 else   ! not UHF
  if(prt_no) then
   allocate(occ_num(nif), coeff(nbf,nif))
   call read_eigenvalues_from_fch(fchname, nif, 'a', occ_num)
   call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)
  else
   allocate(coeff(nbf,nif))
   call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)
  end if
 end if

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

! first we adjust the basis functions in each MO according to the Shell to atom map
 ! 1) split the 'L' into 'S' and 'P', this is to ensure that D comes after L functions
 call split_L_func(k, shell_type, shell2atom_map, length)

 ! 2) sort the shell_type and shell2atom_map by ascending order
 ! MOs will be adjusted accordingly
 call sort_shell_and_mo(length, shell_type, shell2atom_map, nbf, nif, coeff)
! adjust done

! then we adjust the basis functions in each MO according to the type of basis functions
 k = length  ! update k
 allocate(d_mark(k), f_mark(k), g_mark(k), h_mark(k), i_mark(k))

 ! adjust the order of d, f, etc. functions
 allocate(idx(nbf))
 forall(i = 1:nbf) idx(i) = i ! initialization
 allocate(coeff0(nbf,nif), source=coeff)
 if(sph) then
  call read_mark_from_shltyp_sph(k, shell_type, n5dmark, n7fmark, n9gmark, &
                 n11hmark, n13imark, d_mark, f_mark, g_mark, h_mark, i_mark)
  call fch2inporb_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, d_mark, &
                              f_mark, g_mark, h_mark, nbf, idx)
  forall(i=1:nif, j=1:nbf) coeff(j,i) = coeff0(idx(j),i)
 else
  allocate(norm(nbf), source=1d0)
  call read_mark_from_shltyp_cart(k, shell_type, n6dmark, n10fmark, n15gmark, &
                    n21hmark, n28imark, d_mark, f_mark, g_mark, h_mark, i_mark)
  call fch2inporb_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, &
                           d_mark, f_mark, g_mark, h_mark, nbf, idx, norm)
  forall(i=1:nif, j=1:nbf) coeff(j,i) = coeff0(idx(j),i)*norm(j)
  deallocate(norm)
 end if
 deallocate(d_mark, f_mark, g_mark, h_mark, i_mark, coeff0, idx)
! adjustment finished

! move the 2nd, 3rd, ... Zeta basis functions forward
 i = 0
 nbf = 0
 do while(i < k)
  i = i + 1
  j = shell2atom_map(i)
  m = shell_type(i)
  nbf1 = nbf
  select case(m)
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   nbf = nbf + 3
  case(-1)   ! SP or L
   nbf = nbf + 4
  case(-2)   ! 5D
   nbf = nbf + 5
  case( 2)   ! 6D
   nbf = nbf + 6
  case(-3)   ! 7F
   nbf = nbf + 7
  case( 3)   ! 10F
   nbf = nbf + 10
  case(-4)   ! 9G
   nbf = nbf + 9
  case( 4)   ! 15G
   nbf = nbf + 15
  case(-5)   ! 11H
   nbf = nbf + 11
  case( 5)   ! 21H
   nbf = nbf + 21
  end select
  if(m == 0) cycle

  length = 1
  do while(i < k)
   i = i + 1
   if(shell_type(i) /= m) exit
   if(shell2atom_map(i) /= j) exit
   length = length + 1
  end do
  if(i < k) i = i - 1
  if(length > 1) then
   call zeta_mv_forwd(nbf1, m, length, nbf0, nif, coeff)
   nbf = nbf1 + length*(nbf-nbf1)
  end if
 end do
 deallocate(shell_type, shell2atom_map)
! move done

! print MOs into INPORB
 i = SCAN(fchname, '.', BACK=.true.)
 orbfile = fchname(1:i-1)//'.INPORB'
 open(newunit=orbid,file=TRIM(orbfile),status='replace')

 write(orbid,'(A,/,A)') '#INPORB 2.2', '#INFO'

 if(uhf) then
  write(orbid,'(A)') '* UHF orbitals'
  write(orbid,'(A)') '       1       1       4'
  nif = nif/2   ! change to original size
 else
  if(prt_no) then
   write(orbid,'(A)') '* natural orbitals'
   write(orbid,'(A)') '       0       1       0'
  else
   write(orbid,'(A)') '* SCF orbitals'
   write(orbid,'(A)') '       0       1       2'
  end if
 end if

 write(orbid,'(2X,I6,/,2X,I6)') nbf0, nif

 hostname = ' '
 data_string = ' '
 i = getpid()
 j = hostnm(hostname)
 call fdate(data_string)
 write(orbid,'(A,I6,A)') '*BC:HOST '//TRIM(hostname)//' PID', i, ' DATE '//&
                         TRIM(data_string)//' Generated by MOKIT'

 write(orbid,'(A)') '#ORB'
 do i = 1, nif, 1
  write(orbid,'(A14,I5)') '* ORBITAL    1', i
  write(orbid,'(5(1X,ES21.14))') (coeff(j,i),j=1,nbf0)
 end do

 if(uhf) then
  write(orbid,'(A)') '#UORB'
  do i = 1, nif, 1
   write(orbid,'(A14,I5)') '* ORBITAL    1', i
   write(orbid,'(5(1X,ES21.14))') (coeff(j,i+nif),j=1,nbf0)
  end do
 end if
 deallocate(coeff)

 if(uhf) then
  allocate(occ_num(2*nif), source=0d0)
  occ_num(1:nalpha) = 1d0
  occ_num(nif+1:nif+nbeta) = 1d0
  write(orbid,'(A,/,A)') '#OCC', '* OCCUPATION NUMBERS'
  write(orbid,'(5(1X,ES21.14))') (occ_num(i),i=1,nif)
  write(orbid,'(A,/,A)') '#UOCC', '* Beta OCCUPATION NUMBERS'
  write(orbid,'(5(1X,ES21.14))') (occ_num(i),i=nif+1,2*nif)
 else
  if(prt_no) then
   write(orbid,'(A,/,A)') '#OCC', '* OCCUPATION NUMBERS'
   write(orbid,'(5(1X,ES21.14))') (occ_num(i),i=1,nif)
  else
   allocate(occ_num(nif), source=0d0)
   if(nbeta < nalpha) then
    occ_num(1:nbeta) = 2d0
    occ_num(nbeta+1:nalpha) = 1d0
   else if(nbeta == nalpha) then
    occ_num(1:nalpha) = 2d0
   end if
   write(orbid,'(A,/,A)') '#OCC', '* OCCUPATION NUMBERS'
   write(orbid,'(5(1X,ES21.14))') (occ_num(i),i=1,nif)
  end if
 end if

 write(orbid,'(A,/,A)') '#OCHR', '* OCCUPATION NUMBERS (HUMAN-READABLE)'
 write(orbid,'(10(1X,F7.4))') (occ_num(i),i=1,nif)
 if(uhf) then
  write(orbid,'(A,/,A)') '#UOCHR', '* Beta OCCUPATION NUMBERS (HUMAN-READABLE)'
  write(orbid,'(10(1X,F7.4))') (occ_num(i),i=nif+1,2*nif)
 end if

 deallocate(occ_num)
 close(orbid)
! print done
end subroutine fch2inporb

! move the 2nd, 3rd, ... Zeta basis functions forward
subroutine zeta_mv_forwd(i0, shell_type, length, nbf, nif, coeff2)
 implicit none
 integer i, j, k
 integer, intent(in) :: i0, shell_type, length, nbf, nif
 integer, parameter :: num0(-5:5) = [11, 9, 7, 5, 0, 0, 3, 6, 10, 15, 21]
 !                                   11H 9G 7F 5D L  S 3P 6D 10F 15G 21H
 real(kind=8), intent(inout) :: coeff2(nbf,nif)
 real(kind=8), allocatable :: coeff(:,:)

 if(length == 1) return

 if(shell_type==0 .or. shell_type==-1) then
  write(6,'(/,A)') 'ERROR in subroutine zeta_mv_forwd: this element of shell_ty&
                   &pe is 0 or -1.'
  write(6,'(2(A,I0))') 'shell_type=', shell_type, ', length=', length
  stop
 end if

 coeff = coeff2
 k = num0(shell_type)
 do i = 1, k, 1
  do j = 1, length, 1
   coeff(i0+j+(i-1)*length,:) = coeff2(i0+i+(j-1)*k,:)
  end do ! for j
 end do ! for i

 coeff2 = coeff
 deallocate(coeff)
end subroutine zeta_mv_forwd

subroutine add_dftname2molcas_inp(inpname, dftname)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=4) :: str4
 character(len=7) :: str7
 character(len=30), intent(in) :: dftname
 character(len=44), parameter :: error_warn = 'ERROR in subroutine add_dftname2&
                                              &molcas_inp: '
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 call find_specified_suffix(inpname, '.input', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  if(LEN_TRIM(buf) > 0) then
   buf = ADJUSTL(buf)
   str7 = buf(1:7)
   call upper(str7)
   if(str7 == "&SEWARD") exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&SEWARD" cannot be located in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(fid1,'(A)') 'Grid input'
 write(fid1,'(A)') ' grid= ultrafine'
 write(fid1,'(A)') 'End of grid input'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  if(LEN_TRIM(buf) > 0) then
   buf = ADJUSTL(buf)
   str4 = buf(1:4)
   call upper(str4)
   if(str4 == "&SCF") exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"&SCF" cannot be located in'
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if
 write(fid1,'(A)') 'KSDFT= '//TRIM(dftname)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid1)
 close(fid,status='delete')
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine add_dftname2molcas_inp

