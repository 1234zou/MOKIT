! written by jxzou at 20191030: transform basis sets in GAMESS format to (Open)Molcas format
! updated by jxzou at 20200326: move lower angular momentum (occurred 2nd time) forward
!   This special case happens at, e.g. def2TZVP for Co. The MOs moved forward have been
!   considered in fch2inporb. Here the basis sets should be adjusted correspondingly.
! updated by jxzou at 20200402: ECP/PP supported
! updated by jxzou at 20201110: dealing with the extra primitive duplicate, e.g.
!   cc-pVDZ for Co. The original way is to both extend column and row. Now the
!   trick is to only extend column, thus will not generate the primitive duplicate.
!   The primitive duplicate is not allowed in relativistic calculations.
! updated by jxzou at 20201209: move some subroutines to read_gms_inp.f90
! updated by jxzou at 20210207: add contraction strings into '.....'

! Note: Currently isotopes are not tested.
program main
 implicit none
 integer :: i
 character(len=4) :: str
 character(len=240) :: fname
 ! fname: input file contains basis sets and Cartesian coordinates in GAMESS format
 logical :: spherical

 i = iargc()
 if(i<1 .or. i>2) then
  write(6,'(/,A,/)') 'Example1: bas_gms2molcas a.inp (generate an a.input file)'
  write(6,'(A,/)') "Example2: bas_gms2molcas a.inp -sph (without 'Cartesian all')"
  stop
 end if

 str = ' '
 fname = ' '
 spherical = .false.
 call getarg(1, fname)
 call require_file_exist(fname)

 if(i == 2) then
  call getarg(2,str)

  if(str == '-sph') then
   spherical = .true.
  else
   write(6,'(A)') 'ERROR in subroutine bas_gms2molcas: wrong command line arguments!'
   write(6,'(A)') "The 2nd argument can only be '-sph'. But got '"//str//"'"
   stop
  end if
 end if

 call bas_gms2molcas(fname, spherical)
 stop
end program main

! Transform the basis sets in GAMESS format to those in (Open)Molcas format
subroutine bas_gms2molcas(fort7, spherical)
 use pg, only: natom, ram, ntimes, elem, coor, prim_gau, all_ecp, ecp_exist
 implicit none
 integer :: i, nline, rc, rel, charge, mult, fid1, fid2
 character(len=240), intent(in) :: fort7
 character(len=240) :: buf, input
 ! input is the (Open)Molcas .input file
 character(len=1) :: stype
 character(len=21) :: str1, str2
 logical, intent(in) :: spherical
 logical :: uhf, ghf, X2C

 ! initialization
 buf = ' '
 input = ' '

 i = index(fort7, '.', back=.true.)
 input = fort7(1:i-1)//'.input'

 call read_natom_from_gms_inp(fort7, natom)
 allocate(ram(natom), elem(natom), coor(3,natom), ntimes(natom))
 call read_elem_nuc_coor_from_gms_inp(fort7, natom, elem, ram, coor)
 ! ram cannot be deallocated here since subroutine prt_prim_gau will use it

 call calc_ntimes(natom, elem, ntimes)
 call read_charge_and_mult_from_gms_inp(fort7, charge, mult, uhf, ghf, ecp_exist)
 call read_all_ecp_from_gms_inp(fort7)

 ! find the $DATA section
 open(newunit=fid1,file=TRIM(fort7),status='old',position='rewind')
 do while(.true.)
  read(fid1,'(A)',iostat=rc) buf
  if(rc /= 0) exit
  if(buf(2:2) == '$') then
   call upper(buf(3:6))
   if(buf(3:6) == 'DATA') exit
  end if
 end do ! for while

 if(rc /= 0) then
  write(6,'(A)') 'ERROR in subroutine bas_gms2molcas: No $DATA section found&
                   & in file '//TRIM(fort7)//'.'
  close(fid1)
  stop
 end if

 ! skip 2 lines: the Title line and the Point Group line
 read(fid1,'(A)') buf
 read(fid1,'(A)') buf

 open(newunit=fid2,file=TRIM(input),status='replace')
 write(fid2,'(A)') "&GATEWAY"
! write(fid2,'(A)') 'Expert'

 ! initialization: clear all primitive gaussians
 call clear_prim_gau()

 do i = 1, natom, 1
  read(fid1,'(A)',iostat=rc) buf
  if(rc /= 0) exit
  ! 'buf' contains the element, ram and coordinates
  write(fid2,'(A)') 'Basis set'

  ! deal with primitive gaussians
  do while(.true.)
   read(fid1,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit

   read(buf,*) stype, nline
   call read_prim_gau(stype, nline, fid1)
  end do ! for while

  if(all_ecp(i)%ecp) then
   write(fid2,'(A)',advance='no') TRIM(elem(i))//'.ECP..'
  else
   write(fid2,'(A)',advance='no') TRIM(elem(i))//'...'
  end if
  call gen_contracted_string(prim_gau(:)%nline,prim_gau(:)%ncol,str1,str2)
  write(fid2,'(A)') TRIM(str1)//'.'//TRIM(str2)//'.   / inline'

  ! print basis sets and ECP/PP (if any) of this atom in Molcas format
  call prt_prim_gau(i, fid2)
  write(fid2,'(A,I0,3(2X,F15.8),A)') TRIM(elem(i)), ntimes(i), &
                                     coor(1:3,i),'   Angstrom'
  if(.not. spherical) write(fid2,'(A)') 'Cartesian all'
  write(fid2,'(A)') 'End of basis set'

  ! clear all primitive gaussians for next cycle
  call clear_prim_gau()
 end do ! for i

 close(fid1)
 ! now ram can be deallocated
 deallocate(ram, elem, ntimes, coor, all_ecp)

 if(rc /= 0) then
  write(6,'(A)') "ERROR in subroutine bas_gms2molcas: it seems the '$DATA'&
                   & has no corresponding '$END'."
  write(6,'(A)') 'Incomplete file '//TRIM(fort7)
  close(fid2,status='delete')
  stop
 end if

 call check_X2C_in_gms_inp(fort7, X2C)
 if(X2C) write(fid2,'(A)') 'RX2C'
 write(fid2,'(/,A)') "&SEWARD"

 call check_DKH_in_gms_inp(fort7, rel)
 select case(rel)
 case(-2) ! nothing
 case(-1) ! RESC
  write(6,'(A)') 'ERROR in subroutine bas_gms2molcas: RESC keywords detected.'
  write(6,'(A)') 'But RESC is not supported in (Open)Molcas.'
  stop
 case(0,1,2,4)  ! DKH0/1/2/4
  if(.not. X2C) write(fid2,'(A,I2.2,A)') 'Relativistic = R',rel,'O'
  !if(.not. X2C) write(fid2,'(A,I2.2,A)') 'R',rel,'O'
 case default
  write(6,'(A)') 'ERROR in subroutine bas_gms2molcas: rel out of range!'
  write(6,'(A,I0)') 'rel=', rel
  close(fid2,status='delete')
  stop
 end select

 write(fid2,'(/,A)') "&SCF"
 if(uhf) write(fid2,'(A)') 'UHF'
 write(fid2,'(A,I0)') 'Charge = ', charge
 write(fid2,'(A,I0)') 'Spin = ', mult
 i = INDEX(fort7, '.', back=.true.)
 input = fort7(1:i-1)//'.INPORB'
 write(fid2,'(A)') 'FILEORB = '//TRIM(input)

 if((.not.uhf) .and. (mult/=1)) then ! ROHF
  write(fid2,'(A)') "* ROHF is not directly supported in (Open)Molcas. You&
                    & need to write &RASSCF module to realize ROHF."
 end if

 close(fid2)
 return
end subroutine bas_gms2molcas

