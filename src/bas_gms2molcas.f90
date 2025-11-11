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
   write(6,'(/,A)') 'ERROR in subroutine bas_gms2molcas: wrong command line arg&
                    &uments!'
   write(6,'(A)') "The 2nd argument can only be '-sph'. But got '"//str//"'"
   stop
  end if
 end if

 call bas_gms2molcas(fname, spherical)
end program main

! Transform the basis sets in GAMESS format to those in (Open)Molcas format
subroutine bas_gms2molcas(fort7, spherical)
 use pg, only: natom, nuc, ntimes, elem, coor, prim_gau, all_ecp, ecp_exist, &
  ghost
 implicit none
 integer :: i, na, nb, nline, rc, rel, charge, mult, isph, fid1, fid2
 character(len=240), intent(in) :: fort7
 character(len=240) :: buf, input
 ! input is the (Open)Molcas .input file
 character(len=1) :: stype
 character(len=21) :: str1, str2
 logical, intent(in) :: spherical
 logical :: uhf, ghf, X2C

 buf = ' '; input = ' ' ! initialization

 i = INDEX(fort7, '.', back=.true.)
 input = fort7(1:i-1)//'.input'

 call read_natom_from_gms_inp(fort7, natom)
 allocate(nuc(natom), elem(natom), coor(3,natom), ntimes(natom), ghost(natom))
 call read_elem_nuc_coor_from_gms_inp(fort7, natom, elem, nuc, coor, ghost)
 ! nuc cannot be deallocated here since subroutine prt_prim_gau will use it

 ! deal with ghost atoms
 do i = 1, natom, 1
  if(ghost(i)) then
   elem(i) = 'X '
  else
   if(elem(i) == 'Bq') then
    ghost(i) = .true.; elem(i) = 'X '
   end if
  end if
 end do ! for i

 call calc_ntimes(natom, elem, ntimes)
 call read_charge_mult_isph_from_gms_inp(fort7, charge, mult, isph, uhf, ghf, &
                                         ecp_exist)
 call read_all_ecp_from_gms_inp(fort7)

 ! find the $DATA section
 call goto_data_section_in_gms_inp(fort7, fid1)
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
  ! 'buf' contains the element, nuc and coordinates
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

  ! print basis sets and ECP/PP (if any) of this atom into .input file
  call prt_prim_gau2molcas_inp(i, fid2)
  write(fid2,'(A,I0,3(2X,F15.8),A)') TRIM(elem(i)), ntimes(i), &
                                     coor(1:3,i),'   Angstrom'
  if(.not. spherical) write(fid2,'(A)') 'Cartesian all'
  write(fid2,'(A)') 'End of basis set'

  ! clear all primitive gaussians for next cycle
  call clear_prim_gau()
 end do ! for i

 close(fid1)
 ! now nuc can be deallocated
 deallocate(nuc, elem, ntimes, coor, all_ecp, ghost)

 if(rc /= 0) then
  write(6,'(/,A)') "ERROR in subroutine bas_gms2molcas: it seems the '$DATA' ha&
                   &s no corresponding '$END'."
  write(6,'(A)') 'Incomplete file '//TRIM(fort7)
  close(fid2,status='delete')
  stop
 end if

 ! RICD becomes default when MOLCAS_NEW_DEFAULTS=YES. To control RI behavior in
 ! MOKIT, we have to write noCD here. In utilities like automr/autosr, noCD will
 ! be checked and deleted if RI is required.

 call check_X2C_in_gms_inp(fort7, X2C)
 if(X2C) write(fid2,'(A)') 'RX2C'
 write(fid2,'(A)') 'noCD'

 write(fid2,'(/,A)') "&SEWARD"
 call check_DKH_in_gms_inp(fort7, rel)

 select case(rel)
 case(-2) ! nothing
 case(-1) ! RESC
  write(6,'(/,A)') 'ERROR in subroutine bas_gms2molcas: RESC keywords detected.'
  write(6,'(A)') 'But RESC is not supported in (Open)Molcas.'
  stop
 case(0,1,2,4)  ! DKH0/1/2/4
  if(.not. X2C) write(fid2,'(A,I2.2,A)') 'Relativistic= R',rel,'O'
  !if(.not. X2C) write(fid2,'(A,I2.2,A)') 'R',rel,'O'
 case default
  write(6,'(/,A)') 'ERROR in subroutine bas_gms2molcas: rel out of range!'
  write(6,'(A,I0)') 'rel=', rel
  close(fid2,status='delete')
  stop
 end select

 if((.not.uhf) .and. (mult/=1)) then ! ROHF
  write(fid2,'(/,A)') "&RASSCF"
  write(fid2,'(A,I0)') 'Charge= ', charge
  write(fid2,'(A,I0)') 'Spin= ', mult
  write(fid2,'(A)') 'CIMX= 200'
  write(fid2,'(A)') 'Tight= 5d-8 5d-6'
  write(fid2,'(A,I0)') 'nActEl= ', mult-1
  write(fid2,'(A,I0)') 'RAS2= ', mult-1
  write(fid2,'(A)') 'OutOrbitals= Canonical'
 else
  write(fid2,'(/,A)') "&SCF"
  if(uhf) then ! UHF
   write(fid2,'(A)') 'UHF'
   call read_na_and_nb_from_gms_inp(fort7, na, nb)
   write(fid2,'(A,/,I0,/,I0)') 'Occupied', na, nb
   ! For OpenMolcas>=23.02, default initial occupation number scheme is changed,
   ! we cannot use (and no need to use) the "Fermi aufbau" since we have converged
   ! orbitals
  else         ! RHF
   write(fid2,'(A,I0)') 'Charge= ', charge
   write(fid2,'(A,I0)') 'Spin= ', mult
  end if
 end if

 i = INDEX(fort7, '.', back=.true.)
 input = fort7(1:i-1)//'.INPORB'
 write(fid2,'(A)') 'FILEORB= '//TRIM(input)

 close(fid2)
end subroutine bas_gms2molcas

