! written by jxzou at 20210104: copied from file bas_gms2molcas.f90 and modified
! updated by jxzou at 20210407: remove '-uhf', add automatic determination

! Note: Currently isotopes are not tested.
program main
 use pg, only: iout
 implicit none
 integer :: i
 character(len=240) :: fname
 ! fname: input file contains basis sets and Cartesian coordinates in GAMESS format

 i = iargc()
 if(i<1 .or. i>2) then
  write(iout,'(/,A)') ' ERROR in subroutine bas_gms2bdf: wrong command line arguments!'
  write(iout,'(A,/)') ' Example : bas_gms2bdf a.inp (generate an a_bdf.inp file)'
  stop
 end if

 fname = ' '
 call getarg(1, fname)
 call require_file_exist(fname)
 call bas_gms2bdf(fname)
 stop
end program main

! Transform the basis sets in GAMESS format to those in BDF format
subroutine bas_gms2bdf(fort7)
 use pg, only: iout, natom, ram, ntimes, coor, elem, all_ecp, ecp_exist
 implicit none
 integer :: i, k, nline, rc, rel, nbf, nif
 integer :: fid1, fid2
 integer :: charge, mult
 character(len=7) :: str
 character(len=240), intent(in) :: fort7
 character(len=240) :: buf, input, basfile
 ! input is the BDF input file
 ! if you do not like the suffix .bdf, you can change it into .inp
 character(len=1) :: stype
 logical :: X2C, uhf, lin

 ! initialization
 buf = ' '
 input = ' '
 basfile = ' '

 call read_nbf_and_nif_from_gms_inp(fort7, nbf, nif)
 if(nbf > nif) then
  lin = .true.  ! linear dependent
 else
  lin = .false. ! no linear dependence
 end if

 k = index(fort7, '.', back=.true.)
 input = fort7(1:k-1)//'_bdf.inp'
 call read_natom_from_gms_inp(fort7, natom)
 allocate(elem(natom), ram(natom), coor(3,natom), ntimes(natom))
 call read_elem_nuc_coor_from_gms_inp(fort7, natom, elem, ram, coor)
 ! ram cannot be deallocated here since subroutine prt_prim_gau_bdf will use it

 call calc_ntimes(natom, elem, ntimes)
 call read_charge_and_mult_from_gms_inp(fort7, charge, mult, uhf, ecp_exist)
 call read_all_ecp_from_gms_inp(fort7)

 if(ecp_exist) then
  basfile = 'ECP.'//fort7(1:k-1)//'.BAS'
 else
  basfile = fort7(1:k-1)//'.BAS'
 end if
 ! read ECP/PP done

 call upper(basfile)
 open(newunit=fid2,file=TRIM(input),status='replace')
 write(fid2,'(A)') '$COMPASS'
 write(fid2,'(A)') 'Title'
 write(fid2,'(A)') 'auto-generated file by the bas_gms2bdf utility of MOKIT'
 write(fid2,'(A)') 'Geometry'

 call check_X2C_in_gms_inp(fort7, X2C)
 call check_DKH_in_gms_inp(fort7, rel)

 do i = 1, natom, 1
  str = ' '
  write(str,'(A,I0)') TRIM(elem(i)), ntimes(i)
  write(fid2,'(A7,1X,3F18.8)') str, coor(1:3,i)
 end do ! for i
 deallocate(coor)

 write(fid2,'(A)') 'End Geometry'
 write(fid2,'(A)') 'Basis'
 write(fid2,'(A)') TRIM(basfile)
 write(fid2,'(A)') 'Nosymm'
 write(fid2,'(A)') '$END'
 write(fid2,'(/,A)') '$XUANYUAN'

 if(rel>-1 .or. X2C) then
  write(fid2,'(A)') 'Scalar'

  if(rel>-1 .and. (.not.X2C)) then
   write(iout,'(A)') 'Warning in subroutine bas_gms2bdf: BDF program does&
                    & not support DKH Hamiltonian.'
   write(iout,'(A)') "Spin-free X2C keyword 'Scalar' are written in $XUANYUAN&
                    & instead."
  end if
 end if

 write(fid2,'(A)') '$END'
 write(fid2,'(/,A)') '$SCF'
 if(lin) then
  write(fid2,'(A)') 'CheckLin'
  write(fid2,'(A)') 'TolLin'
  write(fid2,'(A)') '1.D-6'
 end if

 if(uhf) then
  write(fid2,'(A3)') 'UHF'
 else
  if(mult /= 1) then
   write(fid2,'(A4)') 'ROHF'
  else
   write(fid2,'(A3)') 'RHF'
  end if
 end if

 write(fid2,'(A,/,I0)') 'Charge', charge
 write(fid2,'(A,/,I0)') 'Spin', mult
 write(fid2,'(A)') 'Guess'
 write(fid2,'(A)') 'read'
 write(fid2,'(A)') '$END'
 close(fid2)

 ! open the .inp/.dat file and find the $DATA section
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
  write(iout,'(A)') 'ERROR in subroutine bas_gms2bdf: No $DATA section found&
                   & in file '//TRIM(fort7)//'.'
  close(fid1)
  stop
 end if

 ! skip 2 lines: the Title line and the Point Group line
 read(fid1,'(A)') buf
 read(fid1,'(A)') buf

 ! initialization: clear all primitive gaussians
 call clear_prim_gau()
 open(newunit=fid2,file=TRIM(basfile),status='replace')

 ! read the element, relative atomic mass (ram) and coordinates
 do i = 1, natom, 1
  read(fid1,'(A)',iostat=rc) buf
  ! 'buf' contains the element, ram and coordinates
  if(rc /= 0) exit

  ! deal with primitive gaussians
  do while(.true.)
   read(fid1,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit

   read(buf,*) stype, nline
   call read_prim_gau(stype, nline, fid1)
  end do ! for while

  ! print basis sets and ECP/PP (if any) of this atom in BDF format
  call prt_prim_gau_bdf(i, fid2)

  ! clear all primitive gaussians for next cycle
  call clear_prim_gau()
 end do ! for i

 deallocate(ram, ntimes, elem, all_ecp)
 close(fid1)

 if(rc /= 0) then
  write(iout,'(A)') "ERROR in subroutine bas_gms2bdf: it seems the '$DATA'&
                   & has no corresponding '$END'."
  write(iout,'(A)') 'Incomplete file '//TRIM(fort7)
  close(fid2,status='delete')
  stop
 end if

 close(fid2)
 return
end subroutine bas_gms2bdf

! print primitive gaussians
subroutine prt_prim_gau_bdf(iatom, fid)
 use pg, only: prim_gau, ram, highest, all_ecp, elem, ntimes, ecp_exist
 implicit none
 integer :: i, j, k, m, n, nline, ncol
 integer, intent(in) :: iatom, fid
 integer, allocatable :: list(:)
 character(len=1), parameter :: am(0:5) = ['s','p','d','f','g','h']

 call get_highest_am()
 write(fid,'(A)') '****'
 write(fid,'(A,I0,3X,I0,3X,I1)') TRIM(elem(iatom)),ntimes(iatom),ram(iatom),highest

 do i = 1, 7, 1
  if(.not. allocated(prim_gau(i)%coeff)) cycle
  write(fid,'(A)',advance='no') prim_gau(i)%stype
  nline = prim_gau(i)%nline
  ncol = prim_gau(i)%ncol
  write(fid,'(2(1X,I4))') nline, ncol-1
  do j = 1, nline, 1
   write(fid,'(3X,E16.9)') prim_gau(i)%coeff(j,1)
  end do ! for j
  do j = 1, nline, 1
   write(fid,'(10(E16.9,3X))') (prim_gau(i)%coeff(j,k), k=2,ncol)
  end do ! for j
 end do ! for i

 if(.not. ecp_exist) return

 if(all_ecp(iatom)%ecp) then
  write(fid,'(A)') 'ECP'
  m = all_ecp(iatom)%highest
  write(fid,'(A,2(I0,1X),I0)') TRIM(elem(iatom)),ntimes(iatom),all_ecp(iatom)%core_e,m
  allocate(list(m+1))
  list(1) = m
  forall(i=2:m+1) list(i) = i-2
  do i = 1, m+1, 1
   n = all_ecp(iatom)%potential(i)%n
   if(i == 1) then
    write(fid,'(A,I0)') am(list(i))//' potential ', n
   else ! i > 1
    write(fid,'(A,I0)') am(list(i))//' potential ', n
   end if

   do j = 1, n, 1
    write(fid,'(I0,2(1X,F16.8))') all_ecp(iatom)%potential(i)%col2(j), all_ecp(iatom)%potential(i)%col3(j), &
                                  all_ecp(iatom)%potential(i)%col1(j)
   end do ! for j
  end do ! for i
 end if

 return
end subroutine prt_prim_gau_bdf

