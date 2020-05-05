! written by jxzou at 20200504: read/write basis, MOs or eigenvalues from/to given files

! read nalpha and nbeta from .fch(k) file
subroutine read_na_and_nb_from_fch(fchname, na, nb)
 implicit none
 integer :: fid
 integer, intent(out) :: na, nb ! number of alpha/beta electrons
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:15) == 'Number of alpha') exit
 end do

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, na
 read(fid,'(A49,2X,I10)') buf, nb
 close(fid)
 return
end subroutine read_na_and_nb_from_fch

! read nbf and nif from .fch(k) file
subroutine read_nbf_and_nif_from_fch(fchname, nbf, nif)
 implicit none
 integer :: fid
 integer, intent(out) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:17) == 'Number of basis f') exit
 end do

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, nbf
 read(fid,'(A49,2X,I10)') buf, nif
 close(fid)
 return
end subroutine read_nbf_and_nif_from_fch

! read nbf and nif from .Orb file of MOLCAS/OpenMOLCAS
subroutine read_nbf_and_nif_from_orb(orbname, nbf, nif)
 implicit none
 integer :: fid
 integer, intent(out) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: orbname

 open(newunit=fid,file=TRIM(orbname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:5) == '#INFO') exit
 end do

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 read(fid,*) nbf
 read(fid,*) nif
 close(fid)
 return
end subroutine read_nbf_and_nif_from_orb

! read Alpha/Beta MOs from a given .fch(k) file
subroutine read_mo_from_fch(fchname, nbf, nif, ab, mo)
 implicit none
 integer :: i, fid, ncoeff
 integer, intent(in) :: nbf, nif
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: mo(nbf,nif)
 real(kind=8), allocatable :: coeff(:)
 character(len=1), intent(in) :: ab
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 character(len=8) :: key
 character(len=8), parameter :: key1 = 'Alpha MO'
 character(len=7), parameter :: key2 = 'Beta MO'

 key = key1
 if(ab/='a' .and. ab/='A') key = key2//' '

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == key) exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_mo_from_fch: no '"//key//"' found&
                   & in file "//TRIM(fchname)//'.'
  stop
 end if

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, ncoeff
 if(ncoeff /= nbf*nif) then
  write(iout,'(A)') 'ERROR in subroutine read_mo_from_fch: ncoeff /= nbf*nif.'
  write(iout,'(A)') 'Inconsistency found between input nbf,nif and those&
                   & in file '//TRIM(fchname)
  stop
 end if

 allocate(coeff(ncoeff), source=0.0d0)
 read(fid,'(5(1X,ES15.8))') (coeff(i),i=1,ncoeff)
 mo = RESHAPE(coeff,(/nbf,nif/))

 deallocate(coeff)
 close(fid)
 return
end subroutine read_mo_from_fch

! read Alpha/Beta MOs from .Orb file of MOLCAS/OpenMolcas
subroutine read_mo_from_orb(orbname, nbf, nif, ab, mo)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nbf, nif
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: mo(nbf,nif)
 character(len=1), intent(in) :: ab
 character(len=240) :: buf
 character(len=240), intent(in) :: orbname
 character(len=4) :: key
 character(len=4), parameter :: key1 = '#ORB', key2 = '#UOR'

 key = key1
 if(ab/='a' .and. ab/='A') key = key2

 open(newunit=fid,file=TRIM(orbname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == key) exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_mo_from_orb: no '"//key//"' found&
                   & in file "//TRIM(orbname)//'.'
  stop
 end if

 mo = 0.0d0
 do i = 1, nif, 1
  read(fid,'(A)') buf
  read(fid,'(5(1X,ES21.14))') (mo(j,i),j=1,nbf)
 end do ! for i

 close(fid)
 return
end subroutine read_mo_from_orb

! read Alpha/Beta eigenvalues in a given .fch(k) file
! Note: the Alpha/Beta Orbital Energies in .fch(k) file can be either energy levels
!       or NOONs, depending on the job type
subroutine read_eigenvalues_from_fch(fchname, nif, ab, noon)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nif
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: noon(nif)
 character(len=1), intent(in) :: ab
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 character(len=8) :: key
 character(len=8), parameter :: key1 = 'Alpha Or'
 character(len=7), parameter :: key2 = 'Beta Or'

 key = key1
 if(ab/='a' .and. ab/='A') key = key2//'b'

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == key) exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_eigenvalues_from_fch: no '"//key//"' found&
                   & in file "//TRIM(fchname)//'.'
  stop
 end if

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, i
 if(i /= nif) then
  write(iout,'(A)') 'ERROR in subroutine read_eigenvalues_from_fch: i /= nif.'
  write(iout,'(A)') 'Inconsistency found between input nif and that in file '//TRIM(fchname)
  stop
 end if

 noon = 0.0d0
 read(fid,'(5(1X,ES15.8))') (noon(i),i=1,nif)
 close(fid)
 return
end subroutine read_eigenvalues_from_fch

! read Alpha/Beta occupation numbers from a given .orb file of MOLCAS/OpenMolcas
subroutine read_on_from_orb(orbname, nif, ab, on)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nif
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: on(nif)
 character(len=1), intent(in) :: ab
 character(len=240) :: buf
 character(len=240), intent(in) :: orbname
 character(len=5) :: key
 character(len=4), parameter :: key1 = '#OCC'
 character(len=5), parameter :: key2 = '#UOCC'

 key = key1//' '
 if(ab/='a' .and. ab/='A') key = key2

 open(newunit=fid,file=TRIM(orbname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(index(buf,TRIM(key)) /= 0) exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_on_from_orb: no '"//TRIM(key)&
                   & //"' found in file "//TRIM(orbname)//'.'
  stop
 end if

 on = 0.0d0
 read(fid,'(A)') buf
 read(fid,'(5(1X,ES21.14))') (on(i),i=1,nif)

 close(fid)
 return
end subroutine read_on_from_orb

! read the array size of shell_type and shell_to_atom_map from a given .fch(k) file
subroutine read_ncontr_from_fch(fchname, ncontr)
 implicit none
 integer :: i, fid
 integer, intent(out) :: ncontr
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 open(unit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:18) == 'Number of contract') exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_ncontr_from_fch: missing&
                   & 'Number of contract' section in file "//TRIM(fchname)
  close(fid)
  return
 end if

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, ncontr
 close(fid)
 return
end subroutine read_ncontr_from_fch

! read shell_type and shell_to_atom_map from a given .fch(k) file
subroutine read_shltyp_and_shl2atm_from_fch(fchname, k, shltyp, shl2atm, sph)
 implicit none
 integer :: i, fid
 integer, intent(in) :: k
 integer, intent(out) :: shltyp(k), shl2atm(k)
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 logical, intent(out) :: sph

 open(unit=fid,file=TRIM(fchname),status='old',position='rewind')

 ! find and read Shell types
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == 'Shell types') exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_shltyp_and_shl2atm_from_fch:&
                   & missing 'Shell types' section in file "//TRIM(fchname)
  close(fid)
  return
 end if

 shltyp = 0
 read(fid,'(6(6X,I6))') (shltyp(i),i=1,k)
 ! read Shell types done

 if( ANY(shltyp>1) ) then ! whether spheical/Cartesian
  sph = .false.
 else
  sph = .true.
 end if

 ! find and read Shell to atom map
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:13) == 'Shell to atom') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_shltyp_and_shl2atm_from_fch:&
                   & missing 'Shell to atom map' section in file "//TRIM(fchname)
  close(fid)
  return
 end if

 shl2atm = 0
 read(fid,'(6(6X,I6))') (shl2atm(i),i=1,k)
 close(fid)
 return
end subroutine read_shltyp_and_shl2atm_from_fch

! write/print eigenvalues/occupation numbers into a .fch(k) file
subroutine write_eigenvalues_to_fch(fchname, nif, ab, on, replace)
 implicit none
 integer :: i, fid1, fid2, RENAME
 integer, intent(in) :: nif
 integer, parameter :: iout = 6
 real(kind=8), intent(in) :: on(nif)
 character(len=1), intent(in) :: ab
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
 character(len=8) :: key
 character(len=8), parameter :: key1 = 'Alpha Or'
 character(len=7), parameter :: key2 = 'Beta Or'
 logical, intent(in) :: replace

 key = key1
 if(ab/='a' .and. ab/='A') key = key2//'b'

 open(newunit=fid1,file=TRIM(fchname),status='old',position='rewind')
 i = index(fchname,'.fch',back=.true.)
 fchname1 = fchname(1:i-1)//'_D.fch'
 open(newunit=fid2,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(buf(1:8) == key) exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine write_eigenvalues_to_fch: no '"//&
                    key//"' found in file "//TRIM(fchname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(5(1X,ES15.8))') (on(i),i=1,nif)

 ! skip the Alpha/Beta Orbital Energies in fname1
 do while(.true.)
  read(fid1,'(A)') buf
  if(index(buf,'=') /= 0) exit
 end do
 BACKSPACE(fid1)

 ! copy remaining content in file fname1
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do
 close(fid2)

 if(replace) then
  close(fid1,status='delete')
  i = RENAME(TRIM(fchname1), TRIM(fchname))
 else
  close(fid1)
 end if

 return
end subroutine write_eigenvalues_to_fch

! read Alpha/Beta occupation numbers from a given .orb file of MOLCAS/OpenMolcas
subroutine write_on_to_orb(orbname, nif, ab, on, replace)
 implicit none
 integer :: i, nline, fid1, fid2, RENAME
 integer, intent(in) :: nif
 integer, parameter :: iout = 6
 real(kind=8), intent(in) :: on(nif)
 character(len=1), intent(in) :: ab
 character(len=240) :: buf, orbname1
 character(len=240), intent(in) :: orbname
 character(len=6) :: key
 character(len=4), parameter :: key1 = '#OCC'
 character(len=5), parameter :: key2 = '#UOCC', key3 = '#OCHR'
 character(len=6), parameter :: key4 = '#UOCHR'
 logical, intent(in) :: replace

 open(newunit=fid1,file=TRIM(orbname),status='old',position='rewind')
 i = index(orbname,'.',back=.true.)
 orbname1 = orbname(1:i-1)//'_D.Orb'
 open(newunit=fid2,file=TRIM(orbname1),status='replace')

 key = key1//'  '
 if(ab/='a' .and. ab/='A') key = key2//' '

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(index(buf,TRIM(key)) /= 0) exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_on_from_orb: no '"//TRIM(key)&
                   & //"' found in file "//TRIM(orbname)//'.'
  stop
 end if

 read(fid1,'(A)') buf
 write(fid2,'(A)') TRIM(buf)
 write(fid2,'(5(1X,ES21.14))') (on(i),i=1,nif)

 nline = nif/5
 if(nif-5*nline > 0) nline = nline + 1
 ! skip the #OCC/#UOCC data in original file
 do i = 1, nline, 1 
  read(fid1,'(A)') buf
 end do ! for i

 key = key3//'  '
 if(ab/='a' .and. ab/='A') key = key4//' '

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(index(buf,TRIM(key)) /= 0) exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_on_from_orb: no '"//TRIM(key)&
                   & //"' found in file "//TRIM(orbname)//'.'
  stop
 end if

 read(fid1,'(A)') buf
 write(fid2,'(A)') TRIM(buf)
 write(fid2,'(10(1X,F7.4))') (on(i),i=1,nif)

 nline = nif/10
 if(nif-10*nline > 0) nline = nline + 1
 ! skip the #OCHR/#UOCHR data in original file
 do i = 1, nline, 1 
  read(fid1,'(A)') buf
 end do ! for i

 ! copy remaining content
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while
 close(fid2)

 if(replace) then
  close(fid1,status='delete')
  i = RENAME(TRIM(orbname1), TRIM(orbname))
 else
  close(fid1)
 end if

 return
end subroutine write_on_to_orb

