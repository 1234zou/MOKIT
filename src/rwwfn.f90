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

subroutine read_mo_from_chk_txt(txtname, nbf, nif, ab, mo)
 implicit none
 integer :: i, k, chkid
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(out) :: mo(nbf,nif)
 character(len=1), intent(in) :: ab
 character(len=240) :: buf
 character(len=240), intent(in) :: txtname
 character(len=8) :: key
 character(len=8), parameter :: key1 = 'Alpha MO'
 character(len=7), parameter :: key2 = 'Beta MO'

 key = key1
 if(ab/='a' .and. ab/='A') then
  key = key2//' '
 end if

 open(unit=chkid,file=TRIM(txtname),status='old',position='rewind')
 do while(.true.)
  read(chkid,'(A)') buf
  if(buf(7:14) == key) exit
 end do
 BACKSPACE(chkid)

 mo = 0.0d0
 do i = 1, nif, 1
  read(chkid,'(A)') buf
  read(chkid,'(4D20.12)') (mo(k,i), k=1,nbf)
 end do

 close(chkid)
 return
end subroutine read_mo_from_chk_txt

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

! read AO-basis overlap matrix from an OpenMolcas output file
subroutine read_ovlp_from_gau_log(logname, nbf, S)
 implicit none
 integer :: i, j, k, fid, nline
 integer, intent(in) :: nbf
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: S(nbf,nbf)
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 S = 0.0d0

 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:16) == '*** Overlap ***') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_ovlp_from_gau_log: no '***&
                   & Overlap ***' found in file "//TRIM(logname)
  stop
 end if

 k = 1
 do while(k <= nbf)
  read(fid,'(A)') buf
  write(6,'(A)') TRIM(buf)//' '

  do i = k, nbf, 1
   read(fid,'(A)') buf
   buf(1:7) = ' '
   buf = ADJUSTL(buf)
   j = i - k + 1
   if(j > 4) j = 5
   read(buf,*) S(k:k+j-1,i)
   write(6,'(I7,1X,D19.12,4D20.12)') i, S(k:k+j-1,i)
  end do ! for i
  k = k + 5
 end do ! for while

 close(fid)

 do i = 1, nbf-1, 1
  do j = i+1, nbf, 1
   S(j,i) = S(i,j)
  end do ! for j
 end do ! for i

 return
end subroutine read_ovlp_from_gau_log

! read AO-basis overlap matrix from an OpenMolcas output file
! Note: add
!          print
!          1
!          112 10
! in &SEWARD section enables printing one-electron integrals
subroutine read_ovlp_from_molcas_out(outname, nbf, S)
 implicit none
 integer :: i, j, k, nline, fid
 integer, intent(in) :: nbf
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: S(nbf,nbf)

 S = 0.0d0

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:30) == 'SO Integrals of type Mltpl  0') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_ovlp_from_molcas_out: no 'SO&
                   & Integrals of type Mltpl  0' found in file "//TRIM(outname)
  stop
 end if

 do i = 1, 3
  read(fid,'(A)') buf
 end do ! for i

 i = index(buf, 'x')
 read(buf(i+1:),*) j
 if(j /= nbf) then
  write(iout,'(A)') 'ERROR in subroutine read_ovlp_from_molcas_out: inconsistent&
   & nbf(i.e. number of basis functions) in orbital file and overlap file.'
   write(iout,'(2(A,I5))') 'j=', j, ', nbf=', nbf
  stop
 end if

 read(fid,'(A)') buf
 do i = 1, nbf, 1
  nline = i/5

  do j = 1, nline, 1
   read(fid,*) S(5*j-4:5*j,i)
  end do ! for k
  k = i - 5*nline

  if(k > 0) read(fid,*) S(5*nline+1:i,i)
 end do ! for i

 close(fid)

 do i = 1, nbf-1, 1
  do j = i+1, nbf, 1
   S(j,i) = S(i,j)
  end do ! for j
 end do ! for i

 return
end subroutine read_ovlp_from_molcas_out

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

! write Alpha/Beta MOs from a given .fch(k) file
subroutine write_mo_into_fch(fchname, nbf, nif, ab, mo)
 implicit none
 integer :: i, fid1, fid2, ncoeff, RENAME
 integer, intent(in) :: nbf, nif
 integer, parameter :: iout = 6
 real(kind=8), intent(in) :: mo(nbf,nif)
 real(kind=8), allocatable :: coeff(:)
 character(len=1), intent(in) :: ab
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
 character(len=8) :: key
 character(len=8), parameter :: key1 = 'Alpha MO'
 character(len=7), parameter :: key2 = 'Beta MO'
 logical :: alive

 inquire(file=TRIM(fchname),exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine write_mo_into_fch: file '//TRIM(fchname)&
                  //' does not exist.'
  stop
 end if

 key = key1
 if(ab/='a' .and. ab/='A') key = key2//' '
 fchname1 = TRIM(fchname)//'.tmp'

 open(newunit=fid1,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(buf(1:8) == key) exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_mo_from_fch: no '"//key//"' found&
                   & in file "//TRIM(fchname)//'.'
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 BACKSPACE(fid1)
 read(fid1,'(A49,2X,I10)') buf, ncoeff
 if(ncoeff /= nbf*nif) then
  write(iout,'(A)') 'ERROR in subroutine read_mo_from_fch: ncoeff /= nbf*nif.'
  write(iout,'(A)') 'Inconsistency found between input nbf,nif and those&
                   & in file '//TRIM(fchname)//'.'
  stop
 end if

 allocate(coeff(ncoeff), source=0.0d0)
         coeff = RESHAPE(mo, (/ncoeff/))
 write(fid2,'(5(1X,ES15.8))') (coeff(i),i=1,ncoeff)
 deallocate(coeff)

 do while(.true.) ! skip MOs in fchname
  read(fid1,'(A)') buf
  if(index(buf,'=') /= 0) exit
 end do ! for while

 BACKSPACE(fid1)
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(fchname1, fchname)
 return
end subroutine write_mo_into_fch

! determine whether sperical harmonic or Cartesian fucntions are used in .fch(k) file
subroutine determine_sph_or_cart(fchname, cart)
 implicit none
 integer :: i, k, fid
 integer, allocatable :: shltyp(:)
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 logical, intent(out) :: cart

 open(unit=fid,file=TRIM(fchname),status='old',position='rewind')

 ! find and read Shell types
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == 'Shell types') exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine determine_sph_or_cart: missing&
                   & 'Shell types' section in file "//TRIM(fchname)
  close(fid)
  return
 end if

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, k
 allocate(shltyp(k), source=0)
 read(fid,'(6(6X,I6))') (shltyp(i),i=1,k)
 ! read Shell types done

 close(fid)

 cart = .false.
 if( ANY(shltyp>1) ) cart = .true. ! whether spheical/Cartesian

 deallocate(shltyp)
 return
end subroutine determine_sph_or_cart

! read 4 variables (npair, nbf, nif, lin_dep) from uno.out
subroutine read_npair_from_uno_out(nbf, nif, ndb, npair, nopen, lin_dep)
 implicit none
 integer :: i, fid, idx(3), nvir
 integer, parameter :: iout = 6
 integer, intent(out) :: nbf, nif, ndb, npair, nopen
 character(len=240) :: buf = ' '
 logical, intent(out) :: lin_dep

 open(newunit=fid,file='uno.out',status='old',position='rewind')
 read(fid,'(A)') buf
 i = index(buf,'=')
 read(buf(i+1:),*) nbf

 read(fid,'(A)') buf
 i = index(buf,'=')
 read(buf(i+1:),*) nif

 if(nbf > nif) then
  lin_dep = .true.
 else if(nbf < nif) then
  write(iout,'(A)') 'ERROR in subroutine read_npair_from_uno_out: nbf<nif.'
  write(iout,'(A)') 'This is impossible. Please check why.'
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3) == 'ndb') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_npair_from_uno_out: 'ndb' not found."
  write(iout,'(A)') "The file 'uno.out' may be incomplete."
  close(fid)
  stop
 end if

 i = index(buf,'=')
 read(buf(i+1:),*) ndb

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3) == 'idx') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_npair_from_uno_out: 'idx' not found."
  write(iout,'(A)') "The file 'uno.out' may be incomplete."
  close(fid)
  stop
 end if

 idx = 0
 i = index(buf,'=')
 read(buf(i+1:),*) idx(1:3)
 close(fid,status='delete')

 npair = (idx(2) - idx(1) - idx(3))/2
 nvir = nif - ndb - 2*npair - idx(3)
 nopen = idx(3)
 write(iout,'(A,I5,4X,A,I5)') 'nbf =', nbf, 'nif =', nif
 write(iout,'(4(A,I5,4X))') 'doubly_occ=', idx(1)-1, 'npair=', npair, 'nopen=',&
                            idx(3), 'nvir=', nvir
 return
end subroutine read_npair_from_uno_out

subroutine read_gvb_energy_from_gms(gmsname, e)
 implicit none
 integer :: i, j, fid
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: e
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname

 e = 0.0d0
 open(newunit=fid,file=TRIM(gmsname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:10) == 'FINAL GVB') exit
 end do

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_gvb_energy_from_gms: no GVB energy&
                   & found in file '//TRIM(gmsname)
  close(fid)
  stop
 end if
 close(fid)

 i = index(buf,'IS'); j = index(buf,'AFTER')
 read(buf(i+2:j-1),*) e

 if(DABS(e) < 1.0d-4) then
  write(iout,'(A)') 'ERROR in subroutine read_gvb_energy_from_gms: GVB computation&
                   & does not converge.'
  stop
 end if

 write(iout,'(/,A,F18.8,1X,A4)') 'E(GVB) = ', e, 'a.u.'
 return
end subroutine read_gvb_energy_from_gms

! read CASCI/CASSCF energy from a Gaussian/PySCF/GAMESS/OpenMolcas/ORCA output file
subroutine read_cas_energy_from_output(cas_prog, outname, e, scf, spin, dmrg)
 implicit none
 integer, intent(in) :: spin
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: e(2)
 character(len=10), intent(in) :: cas_prog
 character(len=10), intent(in) :: outname
 logical, intent(in) :: scf, dmrg

 select case(TRIM(cas_prog))
 case('gaussian')
  call read_cas_energy_from_gaulog(outname, e, scf)
 case('gamess')
  call read_cas_energy_from_gmsgms(outname, e, scf, spin)
 case('pyscf')
  call read_cas_energy_from_pyout(outname, e, scf, spin, dmrg)
 case('openmolcas')
  call read_cas_energy_from_molcas_out(outname, e, scf)
 case('orca')
  
 case default
  write(iout,'(A)') 'ERROR in subroutine read_cas_energy_from_output: cas_prog&
                   & cannot be identified.'
  write(iout,'(A)') 'cas_prog='//TRIM(cas_prog)
  stop
 end select

 return
end subroutine read_cas_energy_from_output

! read CASCI/CASSCF energy from a Gaussian .log/.out file
subroutine read_cas_energy_from_gaulog(outname, e, scf)
 implicit none
 integer :: i, fid
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(13:22) == 'EIGENVALUE') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_cas_energy_from_gaulog: no&
                   & 'EIGENVALUE' found in file "//TRIM(outname)
  close(fid)
  stop
 end if
 close(fid)

 if(scf) then
  read(buf(23:),*) e(2) ! CASSCF
 else
  read(buf(23:),*) e(1) ! CASCI
 end if

 if(scf) then ! read CASCI energy
  open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)') buf
   if(i /= 0) exit
   if(buf(2:8) == 'ITN=  1') exit
  end do ! for while
  close(fid)
  i = index(buf,'E=')
  read(buf(i+2:),*) e(1)
 end if

 return
end subroutine read_cas_energy_from_gaulog

! read CASCI/CASSCF energy from a PySCF output file
subroutine read_cas_energy_from_pyout(outname, e, scf, spin, dmrg)
 implicit none
 integer :: i, j, k, fid
 integer, parameter :: iout = 6
 integer, intent(in) :: spin ! na - nb
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8) :: s_square = 0.0d0, expect = 0.0d0
 real(kind=8), intent(out) :: e(2)
 logical, intent(in) :: scf, dmrg

 e = 0.0d0
 i = 0; j = 0; k = 0
 expect = DBLE(spin)/2.0d0
 expect = expect*(expect + 1.0d0)

 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 if(scf) then ! CASSCF/DMRG-CASSCF
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)',iostat=k) buf
   if(k /= 0) exit
   if(buf(1:23) == '1-step CASSCF converged') then
    i = 1; exit
   end if
   if(buf(1:27) == '1-step CASSCF not converged') then
    j = 1; exit
   end if
  end do ! for while

 else ! CASCI/DMRG-CASCI

  if(.not. dmrg) then
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)',iostat=k) buf
    if(k /= 0) exit
    if(buf(1:15) == 'CASCI converged') then
     i = 1; exit
    end if
    if(buf(1:19) == 'CASCI not converged') then
     j = 1; exit
    end if
   end do ! for while
  else
   BACKSPACE(fid)
  end if
 end if

 if(k /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_cas_energy_from_out:&
                   & incomplete file '//TRIM(outname)//'.'
  stop
 else ! k = 0
  if(j /= 0) then
   write(iout,'(A)') 'ERROR in subroutine read_cas_energy_from_out:&
                    & CASCI or CASSCF not converged.'
   stop
  end if
 end if

 ! read CASCI/CASSCF energy in CASCI/CASSCF job
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == 'CASCI E') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_cas_energy_from_out:&
                   & 'CASCI E' not found in "//TRIM(outname)
  stop
 end if

 i = index(buf, '=')
 if(scf) then
  read(buf(i+1:),*) e(2)
 else
  read(buf(i+1:),*) e(1)
 end if

 i = index(buf, '=', back=.true.)
 read(buf(i+1:),*) s_square
 if( DABS(expect - s_square) > 1.0D-3) then
  if(scf) then
   write(iout,'(A)') 'ERROR in subroutine read_cas_energy_from_pyout: CASSCF&
                    & <S**2> deviates too much from expectation value.'
  else
   write(iout,'(A)') 'ERROR in subroutine read_cas_energy_from_pyout: CASCI&
                    & <S**2> deviates too much from expectation value.'
  end if
  write(iout,'(2(A,F10.6))') 'expectation = ', expect, ', s_square=', s_square
  stop
 end if

 ! Note: in a CASSCF job, there is also a CASCI energy, read it
 if(scf) then
  rewind(fid)
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(1:9) == 'CASCI E =') exit
  end do ! for while
  close(fid)
  read(buf(10:),*) e(1)

  i = index(buf, '=', back=.true.)
  read(buf(i+1:),*) s_square
  if( DABS(expect - s_square) > 1.0D-3) then
   write(iout,'(A)') 'ERROR in subroutine read_cas_energy_from_pyout: in this&
                    & CASSCF job, the 0-th step, i.e., the CASCI'
   write(iout,'(A)') '<S**2> deviates too much from the expectation value.'
   write(iout,'(2(A,F10.6))') 'expectation = ', expect, ', s_square=', s_square
   stop
  end if
 end if

 return
end subroutine read_cas_energy_from_pyout

! read CASCI/CASSCF energy from the GAMESS output file
subroutine read_cas_energy_from_gmsgms(outname, e, scf, spin)
 implicit none
 integer :: i, fid
 integer, intent(in) :: spin
 integer, parameter :: iout = 6
 real(kind=8) :: s_square, expect
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 expect = DBLE(spin)/2.0d0
 expect = expect*(expect + 1.0d0)

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 if(scf) then  ! CASSCF job
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:24) == 'THE DENSITIES ARE STATE') exit
  end do ! for while
 
  if(i /= 0) then
   write(iout,'(A)') "ERROR in subroutine read_cas_energy_from_gmsgms: no&
                   & 'THE DENSITIES ARE STATE' found in file "//TRIM(outname)
   stop
  end if

  read(fid,'(A)') buf
  i = index(buf,'ENERGY=')
  read(buf(i+7:),*) e(1)   ! CASCI energy in the CASSCF job
  i = index(buf,'=',back=.true.)
  read(buf(i+1:),*) s_square
  if( DABS(expect - s_square) > 1.0D-2) then
   write(iout,'(A)') 'ERROR in subroutine read_cas_energy_from_gmsgms: in this&
                    & CASSCF job, the 0-th step, i.e., the CASCI'
   write(iout,'(A)') '<S**2> deviates too much from the expectation value.'
   write(iout,'(2(A,F10.6))') 'expectation = ', expect, ', s_square=', s_square
   stop
  end if

  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:13) == 'STATE   1  E') exit
  end do ! for while
  i = index(buf,'ENERGY=')
  read(buf(i+7:),*) e(2)   ! CASSCF energy
  i = index(buf,'S=')
  read(buf(i+2:),*) s_square
  if( DABS(expect - s_square) > 1.0D-2) then
   write(iout,'(A)') 'ERROR in subroutine read_cas_energy_from_gmsgms: CASSCF&
                    & <S**2> deviates too much from the expectation value.'
   write(iout,'(2(A,F10.6))') 'expectation = ', expect, ', s_square=', s_square
   stop
  end if

 else          ! CASCI job
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:20) == 'DENSITY MATRIX WILL') exit
  end do ! for while
 
  if(i /= 0) then
   write(iout,'(A)') "ERROR in subroutine read_cas_energy_from_gmsgms: no&
                   & 'DENSITY MATRIX' found in file "//TRIM(outname)
   stop
  end if
 
  i = index(buf,'=', back=.true.)
  read(buf(i+1:),*) s_square
  if( DABS(expect - s_square) > 1.0D-2) then
   write(iout,'(A)') 'ERROR in subroutine read_cas_energy_from_gmsgms: CASCI&
                    & <S**2> deviates too much from the expectation value.'
   write(iout,'(2(A,F10.6))') 'expectation = ', expect, ', s_square=', s_square
   stop
  end if

  read(fid,'(A)') buf
  read(fid,'(A)') buf
  i = index(buf,'=', back=.true.)
  read(buf(i+1:),*) e(1)
 end if

 close(fid)
 return
end subroutine read_cas_energy_from_gmsgms

! read CASCI/CASSCF energy from a given OpenMolcas/Molcas output file
subroutine read_cas_energy_from_molcas_out(outname, e, scf)
 implicit none
 integer :: i, j, fid
 integer, parameter :: iout = 6
 real(kind=8) :: add
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(7:27) == 'RASSCF root number  1') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_cas_energy_from_molcas_out: no&
                   & '' found in file "//TRIM(outname)//'.'
  stop
 end if

 e = 0.0d0; add = 0.0d0
 i = index(buf,':',back=.true.)
 if(scf) then   ! CASSCF
  read(buf(i+1:),*) e(2)
  rewind(fid)   ! read CASCI energy in CASSCF job
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:21) == 'Nr of preliminary CI') exit
  end do ! for while

  read(fid,'(A)') buf
  if(index(buf,'No convergence') /= 0) then
   write(iout,'(/,A)') 'Warning in subroutine read_cas_energy_from_molcas_out:'
   write(iout,'(A)') 'The 0-th step in CASSCF, i.e. the CASCI (in the CASSCF)&
                    & iterative diagonalization fails to converge.'
   write(iout,'(A)') 'This is a defect of OpenMolcas when doing CASSCF. If you want a'
   write(iout,'(A)') 'correct CASCI energy, please run a single CASCI job.'
   write(iout,'(A)') 'This may or may not affect the final CASSCF result, so continue.'
   read(fid,'(A)') buf
  end if
  if(index(buf,'Total energies') /= 0) then
   i = index(buf,'Add'); j = index(buf,'au')
   read(buf(i+3:j-1),*) add
   read(fid,'(A)') buf
  end if

  read(buf,*) i, j, i, j, e(1)
  e(1) = e(1) + add

 else           ! CASCI
  read(buf(i+1:),*) e(1)
 end if

 close(fid)
 return
end subroutine read_cas_energy_from_molcas_out

! read NEVPT2 energy from PySCF output file
subroutine read_mrpt2_energy_from_pyscf_out(outname, e)
 implicit none
 integer :: i, fid
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: e

 e = 0.0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:13) == 'Nevpt2 Energy') exit
 end do ! for while
 close(fid)

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_mrpt2_energy_from_pyscf_out:'
  write(iout,'(A)') 'No NEVPT2 energy found in file '//TRIM(outname)
  stop
 end if

 i = index(buf, '=')
 read(buf(i+1:),*) e
 write(iout,'(/,A,F18.8,1X,A4)') 'E_corr(SC-NEVPT2) = ', e, 'a.u.'
 return
end subroutine read_mrpt2_energy_from_pyscf_out

! read CASTP2 energy from OpenMolcas output file
subroutine read_mrpt2_energy_from_molcas_out(outname, e)
 implicit none
 integer :: i, fid
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: e

 e = 0.0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(index(buf,'Total CASPT2 energies:') /= 0) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_mrpt2_energy_from_molcas_out:'
  write(iout,'(A)') 'No CASPT2 energy found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 close(fid)

 i = index(buf, ':', back=.true.)
 read(buf(i+1:),*) e
 return
end subroutine read_mrpt2_energy_from_molcas_out

! find npair0: the number of active pairs (|C2| > 0.1)
! (assuming that the pair coefficients haven been sorted)
subroutine find_npair0_from_dat(datname, npair, npair0)
 implicit none
 integer :: i, k, datid
 integer, intent(in) :: npair
 integer, intent(out) :: npair0
 integer, parameter :: iout = 6
 real(kind=8), allocatable :: pair_coeff(:,:)
 character(len=240) :: buf
 character(len=240), intent(in) :: datname

 open(newunit=datid,file=TRIM(datname),status='old',position='rewind')

 ! find pair coefficients
 do while(.true.)
  read(datid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(index(buf,'CICOEF(') /= 0) exit
 end do
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine find_npair0_from_dat: no 'CICOEF(' found!"
  write(iout,'(A)') 'The input file '//TRIM(datname)//' is not complete!'
  close(datid)
  stop
 end if

 ! read pair coefficients
 BACKSPACE(datid)
 allocate(pair_coeff(2,npair), source=0.0d0)
 do i = 1, npair, 1
  read(datid,'(A)') buf
  k = index(buf,'=')
  read(buf(k+1:),*) pair_coeff(1,i)
  k = index(buf,',')
  read(buf(k+1:),*) pair_coeff(2,i)
 end do
 ! pair coefficients read done

 close(datid)

 npair0 = COUNT(pair_coeff(2,:) <= -0.1d0)
 deallocate(pair_coeff)
 return
end subroutine find_npair0_from_dat

! find the number of active UNO pairs from a given .fch(k) file
! Note that UNO are in pairs naturally, so npair0 from occupied space
!  must be equal to that from unoccupied space
subroutine find_npair0_from_fch(fchname, nopen, npair0)
 implicit none
 integer :: i, fid, nif
 integer, intent(in) :: nopen
 integer, intent(out) :: npair0
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 real(kind=8), parameter :: no_thres = 0.02d0
 real(kind=8), allocatable :: noon(:)

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == 'Alpha O') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine find_npair0_from_fch: keyword&
                  & 'Alpha O' not found in file "//TRIM(fchname)
  stop
 end if

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, nif
 allocate(noon(nif), source=0.0d0)
 read(fid,'(5(1X,ES15.8))') (noon(i),i=1,nif)
 close(fid)

 npair0 = 0
 do i = 1, nif, 1
  if(noon(i)>no_thres .and. noon(i)<(2.0d0-no_thres)) npair0 = npair0 + 1
 end do
 deallocate(noon)

 if(MOD(npair0-nopen,2) /= 0) then
  write(iout,'(A)') 'ERROR in subroutine find_npair0_from_fch: npair0 - nopen&
                   & is not an even integer.'
  write(iout,'(A)') "This is probably because UNO occupation numbers in 'Alpha O'&
                   & are probably incorrect."
  stop
 end if

 npair0 = (npair0 - nopen)/2
 return
end subroutine find_npair0_from_fch

! read variables nbf, nif, ndb, etc from a .fch(k) file containing NOs and NOONs
subroutine read_no_info_from_fch(fchname, nbf, nif, ndb, nopen, nacta, nactb, nacto, nacte)
 implicit none
 integer :: i, na, nb
 integer, intent(out) :: nbf, nif, ndb, nopen, nacta, nactb, nacto, nacte
 integer, parameter :: iout = 6
 real(kind=8), parameter :: no_thres = 0.02d0
 real(kind=8), allocatable :: noon(:)
 character(len=240), intent(in) :: fchname

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 call read_na_and_nb_from_fch(fchname, na, nb)
 nopen = na - nb

 allocate(noon(nif), source=0.0d0)
 call read_eigenvalues_from_fch(fchname, nif, 'a', noon)
 if( ANY(noon<0.0d0) ) then
  write(iout,'(A)') 'ERROR in subroutine read_no_info_from_fch: there exists&
                   & negative occupation number(s), this is not possible.'
  write(iout,'(A)') 'Do you mistake the energy levels for occupation numbers?'
  stop
 end if

 nacto = 0; nacta = 0; nactb = 0

 do i = 1, nif, 1
  if(noon(i)>no_thres .and. noon(i)<(2.0d0-no_thres)) then
   nacto = nacto + 1
   if(i <= nb) then
    nacta = nacta + 1
    nactb = nactb + 1
   else if(i <=na) then
    nacta = nacta + 1
   end if
  end if
 end do ! for i

 deallocate(noon)
 ndb = na - nacta
 nacte = nacta + nactb
 return
end subroutine read_no_info_from_fch

