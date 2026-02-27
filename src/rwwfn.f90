! written by jxzou at 20200504: read/write basis, MOs or eigenvalues from/to given files
! updated by jxzou at 20201213: add subroutines for reading occupation numbers
! updated by jxzou at 20201213: fix bug: read CASCI NOONs of ORCA
! updated by jxzou at 20201214: add read MO subroutines of Molpro
! updated by jxzou at 20210128: add read int1e subroutines of Gaussian log

! modify the IROHF value in a given .fch(k) file
subroutine modify_irohf_in_fch(fchname, k)
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: k
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname

 buf = ' '
 call find_specified_suffix(fchname, '.fch', i)
 fchname1 = fchname(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5) == 'IROHF') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  close(fid)
  close(fid1,status='delete')
  return
 end if

 write(fid1,'(A5,38X,A1,16X,I1)') 'IROHF', 'I', k

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
end subroutine modify_irohf_in_fch

! modify the charge and spin multiplicity in a specified .fch(k) file
subroutine modify_charge_and_mult_in_fch(fchname, charge, mult)
 implicit none
 integer :: i, j, charge0, mult0, ne, na, fid, fid1, RENAME
 integer, intent(in) :: charge, mult
!f2py intent(in) :: charge, mult
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 character(len=240) :: buf, fchname1

 call read_charge_and_mult_from_fch(fchname, charge0, mult0)
 if(charge0==charge .and. mult0==mult) return

 i = MOD(IABS(charge0-charge), 2)
 j = MOD(IABS(mult0-mult), 2)
 if(i /= j) then
  write(6,'(/,A)') 'ERROR in subroutine modify_charge_and_mult_in_fch: the spec&
                   &ified charge and'
  write(6,'(A)') 'mult are inconsistent.'
  write(6,'(2(A,I0))') 'charge=', charge, ', mult=', mult
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 call find_specified_suffix(fchname, '.fch', i)
 fchname1 = fchname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:6) == 'Charge') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 read(fid,'(A)') buf              ! Multiplicity
 read(fid,'(A49,2X,I10)') buf, ne ! Number of electrons
 ne = ne + charge0 - charge

 if(mult > ne+1) then
  write(6,'(/,A)') 'ERROR in subroutine modify_charge_and_mult_in_fch: the spec&
                   &ified mult is'
  write(6,'(A)') 'too large.'
  write(6,'(2(A,I0))') 'charge=', charge, ', mult=', mult
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 na = (ne + mult - 1)/2
 write(fid1,'(A,I17)') 'Charge                                     I',charge
 write(fid1,'(A,I17)') 'Multiplicity                               I',mult
 write(fid1,'(A,I17)') 'Number of electrons                        I',ne
 write(fid1,'(A,I17)') 'Number of alpha electrons                  I',na
 write(fid1,'(A,I17)') 'Number of beta electrons                   I',ne-na

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:16) == 'Number of beta e') exit
 end do ! for while

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))

 ! GaussView checks some information
 if(mult0>1 .and. mult==1) then
  call del_dm_in_fch(fchname, 2)
  call modify_irohf_in_fch(fchname, 0)
 end if
end subroutine modify_charge_and_mult_in_fch

! read nalpha and nbeta from .fch(k) file
subroutine read_na_and_nb_from_fch(fchname, na, nb)
 implicit none
 integer :: i, fid
 integer, intent(out) :: na, nb
!f2py intent(out) :: na, nb
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 na = 0; nb = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:15) == 'Number of alpha') exit
 end do

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_na_and_nb_from_fch: no 'Number of &
                   &alpha' found in"
  write(6,'(A)') 'file '//TRIM(fchname)
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, na
 read(fid,'(A49,2X,I10)') buf, nb
 close(fid)
end subroutine read_na_and_nb_from_fch

! read Alpha/Beta MOs from a given .fch(k) file
subroutine read_mo_from_fch(fchname, nbf, nif, ab, mo)
 implicit none
 integer :: i, j, fid, ncoeff
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(out) :: mo(nbf,nif)
!f2py intent(out) :: mo
!f2py depend(nbf,nif) :: mo
 real(kind=8), allocatable :: coeff(:)
 character(len=1), intent(in) :: ab
!f2py intent(in) :: ab
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 character(len=8) :: key
 character(len=8), parameter :: key1 = 'Alpha MO'
 character(len=7), parameter :: key2 = 'Beta MO'

 mo = 0d0; key = key1
 if(ab/='a' .and. ab/='A') key = key2//' '

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == key) exit
 end do

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mo_from_fch: no '"//key//"' found &
                   &in file "//TRIM(fchname)
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, ncoeff
 if(ncoeff /= nbf*nif) then
  write(6,'(/,A)') 'ERROR in subroutine read_mo_from_fch: ncoeff /= nbf*nif.'
  write(6,'(A)') 'Inconsistency found between input nbf,nif and those in file&
                 & '//TRIM(fchname)
  write(6,'(A,I10,2I5)') 'ncoeff, nbf, nif=', ncoeff, nbf, nif
  close(fid)
  stop
 end if

 allocate(coeff(ncoeff), source=0d0)
 read(fid,'(5(1X,ES15.8))',iostat=j) (coeff(i),i=1,ncoeff)

 if(j /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mo_from_fch: failed to read MOs.'
  write(6,'(A)') 'This file seems problematic: '//TRIM(fchname)
  close(fid)
  stop
 end if

 mo = RESHAPE(coeff, (/nbf,nif/))
 deallocate(coeff)
 close(fid)
end subroutine read_mo_from_fch

! read Alpha/Beta MOs from the output file generated by Gaussian utility chkchk
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

 call open_file(txtname, .true., chkid)
 do while(.true.)
  read(chkid,'(A)') buf
  if(buf(7:14) == key) exit
 end do
 BACKSPACE(chkid)

 mo = 0d0
 do i = 1, nif, 1
  read(chkid,'(A)') buf
  read(chkid,'(4D20.12)') (mo(k,i), k=1,nbf)
 end do

 close(chkid)
end subroutine read_mo_from_chk_txt

! read Alpha/Beta MOs from .Orb file of (Open)Molcas
subroutine read_mo_from_orb(orbname, nbf, nif, ab, mo)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(out) :: mo(nbf,nif)
 character(len=1), intent(in) :: ab
 character(len=240) :: buf
 character(len=240), intent(in) :: orbname
 character(len=4) :: key
 character(len=4), parameter :: key1 = '#ORB', key2 = '#UOR'

 key = key1
 if(ab/='a' .and. ab/='A') key = key2

 call open_file(orbname, .true., fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == key) exit
 end do

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_mo_from_orb: no '"//key//"' found&
                 & in file "//TRIM(orbname)//'.'
  stop
 end if

 mo = 0d0
 do i = 1, nif, 1
  read(fid,'(A)') buf
  read(fid,'(5(1X,ES21.14))') (mo(j,i),j=1,nbf)
!  read(fid,'(4E18.12)') (mo(j,i),j=1,nbf) ! MOLCAS 8.0
 end do ! for i

 close(fid)
end subroutine read_mo_from_orb

! read Alpha/Beta MO coefficients from Molpro .xml file
subroutine read_mo_from_xml(xmlname, nbf, nif, ab, mo)
 implicit none
 integer :: i, j, k, fid, nline
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(out) :: mo(nbf,nif)
 character(len=240) :: buf
 character(len=1), intent(in) :: ab
 character(len=240), intent(in) :: xmlname

 mo = 0d0
 open(newunit=fid,file=TRIM(xmlname),status='old',position='rewind')

 if(ab=='a' .or. ab=='A') then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(INDEX(buf,'type="ALPH')>0 .or. INDEX(buf,'type="CANO')>0 .or. &
      INDEX(buf,'type="NATU')>0 .or. buf(65:74)=='type="RAW"') exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_mo_from_xml: none of 'ALPH'/'CANO&
                    &'/'NATU' is"
   write(6,'(A)') 'found in file '//TRIM(xmlname)
   close(fid)
   stop
  end if

 else ! UHF

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(INDEX(buf,'type="BETA') > 0) exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine read_mo_from_xml: no type="BETA found &
                    &in file '//TRIM(xmlname)
   close(fid)
   stop
  end if
 end if

 read(fid,'(A)') buf
 nline = nif/10
 if(nif-10*nline > 0) nline = nline + 1

 do i = 1, nif, 1
  do while(.true.)
   read(fid,'(A)') buf
   k = LEN_TRIM(buf)
   if(buf(k:k) == '>') exit
  end do ! for while

  do j = 1, nline, 1
   k = min(10*j,nif)
   read(fid,*) mo(10*j-9:k,i)
  end do ! for j
  read(fid,'(A)') buf
 end do ! for i

 close(fid)
end subroutine read_mo_from_xml

! read Alpha/Beta MOs from orbital file of BDF
subroutine read_mo_from_bdf_orb(orbname, nbf, nif, ab, mo)
 implicit none
 integer :: i, j, k, nline, fid
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(out) :: mo(nbf,nif)
 character(len=1), intent(in) :: ab
 character(len=240) :: buf
 character(len=240), intent(in) :: orbname
 character(len=4) :: key
 character(len=4), parameter :: key1 = 'ALPH', key2 = 'BETA'

 if(ab=='a' .or. ab=='A') then
  key = key1
 else if(ab=='b' .or. ab=='B') then
  key = key2
 else
  write(6,'(/,A)') 'ERROR in subroutine read_mo_from_bdf_orb: invalid ab.'
  write(6,'(A)') 'ab = '//ab
  stop
 end if

 call open_file(orbname, .true., fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(24:27) == key) exit
 end do

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mo_from_bdf_orb: no '"//key//"'&
                   & found in file "//TRIM(orbname)
  close(fid)
  stop
 end if

 mo = 0d0
 nline = nbf/5
 if(nbf-5*nline > 0) nline = nline + 1

 do i = 1, nif, 1
  do j = 1, nline, 1
   k = min(5*j,nbf)
   read(fid,*) mo(5*j-4:k,i)
  end do ! for j
 end do ! for i

 close(fid)
end subroutine read_mo_from_bdf_orb

! read Alpha MOs from a DALTON.MOPUN format file
subroutine read_mo_from_dalton_mopun(orbname, nbf, nif, coeff)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: orbname
 real(kind=8), intent(out) :: coeff(nbf,nif)

 coeff = 0d0
 call open_file(orbname, .true., fid)
 read(fid,'(A)') buf

 do i = 1, nif, 1
  read(fid,'(4F18.14)',iostat=k) coeff(:,i)
  if(k /= 0) exit
 end do ! for i

 close(fid)
 if(k /= 0) then
  write(6,'(/,A)') 'ERROR in subrouine read_mo_from_dalton_mopun: failed&
                   & to read MOs from'
  write(6,'(A)') 'file '//TRIM(orbname)//'.'
  write(6,'(4(A,I0))') 'nbf=',nbf,',nif=',nif,',i=',i,',j=',j
  stop
 end if
end subroutine read_mo_from_dalton_mopun

! read Alpha MOs from a Turbomole mos file
subroutine read_mo_from_mos(fname, nbf, nif, coeff)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(out) :: coeff(nbf,nif)
!f2py intent(out) :: coeff
!f2py depend(nbf,nif) :: coeff
 character(len=240) :: buf
 character(len=240), intent(in) :: fname
!f2py intent(in) :: fname

 coeff = 0d0
 open(newunit=fid,file=TRIM(fname),status='old',position='rewind')
 read(fid,'(A)') buf

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:1) /= '#') exit
 end do ! for while
 BACKSPACE(fid)

 do i = 1, nif, 1
  read(fid,'(A)') buf
  read(fid,'(4D20.14)') coeff(:,i)
 end do ! for i

 close(fid)
end subroutine read_mo_from_mos

! read AO-basis 1-e integral matrix from a Gaussian output file
! stype: overlap, kinetic, potential, core
subroutine read_int1e_from_gau_log(logname, itype, nbf, mat)
 implicit none
 integer :: i, j, k, m, n, p, fid
 integer, intent(in) :: itype, nbf
 real(kind=8), intent(out) :: mat(nbf,nbf)
 character(len=7), parameter :: key(4) = &
  ['Overlap', 'Kinetic', '* Poten', '** Core']
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 if(itype<1 .or. itype>4) then
  write(6,'(A,I0)') 'ERROR in subroutine read_int1e_from_gau_log: invalid&
                    & itype = ', itype
  write(6,'(A)') 'Allowed values are 1/2/3/4 for Overlap/Kinetic/Potential/&
                 &Core Hamiltonian.'
  stop
 end if

 mat = 0d0
 call open_file(logname, .true., fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(6:12) == key(itype)) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_int1e_from_gau_log: no key '"&
                 & //key(itype)//"' found in file "//TRIM(logname)
  close(fid)
  stop
 end if

 n = nbf/5
 if(nbf-5*n > 0) n = n + 1

 do i = 1, n, 1
  read(fid,'(A)') buf
  k = 5*i - 4

  do j = k, nbf, 1
   m = min(4, j-k)
   read(fid,*) p, mat(k:k+m,j)
  end do ! for j
 end do ! for i

 close(fid)

 do i = 1, nbf-1, 1
  do j = i+1, nbf, 1
   mat(j,i) = mat(i,j)
  end do ! for j
 end do ! for i
end subroutine read_int1e_from_gau_log

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
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: S(nbf,nbf)

 S = 0d0
 call open_file(outname, .true., fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:30) == 'SO Integrals of type Mltpl  0') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_ovlp_from_molcas_out: no 'SO Integ&
                   &rals of type Mltpl  0'"
  write(6,'(A)') 'found in file '//TRIM(outname)
  stop
 end if

 do i = 1, 3
  read(fid,'(A)') buf
 end do ! for i

 i = INDEX(buf, 'x')
 read(buf(i+1:),*) j
 if(j /= nbf) then
  write(6,'(/,A)') 'ERROR in subroutine read_ovlp_from_molcas_out: inconsistent&
                   & nbf in orbital'
  write(6,'(A)') 'file and overlap file.'
  write(6,'(2(A,I5))') 'j=', j, ', nbf=', nbf
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
 call symmetrize_dmat(nbf, S)
end subroutine read_ovlp_from_molcas_out

! write/print eigenvalues/occupation numbers into a .fch(k) file
subroutine write_eigenvalues_to_fch(fchname, nif, ab, on, replace)
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: nif
!f2py intent(in) :: nif
 real(kind=8), intent(in) :: on(nif)
!f2py intent(in) :: on
!f2py depend(nif) :: on
 character(len=1), intent(in) :: ab
!f2py intent(in) :: ab
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 character(len=8) :: key
 character(len=8), parameter :: key1 = 'Alpha Or'
 character(len=7), parameter :: key2 = 'Beta Or'
 logical, intent(in) :: replace
!f2py intent(in) :: replace

 key = key1
 if(ab/='a' .and. ab/='A') key = key2//'b'

 call find_specified_suffix(fchname, '.fch', i)
 fchname1 = fchname(1:i-1)//'_D.fch'
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  if(buf(1:8) == key) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine write_eigenvalues_to_fch: no '"//&
                    key//"' found in file "//TRIM(fchname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(fid1,'(5(1X,ES15.8))') (on(i),i=1,nif)

 ! skip the Alpha/Beta Orbital Energies in fname1
 do while(.true.)
  read(fid,'(A)') buf
  if(INDEX(buf,'=') > 0) exit
 end do
 BACKSPACE(fid)

 ! copy remaining content in file fname1
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do
 close(fid1)

 if(replace) then
  close(fid,status='delete')
  i = RENAME(TRIM(fchname1), TRIM(fchname))
 else
  close(fid)
 end if
end subroutine write_eigenvalues_to_fch

! read Alpha/Beta occupation numbers from a given .orb file of MOLCAS/OpenMolcas
subroutine write_on_to_orb(orbname, nif, ab, on, replace)
 implicit none
 integer :: i, nline, fid1, fid2, RENAME
 integer, intent(in) :: nif
 real(kind=8), intent(in) :: on(nif)
 character(len=1), intent(in) :: ab
 character(len=240) :: buf, orbname1
 character(len=240), intent(in) :: orbname
 character(len=6) :: key
 character(len=4), parameter :: key1 = '#OCC'
 character(len=5), parameter :: key2 = '#UOCC', key3 = '#OCHR'
 character(len=6), parameter :: key4 = '#UOCHR'
 logical, intent(in) :: replace

 call open_file(orbname, .true., fid1)
 i = INDEX(orbname,'.',back=.true.)
 orbname1 = orbname(1:i-1)//'_D.Orb'
 open(newunit=fid2,file=TRIM(orbname1),status='replace')

 key = key1//'  '
 if(ab/='a' .and. ab/='A') key = key2//' '

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(INDEX(buf,TRIM(key)) /= 0) exit
 end do

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_on_from_orb: no '"//TRIM(key)//&
                  &"' found in"
  write(6,'(A)') 'file '//TRIM(orbname)
  stop
 end if

 read(fid1,'(A)') buf
 write(fid2,'(A)') TRIM(buf)
 write(fid2,'(5(1X,ES21.14))') (on(i),i=1,nif)

 nline = nif/5
 if(nif > 5*nline) nline = nline + 1
 ! skip the #OCC/#UOCC data in original file
 do i = 1, nline, 1 
  read(fid1,'(A)') buf
 end do ! for i

 key = key3//' '
 if(ab/='a' .and. ab/='A') key = key4

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(INDEX(buf,TRIM(key)) /= 0) exit
 end do

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_on_from_orb: no '"//TRIM(key)//&
                  &"' found in"
  write(6,'(A)') 'file '//TRIM(orbname)
  stop
 end if

 read(fid1,'(A)') buf
 write(fid2,'(A)') TRIM(buf)
 write(fid2,'(10(1X,F7.4))') (on(i),i=1,nif)

 nline = nif/10
 if(nif > 10*nline) nline = nline + 1
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
end subroutine write_on_to_orb

! write Alpha/Beta MOs into a given .fch(k) file
subroutine write_mo_into_fch(fchname, nbf, nif, ab, mo)
 implicit none
 integer :: i, fid1, fid2, ncoeff, RENAME
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: mo(nbf,nif)
!f2py intent(in) :: mo
!f2py depend(nbf,nif) :: mo
 real(kind=8), allocatable :: coeff(:)
 character(len=1), intent(in) :: ab
!f2py intent(in) :: ab
 character(len=240) :: buf, fchname1
 character(len=8) :: key
 character(len=8), parameter :: key1 = 'Alpha MO'
 character(len=7), parameter :: key2 = 'Beta MO'
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 key = key1
 if(ab/='a' .and. ab/='A') key = key2//' '

 call find_specified_suffix(fchname, '.fch', i)
 fchname1 = fchname(1:i-1)//'.t'
 open(newunit=fid1,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(buf(1:8) == key) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine write_mo_into_fch: no '"//key//"' found&
                   & in file "//TRIM(fchname)//'.'
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 BACKSPACE(fid1)
 read(fid1,'(A49,2X,I10)') buf, ncoeff
 if(ncoeff /= nbf*nif) then
  write(6,'(/,A)') 'ERROR in subroutine write_mo_into_fch: ncoeff /= nbf*nif.'
  write(6,'(A)') 'Inconsistency found between input nbf,nif and those in file&
                 & '//TRIM(fchname)//'.'
  write(6,'(A,I10,2I5)') 'ncoeff,nbf,nif=', ncoeff,nbf,nif
  stop
 end if

 allocate(coeff(ncoeff), source=0d0)
 coeff = RESHAPE(mo, (/ncoeff/))
 write(fid2,'(5(1X,ES15.8))') (coeff(i),i=1,ncoeff)
 deallocate(coeff)

 do while(.true.) ! skip MOs in fchname
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'=') > 0) exit
 end do ! for while

 if(i == 0) then
  BACKSPACE(fid1)
  do while(.true.)
   read(fid1,'(A)',iostat=i) buf
   if(i /= 0) exit
   write(fid2,'(A)') TRIM(buf)
  end do ! for while
 end if

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(fchname1, fchname)
end subroutine write_mo_into_fch

! write Alpha/Beta MOs into a PSI4 Matrix file
subroutine write_mo_into_psi_mat(matfile, nbf, nif, mo)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nbf, nif
 character(len=240), intent(in) :: matfile
 real(kind=8), intent(in) :: mo(nbf,nif)

 open(newunit=fid,file=TRIM(matfile),status='replace')
 write(fid,'(A)') 'MO coefficients (C)'
 write(fid,'(A)') 'symmetry 0'
 write(fid,'(I0)') nbf*nif

 do i = 1, nbf, 1
  do j = 1, nif, 1
   write(fid,'(A,2(1X,I5),1X,ES15.8)') ' 0',i-1,j-1,mo(i,j)
  end do ! for j
 end do ! for i

 close(fid)
end subroutine write_mo_into_psi_mat

! determine whether sperical harmonic or Cartesian functions are used in .fch(k) file
subroutine determine_sph_or_cart(fchname, cart)
 implicit none
 integer :: i, k, fid
 integer, allocatable :: shltyp(:)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 logical, intent(out) :: cart

 call open_file(fchname, .true., fid)
 ! find and read Shell types
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == 'Shell types') exit
 end do

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine determine_sph_or_cart: missing&
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
end subroutine determine_sph_or_cart

! read variables nbf, nif, ndb, etc from a .fch(k) file containing NOs and NOONs
subroutine read_no_info_from_fch(fchname, on_thres, nbf, nif, ndb, nopen, nacta,&
                                 nactb, nacto, nacte)
 implicit none
 integer :: i, na, nb
 integer, intent(out) :: nbf, nif, ndb, nopen, nacta, nactb, nacto, nacte
 real(kind=8), intent(in) :: on_thres
 real(kind=8), allocatable :: noon(:)
 character(len=240), intent(in) :: fchname

 nacto = 0; nacta = 0; nactb = 0
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 call read_na_and_nb_from_fch(fchname, na, nb)
 nopen = na - nb

 allocate(noon(nif))
 call read_eigenvalues_from_fch(fchname, nif, 'a', noon)
 if( ANY(noon < -1d-2) ) then
  write(6,'(/,A)') 'ERROR in subroutine read_no_info_from_fch: there exists neg&
                   &ative occupation'
  write(6,'(A)') 'number(s), this is not possible. Do you mistake the energy le&
                 &vels for occupation'
  write(6,'(A)') 'numbers? Or do you use relaxed density of MP2/CI/CC/TD- metho&
                 &ds?'
  stop
 end if

 if(on_thres<0d0 .or. on_thres>1d0) then
  write(6,'(/,A)') 'ERROR in subroutine read_no_info_from_fch: input on_thres i&
                   &s invalid.'
  write(6,'(A)') '0.0 < on_thres < 1.0 is required.'
  write(6,'(A,F12.6)') 'Current on_thres = ', on_thres
  stop
 end if

 do i = 1, nif, 1
  if(noon(i)>on_thres .and. noon(i)<(2d0-on_thres)) then
   nacto = nacto + 1
   if(i <= nb) then
    nacta = nacta + 1
    nactb = nactb + 1
   else if(i <= na) then
    nacta = nacta + 1
   end if
  else if(on_thres-noon(i)>0d0 .and. on_thres-noon(i)<1d-3) then
   write(6,'(/,A)') REPEAT('-',79)
   write(6,'(A)') 'Warning from subroutine read_no_info_from_fch: some occupati&
                  &on number is very'
   write(6,'(A)') 'close to ON_thres. You may consider enlarge the active space&
                  & size slightly.'
   write(6,'(A,I0)') 'Check No. orbital: ', i
   write(6,'(A)') REPEAT('-',79)
  end if
 end do ! for i

 deallocate(noon)
 ndb = na - nacta
 nacte = nacta + nactb
end subroutine read_no_info_from_fch

! read various AO density matrix from a .fch(k) file
subroutine read_dm_from_fch(fchname, itype, nbf, dm)
 implicit none
 integer :: i, j, fid, ncoeff1, ncoeff2
 integer, intent(in) :: itype, nbf
!f2py intent(in) :: itype, nbf
 real(kind=8), intent(out) :: dm(nbf,nbf)
!f2py intent(out) :: dm
!f2py depend(nbf) :: dm
 character(len=11), parameter :: key(12) = ['Total SCF D', 'Spin SCF De',&
   'Total CI De', 'Spin CI Den', 'Total MP2 D', 'Spin MP2 De',&
   'Total CC De', 'Spin CC Den', 'Total CI Rh', 'Spin CI Rho',&
   'Total 2nd O', 'Spin 2nd Or']
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 dm = 0d0
 if(itype<1 .or. itype>12) then
  write(6,'(/,A,I0)') 'ERROR in subroutine read_dm_from_fch: invalid itype=',itype
  write(6,'(A)') 'Allowed values are 1~12.'
  stop
 end if

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == key(itype)) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_dm_from_fch: no key '"//key(itype)&
                 //"' found in file "
  write(6,'(A)') TRIM(fchname)
  close(fid)
  stop
 end if

 ncoeff2 = (nbf*nbf + nbf)/2
 read(buf(50:),*) ncoeff1
 if(ncoeff1 /= ncoeff2) then
  write(6,'(/,A)') 'ERROR in subroutine read_dm_from_fch: inconsistent dimensio&
                   &n between the .fch(k)'
  write(6,'(A)') 'file and input nbf!'
  write(6,'(2(A,I0))') 'nbf=', nbf, ', ncoeff1=', ncoeff1
  close(fid)
  stop
 end if

 read(fid,'(5(1X,ES15.8))') ((dm(j,i),j=1,i),i=1,nbf)
 close(fid)
 call symmetrize_dmat(nbf, dm)
end subroutine read_dm_from_fch

! Write 'Total SCF Density' or 'Spin SCF Density' into a .fch(k) file.
! Note: (1) the dm(j,i) j<=i will be used; (2) the density matrix elements in
!  the array dm must be in Gaussian angular momentum order.
subroutine write_dm_into_fch(fchname, total, nbf, dm)
 implicit none
 integer :: i, j, k, fid, fid1, RENAME
 integer, intent(in) :: nbf
!f2py intent(in) :: nbf
 real(kind=8), intent(in) :: dm(nbf,nbf)
!f2py intent(in) :: dm
!f2py depend(nbf) :: dm
 character(len=11) :: key
 character(len=11), parameter :: key1 = 'Total SCF D'
 character(len=11), parameter :: key2 = 'Spin SCF De'
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical, intent(in) :: total
!f2py intent(in) :: total

 key = key1
 if(.not. total) key = key2

 i = INDEX(fchname, '.fch', back=.true.)
 fchname1 = fchname(1:i-1)//'.t'
 call open_file(fchname, .true., fid)
 open(newunit=fid1,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  if(buf(1:11) == key) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine write_dm_into_fch: no key '"//key//"' f&
                   &ound in file "//TRIM(fchname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(fid1,'(5(1X,ES15.8))') ((dm(j,i),j=1,i),i=1,nbf)

 ! skip density in the original .fch file
 k = nbf*(nbf+1)/2
 j = k/5
 if(k - 5*j > 0) j = j + 1
 do i = 1, j, 1
  read(fid,'(A)') buf
 end do ! for i

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
end subroutine write_dm_into_fch

subroutine copy_dm_between_fch(fchname1, fchname2, itype, total)
 implicit none
 integer :: nbf, nif
 integer, intent(in) :: itype
 real(kind=8), allocatable :: dm(:,:)
 character(len=240), intent(in) :: fchname1, fchname2
 logical, intent(in) :: total

 call read_nbf_and_nif_from_fch(fchname1, nbf, nif)
 allocate(dm(nbf,nbf))
 call read_dm_from_fch(fchname1, itype, nbf, dm)
 call write_dm_into_fch(fchname2, total, nbf, dm)
 deallocate(dm)
end subroutine copy_dm_between_fch

! detect whether there exists 'Spin SCF Density' in a given .fch file
subroutine detect_spin_scf_density_in_fch(fchname, alive)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 logical, intent(out) :: alive

 alive = .false.
 call open_file(fchname, .true., fid)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit

  if(buf(1:16) == 'Spin SCF Density') then
   alive = .true.
   close(fid)
   return
  end if
 end do ! for while

 close(fid)
end subroutine detect_spin_scf_density_in_fch

! add a given density string into a .fch(k) file
subroutine add_density_str_into_fch(fchname, itype)
 implicit none
 integer :: i, j, k, nbf, fid, fid1, RENAME
 integer, intent(in) :: itype ! 1/2 for 'Total SCF Density'/'Spin SCF Density'
!f2py intent(in) :: itype
 character(len=17), parameter :: str(2) = ['Total SCF Density','Spin SCF Density ']
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 if(itype<1 .or. itype>2) then
  write(6,'(/,A,I0)') 'ERROR in subroutine add_density_str_into_fch: invalid it&
                      &ype=', itype
  stop
 end if

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:17) == str(itype)) exit
  if(buf(1:15) == 'Number of basis') read(buf(50:),*) nbf
 end do ! for while

 if(i /= 0) then
  rewind(fid)
  call find_specified_suffix(fchname, '.fch', i)
  fchname1 = fchname(1:i-1)//'.t'
  open(newunit=fid1,file=TRIM(fchname1),status='replace')

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   write(fid1,'(A)') TRIM(buf)
   if(itype == 1) then
    if(buf(1:7) == 'Beta MO') exit
   else
    if(buf(1:7) == 'Total S') exit
   end if
  end do ! for while

  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine add_density_str_into_fch: key not foun&
                    &d in file '//TRIM(fchname)
   write(6,'(A,I0)') 'itype=', itype
   close(fid)
   close(fid1,status='delete')
   stop
  end if

  do while(.true.)
   read(fid,'(A)',iostat=j) buf ! use j here, not i or k
   if(j /= 0) exit
   if(buf(49:49) == '=') exit
   write(fid1,'(A)') TRIM(buf)
  end do ! for while

  k = nbf*(nbf+1)/2
  write(fid1,'(A,I12)') str(itype)//'                          R   N=', k
  write(fid1,'(5(1X,ES15.8))') (0d0,i=1,k)
  if(j == 0) write(fid1,'(A)') TRIM(buf)

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   write(fid1,'(A)') TRIM(buf)
  end do ! for while

  close(fid,status='delete')
  close(fid1)
  i = RENAME(TRIM(fchname1), TRIM(fchname))

 else
  close(fid)
 end if
end subroutine add_density_str_into_fch

! update 'Total SCF Density' (and 'Spin SCF Density') using Alpha(and Beta) MOs
! in a Gaussian .fch(k) file
subroutine update_density_using_mo_in_fch(fchname)
 implicit none
 integer :: i, fid, nbf, nif, na, nb
 real(kind=8), allocatable :: alpha_coeff(:,:), beta_coeff(:,:), dm_a(:,:), &
  dm_b(:,:), tot_dm(:,:), spin_dm(:,:), occ(:)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical :: uhf

 uhf = .false.
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == 'Beta O') then
   uhf = .true.
   exit
  end if
  if(buf(1:11) == 'Total SCF D') exit
 end do ! for while
 close(fid)

 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(tot_dm(nbf,nbf), source=0d0)
 allocate(alpha_coeff(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', alpha_coeff)
 call add_density_str_into_fch(fchname, 1)
 allocate(occ(nif), source=0d0)

 if(uhf) then
  if(na > 0) occ(1:na) = 1d0
  allocate(dm_a(nbf,nbf))
  call calc_dm_using_mo_and_on(nbf, nif, alpha_coeff, occ, dm_a)
  deallocate(alpha_coeff)

  allocate(dm_b(nbf,nbf))
  allocate(beta_coeff(nbf,nif))
  call read_mo_from_fch(fchname, nbf, nif, 'b', beta_coeff)
  if(na > nb) occ(nb+1:na) = 0d0
  call calc_dm_using_mo_and_on(nbf, nif, beta_coeff, occ, dm_b)
  deallocate(beta_coeff)

  tot_dm = dm_a + dm_b
  allocate(spin_dm(nbf,nbf), source=dm_a-dm_b)
  deallocate(dm_a, dm_b)
  call add_density_str_into_fch(fchname, 2)
  call write_dm_into_fch(fchname, .false., nbf, spin_dm)
  deallocate(spin_dm)

 else ! R(O)HF
  if(nb > 0) occ(1:nb) = 2d0
  if(na > nb) occ(nb+1:na) = 1d0
  call calc_dm_using_mo_and_on(nbf, nif, alpha_coeff, occ, tot_dm)
 end if

 deallocate(occ)
 call write_dm_into_fch(fchname, .true., nbf, tot_dm)
 deallocate(tot_dm)
end subroutine update_density_using_mo_in_fch

! update 'Total SCF Density' using natural orbitals and occupation numbers
subroutine update_density_using_no_and_on(fchname)
 implicit none
 integer :: nbf, nif
 real(kind=8), allocatable :: noon(:), coeff(:,:), dm(:,:)
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(noon(nif))
 call read_eigenvalues_from_fch(fchname, nif, 'a', noon)

 allocate(coeff(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)

 allocate(dm(nbf,nbf))
 call calc_dm_using_mo_and_on(nbf, nif, coeff, noon, dm)

 deallocate(coeff, noon)
 call write_dm_into_fch(fchname, .true., nbf, dm)
 deallocate(dm)
end subroutine update_density_using_no_and_on

! read Total/Alpha/Beta/Transition Density Matrix from Gaussian output
! If some density matrix occurs more than one time, only the first will be read
subroutine read_density_from_gau_log(logname, itype, nbf, dm)
 implicit none
 integer :: i, j, k, m, n, fid
 integer, intent(in) :: itype, nbf
 ! itype: 1/2/3 for Total/Alpha/Beta density matrix
 real(kind=8), intent(out) :: dm(nbf,nbf)
 character(len=7), parameter :: key(3) = ['Density','Alpha D','Beta De']
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 if(itype<1 .or. itype>3) then
  write(6,'(/,A,I0)') 'ERROR in subroutine read_density_from_gau_log: invalid i&
                      &type = ', itype
  write(6,'(A)') 'Allowed values are 1/2/3 for Total/Alpha/Beta density.'
  stop
 end if

 dm = 0d0
 call open_file(logname, .true., fid)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(6:12) == key(itype)) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_density_from_gau_log: no key '"&
                   &//key(itype)//"'"
  write(6,'(A)') 'found in file '//TRIM(logname)
  close(fid)
  stop
 end if

 n = nbf/5
 if(nbf-5*n > 0) n = n + 1

 do i = 1, n, 1
  read(fid,'(A)') buf
  k = 5*i - 4

  do j = k, nbf, 1
   m = min(4, j-k)
   read(fid,'(21X,5(2X,F8.5))') dm(k:k+m,j)
  end do ! for j
 end do ! for i

 close(fid)

 ! symmetrize the density matrix
 call symmetrize_dmat(nbf, dm)
end subroutine read_density_from_gau_log

! check whether UHF wave function in .fch(k) file is equivalent to RHF
! (by checking the max and average difference between alpha/beta MO coefficients)
subroutine check_if_uhf_equal_rhf(fchname, eq)
 implicit none
 integer :: nbf, nif
 real(kind=8) :: max_diff, ave_diff
 real(kind=8), parameter :: zero1 = 1D-2, zero2 = 1D-4
 real(kind=8), allocatable :: alpha_coeff(:,:), beta_coeff(:,:), diff(:,:)
 character(len=240), intent(in) :: fchname
 logical, intent(out) :: eq

 eq = .false.
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(alpha_coeff(nbf,nif), beta_coeff(nbf,nif), diff(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a',alpha_coeff)
 call read_mo_from_fch(fchname, nbf, nif, 'b', beta_coeff)

 diff = DABS(alpha_coeff - beta_coeff)
 max_diff = MAXVAL(diff)
 ave_diff = SUM(diff)/DBLE(nbf*nif)
 write(6,'(2(A,F12.6))') 'max_diff =', max_diff, ', ave_diff =', ave_diff
 if(max_diff<zero1 .and. ave_diff<zero2) eq = .true.
 deallocate(alpha_coeff, beta_coeff, diff)
end subroutine check_if_uhf_equal_rhf

subroutine copy_ev_and_mo_between_fch(fchname1, fchname2, ab)
 implicit none
 integer :: nbf, nif, nbf2, nif2
 real(kind=8), allocatable :: ev(:), mo(:,:)
 character(len=3), intent(in) :: ab
 character(len=240), intent(in) :: fchname1, fchname2

 call read_nbf_and_nif_from_fch(fchname1, nbf, nif)
 call read_nbf_and_nif_from_fch(fchname2, nbf2, nif2)

 if(nbf/=nbf2 .or. nif/=nif2) then
  write(6,'(/,A)') 'ERROR in subroutine copy_ev_and_mo_between_fch: inconsisten&
                   &t nbf and/or nif'
  write(6,'(A)') 'in two files:'
  write(6,'(A)') TRIM(fchname1)
  write(6,'(A)') TRIM(fchname2)
  stop
 end if

 allocate(ev(nif))
 call read_eigenvalues_from_fch(fchname1, nif, ab(2:2), ev)
 call write_eigenvalues_to_fch(fchname2, nif, ab(3:3), ev, .true.)
 deallocate(ev)

 allocate(mo(nbf,nif))
 call read_mo_from_fch(fchname1, nbf, nif, ab(2:2), mo)
 call write_mo_into_fch(fchname2, nbf, nif, ab(3:3), mo)
 deallocate(mo)
end subroutine copy_ev_and_mo_between_fch

! copy alpha(and beta) orbital energies, alpha(and beta) MOs, Total SCF Density,
! and Spin SCF Density from a .fch(k) file into another one
subroutine copy_orb_and_den_in_fch(fchname1, fchname2, deleted)
 implicit none
 integer :: i, fid, nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname1, fchname2
!f2py intent(in) :: fchname1, fchname2
 real(kind=8), allocatable :: ev(:), mo(:,:), dm(:,:)
 logical :: uhf
 logical, intent(in) :: deleted
!f2py intent(in) :: deleted

 uhf = .false.
 open(newunit=fid,file=TRIM(fchname1),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == 'Beta') then
   uhf = .true.; exit
  end if
 end do ! for while
 close(fid)

 call read_nbf_and_nif_from_fch(fchname1, nbf, nif)
 allocate(ev(nif))
 call read_eigenvalues_from_fch(fchname1, nif, 'a', ev)
 call write_eigenvalues_to_fch(fchname2, nif, 'a', ev, .true.)
 if(uhf) then
  call read_eigenvalues_from_fch(fchname1, nif, 'b', ev)
  call write_eigenvalues_to_fch(fchname2, nif, 'b', ev, .true.)
 end if
 deallocate(ev)

 allocate(mo(nbf,nif))
 call read_mo_from_fch(fchname1, nbf, nif, 'a', mo)
 call write_mo_into_fch(fchname2, nbf, nif, 'a', mo)
 if(uhf) then
  call read_mo_from_fch(fchname1, nbf, nif, 'b', mo)
  call write_mo_into_fch(fchname2, nbf, nif, 'b', mo)
 end if
 deallocate(mo)

 call add_density_str_into_fch(fchname2, 1)
 allocate(dm(nbf,nbf))
 call read_dm_from_fch(fchname1, 1, nbf, dm)
 call write_dm_into_fch(fchname2, .true., nbf, dm)
 if(uhf) then
  call add_density_str_into_fch(fchname2, 2)
  call read_dm_from_fch(fchname1, 2, nbf, dm)
  call write_dm_into_fch(fchname2, .false., nbf, dm)
 end if
 deallocate(dm)

 if(deleted) then
  open(newunit=fid,file=TRIM(fchname1),status='old')
  close(fid,status='delete')
 end if
end subroutine copy_orb_and_den_in_fch

! read AO-basis overlap integral from file .47
subroutine read_ao_ovlp_from_47(file47, nbf, S)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nbf
!f2py intent(in) :: nbf
 real(kind=8), intent(out) :: S(nbf,nbf)
!f2py intent(out) :: S
!f2py depend(nbf) :: S
 character(len=240) :: buf
 character(len=240), intent(in) :: file47
!f2py intent(in) :: file47

 S = 0d0
 open(newunit=fid,file=TRIM(file47),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:9) == '$OVERLAP') exit
 end do ! for while

 if(i /= 0) then
  write(fid,'(/,A)') 'ERROR in subroutine read_ao_ovlp_from_47: failed to read &
                     &AO overlap'
  write(fid,'(A)') 'from file '//TRIM(file47)
  close(fid)
  stop
 end if

 read(fid,'(2X,5E15.7)') ((S(j,i),j=1,i),i=1,nbf)
 close(fid)
 call symmetrize_dmat(nbf, S)
end subroutine read_ao_ovlp_from_47

! 1) compute 1e expectation values of MOs in file mo_fch, using information from
! file no_fch (which usually includes NOs)
! 2) sort paired MOs in mo_fch by 1e expectation values
subroutine get_1e_exp_and_sort_pair(mo_fch, no_fch, npair)
 implicit none
 integer :: i, k, na, nb, ndb, nopen, ncore, nbf, nif, SYSTEM
! here ndb is the number of doubly occupied orbitals in R(O)HF, ndb = nb
 integer, intent(in) :: npair
!f2py intent(in) :: npair
 real(kind=8) :: rtmp
 real(kind=8), allocatable :: coeff(:), mo(:,:), noon(:)
 character(len=240) :: dname
 character(len=240), intent(in) :: mo_fch, no_fch
!f2py intent(in) :: mo_fch, no_fch
 logical :: nat_orb, changed

 i = SYSTEM('solve_ON_matrix '//TRIM(mo_fch)//' '//TRIM(no_fch))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine get_1e_exp_and_sort_pair: failed to cal&
                   &l utility'
  write(6,'(A)') 'solve_ON_matrix. Did you forget to compile this utility?'
  write(6,'(A)') 'mo_fch='//TRIM(mo_fch)
  write(6,'(A)') 'no_fch='//TRIM(no_fch)
  stop
 end if

 call read_na_and_nb_from_fch(mo_fch, na, nb)
 nopen = na - nb; ndb = nb

 call find_specified_suffix(mo_fch, '.fch', i)
 dname = mo_fch(1:i-1)//'_D.txt'
 call require_file_exist(dname)
 call delete_file(TRIM(dname))

 call read_nbf_and_nif_from_fch(no_fch, nbf, nif)
 allocate(noon(nif))
 call read_eigenvalues_from_fch(no_fch, nif, 'a', noon)
 nat_orb = .false.
 if(ALL(noon<2.2d0) .and. ALL(noon>-0.04d0)) nat_orb = .true.
 deallocate(noon)
 write(6,'(/,A,L1,A,3I5)') 'In subroutine get_1e_exp_and_sort_pair: nat_orb=',&
                            nat_orb,', idx=',ndb, nopen, npair

 allocate(coeff(nbf), mo(nbf,nif), noon(nif))
 call read_mo_from_fch(mo_fch, nbf, nif, 'a', mo)
 call read_eigenvalues_from_fch(mo_fch, nif, 'a', noon)

 ! reorder paired MOs as the ascending/descending order of diagonal elements
 ! For NOs, use descending order; for MOs, use ascending order
 k = 2*ndb + nopen + 1
 ncore = ndb - npair

 if(npair > 1) then
  if(nat_orb) then ! assuming natural orbitals
   do while(.true.)
    changed = .false.
    do i = ncore+1, ndb-1, 1
     if(noon(i) < noon(i+1)) then
      rtmp = noon(i); noon(i) = noon(i+1); noon(i+1) = rtmp
      rtmp = noon(k-i); noon(k-i) = noon(k-i-1); noon(k-i-1) = rtmp
      coeff = mo(:,i); mo(:,i) = mo(:,i+1); mo(:,i+1) = coeff
      coeff = mo(:,k-i); mo(:,k-i) = mo(:,k-i-1); mo(:,k-i-1) = coeff
      changed = .true.
     end if
    end do ! for i
    if(.not. changed) exit
   end do ! for while

  else ! assuming MOs with <i|Fock|i>
   do while(.true.)   ! sort paired MOs
    changed = .false.
    do i = ncore+1, ndb-1, 1
     if(noon(i) > noon(i+1)) then
      rtmp = noon(i); noon(i) = noon(i+1); noon(i+1) = rtmp
      rtmp = noon(k-i); noon(k-i) = noon(k-i-1); noon(k-i-1) = rtmp
      coeff = mo(:,i); mo(:,i) = mo(:,i+1); mo(:,i+1) = coeff
      coeff = mo(:,k-i); mo(:,k-i) = mo(:,k-i-1); mo(:,k-i-1) = coeff
      changed = .true.
     end if
    end do ! for i
    if(.not. changed) exit
   end do ! for while

   do while(.true.)   ! sort core MOs
    changed = .false.
    do i = 1, ncore-1, 1
     if(noon(i) > noon(i+1)) then
      rtmp = noon(i); noon(i) = noon(i+1); noon(i+1) = rtmp
      coeff = mo(:,i); mo(:,i) = mo(:,i+1); mo(:,i+1) = coeff
      changed = .true.
     end if
    end do ! for i
    if(.not. changed) exit
   end do ! for while
  end if
 end if

 deallocate(coeff)
 call write_mo_into_fch(mo_fch, nbf, nif, 'a', mo)
 deallocate(mo)
 call write_eigenvalues_to_fch(mo_fch, nif, 'a', noon, .true.)
 deallocate(noon)
end subroutine get_1e_exp_and_sort_pair

! sort NOs by ascending order of NOONs in a .fch file
subroutine sort_no_by_noon(fchname, i1, i2)
 implicit none
 integer :: i, j, nbf, nif
 integer, intent(in) :: i1, i2 ! initial/final index of active orbitals
!f2py intent(in) :: i1, i2
 real(kind=8) :: r1, r2
 real(kind=8), allocatable :: noon(:), coeff(:,:), rtmp(:)
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 if(i1 == i2) return
 if(i1<1 .or. i2<1 .or. i1>i2) then
  write(6,'(/,A)') 'ERROR in subroutine sort_no_by_noon: wrong orbital indices!'
  write(6,'(2(A,I0))') 'i1=', i1, ', i2=', i2
  stop
 end if
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(noon(nif))
 call read_eigenvalues_from_fch(fchname, nif, 'a', noon)
 allocate(coeff(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)
 allocate(rtmp(nbf))

 do i = i1, i2-1, 1
  r1 = noon(i)
  do j = i+1, i2, 1
   r2 = noon(j)
   if(r2 > r1) then
    r1 = r2; noon(j) = noon(i); noon(i) = r1
    rtmp = coeff(:,j); coeff(:,j) = coeff(:,i); coeff(:,i) = rtmp
   end if
  end do ! for j
 end do ! for i

 deallocate(rtmp)
 call write_mo_into_fch(fchname, nbf, nif, 'a', coeff)
 deallocate(coeff)
 call write_eigenvalues_to_fch(fchname, nif, 'a', noon, .true.)
 deallocate(noon)
end subroutine sort_no_by_noon

! Get the index of core-valence orbital separation. The input file fchname must
! include orbital energies (or Fock expectation values) in 'Alpha Or'.
subroutine get_core_valence_sep_idx(fchname, idx)
 implicit none
 integer :: nbf, nif, na, nb
 integer, intent(out) :: idx
!f2py intent(out) :: idx
 real(kind=8), parameter :: thres = 0.3d0
 real(kind=8), allocatable :: ev(:)
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(ev(nif))
 call read_eigenvalues_from_fch(fchname, nif, 'a', ev)
 call read_na_and_nb_from_fch(fchname, na, nb)

 do idx = nb, 2, -1
  if(DABS(ev(idx) - ev(idx-1)) > thres) exit
 end do ! for i

 deallocate(ev)
end subroutine get_core_valence_sep_idx

! Paring each singly occupied orbital with a virtual one, for a high spin GVB
!  wave function. These new pairs would be put after all normal GVB pairs.
! Such type of .fch file may be used for high-spin GVB-BCCC calculations.
subroutine pairing_open_with_vir(fchname)
 implicit none
 integer :: i, j, k, npair, nopen, ndb, na, nb, nbf, nif
 real(kind=8), allocatable :: mo(:,:), new_mo(:,:), ev(:)
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 i = INDEX(fchname, 'gvb', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine pairing_open_with_vir: 'gvb' characters&
                  & not found in"
  write(6,'(A)') 'filename '//TRIM(fchname)
  write(6,'(A)') 'These characters are used to identify the number of pairs.'
  write(6,'(/,A)') "Example: pairing_open_with_vir('h2o_uhf_uno_asrot2gvb3.fch')"
  write(6,'(/,A)') 'Note: do not use h2o_uhf_uno_asrot2gvb3_s.fch'
  stop
 end if

 k = LEN_TRIM(fchname)
 read(fchname(i+3:k-4),*,iostat=j) npair
 if(j /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine pairing_open_with_vir: failed to read n&
                   &pair from'
  write(6,'(A)') 'filename '//TRIM(fchname)
  stop
 end if

 call read_na_and_nb_from_fch(fchname, na, nb)
 nopen = na - nb
 if(nopen == 0) then
  write(6,'(/,A)') 'Warning from subroutine pairing_open_with_vir: this is a si&
                   &nglet .fch file.'
  write(6,'(A)') 'No singly occupied orbitals. Nothing to do.'
  return
 end if
 ndb = na - nopen - npair

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(ev(nif), source=0d0)

 if(ndb > 0) ev(1:ndb) = 2d0 ! doubly occupied MOs
 forall(i = 1:npair)         ! normal GVB pairs
  ev(ndb+2*i-1) = 1.8d0
  ev(ndb+2*i) = 0.2d0
 end forall
 ! 1.8/0.2 are not real occupation numbers, just to remind the user that where
 ! are the normal pairs

 k = ndb + 2*npair
 forall(i = 1:nopen) ev(k+2*i-1) = 1d0 ! singly occupied orbitals

 call write_eigenvalues_to_fch(fchname, nif, 'a', ev, .true.)
 deallocate(ev)

 allocate(mo(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo)
 allocate(new_mo(nbf,nif), source=mo)

 ! normal GVB pairs
 j = ndb + 1
 k = ndb + 2*npair
 new_mo(:,j:k) = mo(:,j+nopen:k+nopen)

 ! new GVB pairs (socc-vir)
 forall(i = 1:nopen)
  new_mo(:,k+2*i-1) = mo(:,ndb+i)
  new_mo(:,k+2*i) = mo(:,na+npair+i)
 end forall

 deallocate(mo)
 call write_mo_into_fch(fchname, nbf, nif, 'a', new_mo)
 deallocate(new_mo)
end subroutine pairing_open_with_vir

! Reorder MOs in the file gvbN_s.fch. A new file gvbN_new.fch would be generated.
! The order of MOs in gvbN_s.fch is expected to be
!  docc-bonding1-bonding2-...-socc-...-antibonding2-antibonding1-vir.
! The order of MOs in gvbN_new.fch would be
!  docc-bonding1-antibonding1-bonding2-antibonding2-...-socc-vir.
subroutine reorder2dbabasv(fchname)
 implicit none
 integer :: i, j, k, nbf, nif, ndb, npair, nopen, na, nb
 real(kind=8), allocatable :: mo(:,:), new_mo(:,:), ev(:), new_ev(:)
 character(len=240) :: new_fch
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 i = LEN_TRIM(fchname)
 if(fchname(i-5:i) /= '_s.fch') then
  write(6,'(/,A)') "ERROR in subroutine reorder2dbabasv: '_s.fch' suffix is not&
                   & found in"
  write(6,'(A)') 'filename '//TRIM(fchname)
  write(6,'(/,A,/)') 'Example: ben_triplet_uhf_uno_asrot2gvb2_s.fch'
  stop
 end if
 new_fch = fchname(1:i-5)//'new.fch'

 k = INDEX(fchname, 'gvb', back=.true.)
 if(k == 0) k = INDEX(fchname, 'GVB', back=.true.)
 if(k == 0) then
  write(6,'(/,A)') "ERROR in subroutine reorder2dbabasv: 'gvb'/'GVB' key not fo&
                   &und in filename "//TRIM(fchname)
  write(6,'(A)') 'Example 1: ben_triplet_uhf_uno_asrot2gvb2_s.fch'
  write(6,'(A)') 'Example 2: ben_triplet_FcGVB14_s.fch'
  stop
 end if

 read(fchname(k+3:i-6),*,iostat=j) npair
 if(j /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine reorder2dbabasv: failed to read npair &
                   &from filename '//TRIM(fchname)
  write(6,'(A)') 'Example: ben_triplet_uhf_uno_asrot2gvb2_s.fch'
  stop
 end if

 call read_na_and_nb_from_fch(fchname, na, nb)
 nopen = na - nb
 if(nopen == 0) then
  write(6,'(/,A)') 'Warning from subroutine reorder2dbabasv: this is a singlet &
                   &.fch file.'
  write(6,'(A)') 'No singly occupied orbitals. Only pair orbitals will be reord&
                 &ered.'
 end if
 ndb = na - nopen - npair

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(mo(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo)
 allocate(ev(nif))
 call read_eigenvalues_from_fch(fchname, nif, 'a', ev)

 allocate(new_ev(nif), new_mo(nbf,nif))
 new_ev = 0d0
 new_mo(:,na+npair+1:nif) = mo(:,na+npair+1:nif) ! virtual orbitals

 if(ndb > 0) then     ! doubly occupied orbitals
  new_ev(1:ndb) = 2d0
  new_mo(:,1:ndb) = mo(:,1:ndb)
 end if

 do i = 1, npair, 1   ! GVB pair orbitals
  new_ev(ndb+2*i-1) = ev(ndb+i)
  new_ev(ndb+2*i) = ev(na+npair-i+1)
  new_mo(:,ndb+2*i-1) = mo(:,ndb+i)
  new_mo(:,ndb+2*i) = mo(:,na+npair-i+1)
 end do ! for i
 deallocate(ev)

 if(nopen > 0) then   ! singly occupied orbitals
  k = ndb + 2*npair
  new_ev(k+1:k+nopen) = 1d0
  new_mo(:,k+1:k+nopen) = mo(:,nb+1:na)
 end if

 deallocate(mo)
 call copy_file(fchname, new_fch, .false.)
 call write_mo_into_fch(new_fch, nbf, nif, 'a', new_mo)
 deallocate(new_mo)
 call write_eigenvalues_to_fch(new_fch, nif, 'a', new_ev, .true.)
 deallocate(new_ev)
end subroutine reorder2dbabasv

! Convert an R(O)HF-type .fch(k) file into a UHF-type one. If ibrosym > 0,
! alpha/beta spatial symmetry will be broken via mixing ibrosym pair of spin
! orbitals.
subroutine fch_r2u(fchname, ibrosym)
 implicit none
 integer :: i, j, k, nbf, nif, ncoeff, na, nb, nline, len_dm, ibrok, fid, fid1
 integer, intent(in) :: ibrosym
!f2py intent(in) :: ibrosym
 real(kind=8), parameter :: fac = 0.5d0*DSQRT(2d0)
 real(kind=8), allocatable :: e_a(:), tmp_mo(:,:), occ_a(:), occ_b(:), &
  mo_a(:,:), mo_b(:,:), spin_dm(:,:)
 character(len=29), parameter :: error_str = 'ERROR in subroutine fch_r2u: '
 character(len=240) :: buf, uhf_fch
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call find_specified_suffix(fchname, '.fch', i)
 uhf_fch = fchname(1:i-1)//'_u.fch'

 call read_na_and_nb_from_fch(fchname, na, nb)
 if(ibrosym > nb) then
  write(6,'(/,A)') error_str//'input ibrosym is nonsense.'
  write(6,'(2(A,I0))') 'ibrosym=', ibrosym, ', nb=', nb
  stop
 end if

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(e_a(nif), mo_a(nbf,nif))
 call read_eigenvalues_from_fch(fchname, nif, 'a', e_a)
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo_a)
 !write(6,'(A,I0)') 'ibrosym=', ibrosym

 if(ibrosym > 0) then
  i = MINVAL([na, nb, nif-na, nif-nb])
  if(ibrosym > i) then
   write(6,'(/,A)') 'Remark from subroutine fch_r2u: ibrosym is larger than max&
                    &imum pair for mixing'
   write(6,'(A,I0)') 'spin orbitals. Automatically adjust ibrosym to ', i
   ibrok = i
  else
   ibrok = ibrosym
  end if
  allocate(mo_b(nbf,nif), source=mo_a)
  j = ibrok; k = 2*j
  allocate(tmp_mo(nbf,k))

  do i = 1, k, 1
   tmp_mo(:,i) = mo_a(:,na+i-j)
  end do ! for i
  do i = 1, j, 1
   mo_a(:,na-j+i) = fac*(tmp_mo(:,i) + tmp_mo(:,k-i+1))
   mo_a(:,na+i) = fac*(tmp_mo(:,j-i+1) - tmp_mo(:,k+i-j))
  end do ! for i

  do i = 1, k, 1
   tmp_mo(:,i) = mo_b(:,nb+i-j)
  end do ! for i
  do i = 1, j, 1
   mo_b(:,nb-j+i) = fac*(tmp_mo(:,i) - tmp_mo(:,k-i+1))
   mo_b(:,nb+i) = fac*(tmp_mo(:,j-i+1) + tmp_mo(:,k+i-j))
  end do ! for i
  deallocate(tmp_mo)
 end if

 allocate(occ_a(nif), source=0d0)
 occ_a(1:na) = 1d0
 allocate(occ_b(nif), source=0d0)
 if(nb > 0) occ_a(1:nb) = 1d0
 allocate(spin_dm(nbf,nbf))
 call calc_spin_dm_using_mo_and_on(nbf, nif, occ_a, occ_b, mo_a, mo_b, spin_dm)
 deallocate(occ_a, occ_b)

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(uhf_fch),status='replace')
 read(fid,'(A)') buf
 write(fid1,'(A)') TRIM(buf)

 read(fid,'(A)') buf
 k = LEN_TRIM(buf)

 if(k > 5) then
  i = INDEX(buf, ' ')
  do j = i+1, k, 1
   if(buf(j:j) /= ' ') exit
  end do ! for j
  if(j < k) then
   i = j - 1
   if(buf(i+1:i+2) == 'RO') then
    buf = buf(1:i)//'U'//TRIM(buf(i+3:))
   else if(buf(i+1:i+1) == 'R') then
    buf = buf(1:i)//'U'//TRIM(buf(i+2:))
   end if
  end if
 end if
 write(fid1,'(A)') TRIM(buf)

 ! set the 1st element of the integer array ILSW as 1
 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)',iostat=i) TRIM(buf)
  if(i /= 0) exit
  if(buf(1:4) == 'ILSW') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_str//"'ILSW' not found in file "//TRIM(fchname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if
 read(fid,'(A)') buf
 buf(12:12) = '1'
 write(fid1,'(A)') TRIM(buf)

 i = 0
 ! set IOpCl as 1
 do while(.true.)
  read(fid,'(A)') buf
  select case(buf(1:5))
  case('IOpCl')
   i = 1; exit
  case('IROHF')
   i = 2; exit
  case('Alpha')
   i = 3; exit
  end select
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 write(fid1,'(A,38X,A,16X,A)') 'IOpCl','I','1'
 write(fid1,'(A,38X,A,16X,A)') 'IROHF','I','0'

 select case(i)
 case(1)
  read(fid,'(A)') buf
  if(buf(1:5) /= 'IROHF') then
   BACKSPACE(fid)
   BACKSPACE(fid)
  end if
 case(2) ! do nothing
 case(3)
  BACKSPACE(fid)
  BACKSPACE(fid)
 case default
  write(6,'(/,A,I0)') error_str//'invalid case i=', i
  write(6,'(A)') 'fchname='//TRIM(fchname)
  close(fid)
  close(fid1,status='delete')
  stop
 end select

 ! find Alpha Orbital Energies
 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(buf(1:7) == 'Alpha O') exit
 end do ! for while

 nline = nif/5
 if(nif > nline*5) nline = nline + 1

 do i = 1, nline, 1
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
 end do ! for i

 write(fid1,'(A,22X,A,I12)') 'Beta Orbital Energies','R   N=', nif
 write(fid1,'(5(1X,ES15.8))') e_a
 deallocate(e_a)

 ! find Alpha MO coefficients
 read(fid,'(A)') buf
 write(fid1,'(A)') TRIM(buf)
 if(buf(1:7) /= 'Alpha M') then
  write(6,'(/,A)') error_str//'"Alpha M" is not found at the expected position.'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 ncoeff = nbf*nif
 nline = ncoeff/5
 if(ncoeff > nline*5) nline = nline + 1

 if(ibrok > 0) then
  write(fid1,'(5(1X,ES15.8))') mo_a
  deallocate(mo_a)
  write(fid1,'(A,23X,A,I12)') 'Beta MO coefficients','R   N=', ncoeff
  write(fid1,'(5(1X,ES15.8))') mo_b
  deallocate(mo_b)
  do i = 1, nline, 1 ! skip Alpha MO coefficients
   read(fid,'(A)') buf
  end do ! for i
 else
  do i = 1, nline, 1 ! simply copy Alpha MO coefficients
   read(fid,'(A)') buf
   write(fid1,'(A)') TRIM(buf)
  end do ! for i
  write(fid1,'(A,23X,A,I12)') 'Beta MO coefficients','R   N=', ncoeff
  write(fid1,'(5(1X,ES15.8))') mo_a
  deallocate(mo_a)
 end if

 ! For g09, `Spin SCF Density` section is required otherwise the Gaussian unfchk
 ! utility will signal errors. So we try to add this section.
 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)',iostat=i) TRIM(buf)
  if(i /= 0) exit
  if(buf(1:11) == 'Total SCF D') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_str//'"Total SCF D" not found in file '//TRIM(fchname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 len_dm = nbf*(nbf+1)/2
 nline = len_dm/5
 if(5*nline < len_dm) nline = nline + 1
 ! skip "Total SCF D" section
 do i = 1, nline, 1
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
 end do ! for i
 read(fid,'(A)',iostat=k) buf

 if(k == 0) then
  if(buf(1:10) == 'Spin SCF D') then
   write(fid1,'(A)') TRIM(buf)
   write(fid1,'(5(1X,ES15.8))') ((spin_dm(j,i),j=1,i),i=1,nbf)
   do i = 1, nline, 1
    read(fid,'(A)') buf
   end do ! for i
  else ! other strings like 'QEq coupling tensors'
   write(fid1,'(A,I12)') 'Spin SCF Density                           R   N=',len_dm
   write(fid1,'(5(1X,ES15.8))') ((spin_dm(j,i),j=1,i),i=1,nbf)
   write(fid1,'(A)') TRIM(buf)
  end if
 else ! assuming there is nothing after 'Total SCF D'
  close(fid)
  write(fid1,'(A,I12)') 'Spin SCF Density                           R   N=',len_dm
  write(fid1,'(5(1X,ES15.8))') ((spin_dm(j,i),j=1,i),i=1,nbf)
  close(fid1)
  return
 end if

 ! copy remaining content
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid)
 close(fid1)
end subroutine fch_r2u

