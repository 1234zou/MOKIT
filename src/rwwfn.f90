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

! read the total charge and the spin mltiplicity from a given .fch(k) file
subroutine read_charge_and_mult_from_fch(fchname, charge, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: charge, mult
!fp2y intent(out) :: charge, mult
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 charge = 0; mult = 1
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == 'Charge') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_charge_and_mult_from_fch: no 'Char&
                   &ge' found in"
  write(6,'(A)') 'file '//TRIM(fchname)
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, charge
 read(fid,'(A49,2X,I10)') buf, mult
 close(fid)
end subroutine read_charge_and_mult_from_fch

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

! read the total charge and the spin mltiplicity from a given .mkl file
subroutine read_charge_and_mult_from_mkl(mklname, charge, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: charge, mult
!f2py intent(out) :: charge, mult
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname
!f2py intent(in) :: mklname

 charge = 0; mult = 1
 call require_file_exist(mklname)
 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == '$CHAR_M') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_charge_and_mult_from_mkl: no '$CHA&
                   &R_M' found in"
  write(6,'(A)') 'file '//TRIM(mklname)
  close(fid)
  stop
 end if

 read(fid,*) charge, mult
 close(fid)
end subroutine read_charge_and_mult_from_mkl

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

! read nbf and nif from .fch(k) file
subroutine read_nbf_and_nif_from_fch(fchname, nbf, nif)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nbf, nif
!f2py intent(out) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call open_file(fchname, .true., fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:17) == 'Number of basis f') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') "ERROR in subroutine read_nbf_and_nif_from_fch: 'Number of b&
                   &asis f' not found in"
  write(6,'(A)') 'file '//TRIM(fchname)
  stop
 end if
 read(buf(52:),*) nbf

 ! In case that this is a Q-Chem .fch file, let's assume nif is not just below
 ! nbf
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) then
   close(fid)
   write(6,'(/,A)') "ERROR in subroutine read_nbf_and_nif_from_fch: 'Number of &
                    &indepen' not found in"
   write(6,'(A)') 'file '//TRIM(fchname)
   stop
  end if
  if(buf(1:17) == 'Number of indepen') exit
 end do ! for while

 read(buf(52:),*) nif
 close(fid)
end subroutine read_nbf_and_nif_from_fch

! read nbf and nif from .Orb file of MOLCAS/OpenMOLCAS
subroutine read_nbf_and_nif_from_orb(orbname, nbf, nif)
 implicit none
 integer :: fid
 integer, intent(out) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: orbname

 call open_file(orbname, .true., fid)
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:5) == '#INFO') exit
 end do

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 read(fid,*) nbf
 read(fid,*) nif
 close(fid)
end subroutine read_nbf_and_nif_from_orb

! read nbf from a GAMESS .dat file
subroutine read_cart_nbf_from_dat(datname, nbf)
 implicit none
 integer :: i, j, fid
 integer, intent(out) :: nbf
 character(len=240) :: buf
 character(len=240), intent(in) :: datname

 nbf = 0
 open(newunit=fid,file=TRIM(datname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:2) == '$') then
   call upper(buf(3:5))
   if(buf(2:5) == '$VEC') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nbf_from_dat: no '$VEC' found in f&
                   &ile "//TRIM(datname)
  close(fid)
  stop
 end if

 j = 0
 do while(.true.)
  read(fid,'(A)') buf
  read(buf(1:2),*) i
  if(i == 2) exit
  j = j + 1
 end do ! for while

 nbf = (j-1)*5
 BACKSPACE(fid)
 BACKSPACE(fid)
 read(fid,'(A)') buf
 close(fid)

 nbf = nbf + LEN_TRIM(buf(6:))/15
end subroutine read_cart_nbf_from_dat

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

 call open_file(xmlname, .true., fid)

 if(ab=='a' .or. ab=='A') then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(INDEX(buf,"type=""ALPH")/=0 .or. INDEX(buf,"type=""CANO")/=0 .or. &
      INDEX(buf,"type=""NATU")/=0) exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_mo_from_xml: none of 'ALPH', 'CAN&
                    &O' or 'NATU' is"
   write(6,'(A)') 'found in file '//TRIM(xmlname)
   close(fid)
   stop
  end if

 else ! UHF

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(INDEX(buf,"type=""BETA") /= 0) exit
  end do ! for while
  if(i /= 0) then
   write(6,'(A)') "ERROR in subroutine read_mo_from_xml: no type=""BETA&
                  & found in file "//TRIM(xmlname)//'.'
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
 integer :: i, fid1, fid2, RENAME
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

 call open_file(fchname, .true., fid1)
 i = INDEX(fchname,'.fch',back=.true.)
 fchname1 = fchname(1:i-1)//'_D.fch'
 open(newunit=fid2,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(buf(1:8) == key) exit
 end do

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine write_eigenvalues_to_fch: no '"//&
                  key//"' found in file "//TRIM(fchname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(5(1X,ES15.8))') (on(i),i=1,nif)

 ! skip the Alpha/Beta Orbital Energies in fname1
 do while(.true.)
  read(fid1,'(A)') buf
  if(INDEX(buf,'=') /= 0) exit
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
 fchname1 = TRIM(fchname)//'.t'

 call open_file(fchname, .true., fid1)
 open(newunit=fid2,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(buf(1:8) == key) exit
 end do

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

subroutine read_npair_from_uno_out(unofile,nbf,nif,ndb,npair,nopen,lin_dep)
 implicit none
 integer :: i, fid, idx(3), nvir
 integer, intent(out) :: nbf, nif, ndb, npair, nopen
 character(len=240), intent(in) :: unofile
 character(len=240) :: buf
 logical, intent(out) :: lin_dep

 buf = ' '; lin_dep = .false.
 open(newunit=fid,file=TRIM(unofile),status='old',position='rewind')

 read(fid,'(A)') buf
 i = INDEX(buf,'=')
 read(buf(i+1:),*) nbf

 read(fid,'(A)') buf
 i = INDEX(buf,'=')
 read(buf(i+1:),*) nif

 if(nbf > nif) then
  lin_dep = .true.
 else if(nbf < nif) then
  write(6,'(/,A)') 'ERROR in subroutine read_npair_from_uno_out: nbf<nif.'
  write(6,'(A)') 'This is impossible. Please check why.'
  stop
 end if

 close(fid)
 open(newunit=fid,file=TRIM(unofile),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3) == 'ndb') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_npair_from_uno_out: 'ndb' not found."
  write(6,'(A)') 'Please open file '//TRIM(unofile)//' and check.'
  close(fid)
  stop
 end if

 i = INDEX(buf,'=')
 read(buf(i+1:),*) ndb

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3) == 'idx') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_npair_from_uno_out: 'idx' not found."
  write(6,'(A)') 'Please open file '//TRIM(unofile)//' and check.'
  close(fid)
  stop
 end if

 idx = 0
 i = INDEX(buf,'=')
 read(buf(i+1:),*) idx(1:3)
 close(fid,status='delete')

 npair = (idx(2) - idx(1) - idx(3))/2
 nvir = nif - ndb - 2*npair - idx(3)
 nopen = idx(3)
 write(6,'(6(A,I0))') 'nbf=', nbf, ', nif=', nif, ', doubly_occ=', idx(1)-1, &
                      ', npair=', npair, ', nopen=', idx(3), ', nvir=', nvir
end subroutine read_npair_from_uno_out

! read GVB electronic energy from a given GAMESS .gms file
subroutine read_gvb_energy_from_gms(gmsname, e)
 implicit none
 integer :: i, j, fid
 real(kind=8), intent(out) :: e
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname

 e = 0d0
 call open_file(gmsname, .true., fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:10) == 'FINAL GVB') exit
 end do

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_gvb_energy_from_gms: no GVB energ&
                   &y found in'
  write(6,'(A)') 'file '//TRIM(gmsname)
  write(6,'(/,A)') 'You can open this file and check whether the SCF oscillates.'
  write(6,'(A)') 'If yes, reducing the number of processors and re-run may do&
                 & dome help.'
  write(6,'(A)') "If no, check if there is any error message like 'gamess.01.x&
                 & could not be found'."
  write(6,'(A)') 'In the latter case, you should read Section 4.4.10 in MOKIT&
                 & manual.'
  close(fid)
  stop
 end if
 close(fid)

 i = INDEX(buf,'IS'); j = INDEX(buf,'AFTER')
 read(buf(i+2:j-1),*) e

 if(DABS(e) < 1d-5) then
  write(6,'(/,A)') 'ERROR in subroutine read_gvb_energy_from_gms: it seems tha&
                   &t GVB computation does not'
  write(6,'(A)') 'converge. You can try to reduce the number of processors and&
                 & re-run.'
  stop
 end if
end subroutine read_gvb_energy_from_gms

! read CASCI/CASSCF energy from a Gaussian/PySCF/GAMESS/OpenMolcas/ORCA output file
subroutine read_cas_energy_from_output(cas_prog, outname, e, scf, spin, dmrg,&
                                       ptchg_e, nuc_pt_e)
 implicit none
 integer, intent(in) :: spin
 real(kind=8), intent(in) :: ptchg_e, nuc_pt_e
 real(kind=8), intent(out) :: e(2)
 character(len=10), intent(in) :: cas_prog
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf, dmrg

 select case(TRIM(cas_prog))
 case('gaussian')
  call read_cas_energy_from_gaulog(outname, e, scf)
 case('gamess')
  call read_cas_energy_from_gmsgms(outname, e, scf, spin)
  e(1) = e(1) + ptchg_e + nuc_pt_e
 case('pyscf')
  call read_cas_energy_from_pyout(outname, e, scf, spin, dmrg)
  e = e + ptchg_e
 case('openmolcas')
  call read_cas_energy_from_molcas_out(outname, e, scf)
  e = e + ptchg_e
 case('orca')
  call read_cas_energy_from_orca_out(outname, e, scf)
 case('molpro')
  call read_cas_energy_from_molpro_out(outname, e, scf)
  e = e + ptchg_e
 case('bdf')
  call read_cas_energy_from_bdf_out(outname, e, scf)
  e = e + ptchg_e + nuc_pt_e
 case('psi4')
  call read_cas_energy_from_psi4_out(outname, e, scf)
  e = e + ptchg_e + nuc_pt_e
 case('dalton')
  call read_cas_energy_from_dalton_out(outname, e, scf)
  e = e + ptchg_e
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_cas_energy_from_output: CAS_prog c&
                   &annot be identified.'
  write(6,'(A)') 'CAS_prog='//TRIM(cas_prog)
  stop
 end select
end subroutine read_cas_energy_from_output

! read CASCI/CASSCF energy from a Gaussian .log file
subroutine read_cas_energy_from_gaulog(outname, e, scf)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 call open_file(outname, .false., fid)
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(13:22)=='EIGENVALUE' .or. buf(20:29)=='EIGENVALUE' .or. &
     buf(23:32)=='Eigenvalue') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_cas_energy_from_gaulog: no&
                 & 'EIGENVALUE' or 'Eigenvalue'"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if
 close(fid)

 i = INDEX(buf,'lue')
 if(i == 0) i = INDEX(buf,'LUE')

 if(scf) then
  read(buf(i+3:),*) e(2) ! CASSCF
 else
  read(buf(i+3:),*) e(1) ! CASCI
 end if

 if(scf) then ! read CASCI energy
  call open_file(outname, .true., fid)
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:8) == 'ITN=  1') exit
  end do ! for while
  close(fid)
  i = INDEX(buf,'E=')
  read(buf(i+2:),*) e(1)
 end if
end subroutine read_cas_energy_from_gaulog

! read CASCI/CASSCF energy from a PySCF output file
subroutine read_cas_energy_from_pyout(outname, e, scf, spin, dmrg)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: spin ! na - nb
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8) :: s_square = 0d0, expect = 0d0
 real(kind=8), intent(out) :: e(2)
 real(kind=8), parameter :: max_diff = 1d-3
 logical, intent(in) :: scf, dmrg
 logical :: state_specific = .false.

 e = 0d0; i = 0; j = 0; k = 0
 expect = 0.5d0*DBLE(spin)
 expect = expect*(expect + 1d0)

 if(scf) then ! (DMRG-)CASSCF
  open(newunit=fid,file=TRIM(outname),status='old',position='append')
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
   if(buf(1:3) == 'SSS') state_specific = .true.
  end do ! for while

 else ! (DMRG-)CASCI

  if(dmrg) then
   open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
  else
   open(newunit=fid,file=TRIM(outname),status='old',position='append')
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
  end if
 end if

 if(k /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_cas_energy_from_out: problematic f&
                   &ile '//TRIM(outname)//'.'
  close(fid)
  stop
 else ! k = 0
  if(j/=0 .and. (.not.state_specific)) then
   write(6,'(/,A)') 'ERROR in subroutine read_cas_energy_from_out: CASCI or CAS&
                    &SCF not converged.'
   close(fid)
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
  write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_out: 'CASCI E' not&
                   & found in "//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 if(scf) then
  read(buf(i+1:),*) e(2)
 else
  read(buf(i+1:),*) e(1)
 end if

 i = INDEX(buf, '=', back=.true.)
 read(buf(i+1:),*) s_square

 if(DABS(expect - s_square) > max_diff) then
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Warning from subroutine read_cas_energy_from_pyout: <S**2> de&
                 &viates too much'
  write(6,'(2(A,F11.5))') 'from the expectation value. Expectation=', expect, &
                          ', S_square=', s_square
  write(6,'(A)') 'If this is a ground state calculation, it is probably because&
                 & this CASSCF is'
  write(6,'(A)') 'converged to a wrong spin state. You may try to add the keywo&
                 &rd CrazyWFN in'
  write(6,'(A)') 'mokit{} in .gjf file.'
  write(6,'(A)') 'If this is an excited state calculation where the spin of the&
                 & target excited'
  write(6,'(A)') 'state is different from that of the ground state, you can ign&
                 &ore this warning.'
  write(6,'(A)') REPEAT('-',79)
 end if

 ! Note: in a CASSCF job, there is also a CASCI energy, read it.
 if(scf) then
  rewind(fid)

  do while(.true.)
   read(fid,'(A)') buf
   if(buf(1:9) == 'CASCI E =') exit
  end do ! for while

  close(fid)
  read(buf(10:),*) e(1)
  i = INDEX(buf, '=', back=.true.)
  read(buf(i+1:),*) s_square

  if( DABS(expect - s_square) > max_diff) then
   write(6,'(/,A)') REPEAT('-',79)
   write(6,'(A)') 'Warning in subroutine read_cas_energy_from_pyout: the 0-th s&
                  &tep in this CASSCF'
   write(6,'(A)') 'job, i.e. the CASCI <S**2> deviates too much from the expect&
                  &ation value.'
   write(6,'(2(A,F11.5))') 'Expectation=', expect, ', S_square=', s_square
   write(6,'(A)') 'If this is a ground state calculation, it is probably becaus&
                  &e this CASCI is'
   write(6,'(A)') 'unconverged, or converged to a wrong spin state. If this CAS&
                  &CI energy is'
   write(6,'(A)') 'useless to you, or if the following CASSCF happens to be con&
                  &verged to the'
   write(6,'(A)') 'desired spin, you can ignore this warning. Otherwise, you ma&
                  &y try to add'
   write(6,'(A)') 'the keyword CrazyWFN in mokit{} in .gjf file.'
   write(6,'(A)') 'If this is an excited state calculation where the spin of th&
                  &e target excited'
   write(6,'(A)') 'state is different from that of the ground state, you can ig&
                  &nore this warning.'
   write(6,'(A)') REPEAT('-',79)
  end if
 else
  close(fid)
 end if
end subroutine read_cas_energy_from_pyout

! read CASCI/CASSCF energy from the GAMESS output file
subroutine read_cas_energy_from_gmsgms(outname, e, scf, spin)
 implicit none
 integer :: i, fid
 integer, intent(in) :: spin
 real(kind=8) :: s_square, expect
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 expect = DBLE(spin)/2d0
 expect = expect*(expect + 1d0)
 call open_file(outname, .true., fid)

 if(scf) then  ! CASSCF job
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:24) == 'THE DENSITIES ARE STATE') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(A)') "ERROR in subroutine read_cas_energy_from_gmsgms: no&
                  & 'THE DENSITIES ARE STATE' found in file "//TRIM(outname)
   stop
  end if

  read(fid,'(A)') buf
  i = INDEX(buf,'ENERGY=')
  read(buf(i+7:),*) e(1)   ! CASCI energy in the CASSCF job
  i = INDEX(buf,'=',back=.true.)
  read(buf(i+1:),*) s_square
  s_square = s_square*(s_square+1d0)
  if( DABS(expect - s_square) > 1D-2) then
   write(6,'(A)') 'ERROR in subroutine read_cas_energy_from_gmsgms: in this&
                  & CASSCF job, the 0-th step, i.e., the CASCI'
   write(6,'(A)') '<S**2> deviates too much from the expectation value.'
   write(6,'(2(A,F10.6))') 'expectation = ', expect, ', s_square=', s_square
   stop
  end if

  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:13) == 'STATE   1  E') exit
  end do ! for while
  i = INDEX(buf,'ENERGY=')
  read(buf(i+7:),*) e(2)   ! CASSCF energy
  i = INDEX(buf,'S=')
  read(buf(i+2:),*) s_square
  s_square = s_square*(s_square+1d0)
  if( DABS(expect - s_square) > 1D-2) then
   write(6,'(A)') 'ERROR in subroutine read_cas_energy_from_gmsgms: CASSCF&
                  & <S**2> deviates too much from the expectation value.'
   write(6,'(2(A,F10.6))') 'expectation = ', expect, ', s_square=', s_square
   stop
  end if

 else          ! CASCI job
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:20) == 'DENSITY MATRIX WILL') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(A)') "ERROR in subroutine read_cas_energy_from_gmsgms: no&
                  & 'DENSITY MATRIX' found in file "//TRIM(outname)
   stop
  end if
 
  i = INDEX(buf,'=', back=.true.)
  read(buf(i+1:),*) s_square
  s_square = s_square*(s_square+1d0)
  if( DABS(expect - s_square) > 1D-2) then
   write(6,'(A)') 'ERROR in subroutine read_cas_energy_from_gmsgms: CASCI&
                  & <S**2> deviates too much from the expectation value.'
   write(6,'(2(A,F10.6))') 'expectation = ', expect, ', s_square=', s_square
   stop
  end if

  read(fid,'(A)') buf
  read(fid,'(A)') buf
  i = INDEX(buf,'=', back=.true.)
  read(buf(i+1:),*) e(1)
 end if

 close(fid)
end subroutine read_cas_energy_from_gmsgms

! read CASCI/CASSCF energy from a given OpenMolcas/Molcas output file
subroutine read_cas_energy_from_molcas_out(outname, e, scf)
 implicit none
 integer :: i, j, fid
 real(kind=8) :: add
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 character(len=53), parameter :: error_warn = 'ERROR in subroutine read_cas_en&
                                              &ergy_from_molcas_out: '
 logical, intent(in) :: scf

 call open_file(outname, .true., fid)
 e = 0d0; add = 0d0

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:21) == 'Nr of preliminary CI') exit
 end do ! for while

 read(fid,'(A)') buf
 if(INDEX(buf,'No convergence') > 0) then
  if(scf) then
   write(6,'(/,A)') 'Warning in subroutine read_cas_energy_from_molcas_out:'
   write(6,'(A)') 'The CASCI iterative diagonalization fails to converge. This &
                  &is a defect'
   write(6,'(A)') 'of OpenMolcas when doing CASSCF. If you want a correct CASCI&
                  & energy, please'
   write(6,'(A)') 'run a single CASCI job. This may or may not affect the final&
                  & CASSCF result,'
   write(6,'(A)') 'so the program will continue.'
  else ! CASCI
   write(6,'(/,A)') error_warn
   write(6,'(A)') 'The CASCI iterative diagonalization fails to converge.'
   close(fid)
   stop
  end if
  read(fid,'(A)') buf
 end if

 if(INDEX(buf,'Total energies') > 0) then
  i = INDEX(buf,'Add'); j = INDEX(buf,'au')
  read(buf(i+3:j-1),*) add
  read(fid,'(A)') buf
 end if

 read(buf,*) i, j, i, j, e(1) ! CASCI energy
 e(1) = e(1) + add
 if(.not. scf) then ! CASCI
  close(fid)
  return
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'Convergence after') > 0) exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') error_warn
  write(6,'(A)') "'Convergence after' is not found in file "//TRIM(outname)
  stop
 end if

 read(fid,*) i, j, i, j, e(2) ! CASSCF energy
 e(2) = e(2) + add
 close(fid)
end subroutine read_cas_energy_from_molcas_out

! read CASCI/CASSCF energy from a specified ORCA output file
subroutine read_cas_energy_from_orca_out(outname, e, scf)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(2)
 character(len=51), parameter :: error_warn = &
  'ERROR in subroutine read_cas_energy_from_orca_out: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 e = 0d0

 if(scf) then ! CASSCF
  open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:9) == '   E(CAS)') exit
  end do ! for while

  if(i /= 0) then
   write(6,'(/,A)') error_warn//"'   E(CAS)' not found in"
   write(6,'(A)') 'file '//TRIM(outname)
   close(fid)
   stop
  end if
  i = INDEX(buf, '=')
  read(buf(i+1:),*) e(1)

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:19) == 'Final CASSCF energy') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(/,A)') error_warn//"'Final CASSCF energy' not found in"
   write(6,'(A)') 'file '//TRIM(outname)
   close(fid)
   stop
  end if
  i = INDEX(buf, ':')
  read(buf(i+1:),*) e(2)

 else ! CASCI
  open(newunit=fid,file=TRIM(outname),status='old',position='append')

  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:5) == 'STATE') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(/,A)') error_warn//"'^STATE' not found in file "//TRIM(outname)
   close(fid)
   stop
  end if

  i = INDEX(buf, '=')
  read(buf(i+1:),*) e(1)
 end if

 close(fid)
end subroutine read_cas_energy_from_orca_out

! read CASCI/CASSCF energy from a given Molpro output file
subroutine read_cas_energy_from_molpro_out(outname, e, scf)
 implicit none
 integer :: i, j, k, fid
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  select case(buf(2:8))
  case('ITER. M','ITER  M','ITE  MI')
   exit
  end select
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_molpro_out: 'ITER.&
                   & M' not found in"
  write(6,'(A)') 'file '//TRIM(outname)
  write(6,'(A)') 'Error termination of the Molpro CASCI job.'
  close(fid)
  stop
 end if

 do i = 1, 4
  read(fid,fmt=*,iostat=j) k,k,k,k, e(1) ! CASCI energy
  if(j == 0) exit
 end do ! for i
 !read(fid,'(A)') buf
 !read(fid,'(A)') buf
 !if(LEN_TRIM(buf) == 0) then ! Multipassing in transformation
 ! read(fid,'(A)') buf
 ! read(fid,'(A)') buf
 !end if
 !read(buf,*) k,k,k,k, e(1) ! CASCI energy

 if(scf) then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:9) == '!MCSCF S') exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_molpro_out: '!MCS&
                    &CF S' not found in"
   write(6,'(A)') 'file '//TRIM(outname)
   write(6,'(A)') 'Error termination of the Molpro CASSCF job.'
   close(fid)
   stop
  end if
  i = INDEX(buf, 'Energy')
  read(buf(i+6:),*) e(2) ! CASSCF energy
 end if

 close(fid)
end subroutine read_cas_energy_from_molpro_out

! read CASCI/CASSCF energy from a given BDF output file
! Note: BDF changes output format frequently, one must frequently update this subroutine
subroutine read_cas_energy_from_bdf_out(outname, e, scf)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 e = 0d0
 call open_file(outname, .true., fid)

 if(scf) then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:11) == 'mcscf_eneci') exit
  end do ! for while

  if(i /= 0) then
   write(6,'(A)') "ERROR in subroutine read_cas_energy_from_bdf_out:&
                  & 'mcscf_eneci' not found in file "//TRIM(outname)
   write(6,'(A)') 'Error termination of the BDF CASSCF job.'
   close(fid)
   stop
  end if
  read(buf(15:),*) e(1) ! CASCI energy in CASSCF job
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(3:19) == 'CHECKDATA:MCSCF:M') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_cas_energy_from_bdf_out:&
                 & 'CHECKDATA:MCSCF:M' not found in file "//TRIM(outname)
  write(6,'(A)') 'Error termination of the BDF CASCI/CASSCF job.'
  close(fid)
  stop
 end if

 i = INDEX(buf,':', back=.true.)
 if(scf) then
  read(buf(i+1:),*) e(2) ! CASSCF energy in CASSCF job
 else
  read(buf(i+1:),*) e(1) ! CASCI energy in CASCI job
 end if

 close(fid)
end subroutine read_cas_energy_from_bdf_out

! read CASCI/CASSCF energy from a given PSI4 output file
subroutine read_cas_energy_from_psi4_out(outname, e, scf)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 e = 0d0
 call open_file(outname, .true., fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(scf) then
   if(buf(1:12)=='        Iter' .or. buf(1:15)=='           Iter') exit
  else
   if(buf(5:19) == 'Total CI energy') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_cas_energy_from_psi4_out: no '&
                 &Iter' found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 if(.not. scf) then ! CASCI
  close(fid)
  i = INDEX(buf,'=')
  read(buf(i+1:),*) e(1)
  return
 end if

 read(fid,'(A)') buf
 ! sometimes there will be extra output in PSI4, e.g.
 ! '(sem_iter): H0block_->H0b_diag'...
 if(INDEX(buf,'MCSCF  1:') == 0) then
  do while(.true.)
   read(fid,'(A)') buf
   if(INDEX(buf,'MCSCF  1:') > 0) exit
  end do ! for while
 end if

 i = INDEX(buf,':')
 read(buf(i+1:),*) e(1)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'MCSCF Final E') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_cas_energy_from_psi4_out: no '&
                 &MCSCF Final E' found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf,':')
 read(buf(i+1:),*) e(2)
 close(fid)
end subroutine read_cas_energy_from_psi4_out

! read CASCI/CASSCF energy from a given Dalton output file
subroutine read_cas_energy_from_dalton_out(outname, e, scf)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(2)
 character(len=1) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 e = 0d0
 call open_file(outname, .true., fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:19) == '@ Final CI energies') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_cas_energy_from_dalton_out: no '&
                 &@ Final CI energies' found in"
  write(6,'(A)') 'file '//TRIM(outname)//'.'
  close(fid)
  stop
 end if

 read(fid,*) str, i, e(1) ! CASCI energy
 if(.not. scf) then
  close(fid)
  return
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:19) == '@    Final MCSCF en') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_cas_energy_from_dalton_out: no '&
                 &@    Final MCSCF en' found in"
  write(6,'(A)') 'file '//TRIM(outname)//'.'
  close(fid)
  stop
 end if

 i = INDEX(buf,':')
 read(buf(i+1:),*) e(2)
 close(fid)
end subroutine read_cas_energy_from_dalton_out

! read NEVPT2 energy from PySCF output file
subroutine read_mrpt_energy_from_pyscf_out(outname, troot, ref_e, corr_e)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: troot ! 0 for the ground state, >0 for excited state
 character(len=7), parameter :: str1 = 'CASCI E'
 character(len=15) :: str2 = ' '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: ref_e, corr_e

 ref_e = 0d0; corr_e = 0d0
 call open_file(outname, .false., fid)

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:13) == 'Nevpt2 Energy') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine read_mrpt_energy_from_pyscf_out:'
  write(6,'(A)') 'No NEVPT2 energy found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) corr_e

 if(troot == 0) then
  str2 = str1
  k = 7
 else
  write(str2,'(A,1X,I0)') 'CASCI state', troot
  k = LEN_TRIM(str2)
 end if

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:k) == TRIM(str2)) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_pyscf_out:'
  write(6,'(A)') 'No CASCI energy found in file '//TRIM(outname)
  close(fid)
  stop
 end if
 close(fid)

 i = INDEX(buf, '=')
 read(buf(i+1:),*) ref_e
end subroutine read_mrpt_energy_from_pyscf_out

! read CASTP2 energy from OpenMolcas output file
subroutine read_mrpt_energy_from_molcas_out(outname, itype, ref_e, corr_e)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: itype ! 1/2/3 for SC-NEVPT2/FIC-NEVPT2/CASPT2
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: ref_e, corr_e

 ref_e = 0d0
 corr_e = 0d0
 if(itype<1 .or. itype>3) then
  write(6,'(A,I0)') 'ERROR in subroutine read_mrpt_energy_from_molcas_out:&
                    & invalid itype=', itype
  write(6,'(A)') 'Allowed values are 1/2/3 for SC-NEVPT2/FIC-NEVPT2/CASPT2.'
  stop
 end if

 call open_file(outname, .true., fid)
 select case(itype)
 case(1,2) ! NEVPT2
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:30) == 'Energies of zeroth-order DMRG') exit
  end do ! for while

  if(i /= 0) then
   write(6,'(A)') 'ERROR in subroutine read_mrpt_energy_from_molcas_out:'
   write(6,'(A)') 'No reference energy found in file '//TRIM(outname)
   close(fid)
   stop
  end if

  read(fid,'(A)') buf
  read(fid,'(A)') buf
  i = INDEX(buf, '=')
  read(buf(i+1:),*) ref_e

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:18) == 'state number:   1') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(A)') 'ERROR in subroutine read_mrpt_energy_from_molcas_out:'
   write(6,'(A)') "No 'state number:   1' found in file "//TRIM(outname)
   close(fid)
   stop
  end if

  do j = 1, 15
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:6) == 'Total') exit
  end do ! for j

  if(i/=0 .or. j==16) then
   write(6,'(A)') 'ERROR in subroutine read_mrpt_energy_from_molcas_out:'
   write(6,'(A)') 'No NEVPT2 energy found in file '//TRIM(outname)
   close(fid)
   stop
  end if

  if(itype == 1) then ! SC-NEVPT2
   read(buf(7:),*) rtmp, corr_e
  else                ! FIC-NEVPT2
   read(buf(7:),*) rtmp, rtmp, rtmp, corr_e
  end if

 case(3) ! CASPT2
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(7:23) == 'Reference energy:') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(A)') 'ERROR in subroutine read_mrpt_energy_from_molcas_out:'
   write(6,'(A)') 'No reference energy found in file '//TRIM(outname)
   close(fid)
   stop
  end if

  i = INDEX(buf,':')
  read(buf(i+1:),*) ref_e

  do i = 1, 3
   read(fid,'(A)') buf
  end do ! for i 

  i = INDEX(buf, ':')
  read(buf(i+1:),*) corr_e
  corr_e = corr_e - ref_e

 case default
  write(6,'(A)') 'ERROR in subroutine read_mrpt_energy_from_molcas_out:'
  write(6,'(A,I0)') 'Invalid itype=', itype
  stop
 end select

 close(fid)
end subroutine read_mrpt_energy_from_molcas_out

! read NEVPT2/CASPT2 energy from a Molpro output file
subroutine read_mrpt_energy_from_molpro_out(outname, itype, ref_e, corr_e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: itype ! 1/2/3 for SC-NEVPT2/FIC-NEVPT2/CASPT2
 real(kind=8), intent(out) :: ref_e, corr_e
 character(len=8), parameter :: key(3) = ['Strongly','!NEVPT2 ','!RSPT2 S']
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; corr_e = 0d0
 if(itype<1 .or. itype>3) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_molpro_out: itype&
                   & out of range.'
  write(6,'(A,I0,A)') 'itype=', itype, ', outname='//TRIM(outname)
  stop
 end if

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:20)=='!MCSCF STATE  1.1 E' .or. buf(2:19)=='!MCSCF STATE 1.1 E') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_molpro_out: CASSC&
                   &F energy cannot'
  write(6,'(A)') "be found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf,'ergy')
 read(buf(i+4:),*) ref_e

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:9) == key(itype)) exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_molpro_out:'
  write(6,'(A)') "'"//key(itype)//"' not found in file "//TRIM(outname)
  stop
 end if

 i = INDEX(buf,'ergy')
 read(buf(i+4:),*) corr_e
 corr_e = corr_e - ref_e
end subroutine read_mrpt_energy_from_molpro_out

! read NEVPT2/CASPT2 energy from a ORCA .out file
subroutine read_mrpt_energy_from_orca_out(outname, itype, ref_e, corr_e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: itype ! 1/2/3 for SC-NEVPT2/FIC-NEVPT2/CASPT2
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: ref_e, corr_e

 ref_e = 0d0
 corr_e = 0d0
 call open_file(outname, .true., fid)

 select case(itype)
 case(1,2,3)
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(INDEX(buf,'Total Energy C') /=0) exit
  end do ! for while

  if(i /= 0) then
   write(6,'(A)') "ERROR in subroutine read_mrpt_energy_from_orca_out: no&
                  & 'Total Energy (' found in file "//TRIM(outname)
   close(fid)
   stop
  end if

  i = INDEX(buf, '=')
  read(buf(i+1:),*) corr_e
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  i = INDEX(buf, '=')
  read(buf(i+1:),*) ref_e
 case default
  write(6,'(A,I0)') 'ERROR in subroutine read_mrpt_energy_from_orca_out:&
                    & invalid itype=', itype
  stop
 end select

 close(fid)
end subroutine read_mrpt_energy_from_orca_out

! read MRMP2 energy from GAMESS output file (.gms)
subroutine read_mrpt_energy_from_gms_out(outname, ref_e, corr_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, corr_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0
 corr_e = 0d0
 call open_file(outname, .true., fid)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:16) == 'TOTAL   (MCSCF)') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine read_mrpt_energy_from_gms_out:'
  write(6,'(A)') 'No reference energy found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) ref_e

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(22:42) == '2ND ORDER ENERGY CORR') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine read_mrpt_energy_from_gms_out:'
  write(6,'(A)') 'No MRMP2 energy found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) corr_e
 close(fid)
end subroutine read_mrpt_energy_from_gms_out

! read NEVPT2 energy from a BDF .out file
subroutine read_mrpt_energy_from_bdf_out(outname, itype, ref_e, corr_e, dav_e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: itype   ! 1/2 for SDSPT2/FIC-NEVPT2
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8) :: rtmp(6)
 real(kind=8), intent(out) :: ref_e, corr_e, dav_e
 ! dav_e: Davidson correction energy for size-inconsistency error

 ref_e = 0d0
 corr_e = 0d0
 dav_e = 0d0
 if(.not. (itype==1 .or. itype==2)) then
  write(6,'(A)') 'ERROR in subroutine read_mrpt_energy_from_bdf_out: currently&
                 & only reading SDSPT2/NEVPT2 energy is supported.'
  stop
 end if

 call open_file(outname, .true., fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == 'Print final') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_mrpt_energy_from_bdf_out:&
                 & 'Print final' not found in file "//TRIM(outname)
  write(6,'(A)') 'Error termination of the BDF CASSCF in SDSPT2/NEVPT2 job.'
  close(fid)
  stop
 end if

 do i = 1, 4
  read(fid,'(A)') buf
 end do ! for i
 i = INDEX(buf, '=')
 read(buf(i+1:),*) ref_e ! CASCI/CASSCF energy

 if(itype == 1) then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:11) == 'MRPT2 calc') exit
  end do ! for while
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  read(fid,*) i, rtmp(1:6), dav_e
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:13) == 'TARGET_XIANCI') exit
 end do ! for while
 close(fid)

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_mrpt_energy_from_bdf_out: no&
                 & 'TARGET_XIANCI' found in file "//TRIM(outname)
  stop
 end if

 read(buf(26:),*) corr_e
 if(itype == 1) dav_e = dav_e - corr_e
 corr_e = corr_e - ref_e
end subroutine read_mrpt_energy_from_bdf_out

! whether the Davidson Q CORRECTION can be found in a GAMESS output file
function has_davidson_q(outname) result(alive)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: alive

 alive = .false.
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:18) == 'CALC. OF DAVIDSON') exit
 end do ! for while

 close(fid)
 if(i /= 0) return

 i = INDEX(buf, '=')
 read(buf(i+1:),*) alive
end function has_davidson_q

! read Davidson correction and MRCISD energy from OpenMolcas, ORCA, Gaussian or
! Molpro output file
subroutine read_mrci_energy_from_output(CtrType, mrcisd_prog, outname, ptchg_e,&
                                        nuc_pt_e, davidson_e, e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: CtrType
 real(kind=8) :: casci_e, ref_weight
 real(kind=8), intent(in) :: ptchg_e, nuc_pt_e
 real(kind=8), intent(out) :: davidson_e, e
 character(len=10), intent(in) :: mrcisd_prog
 character(len=240), intent(in) :: outname
 character(len=240) :: buf
 logical, external :: has_davidson_q

 davidson_e = 0d0; e = 0d0; ref_weight = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 select case(TRIM(mrcisd_prog))
 case('pyscf')
  if(CtrType == 3) then ! DMRG-FIC-MRCISD
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(1:7) == 'E(MRCI)') exit
   end do ! for while
   i = INDEX(buf, 'DC =', back=.true.)
   read(buf(i+4:),*) davidson_e
   read(fid,'(A)') buf
   i = INDEX(buf, '=')
   read(buf(i+1:),*) e
  else
   write(6,'(/,A)') 'ERROR in subroutine read_mrci_energy_from_output: invalid &
                    &CtrType.'
   write(6,'(A,I0)') 'Current CtrType=', CtrType
   close(fid)
   stop
  end if
  e = e + ptchg_e
 case('openmolcas')
  if(CtrType == 1) then ! uncontracted MRCISD
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(1:29) == '::    RASSCF root number  1 T') exit
   end do ! for while
   i = INDEX(buf,':',back=.true.)
   read(buf(i+1:),*) e
  else if(CtrType == 2) then ! ic-MRCISD
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(23:31) == 'CI ENERGY') exit
   end do ! for while
   i = INDEX(buf,':',back=.true.)
   read(buf(i+1:),*) e
   read(fid,'(A)') buf
   i = INDEX(buf,':',back=.true.)
   read(buf(i+1:),*)  davidson_e
  end if
  e = e + ptchg_e

 case('orca')
  if(CtrType == 1) then ! uncontracted MRCISD
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(20:27) == 'E(MR-CI)') exit
   end do ! for while
   read(fid,'(A)') buf
   read(fid,'(A)') buf
   read(buf(13:),*) e
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(14:26) == 'residual conv') exit
   end do ! for while
   read(fid,'(A)') buf
   read(buf(67:),*) ref_weight
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(1:17) == 'Computing the ref') exit
   end do ! for while
   BACKSPACE(fid)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   read(buf(22:),*) casci_e
   davidson_e = (1d0 - ref_weight)*(e - casci_e)
  else if(CtrType == 3) then ! FIC-MRCISD
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(1:5) == 'ETOT ') exit
   end do ! for while
   i = INDEX(buf,'...',back=.true.)
   read(buf(i+3:),*) e
   read(fid,'(A)') buf
   read(fid,'(A)') buf
   i = INDEX(buf,'...',back=.true.)
   read(buf(i+3:),*) davidson_e
  end if

 case('gaussian') ! uncontracted MRCISD
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   !               Davidson                           Lanczos
   if(buf(17:32)=='Final Eigenvalue' .or. buf(20:29)=='EIGENVALUE') exit
  end do ! for while
  i = INDEX(buf,'lue',back=.true.)
  if(i == 0) i = INDEX(buf,'LUE',back=.true.)
  read(buf(i+3:),*) e

 case('molpro')
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(INDEX(buf(2:10),'!Total e') > 0) exit
  end do ! for while

  read(buf(17:),*) e
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) read(fid,'(A)') buf
  i = INDEX(buf,':')
  read(buf(i+1:),*) davidson_e
  davidson_e = davidson_e - e
  e = e + ptchg_e

 case('psi4')
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(5:19) == 'Total CI energy') exit

   if(buf(11:23) == 'Psi4: An Open') then
    write(6,'(/,A)') 'ERROR in subroutine read_mrci_energy_from_output:'
    write(6,'(A)') "No 'Total CI energy' found in file "//TRIM(outname)
    close(fid)
    stop
   end if
  end do ! for while

  i = INDEX(buf,'=')
  read(buf(i+1:),*) e
  e = e + ptchg_e + nuc_pt_e

 case('dalton')
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(1:19) == '@ Final CI energies') exit

   if(buf(22:32) == 'Dalton - An') then
    write(6,'(/,A)') 'ERROR in subroutine read_mrci_energy_from_output:'
    write(6,'(A)') "No '@ Final CI energies' found in file "//TRIM(outname)
    close(fid)
    stop
   end if
  end do ! for while

  read(fid,'(A)') buf
  read(buf(7:),*) e
  e = e + ptchg_e

 case('gamess')
  close(fid)

  if(has_davidson_q(outname)) then
   call open_file(outname, .false., fid)
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(1:13) == ' E(MR-CISD) =') exit
   end do ! for while

   read(buf(14:),*) e
   read(fid,'(A)') buf
   read(fid,'(A)') buf
   read(buf(14:),*) ref_weight

   do while(.true.)
    read(fid,'(A)') buf
    if(buf(2:7) == 'E(REF)') exit
   end do ! for while
   i = INDEX(buf, '=')
   read(buf(i+1:),*) casci_e

   ! GAMESS uses renormalized Davidson correction. Here we compute the Davidson
   ! correction
   davidson_e = (1d0 - ref_weight)*(e - casci_e)

  else ! no 'CALC. OF DAVIDSON', i.e. no Davidson Q

   call open_file(outname, .true., fid)
   do while(.true.)
    read(fid,'(A)')  buf
    if(buf(2:14) == 'CI EIGENSTATE') exit
   end do ! for while

   close(fid)
   i = INDEX(buf, '=')
   read(buf(i+1:),*) e
  end if

 case default
  write(6,'(/,A)') 'ERROR in subroutine read_mrci_energy_from_output: invalid m&
                   &rcisd_prog='//TRIM(mrcisd_prog)
  close(fid)
  stop
 end select

 close(fid)
end subroutine read_mrci_energy_from_output

! find npair0: the number of active pairs (|C2| > 0.1)
! (assuming that the pair coefficients haven been sorted)
subroutine find_npair0_from_dat(datname, npair, npair0)
 implicit none
 integer :: i, k, datid
 integer, intent(in) :: npair
!f2py intent(in) :: npair
 integer, intent(out) :: npair0
!f2py intent(out) :: npair0
 real(kind=8) :: rtmp
 real(kind=8), allocatable :: pair_coeff(:,:)
 character(len=240) :: buf
 character(len=240), intent(in) :: datname
!f2py intent(in) :: datname

 if(npair == 0) then
  npair0 = 0
  write(6,'(/,A)') 'Warning in subroutine find_npair0_from_dat: npair=npair0=0.'
  return
 end if

 call open_file(datname, .true., datid)
 ! find pair coefficients
 do while(.true.)
  read(datid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'CICOEF(') /= 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine find_npair0_from_dat: 'CICOEF(' keyword&
                   & not found in"
  write(6,'(A)') 'file '//TRIM(datname)
  close(datid)
  stop
 end if

 ! read pair coefficients
 BACKSPACE(datid)
 allocate(pair_coeff(2,npair), source=0d0)
 do i = 1, npair, 1
  read(datid,'(A)') buf
  k = INDEX(buf,'=')
  read(buf(k+1:),*) pair_coeff(1,i)
  k = INDEX(buf,',')
  read(buf(k+1:),*) pair_coeff(2,i)
 end do ! for i
 ! pair coefficients read done

 close(datid)
 npair0 = COUNT(pair_coeff(2,:) <= -0.1d0)

 do i = 1, npair, 1
  rtmp = pair_coeff(2,i) + 0.1d0
  if(rtmp>0d0 .and. rtmp<2.6d-3) then
   write(6,'(/,A)') REPEAT('-',79)
   write(6,'(A)') 'Warning from subroutine find_npair0_from_dat: some occupatio&
                  &n number is very'
   write(6,'(A)') 'close to ON_thres. You may consider enlarge the active space&
                  & size slightly in'
   write(6,'(A,I0)') 'later CAS calculation. Check No. pair: ', i
   write(6,'(A)') REPEAT('-',79)
  end if
 end do ! for i

 deallocate(pair_coeff)
end subroutine find_npair0_from_dat

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

! check whether pure Cartesian functions
subroutine check_cart_in_fch(fchname, cart)
 implicit none
 integer :: i, k, fid
 integer, allocatable :: shltyp(:)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 character(len=38), parameter :: error_warn='ERROR in subroutine check_cart_in_fch:'
 logical, intent(in) :: cart

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == 'Shell types') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//" missing 'Shell types'"
  write(6,'(A)') 'in file '//TRIM(fchname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=', back=.true.)
 read(buf(i+1:),*) k
 allocate(shltyp(k), source=0)
 read(fid,'(6(6X,I6))') (shltyp(i),i=1,k)
 close(fid)

 if(ANY(shltyp<-1) .and. ANY(shltyp>1)) then
  write(6,'(/,A)') error_warn//' mixed spherical harmonic/Cartesian functions detected.'
  write(6,'(A)') 'You probably used the 6-31G(d) basis set in Gaussian. Its&
                 & default setting is (6D,7F).'
  write(6,'(A)') 'AutoMR can deal only pure spherical harmonic or pure Cartesian&
                 & functions.'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 if(ANY(shltyp<-1) .and. cart) then
  write(6,'(/,A)') error_warn//' Cartesian functions required. But you provided'
  write(6,'(A)') 'a .fch file which uses spherical harmonic functions. Two&
                 & possible solutions:'
  write(6,'(A)') "1) delete keyword 'Cart' in MOKIT{} ; 2) provide another .fch&
                 & file which uses"
  write(6,'(A)') 'pure Cartesian functions. fchname='//TRIM(fchname)
  stop
 end if

 if(ANY(shltyp>1) .and. (.not.cart)) then
  write(6,'(/,A)') error_warn//' spherical harmonic functions default. But you'
  write(6,'(A)') 'provided a .fch file which has Cartesian functions. Two&
                 & possible solutions:'
  write(6,'(A)') "1) add keyword 'Cart' in MOKIT{}; 2) provide another .fch fil&
                 &e which uses pure"
  write(6,'(A)') 'spherical harmonic functions. fchname='//TRIM(fchname)
  stop
 end if

 deallocate(shltyp)
end subroutine check_cart_in_fch

! return whether pure Cartesian or spherical harmonic
subroutine check_sph_in_fch(fchname, sph)
 implicit none
 integer :: i, k, fid
 integer, allocatable :: shltyp(:)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 character(len=37), parameter :: error_warn='ERROR in subroutine check_sph_in_fch:'
 logical, intent(out) :: sph

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == 'Shell types') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//" missing 'Shell types'"
  write(6,'(A)') 'in file '//TRIM(fchname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=', back=.true.)
 read(buf(i+1:),*) k
 allocate(shltyp(k), source=0)
 read(fid,'(6(6X,I6))') (shltyp(i),i=1,k)
 ! read Shell types done
 close(fid)

 if(ANY(shltyp<-1) .and. ANY(shltyp>1)) then
  write(6,'(/,A)') error_warn//' mixed spherical harmonic/Cartesian functions &
                  &detected.'
  write(6,'(A)') 'You probably used the 6-31G(d) basis set in Gaussian. Its&
                 & default setting is (6D,7F).'
  write(6,'(A)') 'Only pure Cartesian or spherical harmonic is allowed'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 else if(ANY(shltyp>1)) then
  sph = .false.
 else
  sph = .true.
 end if
end subroutine check_sph_in_fch

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

 ! symmetrize the density matrix
 call symmetrize_dmat(nbf, dm)
end subroutine read_dm_from_fch

! Write 'Total SCF Density' or 'Spin SCF Density' into a .fch(k) file.
! Note: (1) the dm(j,i) j<=i will be used; (2) the density matrix elements in
!  the array dm must be in Gaussian angular momentum order.
subroutine write_dm_into_fch(fchname, total, nbf, dm)
 implicit none
 integer :: i, j, k, fid, fid1, RENAME
 integer, intent(in) :: nbf
 real(kind=8), intent(in) :: dm(nbf,nbf)
 character(len=11) :: key
 character(len=11), parameter :: key1 = 'Total SCF D'
 character(len=11), parameter :: key2 = 'Spin SCF De'
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: total

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

 call open_file(fchname1, .true., fid)
 uhf = .false.
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
 call open_file(file47, .true., fid)

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
                   &l utility solve_ON_matrix.'
  write(6,'(A)') 'Did you forget to compile this utility?'
  write(6,'(A)') 'mo_fch='//TRIM(mo_fch)
  write(6,'(A)') 'no_fch='//TRIM(no_fch)
  stop
 end if

 call read_na_and_nb_from_fch(mo_fch, na, nb)
 nopen = na - nb; ndb = nb

 i = INDEX(mo_fch, '.fch', back=.true.)
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
 real(kind=8) :: r1, r2
 real(kind=8), allocatable :: noon(:), coeff(:,:), rtmp(:)
 character(len=240), intent(in) :: fchname

 if(i1 == i2) return
 if(i1<1 .or. i2<1 .or. i1>i2) then
  write(6,'(/,A)') 'ERROR in subroutine sort_no_by_noon: wrong orbital indices!'
  write(6,'(2(A,I0))') 'i1=', i1, ', i2=', i2
  stop
 end if
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(noon(nif))
 call read_eigenvalues_from_fch(fchname, nif, 'a', noon)
 allocate(coeff(nbf,nif), rtmp(nbf))
 call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)

 do i = i1, i2-1, 1
  r1 = noon(i); rtmp = coeff(:,i)
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
 integer :: nif, na, nb
 integer, intent(out) :: idx
!f2py intent(out) :: idx
 real(kind=8), parameter :: thres = 0.3d0
 real(kind=8), allocatable :: ev(:)
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call read_nif_from_fch(fchname, nif)
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

! Convert an R(O)HF-type .fch(k) file into a UHF-type one. If
! brokensym=.True., beta HOMO-LUMO will be interchanged.
subroutine fch_r2u(fchname, brokensym)
 implicit none
 integer :: i, j, k, nbf, nif, ncoeff, na, nb, nline, fid, fid1
 real(kind=8), allocatable :: e_a(:), rtmp(:), mo_a(:,:)
 character(len=240) :: buf, uhf_fch
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical, intent(in) :: brokensym
!f2py intent(in) :: brokensym

 call find_specified_suffix(fchname, '.fch', i)
 uhf_fch = fchname(1:i-1)//'_u.fch'

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(e_a(nif), mo_a(nbf,nif))
 call read_eigenvalues_from_fch(fchname, nif, 'a', e_a)
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo_a)
 call read_na_and_nb_from_fch(fchname, na, nb)

 if(brokensym) then
  allocate(rtmp(nbf), source=mo_a(:,nb))
  mo_a(:,nb) = mo_a(:,nb+1)
  mo_a(:,nb+1) = rtmp
  rtmp(1) = e_a(nb)
  e_a(nb) = e_a(nb+1)
  e_a(nb+1) = rtmp(1)
  deallocate(rtmp)
 end if

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
  write(6,'(/,A)') 'ERROR in subroutine fch_r2u: invalid case.'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  close(fid)
  close(fid1)
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
 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(buf(1:7) == 'Alpha M') exit
 end do ! for while

 ncoeff = nbf*nif
 nline = ncoeff/5
 if(ncoeff > nline*5) nline = nline + 1

 do i = 1, nline, 1
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
 end do ! for i

 write(fid1,'(A,23X,A,I12)') 'Beta MO coefficients','R   N=', ncoeff
 write(fid1,'(5(1X,ES15.8))') mo_a
 deallocate(mo_a)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid)
 close(fid1)
end subroutine fch_r2u

