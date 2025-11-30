! read basic information like nbf/nif, mult, cart

! Read nbf and nif from a GAMESS .dat file. The number of basis functions in
! GAMESS .dat is always at Cartesian basis, no matter the actual calculation
! is performed at Cartesian or spherical harmonic basis.
subroutine read_cart_nbf_nif_from_dat(datname, only_nbf, nbf, nif)
 implicit none
 integer :: i, j, k, fid
 integer, intent(out) :: nbf, nif
!f2py intent(out) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: datname
!f2py character(len=240), intent(in) :: datname
 logical, intent(in) :: only_nbf
!f2py intent(in) :: only_nbf

 nbf = 0; nif = 0
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
  write(6,'(/,A)') "ERROR in subroutine read_cart_nbf_nif_from_dat: no '$VEC' f&
                   &ound in file "//TRIM(datname)
  close(fid)
  stop
 end if

 j = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:2) /= ' 1') exit
  j = j + 1
 end do ! for while

 k = j ! backup

 BACKSPACE(fid)
 BACKSPACE(fid)
 read(fid,'(A)') buf
 nbf = (j-1)*5 + (LEN_TRIM(buf)-5)/15

 if(only_nbf) then ! no need to find nif
  close(fid)
  return
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:2) == '$') then
   call upper(buf(3:5))
   if(buf(2:5) == '$END') exit
  end if
  j = j + 1
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cart_nbf_nif_from_dat: no '$END' c&
                   &orresponds to '$VEC'"
  write(6,'(A)') 'in file '//TRIM(datname)
  stop
 end if

 nif = j/k
end subroutine read_cart_nbf_nif_from_dat

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
end subroutine read_nbf_and_nif_from_orb

! read spin multiplicity from a specified GAMESS output file (.gms)
subroutine read_mult_from_gms_gms(gmsname, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: mult
!f2py intent(out) :: mult
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname
!f2py intent(in) :: gmsname

 mult = 1
 open(newunit=fid,file=TRIM(gmsname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:10) == 'SPIN MULT') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mult_from_gms_gms: 'SPIN MULT' not&
                   & found in file "//TRIM(gmsname)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) mult
end subroutine read_mult_from_gms_gms