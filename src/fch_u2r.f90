! written by jxzou at 20180215
! This is a program/subroutine used to transform a Gaussian UHF type .fchk
! into a RHF/ROHF type .fchk file.

! Note: the Spin SCF Density section will be deleted.

! updated by jxzou at 20180328: input file can be a RHF type .fchk file, nothing will be done
! updated by jxzou at 20180424: output file can be a ROHF type .fchk file
! updated by jxozu at 20180818: fix the ROHF bug in Route Section ('uhf' detected in many lines)
! updated by jxozu at 20200202: simplify code
! updated by jxozu at 20200402: rename from fchk_uhf2rhf to fch_u2r

program main
 implicit none
 integer :: i
 character(len=240) :: fchname, newfch

 i = iargc()
 if(.not. (i==1 .or. i==2)) then
  write(6,'(/,A)') 'ERROR in subroutine fch_u2r: .fch(k) file are required.'
  write(6,'(A)')   'Example 1: fch_u2r a.fch       (a_r.fch generated)'
  write(6,'(A,/)') 'Example 2: fch_u2r a.fch b.fch (generate a_r.fch and rename&
                  & it to b.fch)'
  stop
 end if

 call getarg(1, fchname)
 call require_file_exist(fchname)

 newfch = REPEAT(' ',240)
 if(i == 2) call getarg(2, newfch)
 call fch_u2r(fchname, newfch)
end program main

! transform a Gaussian UHF type .fchk into a RHF/ROHF type .fchk file
subroutine fch_u2r(fchname, newfch)
 implicit none
 integer :: i, k, nalpha, nbeta, fchid, fchid1, RENAME
 character(len=240), intent(in) :: fchname, newfch
 character(len=240) :: fchname1, buf
 logical :: rhf, no_remain

 k = 0; buf = ' '; fchname1 = ' '; rhf = .false.; no_remain = .false.

 call find_specified_suffix(fchname, '.fch', k)
 fchname1 = fchname(1:k-1)//'_r.fch'
 open(newunit=fchid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fchid1,file=TRIM(fchname1),status='replace')

 ! Step 1: get nalpha and nbeta
 do while(.true.)
  read(fchid,'(A)') buf
  if(buf(1:25) == 'Number of alpha electrons') exit
 end do
 BACKSPACE(fchid)
 read(fchid,'(A51,I10)') buf, nalpha
 read(fchid,'(A51,I10)') buf, nbeta

 ! Step 2: modify the second line: UHF->RHF or ROHF
 rewind(fchid)
 read(fchid,'(A)') buf
 write(fchid1,'(A)') TRIM(buf)
 read(fchid,'(A)') buf
 if(buf(11:13) == 'UHF') then
  if(nalpha == nbeta) then
   buf(11:13) = 'RHF'
  else
   buf(11:14) = 'ROHF'
  end if
 end if
 write(fchid1,'(A)') TRIM(buf)

 ! Step 3: modify the 'Route' section
 if(nalpha /= nbeta) then
  do k = 1, 14
   read(fchid,'(A)') buf
   write(fchid1,'(A)') TRIM(buf)
   if(buf(1:5)=='Route' .or. buf(1:5) == 'Charg') exit
  end do ! for i

  if(buf(1:5)=='Route') then
   do while(.true.)
    read(fchid,'(A)') buf
    if(buf(1:6) == 'Charge') then
     BACKSPACE(fchid)
     exit
    end if
    k = INDEX(buf,'UHF')
    if(k == 0) k = INDEX(buf,'uhf')
    if(k == 0) k = INDEX(buf,'Uhf')
    if(k /= 0) then
     buf(k+4:240) = buf(k+3:239)
     buf(k:k+3) = 'ROHF'
    end if
    write(fchid1,'(A)') TRIM(buf)
   end do ! for while
  end if
 end if

 ! Step 4: set the 1st value of ILSW as 0
 do while(.true.)
  read(fchid,'(A)') buf
  write(fchid1,'(A)') TRIM(buf)
  if(buf(1:4) == 'ILSW') exit
 end do
 read(fchid,'(A)') buf
 buf(12:12) = '0'
 write(fchid1,'(A)') TRIM(buf)

 ! Step 5: change IOpCl value to 0
 do while(.true.)
  read(fchid,'(A)') buf
  if(buf(1:5) == 'IOpCl') exit
  write(fchid1,'(A)') TRIM(buf)
  if(buf(1:8) == 'Alpha Or') then
   BACKSPACE(fchid)
   exit
  end if
 end do ! for while
 write(fchid1,'(A5,38X,A1,16X,A1)') 'IOpCl','I','0'
 write(fchid1,'(A5,38X,A1,16X,A1)') 'IROHF','I','0'
 read(fchid,'(A)') buf

 ! Step 6: skip the 'Beta Orbital Energies'
 do while(.true.)
  read(fchid,'(A)') buf
  if(buf(1:12) == 'Beta Orbital') exit
  if(buf(1:8) == 'Alpha MO') exit
  write(fchid1,'(A)') TRIM(buf)
 end do ! for while
 if(buf(1:8) == 'Alpha MO') rhf = .true.
 if(.not. rhf) then
  do while(.true.)
   read(fchid,'(A)') buf
   if(buf(1:8) == 'Alpha MO') exit
  end do
 end if
 BACKSPACE(fchid)

 ! Step 7: copy the Alpha MO coefficients and skip the Beta MO coefficients
 do while(.true.)
  read(fchid,'(A)') buf
  if(buf(1:7) == 'Beta MO') exit
  if(buf(1:11) == 'Orthonormal') exit
  if(buf(1:9) == 'Total SCF') exit
  write(fchid1,'(A)') TRIM(buf)
 end do ! for while

 if(.not. rhf) then
  do while(.true.)
   read(fchid,'(A)') buf
   if(buf(1:11) == 'Orthonormal') exit
   if(buf(1:9) == 'Total SCF') exit
  end do
 end if
 BACKSPACE(fchid)

 ! Step 8: skip the 'Spin SCF Density'
 do while(.true.)
  read(fchid,'(A)') buf
  if(buf(1:8) == 'Spin SCF') exit
  if(buf(1:16) == 'Mulliken Charges') exit
  write(fchid1,'(A)') TRIM(buf)
 end do ! for while

 if(.not. rhf) then
  do while(.true.)
   read(fchid,'(A)',iostat=i) buf
   if(i /= 0) then
    no_remain = .true.
    exit
   end if
   if(buf(1:16) == 'Mulliken Charges') exit
  end do ! for while
 end if

 if(.not. no_remain) then
  BACKSPACE(fchid)

  ! copy the remaining content
  do while(.true.)
   read(fchid,'(A)',iostat=k) buf
   if(k /= 0) exit
   write(fchid1,'(A)') TRIM(buf)
  end do ! for while
  ! copy done
 end if

 close(fchid)
 close(fchid1)
 if(LEN_TRIM(newfch) > 0) k = RENAME(TRIM(fchname1), TRIM(newfch))
end subroutine fch_u2r

