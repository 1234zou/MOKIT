! written by jxzou at 20180208
! updated by jxzou at 20200402: renamed from fchk_mo_copy to fch_mo_copy

! This is a program/subroutine used to copy MOs from one .fchk file into another .fchk file
! Usage: ./fch_mo_copy a.fchk b.fchk [-aa]
! '-aa' means copying the Alpha MO from a.fchk into the Alpha MO part in b.fchk,
! '-ab' means copying the Alpha MO from a.fchk into the Beta MO part in b.fchk, and
! '-ba' and '-bb' are in a simila way.

program main
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=3) :: ab
 character(len=240) :: fchname1, fchname2

 i = 0
 ab = '-aa'
 fchname1 = ' '
 fchname2 = ' '

 i = iargc()
 if(i<1 .or. i>3) then
  write(iout,'(/,A)') ' ERROR in subroutine fch_mo_copy: wrong command line arguments!'
  write(iout,'(A)')   ' Example1: ./fch_mo_copy a.fch b.fch     (copy Alpha MOs in a.fch to Alpha MOs in b.fch)'
  write(iout,'(A)')   ' Example2: ./fch_mo_copy a.fch b.fch -aa (eqv. to Example1)'
  write(iout,'(A)')   ' Example3: ./fch_mo_copy a.fch b.fch -ab (copy Alpha MOs in a.fch to Beta  MOs in b.fch)'
  write(iout,'(A)')   ' Example4: ./fch_mo_copy a.fch b.fch -ba (copy Beta  MOs in a.fch to Alpha MOs in b.fch)'
  write(iout,'(A)')   ' Example5: ./fch_mo_copy a.fch b.fch -bb (copy Beta  MOs in a.fch to Beta  MOs in b.fch)'
  stop
 end if

 call getarg(1, fchname1)
 call getarg(2, fchname2)
 if(i == 3) then
  call getarg(3,ab)
  if(ab(1:1) /= '-') then
   write(iout,'(A)') "ERROR in subroutine fch_mo_copy: the 3rd argument must be '-aa'/'-ab'/'-ba'/'-bb'."
   stop
  end if
 end if

 call fch_mo_copy(fchname1, fchname2, ab)
 stop
end program main

! copy Alpha/Beta MOs from a .fch file into Alpha/Beta MOs of another .fch file
subroutine fch_mo_copy(fchname1, fchname2, ab)
 implicit none
 integer :: i, k, remain
 integer :: ncoeff1, ncoeff2
 integer, parameter :: fid1 = 11, fid2 = 12, fid3 = 13
 integer :: RENAME
 character(len=3), intent(in) :: ab
 character(len=13) :: mark
 character(len=13), parameter :: mark1 = 'Alpha MO coef'
 character(len=13), parameter :: mark2 = 'Beta MO coeff'
 character(len=240) :: buffer, fchname3
 character(len=240), intent(in) :: fchname1, fchname2

 k = 0    ! initialization
 remain = 0
 ncoeff1 = 0
 ncoeff2 = 0
 mark = ' '
 buffer = ' '
 fchname3 = ' '

 k = index(fchname2,'.fch')
 fchname3 = fchname2(1:k-1)//'.fchk.copy.tmp'

 k = 0
 if(ab(2:2)  == 'a') then
  mark = mark1
 else
  mark = mark2
 end if
 open(unit=fid1,file=fchname1,status='old',position='rewind')
 do while(.true.)
  read(fid1,'(A)',iostat=k) buffer
  if(k/=0 .or. buffer(1:13)==mark) exit
 end do
 if(k /= 0) then
  write(*,'(A)') 'ERROR in subroutine fch_mo_copy: '//TRIM(fchname1)
  write(*,'(A)') 'No '//mark//' found!'
  close(fid1)
  return
 end if
 call get_int_after_flag(buffer, '=', ncoeff1)

 k = 0
 buffer = ' '
 if(ab(3:3)  == 'a') then
  mark = mark1
 else
  mark = mark2
 end if
 open(unit=fid2,file=fchname2,status='old',position='rewind')
 open(unit=fid3,file=fchname3,status='replace')
 do while(.true.)
  read(fid2,'(A)',iostat=k) buffer
  if(k/=0 .or. buffer(1:13)==mark) exit
  write(fid3,'(A)') TRIM(buffer)
 end do
 if(k /= 0) then
  write(*,'(A)') 'ERROR in subroutine fch_mo_copy: '//TRIM(fchname2)
  write(*,'(A)') 'No '//mark//' found!'
  close(fid1)
  close(fid2)
  close(fid3,status='delete')
  return
 end if
 write(fid3,'(A)') TRIM(buffer)
 call get_int_after_flag(buffer, '=', ncoeff2)

 if(ncoeff1 /= ncoeff2) then
  write(*,'(A)') 'Warning: the total number of MO coefficients in the two .fchk files:'
  write(*,'(A)') TRIM(fchname1)//' '//TRIM(fchname2)//' is not equal!!!'
  write(*,'(A)') "The reasons could possibly be: 1) basis sets in two .fchk files are not identical;"
  write(*,'(A)') "2) linear dependence may exist in these two files."
 end if

 ! copy the MOs in fchk1 into fchk3
 k = ncoeff1/5
 remain = ncoeff1 - 5*k
 do i = 1, k, 1
  read(fid1,'(A)') buffer
  write(fid3,'(A)') TRIM(buffer)
 end do
 if(remain > 0) then
  read(fid1,'(A)') buffer
  write(fid3,'(A)') TRIM(buffer)
 end if
 ! copy done
 close(fid1)

 ! skip the MOs in fchk2
 do while(.true.)
  read(fid2,'(A)',iostat=i) buffer
  if(i /= 0) exit
  k = index(buffer, '=')
  if(k /= 0) exit
 end do
 ! skip done

 if(i == 0) then
  BACKSPACE(fid2)
  ! copy the rest of fchk2 into fchk3
  do while(.true.)
   read(fid2,'(A)',iostat=k) buffer
   if(k /= 0) exit
   write(fid3,'(A)') TRIM(buffer)
  end do
  ! copy done
 end if

 close(fid2,status='delete')
 close(fid3)
 i = RENAME(fchname3,fchname2)
 if(i /= 0) then
  write(*,'(A)') 'ERROR in subroutine fch_mo_copy: RENAME .fchk error!'
  write(*,'(A)') TRIM(fchname2)
 end if
 return
end subroutine fch_mo_copy

subroutine get_int_after_flag(buffer,flag,k)
 implicit none
 integer,intent(out) :: k
 character(len=1),intent(in) :: flag
 character(len=240),intent(inout) :: buffer

 k = index(buffer,flag)
 if(k == 0) then
  write(*,'(A)') 'ERROR in subroutine get_int_after_flag: sth wrong in this line:'
  write(*,'(A)') TRIM(buffer)
  stop
 end if
 buffer(1:k) = ' '
 buffer = ADJUSTL(buffer)
 read(buffer,*) k
 return
end subroutine get_int_after_flag

