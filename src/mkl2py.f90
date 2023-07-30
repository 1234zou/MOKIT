! written by jxzou at 20230728: a wrapper of mkl2fch and bas_fch2py for ORCA ->
! PySCF

program mkl2py
 use util_wrapper, only: mkl2fch_wrap
 implicit none
 integer :: i, SYSTEM, RENAME
 character(len=240) :: mklname, fchname, inpname0, inpname

 i = iargc()
 if(.not. (i==1 .or. i==2)) then
  write(6,'(/,A)') 'ERROR in program mkl2py: wrong command line argument!'
  write(6,'(A)')   'Example 1: mkl2py a.mkl'
  write(6,'(A,/)') 'Example 2: mkl2py a.mkl b.py'
  stop
 end if

 mklname = ' '
 call getarg(1, mklname)
 call require_file_exist(mklname)

 inpname = ' '
 if(i == 2) then
  call getarg(2, inpname)
 else
  call find_specified_suffix(mklname, '.mkl', i)
  inpname = mklname(1:i-1)//'.py'
 end if

 call find_specified_suffix(mklname, '.mkl', i)
 fchname = mklname(1:i-1)//'.fch'

 call mkl2fch_wrap(mklname, fchname, .false.)

 i = SYSTEM('bas_fch2py '//TRIM(fchname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in program mkl2py: failed to call utility bas_fch2py.'
  write(6,'(A)') 'mkl2py is a wrapper of mkl2fch and bas_fch2py, so bas_fch2py &
                 &is required to'
  write(6,'(A)') 'be called successfully.'
  stop
 end if

 call find_specified_suffix(fchname, '.fch', i)
 inpname0 = fchname(1:i-1)//'.py'

 i = RENAME(TRIM(inpname0), TRIM(inpname))
end program mkl2py

