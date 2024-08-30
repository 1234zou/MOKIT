! written by jxzou at 20230728: a wrapper of mkl2fch and bas_fch2py for ORCA ->
! PySCF

program mkl2py
 use util_wrapper, only: mkl2fch_wrap
 implicit none
 integer :: i, k, narg, irel, SYSTEM, RENAME
 character(len=240), allocatable :: str_arg(:)
 character(len=240) :: mklname, fchname, inpname0, inpname

 narg = iargc()
 if(narg<1 .or. narg>3) then
  write(6,'(/,A)') ' ERROR in program mkl2py: wrong command line argument!'
  write(6,'(A)')   ' Example 1: mkl2py a.mkl'
  write(6,'(A)')   ' Example 2: mkl2py a.mkl b.py'
  write(6,'(A)')   ' Example 3: mkl2py a.mkl -sfx2c'
  write(6,'(A,/)') ' Example 4: mkl2py a.mkl b.py -sfx2c'
  stop
 end if

 allocate(str_arg(narg))
 do i = 1, narg, 1
  call getarg(i, str_arg(i))
 end do ! for i

 mklname = str_arg(1)
 call require_file_exist(mklname)
 call find_specified_suffix(mklname, '.mkl', i)
 inpname = mklname(1:i-1)//'.py'
 k = 2; irel = -1

 if(narg > 1) then
  i = LEN_TRIM(str_arg(2))
  if(str_arg(2)(i-3:i) == '.py') then
   inpname = str_arg(2)
   k = 3
  end if
  if(k <= narg) then
   if(TRIM(str_arg(k)) == '-sfx2c') then
    irel = -3
   else
    write(6,'(/,A)') 'ERROR in subroutine mkl2py: wrong command line argument!'
    write(6,'(A)')   'Example: mkl2py a.mkl -sfx2c'
    stop
   end if
  end if
 end if

 deallocate(str_arg)
 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 fchname = mklname(1:i-1)//'.fch'

 call mkl2fch_wrap(mklname=mklname,fchname=fchname,irel=irel)

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

