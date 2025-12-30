! written by jxzou at 20230728: a wrapper of mkl2fch and fch2inp for ORCA ->
! GAMESS

program mkl2inp
 use util_wrapper, only: mkl2fch_wrap, fch2inp_wrap
 implicit none
 integer :: i, k, narg, irel, RENAME
 character(len=240), allocatable :: str_arg(:)
 character(len=240) :: mklname, fchname, inpname0, inpname

 narg = iargc()
 if(narg<1 .or. narg>3) then
  write(6,'(/,A)') ' ERROR in program mkl2inp: wrong command line argument!'
  write(6,'(A)')   ' Example 1: mkl2inp a.mkl'
  write(6,'(A)')   ' Example 2: mkl2inp a.mkl b.inp'
  write(6,'(A)')   ' Example 3: mkl2inp a.mkl -dkh2'
  write(6,'(A,/)') ' Example 4: mkl2inp a.mkl b.inp -dkh2'
  stop
 end if

 allocate(str_arg(narg))
 do i = 1, narg, 1
  call getarg(i, str_arg(i))
 end do ! for i

 mklname = str_arg(1)
 call require_file_exist(mklname)
 call find_specified_suffix(mklname, '.mkl', i)
 inpname = mklname(1:i-1)//'.inp'
 k = 2; irel = -1

 if(narg > 1) then
  i = LEN_TRIM(str_arg(2))
  if(str_arg(2)(i-3:i) == '.inp') then
   inpname = str_arg(2)
   k = 3
  end if
  if(k <= narg) then
   if(TRIM(str_arg(k)) == '-dkh2') then
    irel = 2
   else
    write(6,'(/,A)') 'ERROR in subroutine mkl2inp: wrong command line argument!'
    write(6,'(A)')   'Example: mkl2inp a.mkl -dkh2'
    stop
   end if
  end if
 end if

 deallocate(str_arg)
 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname=mklname,fchname=fchname,irel=irel)
 call fch2inp_wrap(fchname, .false., 0, 0, .false., .false., .false.)
 call delete_file(fchname)
 call find_specified_suffix(fchname, '.fch', i)
 inpname0 = fchname(1:i-1)//'.inp'

 i = RENAME(TRIM(inpname0), TRIM(inpname))
end program mkl2inp

