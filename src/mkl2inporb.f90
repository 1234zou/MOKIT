! written by jxzou at 20230728: a wrapper of mkl2fch and fch2inporb for ORCA ->
! (Open)Molcas

program mkl2inporb
 use util_wrapper, only: mkl2fch_wrap, fch2inporb_wrap
 implicit none
 integer :: i, k, narg, irel
 character(len=240), allocatable :: str_arg(:)
 character(len=240) :: mklname, fchname, inpname

 narg = iargc()
 if(narg<1 .or. narg>3) then
  write(6,'(/,A)') ' ERROR in program mkl2inporb: wrong command line argument!'
  write(6,'(A)')   ' Example 1: mkl2inporb a.mkl'
  write(6,'(A)')   ' Example 2: mkl2inporb a.mkl b.input'
  write(6,'(A)')   ' Example 3: mkl2inporb a.mkl -dkh2'
  write(6,'(A)')   ' Example 4: mkl2inporb a.mkl -sfx2c'
  write(6,'(A,/)') ' Example 5: mkl2inporb a.mkl b.input -sfx2c'
  stop
 end if

 allocate(str_arg(narg))
 do i = 1, narg, 1
  call getarg(i, str_arg(i))
 end do ! for i

 mklname = str_arg(1)
 call require_file_exist(mklname)
 call find_specified_suffix(mklname, '.mkl', i)
 inpname = mklname(1:i-1)//'.input'
 k = 2; irel = -1

 if(narg > 1) then
  i = LEN_TRIM(str_arg(2))
  if(str_arg(2)(i-6:i) == '.input') then
   inpname = str_arg(2)
   k = 3
  end if
  if(k <= narg) then
   select case(TRIM(str_arg(k)))
   case('-sfx2c')
    irel = -3
   case('-dkh2')
    irel = 2
   case default
    write(6,'(/,A)') 'ERROR in subroutine mkl2inporb: wrong command line argument!'
    write(6,'(A)')   'Example: mkl2inporb a.mkl -sfx2c'
    stop
   end select
  end if
 end if

 deallocate(str_arg)
 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname=mklname,fchname=fchname,irel=irel)
 call fch2inporb_wrap(fchname, .false., inpname)
 call delete_file(TRIM(fchname))
end program mkl2inporb

