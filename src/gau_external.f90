! written by jxzou at 20231206: utility called in Gaussian external

program gau_external
 use util_wrapper, only: gbw2mkl, mkl2fch_wrap
 implicit none
 integer :: i, k, fid, natom, SYSTEM
 real(kind=8) :: e
 real(kind=8), allocatable :: grad(:)
 character(len=240) :: gjfname, inpname, hfile, bakfile, numfile, gbwname, &
  mklname, fchname, new_fch, gradname, orca_path
 character(len=720) :: EIn, EOu
 logical :: alive, update_coor

 i = iargc()
 if(i < 3) then
  write(6,'(/,A)') ' ERROR in program gau_external: wrong command line arguments!'
  write(6,'(A,/)') ' Example: gau_external R Gau-17156.EIn Gau-17156.EOu'
  stop
 end if

 call getarg(2, EIn)
 call getarg(3, EOu)

 call get_gjfname_from_EIn(EIn, gjfname)
 i = LEN_TRIM(gjfname)
 inpname = gjfname(1:i-6)//'_o.inp'
 gbwname = gjfname(1:i-6)//'_o.gbw'
 mklname = gjfname(1:i-6)//'_o.mkl'
 fchname = gjfname(1:i-6)//'_o.fch'
 gradname = gjfname(1:i-6)//'_o.engrad'
 numfile = gjfname(1:i-6)//'_o.num'
 bakfile = gjfname(1:i-6)//'_o.bak'
 hfile = gjfname(1:i-6)//'_o.nc'

 open(newunit=fid,file=TRIM(numfile),status='old',position='rewind')
 read(fid,*) k
 BACKSPACE(fid)
 write(fid,'(I0)') k+1
 close(fid)
 write(new_fch,'(A,I0,A)') gjfname(1:i-6)//'_', k, '.fch'

 inquire(file=hfile,exist=alive)
 if(alive) then
  update_coor = .false.
  call delete_file(hfile) ! supposed not exist in the next cycle
 else
  update_coor = .true.
  i = SYSTEM('replace_xyz_in_inp '//TRIM(EIn)//' '//TRIM(inpname)//' -orca')
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in program gau_external: failed to call utility repl&
                    &ace_xyz_in_inp.'
   stop
  end if
 end if

 call get_orca_path(orca_path)
 call submit_orca_job(orca_path, inpname, .false.)
 call read_natom_from_engrad(gradname, natom)
 allocate(grad(3*natom))
 call read_grad_from_engrad(gradname, natom, e, grad)

 inquire(file=TRIM(bakfile),exist=alive)
 if(alive) then   ! update MO coefficients (and coordinates if necessary)
  call gbw2mkl(gbwname)
  call sys_copy_file(fchname, new_fch, .false.)
  call mkl2fch_wrap(mklname, new_fch)
  call write_grad_into_fch(new_fch, natom, grad)
  if(update_coor) call replace_coor_in_fch_by_engrad(gradname, new_fch)
 end if

 call write_EOu(EOu, e, natom, grad)
 deallocate(grad)
end program gau_external

! find the corresponding .gjf filename from a given Gaussian .EIn filename by
! backtracing PIDs
subroutine get_gjfname_from_EIn(EIn, gjfname)
 implicit none
 integer :: i, k, pid
 character(len=2100) :: longbuf
 character(len=700), intent(in) :: EIn
 character(len=240), intent(out) :: gjfname

 k = LEN_TRIM(EIn)
 i = INDEX(EIn(1:k), 'Gau-', back=.true.)
 read(EIn(i+4:k-4),*) pid
 call get_first_command_from_pid(pid, longbuf)
 k = INDEX(longbuf, '.inp', back=.true.)
 if(k == 0) then
  write(6,'(/,A)') "ERROR in subroutine get_gjfname_from_EIn: no '.inp' keyword&
                   & found in Gaussian"
  write(6,'(A)') 'l402 command. longbuf:'
  write(6,'(A)') TRIM(longbuf)
  stop
 end if

 i = INDEX(longbuf(1:k-1), 'Gau-', back=.true.)
 read(longbuf(i+4:k-1),*) pid
 call get_first_command_from_pid(pid, longbuf)
 k = LEN_TRIM(longbuf)
 i = INDEX(longbuf(1:k), ' ', back=.true.)
 gjfname = longbuf(i+1:k)
end subroutine get_gjfname_from_EIn

! get the 1st command line from `ps aux|grep $PID`
subroutine get_first_command_from_pid(pid, longbuf)
 implicit none
 integer :: i, j, k, fid, SYSTEM
 integer, intent(in) :: pid
 character(len=8) :: s_pid
 character(len=10) :: tfile
 character(len=2100), intent(out) :: longbuf

 write(s_pid,'(I0)') pid
 write(tfile,'(I0,A2)') pid, '.t'
 i = SYSTEM("ps aux|grep '"//TRIM(s_pid)//"' >"//TRIM(tfile))

 open(newunit=fid,file=TRIM(tfile),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) longbuf
  if(i /= 0) exit
  j = INDEX(longbuf, ' ')
  read(longbuf(j+1:),*) k
  if(k == pid) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine get_first_command_from_pid: internal er&
                   &ror.'
  write(6,'(A,I0,A)') 'PID=', pid, ', longbuf='
  write(6,'(A)') TRIM(longbuf)
  close(fid)
  stop
 end if

 close(fid,status='delete')
end subroutine get_first_command_from_pid

