! written by jxzou at 20231206: utility called in Gaussian external

program gau_external
 implicit none
 integer :: i, SYSTEM
 character(len=240) :: gjfname, inpname, hfile, orca_path, gradname
 character(len=720) :: EIn, EOu
 logical :: alive

 call getarg(2, EIn)
 call getarg(3, EOu)

 call get_gjfname_from_EIn(EIn, gjfname)
 i = LEN_TRIM(gjfname)
 inpname = gjfname(1:i-6)//'_o.inp'
 gradname = gjfname(1:i-6)//'_o.engrad'
 hfile = gjfname(1:i-6)//'_o.nc'
 ! here the gjfname, inpname and hfile are correlated, for example,
 !  h2o_g.gjf, h2o_o.inp, and h2o_o.nc

 inquire(file=hfile,exist=alive)
 if(alive) then
  call delete_file(hfile) ! supposed not exist in the next cycle
 else
  i = SYSTEM('replace_xyz_in_inp '//TRIM(EIn)//' '//TRIM(inpname)//' -orca')
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in program gau_external: failed to call utility gau_&
                    &external.'
   stop
  end if
 end if

 call get_orca_path(orca_path)
 call submit_orca_job(orca_path, inpname, .false.)
 ! TODO: using initial Hessian calculated from empirical methods and write into EOu
 call engrad2EOu(gradname, EOu)
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

