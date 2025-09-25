! written by jxzou at 20210113: file operations

subroutine require_file_exist(fname)
 implicit none
 character(len=240), intent(in) :: fname
 logical :: alive

 inquire(file=TRIM(fname),exist=alive)
 if(.not. alive) then
  write(6,'(/,A)') 'ERROR in subroutine require_file_exist: file does not exist!'
  write(6,'(A)') 'Filename='//TRIM(fname)
  stop
 end if
end subroutine require_file_exist

subroutine open_file(fname, is_rewind, fid)
 implicit none
 integer :: i
 integer, intent(out) :: fid
 character(len=240), intent(in) :: fname
 logical, intent(in) :: is_rewind

 call require_file_exist(fname)
 if(is_rewind) then
  open(newunit=fid,file=TRIM(fname),status='old',position='rewind',iostat=i)
 else
  open(newunit=fid,file=TRIM(fname),status='old',position='append',iostat=i)
 end if

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine open_file: failed to open file '//&
                    TRIM(fname)
  stop
 end if
end subroutine open_file

! delete the specified file (if not exist, return)
subroutine delete_file(fname)
 implicit none
 integer :: fid
 character(len=*), intent(in) :: fname
 logical :: alive

 inquire(file=TRIM(fname),exist=alive)

 if(alive) then
  inquire(file=TRIM(fname),opened=alive,number=fid)
  if(.not. alive) open(newunit=fid,file=TRIM(fname),status='old')
  close(fid,status='delete')
 end if
end subroutine delete_file

! delete a set of files
subroutine delete_files(n, fname)
 implicit none
 integer :: i
 integer, intent(in) :: n
 character(len=240), intent(in) :: fname(n)

 do i = 1, n, 1
  call delete_file(TRIM(fname(i)))
 end do ! for i
end subroutine delete_files

subroutine delete_files_in_path(path, n, fname)
 implicit none
 integer :: i
 integer, intent(in) :: n
 character(len=240), intent(in) :: path, fname(n)

 do i = 1, n, 1
  call delete_file(TRIM(path)//'/'//TRIM(fname(i)))
 end do ! for i
end subroutine delete_files_in_path

subroutine delete_cfour_tmp_files(nproc)
 implicit none
 integer :: i
 integer, intent(in) :: nproc
 integer, parameter :: nfile = 25
 character(len=7) :: dirname
 character(len=12) :: filelist(nfile)
 data filelist /'BASINFO.DATA','DIPOL','FILES','fort.81','fort.82','fort.83',&
  'GAMLAM','HFNES1','HFNES2','HFNES3','HFNES4','HFNES5','HFNES6','HFNES7',&
  'IIII','JAINDX','JMOLplot','JOBARC','MOABCD','MOINTS','MOL','MOLDEN','ncpu',&
  'NTOTAL','OPTARC'/

 if(nproc > 0) then
  do i = 0, nproc-1, 1
   write(dirname,'(A4,I3.3)') 'rank', i
   call remove_dir(TRIM(dirname))
  end do ! for i
 end if

 do i = 1, nfile, 1
  call delete_file(TRIM(filelist(i)))
 end do ! for i
end subroutine delete_cfour_tmp_files

! copy file fname1 to fname2 (if delete=.True., delete fname1)
subroutine copy_file(fname1, fname2, delete)
 implicit none
 integer :: i, fid1, fid2
 character(len=240) :: buf
 character(len=240), intent(in) :: fname1, fname2
 logical, intent(in) :: delete

 call require_file_exist(fname1)
 if(LEN_TRIM(fname2) == 0) then
  write(6,'(/,A)') 'ERROR in subroutine copy_file: empty filename for fname2.'
  write(6,'(A)') 'For bug tracking: fname1='//TRIM(fname1)
  stop
 end if

 open(newunit=fid1,file=TRIM(fname1),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(fname2),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid2)
 if(delete) then
  close(fid1, status='delete')
 else
  close(fid1)
 end if
end subroutine copy_file

! Copy a file using system commands. If delete=.True., delete fname1.
subroutine sys_copy_file(fname1, fname2, delete)
 implicit none
 integer :: i, SYSTEM
 character(len=*), intent(in) :: fname1, fname2
 logical, intent(in) :: delete

 if(LEN_TRIM(fname1) == 0) then
  write(6,'(/,A)') 'ERROR in subroutine sys_copy_file: empty filename for fname1.'
  write(6,'(A)') 'For bug tracking: fname2='//TRIM(fname2)
  stop
 end if

#ifdef _WIN32
 i = SYSTEM('copy /Y '//TRIM(fname1)//' '//TRIM(fname2)//' > NUL')
#else
 i = SYSTEM('cp '//TRIM(fname1)//' '//TRIM(fname2))
#endif

 if(delete) then
#ifdef _WIN32
 i = SYSTEM('del '//TRIM(fname1))
#else
 i = SYSTEM('rm -f '//TRIM(fname1))
#endif
 end if

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine sys_copy_file: fail to copy file from'
  write(6,'(A)') TRIM(fname1)//' to '//TRIM(fname2)
  stop
 end if
end subroutine sys_copy_file

! copy a file from a path/directory to the current directory
subroutine copy_file_from_a_path(path, fname)
 implicit none
 character(len=1) :: s
 character(len=700) :: buf
 character(len=*), intent(in) :: path, fname

#ifdef _WIN32
 s = '\'
#else
 s = '/'
#endif

 buf = TRIM(path)//s//TRIM(fname)
 call sys_copy_file(TRIM(buf), TRIM(fname), .false.)
end subroutine copy_file_from_a_path

! move a file into a specified directory
subroutine move_file(fname, dirname)
 implicit none
 integer :: i, SYSTEM
 character(len=*), intent(in) :: fname, dirname

#ifdef _WIN32
 i = SYSTEM('move '//TRIM(fname)//' '//TRIM(dirname)//'\')
#else
 i = SYSTEM('mv '//TRIM(fname)//' '//TRIM(dirname)//'/')
#endif

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine move_file: failed to move file to the s&
                   &pecified directory.'
  write(6,'(A)') 'fname='//TRIM(fname)
  write(6,'(A)') 'dirname='//TRIM(dirname)
  stop
 end if
end subroutine move_file

! move several files into a specified directory
subroutine move_files(n, fname, dirname)
 implicit none
 integer :: i
 integer, intent(in) :: n
 character(len=240), intent(in) :: fname(n), dirname

 do i = 1, n, 1
  call move_file(TRIM(fname(i)), TRIM(dirname))
 end do ! for i
end subroutine move_files

! create a directory
subroutine create_dir(dirname)
 implicit none
 integer :: i, SYSTEM
 character(len=*), intent(in) :: dirname

#ifdef _WIN32
 i = SYSTEM('CD '//TRIM(dirname))
 if(i == 0) then
  i = SYSTEM('CD ..')
 else
  i = SYSTEM('MD '//TRIM(dirname))
 end if
#else
 i = SYSTEM('mkdir -p '//TRIM(dirname))
#endif
end subroutine create_dir

! remove/delete a directory. USE WITH CAUTION!!!
subroutine remove_dir(dirname)
 implicit none
 integer :: i, SYSTEM
 character(len=*), intent(in) :: dirname

#ifdef _WIN32
 i = SYSTEM('rd '//dirname)
#else
 ! in case that important directories are deleted, check the dirname firstly
 i = LEN_TRIM(dirname)
 if(i == 1) then
  write(6,'(A)') 'ERROR in subroutine remove_dir: wrong directory name: '//&
                  TRIM(dirname)
  stop
 end if

 if(i > 3) then
  select case(dirname(1:4))
  case('/usr','/bin','/lib','/etc')
   write(6,'(A)') 'ERROR in subroutine remove_dir: dangerous directory name: '&
                   //TRIM(dirname)
   stop
  case default
  end select
 end if

 if(i > 4) then
  if(dirname(1:5) == '/sbin') then
   write(6,'(A)') 'ERROR in subroutine remove_dir: dangerous directory name: '&
                   //TRIM(dirname)
   stop
  end if
 end if

 if(i > 5) then
  if(dirname(1:6) == '/lib64') then
   write(6,'(A)') 'ERROR in subroutine remove_dir: dangerous directory name: '&
                   //TRIM(dirname)
   stop
  end if
 end if

 i = SYSTEM('rm -rf '//dirname)
#endif
end subroutine remove_dir

! check if MOKIT_ROOT present. Set it if not present.
subroutine check_mokit_root()
 implicit none
 integer :: i, SYSTEM
 character(len=240) :: mokit_root

 write(6,'(/,A)') 'Checking MOKIT_ROOT ... '
 mokit_root = ' '
 call getenv('MOKIT_ROOT', mokit_root)

 if (LEN_TRIM(mokit_root) == 0) then
  ! assume we are under conda install
  i = SYSTEM("echo `get_mokit_loc.py` > ~/.mokitrc ")
  !i = SYSTEM("echo `get_mokit_loc.py` ")
  !i = SYSTEM('echo $MOKIT_ROOT')
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine check_mokit_root: '
   write(6,'(A)') '    fail to set MOKIT_ROOT for conda-installed MOKIT.'
   write(6,'(A)') 'If MOKIT is installed via conda, please report this issue.'
   write(6,'(A)') 'Otherwise, it means your MOKIT_ROOT is not properly set.'
   stop
  else
   write(6,'(A)') 'reset MOKIT_ROOT for conda-installed version'
  end if
 endif
end subroutine check_mokit_root

! check whether UHF-type MOs are hold in a given .fch(k) file
subroutine check_uhf_in_fch(fchname, uhf)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical, intent(out) :: uhf
!f2py intent(out) :: uhf

 uhf = .false.
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i < 0) then
   exit ! end-of-file
  else if(i > 0) then
   write(6,'(/,A)') 'ERROR in subroutine check_uhf_in_fch: failed to read file &
                   &'//TRIM(fchname)
   close(fid)
   stop
  end if

  if(buf(1:7) == 'Beta MO') then
   uhf = .true.
   exit
  end if

  select case(buf(1:11))
  case('Orthonormal','Total SCF D','Mulliken Ch')
   exit
  end select
 end do ! for while

 close(fid)
end subroutine check_uhf_in_fch

! Check whether pure Cartesian or spherical harmonic type basis functions are
!  used in a specified .fch(k) file.
! icart=-1: mixed Cartesian and spherical harmonic
! icart= 0: only S/P/L angular momenta, can be viewed as either pure Cartesian
!           or spherical harmonic
! icart= 1: spherical harmonic includes any of 5D, 7F, 9G, ...
! icart= 2: pure Cartesian includes any of 6D, 10F, 15G, ...
subroutine find_icart_in_fch(fchname, allow_mixed, icart)
 implicit none
 integer :: i, k, fid
 integer, intent(out) :: icart
!f2py intent(out) :: icart
 integer, allocatable :: shltyp(:)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical, intent(in) :: allow_mixed
!f2py intent(in) :: allow_mixed

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == 'Shell types') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine find_icart_in_fch: missing 'Shell types&
                   &' in file"
  write(6,'(A)') TRIM(fchname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=', back=.true.)
 read(buf(i+1:),*) k
 allocate(shltyp(k), source=0)
 read(fid,*) shltyp
 close(fid)

 if(ANY(shltyp<-1) .and. ANY(shltyp>1)) then
  icart = -1
  if(.not. allow_mixed) then
   write(6,'(/,A)') 'ERROR in subroutine find_icart_in_fch: mixed Cartesian and&
                    & spherical harmonic'
   write(6,'(A)') '(like 6D 7F) is found. Only spherical hamonic (5D 7F) or pur&
                  &e Cartesian (6D 10F)'
   write(6,'(A)') 'is allowed.'
   stop
  end if
 else if(ALL(shltyp>-2) .and. ALL(shltyp<2)) then
  icart = 0
 else if(ALL(shltyp<2) .and. ANY(shltyp<-1)) then
  icart = 1
 else if(ALL(shltyp>-2) .and. ANY(shltyp>1)) then
  icart = 2
 end if

 deallocate(shltyp)
end subroutine find_icart_in_fch

! check whether UHF is used in a specified CFOUR output file
subroutine check_uhf_in_cfour_out(outname, uhf)
 implicit none
 integer :: i, k, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname
 logical, intent(out) :: uhf
!f2py intent(out) :: uhf

 uhf = .false.
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(9:23) == 'SCF reference f') then
   k = INDEX(buf, ':')
   buf = ADJUSTL(buf(k+1:))
   if(buf(1:3) == 'UHF') uhf = .true.
   exit
  end if
 end do ! for while

 close(fid)
end subroutine check_uhf_in_cfour_out

! export a real(kind=8) array into a plain .txt file
subroutine export_rarray2txt(txtname, label, n, a)
 implicit none
 integer :: i, fid
 integer, intent(in) :: n
!f2py intent(in) :: n
 real(kind=8), intent(in) :: a(n)
!f2py intent(in) :: a
!f2py depend(n) :: a
 character(len=25), parameter :: star = '*************************'
 character(len=*), intent(in) :: label
!f2py intent(in) :: label
 character(len=240), intent(in) :: txtname
!f2py intent(in) :: txtname

 open(newunit=fid,file=TRIM(txtname),status='replace')
 write(fid,'(A,I0)') 'n=', n
 write(fid,'(A)') ' '//star//' '//label//' '//star

 do i = 1, n, 1
  write(fid,'(I5,1X,F11.5)') i, a(i)
 end do ! for i

 close(fid)
end subroutine export_rarray2txt

! export a square matrix into a plain .txt file
subroutine export_mat_into_txt(txtname, n, mat, lower, label)
 implicit none
 integer :: i, j, k, m, na, fid
 integer, intent(in) :: n
 integer, allocatable :: a(:)
 real(kind=8), intent(in) :: mat(n,n)
 character(len=25), parameter :: star = '*************************'
 character(len=25), intent(in) :: label
 character(len=240), intent(in) :: txtname
 logical, intent(in) :: lower ! True: lower triangle

 open(newunit=fid,file=TRIM(txtname),status='replace')
 write(fid,'(A,L1)') 'Lower_Triangle=',lower
 write(fid,'(A,I0)') 'n=', n
 write(fid,'(A)') ' '//star//' '//label//' '//star

 m = n/5
 if(n-5*m > 0) m = m + 1

 if(lower) then ! lower triangle
  do i = 1, m, 1
   if(i < m) then
    na = 5
   else
    na = n - 5*(m-1)
   end if
   allocate(a(na))
   forall(j=1:na) a(j) = 5*i - 5 + j
   write(fid,'(5I14)') a(1:na)
   deallocate(a)

   k = 5*i - 4
   do j = k, n, 1
    write(fid,'(I6,5(1X,ES15.8))') j, mat(k:min(k+4,j),j)
   end do ! for j
  end do ! for i

 else           ! full matrix
  do i = 1, m, 1
   if(i < m) then
    na = 5
   else
    na = n - 5*(m-1)
   end if
   allocate(a(na))
   forall(j=1:na) a(j) = 5*i - 5 + j
   write(fid,'(5I14)') a(1:na)
   deallocate(a)

   k = 5*i - 4
   do j = 1, n, 1
    write(fid,'(I6,5(1X,ES15.8))') j, mat(j,k:k+na-1)
   end do ! for j
  end do ! for i
 end if

 close(fid)
end subroutine export_mat_into_txt

! simplify a .fch(k) file
subroutine simplify_fch(fchname)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 i = INDEX(fchname, '.fch', back=.true.)
 fchname1 = fchname(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:12)=='Mulliken Cha' .or. buf(1:12)=='Number of ex') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine simplify_fch: keywords for termination &
                   &not found in'
  write(6,'(A)') 'file '//TRIM(fchname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
end subroutine simplify_fch

! read nbf and nif from .fch(k) file
subroutine read_nbf_and_nif_from_fch(fchname, nbf, nif)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nbf, nif
!f2py intent(out) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 nbf = 0; nif = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

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

 ! In case that this is a Q-Chem .fch file, let's assume that nif is not just
 ! below nbf
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

! read spin multipliticity from a given .fch(k) file
subroutine read_mult_from_fch(fchname, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: mult
!f2py intent(out) :: mult
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call open_file(fchname, .true., fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == 'Mult') exit
 end do ! for while
 close(fid)

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_mult_from_fch: no 'Mult' found in&
                & file "//TRIM(fchname)
  stop
 end if
 read(buf(50:),*) mult
end subroutine read_mult_from_fch

! read spin multiplicity from a specified ORCA input file
subroutine read_mult_from_orca_inp(inpname, mult)
 implicit none
 integer :: i, j, itype, fid
 integer, intent(out) :: mult
!f2py intent(out) :: mult
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 mult = 1; itype = 1
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:1) == '*') exit       ! * xyz 0 1
  if(buf(1:7) == '%coords') then ! Mult = 1
   itype = 2
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mult_from_orca_inp: spin multiplic&
                   &ity cannot be'
  write(6,'(A)') 'found in file '//TRIM(inpname)
  close(fid)
  stop
 end if

 select case(itype)
 case(1)
  close(fid)
  i = LEN_TRIM(buf)
  j = INDEX(buf(1:i), ' ', back=.true.)
  read(buf(j+1:i),*) mult
 case(2)
  do i = 1, 4
   read(fid,'(A)') buf
   if(INDEX(buf(1:5),'Mult') > 0) exit
  end do ! for i
  close(fid)
  if(i == 5) then
   write(6,'(/,A,I0)') 'ERROR in subroutine read_mult_from_orca_inp: spin multi&
                       &plicity cannot be'
   write(6,'(A)') 'found in file '//TRIM(inpname)
   stop
  end if
  i = INDEX(buf, '=')
  read(buf(i+1:),*) mult
 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine read_mult_from_orca_inp: invalid ity&
                      &pe=', itype
  write(6,'(A)') 'inpname='//TRIM(inpname)
  close(fid)
  stop
 end select
end subroutine read_mult_from_orca_inp

! write MO coefficients into a specified binary file
subroutine write_mo2bin(binfile, nbf, nif, mo)
 implicit none
 integer :: fid
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: mo(nbf,nif)
!f2py intent(in) :: mo
!f2py depend(nbf,nif) :: mo
 character(len=240), intent(in) :: binfile
!f2py intent(in) :: binfile

 open(newunit=fid,file=TRIM(binfile),status='replace',form='unformatted')
 write(unit=fid) mo
 close(fid)
end subroutine write_mo2bin

! read MO coefficients from a specified binary file
subroutine read_mo_from_bin(binfile, nbf, nif, mo)
 implicit none
 integer :: fid
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(out) :: mo(nbf,nif)
!f2py intent(out) :: mo
!f2py depend(nbf,nif) :: mo
 character(len=240), intent(in) :: binfile
!f2py intent(in) :: binfile

 mo = 0d0
 open(newunit=fid,file=TRIM(binfile),status='old',form='unformatted')
 read(unit=fid) mo
 close(fid)
end subroutine read_mo_from_bin

subroutine write_grad_into_fch(fchname, natom, grad)
 implicit none
 integer :: fid
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 real(kind=8), intent(in) :: grad(3*natom)
!f2py intent(in) :: grad
!f2py depend(natom) :: grad
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 open(newunit=fid,file=TRIM(fchname),status='old',position='append')
 write(fid,'(A,4X,I8)') 'Cartesian Gradient                         R   N=', &
                         natom*3
 write(fid,'(5(1X,ES15.8))') grad
 close(fid)
end subroutine write_grad_into_fch

! read the size of a symmetric matrix from .npy file
subroutine read_sym_mat_size_from_npy(npyname, n)
 implicit none
 integer :: i, fid
 integer, intent(out) :: n
!f2py intent(out) :: n
 character(len=128) :: str
 character(len=240), intent(in) :: npyname
!f2py intent(in) :: npyname

 n = 0; str = ' '
 open(newunit=fid,file=TRIM(npyname),status='old',form='unformatted',&
      access='stream')
 read(fid) str
 close(fid)

 i = INDEX(str, '(')
 read(str(i+1:),*) n
end subroutine read_sym_mat_size_from_npy

! read a symmetric matrix from .npy file
subroutine read_sym_mat_from_npy(npyname, n, a)
 implicit none
 integer :: i, fid
 integer, intent(in) :: n
!f2py intent(in) :: n
 real(kind=8), intent(out) :: a(n,n)
!f2py intent(out) :: a
!f2py depend(n) :: a
 character(len=128) :: str
 character(len=240), intent(in) :: npyname
!f2py intent(in) :: npyname

 str = ' '
 open(newunit=fid,file=TRIM(npyname),status='old',form='unformatted',&
      access='stream')
 read(fid) str
 read(fid,iostat=i) a
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_sym_mat_from_npy: failed to read s&
                   &ymmetric matrix from'
  write(6,'(A)') 'the file '//TRIM(npyname)
  stop
 end if
end subroutine read_sym_mat_from_npy

! find and delete the target .pyc file
subroutine find_and_del_pyc(proname, py_ver)
 implicit none
 integer :: i
 character(len=5) :: ver
 character(len=300) :: pycname, py_ver1
 character(len=240), intent(in) :: proname, py_ver
!f2py intent(in) :: proname, py_ver

 if(py_ver(1:6) == 'Python') then
  py_ver1 = ADJUSTL(TRIM(py_ver(7:)))
 else
  py_ver1 = TRIM(py_ver)
 end if

 i = INDEX(py_ver1, '.')
 ver = py_ver1(1:i-1)
 py_ver1 = ADJUSTL(py_ver1(i+1:))
 i = INDEX(py_ver1, '.')
 ver = TRIM(ver)//py_ver1(1:i-1)

 pycname = '__pycache__/'//TRIM(proname)//'.cpython-'//TRIM(ver)//'.pyc'
 call delete_file(TRIM(pycname))
end subroutine find_and_del_pyc

! modify the file xxx_uno.txt
subroutine modify_uno_out(uno_out, ndb, npair, nopen)
 implicit none
 integer :: k, fid, idx(3)
 integer, intent(in) :: ndb, npair, nopen
!f2py intent(in) :: ndb, npair, nopen
 character(len=240) :: buf
 character(len=240), intent(in) :: uno_out
!f2py intent(in) :: uno_out

 open(newunit=fid,file=TRIM(uno_out),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(1:3) == 'ndb') exit
 end do ! for while

 BACKSPACE(fid)
 k = ndb + 1
 idx = [k, k+nopen+2*npair, nopen]
 write(fid,'(A6,I5)') 'ndb  =', ndb
 write(fid,'(A6,I5)') 'nact =', npair+nopen
 write(fid,'(A6,I5)') 'nact0=', npair
 write(fid,'(A6,3I5)')'idx  =', idx
 close(fid)
end subroutine modify_uno_out

! find the target CASCI root in a specified PySCF CASCI output file
subroutine read_target_root_from_pyscf_out(outname, target_root, found)
 implicit none
 integer :: i, fid
 integer, intent(out) :: target_root
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(out) :: found

 found = .false.; target_root = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:12) == 'target_root=') then
   read(buf(13:),*) target_root
   found = .true.
   exit
  end if
 end do ! for while

 close(fid)
end subroutine read_target_root_from_pyscf_out

! Note: the parameter mem must be provided in unit MB.
subroutine gen_gau_opt_gjf(gjfname, mem, charge, mult, natom, elem, coor, &
                           numfreq)
 implicit none
 integer :: i, fid
 integer,intent(in) :: mem, charge, mult, natom
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=2), intent(in) :: elem(natom)
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: numfreq

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A,I0,A)') '%mem=', min(mem,6000), 'MB'
 write(fid,'(A)') '%nprocshared=1'
 write(fid,'(A)',advance='no') '# opt(nomicro,maxcycles=300)'
 if(numfreq) write(fid,'(A)',advance='no') ' freq=numer'
 write(fid,'(A)') " def2SVPP nosymm external='gau_external'"
 ! TODO: automatically switch to UGBS or mixed basis set when there is any
 ! element which is out of range of def2SVPP.

 write(fid,'(/,A)') 'Using Gaussian as the geometry optimizer'
 write(fid,'(/,I0,1X,I0)') charge, mult
 do i = 1, natom, 1
  write(fid,'(A2,3(1X,F18.8))') elem(i), coor(:,i)
 end do ! for i

 write(fid,'(/)',advance='no')
 close(fid)
end subroutine gen_gau_opt_gjf

! delete specified density matrix in a given .fch(k) file
subroutine del_dm_in_fch(fchname, itype)
 implicit none
 integer :: i, ncoeff, nline, fid, fid1, RENAME
 integer, intent(in) :: itype
!f2py intent(in) :: itype
 character(len=11), parameter :: key(12) = ['Total SCF D', 'Spin SCF De',&
   'Total CI De', 'Spin CI Den', 'Total MP2 D', 'Spin MP2 De',&
   'Total CC De', 'Spin CC Den', 'Total CI Rh', 'Spin CI Rho',&
   'Total 2nd O', 'Spin 2nd Or']
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 if(itype<1 .or. itype>12) then
  write(6,'(/,A,I0)') 'ERROR in subroutine del_dm_in_fch: invalid itype=',itype
  write(6,'(A)') 'Allowed values are 1~12.'
  stop
 end if

 i = INDEX(fchname, '.fch', back=.true.)
 fchname1 = fchname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == key(itype)) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then ! the target density does not exist, simply return
  close(fid)
  close(fid1,status='delete')
  return
 end if

 read(buf(50:),*) ncoeff
 nline = ncoeff/5
 if(ncoeff - 5*nline > 0) nline = nline + 1

 do i = 1, nline, 1
  read(fid,'(A)') buf
 end do ! for i

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5)=='ONIOM' .or. buf(1:5)=='ClPar') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
end subroutine del_dm_in_fch

! print ORCA CASPT2/NEVPT2 selected roots
subroutine prt_orca_mrpt_sel_root(fid, iver, iroot, mult, xmult)
 implicit none
 integer :: i
 integer, intent(in) :: fid, iver, iroot, mult, xmult
 integer, allocatable :: int01(:)

 if(iver > 5) then ! >= ORCA 6.0.0
  if(xmult == mult) then
   write(fid,'(A,I0)') '  SelectedRoots[0]=',iroot
  else
   write(fid,'(A,I0)') '  SelectedRoots[1]=',iroot-1
  end if
 else              ! <= ORCA 5.0.4
  if(xmult == mult) then
   allocate(int01(iroot+3), source=0)
   int01(iroot+1) = 1
   write(fid,'(A)',advance='no') '  SelectedRoots[0]='
   do i = 1, iroot+2, 1
    write(fid,'(I0,A)',advance='no') int01(i),','
   end do ! for i
   write(fid,'(I0)') int01(iroot+3)
  else
   allocate(int01(iroot+2), source=0)
   int01(iroot) = 1
   write(fid,'(A)') '  SelectedRoots[0]=0'
   write(fid,'(A)',advance='no') '  SelectedRoots[1]='
   do i = 1, iroot+1, 1
    write(fid,'(I0,A)',advance='no') int01(i),','
   end do ! for i
   write(fid,'(I0)') int01(iroot+2)
  end if
  deallocate(int01)
 end if
end subroutine prt_orca_mrpt_sel_root

! split characters like '12H' into integer 12 and element 'H'
subroutine split_iatom_elem(iatom_elem, iatom, elem)
 implicit none
 integer :: i, j, k
 integer, intent(out) :: iatom
 character(len=2), intent(out) :: elem
 character(len=11), intent(in) :: iatom_elem

 iatom = 0; elem = '  '; k = LEN_TRIM(iatom_elem)

 do i = 1, k, 1
  j = IACHAR(iatom_elem(i:i))
  if(j<48 .or. j>57) exit
 end do ! for i

 if(i == k+1) then
  write(6,'(/,A)') 'ERROR in subroutine split_iatom_elem: no element is found i&
                   &n '//TRIM(iatom_elem)
  stop
 end if

 read(iatom_elem(1:i-1),*) iatom
 read(iatom_elem(i:k),*) elem
end subroutine split_iatom_elem

! find pNMR isotropic shieldings of target atoms in ORCA pNMR output, and
! calculate the average value
subroutine average_pnmr_shield_in_orca_pnmr_out(pnmr_out, natom, atom_list, &
                                                ave_val)
 implicit none
 integer :: i, iatom, fid
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 integer, intent(in) :: atom_list(natom)
!f2py intent(in) :: atom_list
!f2py depend(natom) :: atom_list
 real(kind=8), intent(out) :: ave_val ! in ppm
!f2py intent(out) :: ave_val
 real(kind=8) :: iso_shield
 real(kind=8), allocatable :: pnmr_shielding(:)
 character(len=2) :: elem
 character(len=11) :: iatom_elem
 character(len=240), intent(in) :: pnmr_out
!f2py intent(in) :: pnmr_out
 character(len=240) :: buf
 logical, allocatable :: found(:) ! size natom

 ave_val = 0d0
 buf = pnmr_out ! buf has length declared as 240
 call require_file_exist(buf)
 open(newunit=fid,file=TRIM(pnmr_out),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:22) == 'Paramagnetic shielding') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine average_pnmr_shield_in_orca_pnmr_out: '&
                   &Paramagnetic shielding'"
  write(6,'(A)') 'not located in file '//TRIM(pnmr_out)
  close(fid)
  stop
 end if

 do i = 1, 5
  read(fid,'(A)') buf
 end do

 allocate(found(natom), pnmr_shielding(natom))
 found = .false.; pnmr_shielding = 0d0

 do while(.true.)
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  read(buf,*) iatom_elem, iso_shield
  call split_iatom_elem(iatom_elem, iatom, elem)
  iatom = iatom + 1
  do i = 1, natom, 1
   if(atom_list(i) == iatom) exit
  end do ! for i
  if(i < natom+1) then
   found(i) = .true.
   pnmr_shielding(i) = iso_shield
  end if
  if(ALL(found .eqv. .true.)) exit
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  if(buf(1:5) == '-----') exit
 end do ! for while

 close(fid)
 if(ANY(found .eqv. .false.)) then
  do i = 1, natom, 1
   if(.not. found(i)) exit
  end do ! for i
  write(6,'(/,A)') 'ERROR in subroutine average_pnmr_shield_in_orca_pnmr_out: t&
                   &he pNMR isotropic'
  write(6,'(A,I0,A)') 'shielding of atom label ',atom_list(i),' is not found in&
                      & file '//TRIM(pnmr_out)
  deallocate(found, pnmr_shielding)
  stop
 end if

 deallocate(found)
 if(natom == 1) then
  ave_val = pnmr_shielding(1)
 else
  ave_val = SUM(pnmr_shielding)/DBLE(natom)
 end if
 deallocate(pnmr_shielding)
end subroutine average_pnmr_shield_in_orca_pnmr_out

! find NMR isotropic shieldings of target atoms in ORCA NMR output, and
! calculate the average value
subroutine average_nmr_shield_in_orca_out(outname, natom, atom_list, ave_val)
 implicit none
 integer :: i, iatom, fid
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 integer, intent(in) :: atom_list(natom)
!f2py intent(in) :: atom_list
!f2py depend(natom) :: atom_list
 real(kind=8), intent(out) :: ave_val
!f2py intent(out) :: ave_val
 real(kind=8) :: iso_shield
 real(kind=8), allocatable :: nmr_shielding(:)
 character(len=2) :: elem
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname
 character(len=240) :: buf
 logical, allocatable :: found(:) ! size natom

 ave_val = 0d0
 buf = outname ! buf has length declared as 240
 call require_file_exist(buf)
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:22) == 'CHEMICAL SHIELDING SUM') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine average_nmr_shield_in_orca_out: 'CHEMIC&
                   &AL SHIELDING SUM'"
  write(6,'(A)') 'not located in file '//TRIM(outname)
  close(fid)
  stop
 end if

 allocate(found(natom), nmr_shielding(natom))
 found = .false.; nmr_shielding = 0d0
 do i = 1, 5
  read(fid,'(A)') buf
 end do ! for i

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  read(buf,*) iatom, elem, iso_shield
  iatom = iatom + 1
  do i = 1, natom, 1
   if(atom_list(i) == iatom) exit
  end do ! for i
  if(i < natom+1) then
   found(i) = .true.
   nmr_shielding(i) = iso_shield
  end if
  if(ALL(found .eqv. .true.)) exit
 end do ! for while

 close(fid)
 if(ANY(found .eqv. .false.)) then
  do i = 1, natom, 1
   if(.not. found(i)) exit
  end do ! for i
  write(6,'(/,A)') 'ERROR in subroutine average_nmr_shield_in_orca_out: the NMR&
                   & isotropic sh-'
  write(6,'(A,I0,A)') 'ielding of atom label ',atom_list(i),' is not found in f&
                      &ile '//TRIM(outname)
  deallocate(found, nmr_shielding)
  stop
 end if

 deallocate(found)
 if(natom == 1) then
  ave_val = nmr_shielding(1)
 else
  ave_val = SUM(nmr_shielding)/DBLE(natom)
 end if
 deallocate(nmr_shielding)
end subroutine average_nmr_shield_in_orca_out

! find NMR isotropic shieldings of target atoms in Gaussian NMR output, and
! calculate the average value
subroutine average_nmr_shield_in_gau_log(logname, natom, atom_list, ave_val)
 implicit none
 integer :: i, j, iatom, fid
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 integer, intent(in) :: atom_list(natom)
!f2py intent(in) :: atom_list
!f2py depend(natom) :: atom_list
 real(kind=8), intent(out) :: ave_val
!f2py intent(out) :: ave_val
 real(kind=8), allocatable :: nmr_shielding(:)
 character(len=240), intent(in) :: logname
!f2py intent(in) :: logname
 character(len=240) :: buf
 logical, allocatable :: found(:) ! size natom

 ave_val = 0d0
 buf = logname ! buf has length declared as 240
 call require_file_exist(buf)
 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:13) == 'SCF GIAO Mag') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine average_nmr_shield_in_gau_log: 'SCF GIA&
                   &O Mag' not"
  write(6,'(A)') 'located in file '//TRIM(logname)
  close(fid)
  stop
 end if

 allocate(found(natom), nmr_shielding(natom))
 found = .false.; nmr_shielding = 0d0

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:11)=='g value of' .or. buf(2:11)=='End of Min') exit
  read(buf,*) iatom
  do i = 1, natom, 1
   if(atom_list(i) == iatom) exit
  end do ! for i
  if(i < natom+1) then
   found(i) = .true.
   j = INDEX(buf, '=')
   read(buf(j+1:),*) nmr_shielding(i)
  end if
  if(ALL(found .eqv. .true.)) exit
  do i = 1, 4
   read(fid,'(A)') buf
  end do ! for i
 end do ! for while

 close(fid)
 if(ANY(found .eqv. .false.)) then
  do i = 1, natom, 1
   if(.not. found(i)) exit
  end do ! for i
  write(6,'(/,A)') 'ERROR in subroutine average_nmr_shield_in_gau_log: the NMR &
                   &isotropic sh-'
  write(6,'(A,I0,A)') 'ielding of atom label ',atom_list(i),' is not found in f&
                      &ile '//TRIM(logname)
  deallocate(found, nmr_shielding)
  stop
 end if

 deallocate(found)
 if(natom == 1) then
  ave_val = nmr_shielding(1)
 else
  ave_val = SUM(nmr_shielding)/DBLE(natom)
 end if
 deallocate(nmr_shielding)
end subroutine average_nmr_shield_in_gau_log

! save the previous K matrix of Cayley transformation
subroutine save_cayley_k_old(nmo, ndiis, k_old, binfile)
 implicit none
 integer :: i, j, nfile, fid, RENAME
 integer, intent(in) :: nmo, ndiis
 real(kind=8), intent(in) :: k_old(nmo,nmo)
 character(len=240), intent(in) :: binfile(2*ndiis-1)
 logical :: alive

 nfile = 2*ndiis - 1
 do i = nfile-1, ndiis+1, -1
  inquire(file=TRIM(binfile(i)), exist=alive)
  if(alive) j = RENAME(TRIM(binfile(i)), TRIM(binfile(i+1)))
 end do ! for i

 open(newunit=fid,file=TRIM(binfile(ndiis+1)),status='new',form='unformatted')
 write(unit=fid) ((k_old(j,i),j=i+1,nmo,1),i=1,nmo-1,1)
 close(fid)
end subroutine save_cayley_k_old

! save the difference of current and previous K matrices of Cayley transformation
subroutine save_cayley_k_diff(nmo,ndiis,ijmap, binfile, k_old, cayley_k, k_diis)
 implicit none
 integer :: i, j, n, npair, nfile, fid, RENAME
 integer, intent(in) :: nmo, ndiis
 integer, intent(in) :: ijmap(2,nmo*(nmo-1)/2)
 real(kind=8), intent(in) :: k_old(nmo,nmo), cayley_k(nmo,nmo)
 real(kind=8), intent(inout) :: k_diis(ndiis+1,ndiis+1)
 real(kind=8), allocatable :: k_diff(:), old_diff(:)
 character(len=240), intent(in) :: binfile(2*ndiis-1)
 logical :: alive

 npair = nmo*(nmo-1)/2; nfile = 2*ndiis - 1

 ! renaming files: N-5->N-6, N-4->N-5, N-3->N-4, N-2->N-3, N-1->N-2, N->N-1
 do i = ndiis-1, 1, -1
  inquire(file=TRIM(binfile(i)), exist=alive)
  if(alive) j = RENAME(TRIM(binfile(i)), TRIM(binfile(i+1)))
 end do ! for i

 do i = 1, ndiis-1, 1
  do j = 1, ndiis-1, 1
   k_diis(j,i) = k_diis(j+1,i+1)
  end do ! for j
 end do ! for i

 allocate(k_diff(npair))
!$omp parallel do schedule(dynamic) default(shared) private(i,j,n)
 do n = 1, npair, 1
  i = ijmap(1,n); j = ijmap(2,n)
  k_diff(n) = cayley_k(j,i) - k_old(j,i)
 end do ! for i
!$omp end parallel do

 do i = 1, ndiis-1, 1
  inquire(file=TRIM(binfile(ndiis+1-i)),exist=alive)
  if(.not. alive) cycle
  allocate(old_diff(npair), source=0d0)
  open(newunit=fid,file=TRIM(binfile(ndiis+1-i)),status='old',form='unformatted')
  read(unit=fid) old_diff
  close(fid)
  k_diis(i,ndiis) = DOT_PRODUCT(k_diff, old_diff)
  deallocate(old_diff)
 end do ! for k

 k_diis(ndiis,ndiis) = DOT_PRODUCT(k_diff, k_diff)
 k_diis(ndiis,1:ndiis-1) = k_diis(1:ndiis-1,ndiis)

 open(newunit=fid,file=TRIM(binfile(1)),status='new',form='unformatted')
 write(unit=fid) k_diff
 deallocate(k_diff)
 close(fid)
end subroutine save_cayley_k_diff

