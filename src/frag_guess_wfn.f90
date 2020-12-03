! written by jxzou at 20201121: perform fragment-guess wavefunction calculations

! The fragment-guess wavefunction calculation in Gaussian has the following
! shortcomings:
! (1) cannot ensure the wavefunction stability of each fragment. This often
!     makes the calculations of transition-metal molecules take more SCF cycles,
!     since the initial guess is not as good as we expect.
! (2) do not support mixed methods, e.g. cannot use RHF/RDFT for some fragments
!     and UHF/UDFT for other fragments
! This utility aims to overcome these shortcomings.

module frag_info
 implicit none
 integer :: nfrag   ! number of fragments

 type :: frag
  integer :: charge, mult, natom
  integer, allocatable :: atm_map(:) ! map to parent system
  real(kind=8), allocatable :: coor(:,:)
  character(len=2), allocatable :: elem(:)
 end type frag

 type(frag), allocatable :: frags(:)
end module frag_info

program main
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=3) :: gver
 character(len=240) :: gjfname
 character(len=64), parameter :: error_warn = &
  'ERROR in subroutine frag_comb_wfn: wrong command line arguments!'

 i = iargc()
 if(i /= 2) then
  write(iout,'(/,A)') error_warn
  write(iout,'(A,/)') "Example: frag_comb_wfn g16 water_dimer.gjf &"
  stop
 end if

 gver = ' '; gjfname = ' '
 call getarg(1, gver)
 if(.not. (gver=='g03' .or. gver=='g09' .or. gver=='g16')) then
  write(iout,'(/,A)') error_warn
  write(iout,'(A)') 'Only g03/g09/g16 is allowed.'
  stop
 end if

 call getarg(2, gjfname)
 if(index(gjfname,'.gjf') == 0) then
  write(iout,'(/,A)') error_warn
  write(iout,'(A)') "'.gjf' suffix not found in filename "//TRIM(gjfname)
  stop
 end if

 call frag_guess_wfn(gver, gjfname)
 stop
end program main 

! perform fragment-guess wavefunction calculations, one by one fragment
subroutine frag_guess_wfn(gver, gjfname)
 use frag_info, only: nfrag, frags
 implicit none
 integer :: i, j, fid
 integer :: charge, mult, natom
 integer, parameter :: iout = 6
 integer, allocatable :: cm(:)
 real(kind=8), allocatable :: coor(:,:)
 character(len=2), allocatable :: elem(:)
 character(len=3), intent(in) :: gver
 character(len=8) :: method
 character(len=17) :: basis ! 17 is for longest 6-311++G(3df,2pd)
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
 character(len=1200) :: longbuf

 buf = ' '; longbuf = ' '
 call read_natom_from_gjf(gjfname, natom) ! read the total number of atoms

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:1) == '#') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine frag_guess_wfn: incomplete file '//TRIM(gjfname)
  write(iout,'(A)') "Failed to locate the Route Section '#' in .gjf file."
  close(fid)
  stop
 end if

 longbuf = TRIM(buf)
 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  longbuf = TRIM(longbuf)//TRIM(buf)
 end do ! for while

 call lower(longbuf)
 call read_nfrag_from_buf(longbuf, nfrag)
 call read_method_and_basis_from_buf(longbuf, method, basis)
 write(iout,'(A,I0)') 'nfrag=', nfrag
 write(iout,'(A)') 'method='//TRIM(method)//', basis='//TRIM(basis)

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
 end do ! for while
 allocate(cm(2+2*nfrag),source=0)
 read(fid,*,iostat=i) cm

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine frag_guess_wfn: incomplete charges and&
                   & multiplicities in file '//TRIM(gjfname)
  
  stop
 end if

 charge = cm(1); mult = cm(2)
 allocate(frags(nfrag))
 do i = 1, nfrag, 1
  frags(i)%charge = cm(2*i+1); frags(i)%mult = cm(2*i+2)
 end do ! for i
 deallocate(cm)

 allocate(coor(3,natom), source=0d0)
 allocate(elem(natom), source='  ')
 do i = 1, natom, 1
  read(fid,'(A)',iostat=j) buf
  if(j/=0 .or. LEN_TRIM(buf)==0) exit

 end do ! for i

 deallocate(frags)
 return
end subroutine frag_guess_wfn

! read the method and basis set from a string
! Note: please call subroutine lower before calling this subroutine,
!       in order to transform all letters to lower case
subroutine read_method_and_basis_from_buf(buf, method, basis)
 implicit none
 integer :: i, j
 integer, parameter :: iout = 6
 character(len=1200), intent(in) :: buf
 character(len=8), intent(out) :: method
 character(len=17), intent(out) :: basis

 j = index(buf, '/')
 if(j == 0) then
  write(iout,'(A)') "ERROR in subroutine read_method_and_basis_from_buf: no&
                   & '/' symbol found in '#' line. Fail to identify the method."
  stop
 end if

 i = index(buf(1:j-1), ' ', back=.true.)
 method = ' '
 method = buf(i+1:j-1)

 if(method(1:1) == 'r') then
  write(iout,'(A)') 'ERROR in subroutine read_method_and_basis_from_buf:&
                   & restricted method '//TRIM(method)//' detected.'
  write(iout,'(A)') 'Only unrestricted method is supported currently.'
  stop
 end if

 i = index(buf(j+1:), ' ')
 basis = ' '
 basis = buf(j+1:j+i-1)
 return
end subroutine read_method_and_basis_from_buf

! read nfrag (the number of fragments) from a string
! Note: please call subroutine lower before calling this subroutine,
!       in order to transform all letters to lower case
subroutine read_nfrag_from_buf(buf, nfrag)
 implicit none
 integer :: i, j, k
 integer, intent(out) :: nfrag
 integer, parameter :: iout = 6
 character(len=1200), intent(in) :: buf

 i = index(buf, 'guess')
 if(i == 0) then
  write(iout,'(A)') "ERROR in subroutine read_nfrag_from_buf: keyword 'guess'&
                   & not found in buf."
  write(iout,'(A)') 'buf = '//TRIM(buf)
  stop
 end if

 i = i + 4
 nfrag = 0
 j = i+1 + index(buf(i+2:), '=')

 select case(buf(i+1:i+1))
 case('=')
  k = i+1 + index(buf(i+2:), ' ')
 case('(')
  k = i+1 + index(buf(i+2:), ')')
 case default
  write(iout,'(A)') 'ERROR in subroutine read_nfrag_from_buf: syntax error.'
  write(iout,'(A)') 'buf = '//TRIM(buf)
  stop
 end select
 read(buf(j+1:k-1),*) nfrag

 if(nfrag == 1) then
  write(iout,'(A)') 'ERROR in subroutine read_nfrag_from_buf: nfrag = 1.'
  write(iout,'(A)') 'Only one fragment. No need to use fragment-guess.'
  stop
 else if(nfrag < 1) then
  write(iout,'(A,I0)') 'ERROR in subroutine read_nfrag_from_buf: invalid&
                      & nfrag = ', nfrag
  stop
 end if

 return
end subroutine read_nfrag_from_buf

! turn letters in buf into lower case, except those in symbol ''
subroutine lower(buf)
 implicit none
 integer :: i, j, k1, k2
 character(len=240), intent(inout) :: buf

 k1 = index(buf,"'")
 k2 = index(buf,"'",back=.true.)

 do i = 1, LEN_TRIM(buf), 1
  if(i>k1 .and. i<k2) cycle

  j = IACHAR(buf(i:i))
  if(j>64 .and. j<91) buf(i:i) = ACHAR(j+32)
 end do ! for i

 return
end subroutine lower

