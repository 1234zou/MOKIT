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
  integer :: charge = 0
  integer :: mult = 0
  integer :: natom = 0
  integer, allocatable :: atm_map(:) ! map to parent system
  real(kind=8) :: e = 0.0d0          ! electronic energy
  real(kind=8) :: ssquare = 0d0      ! S(S+1)
  real(kind=8), allocatable :: coor(:,:)
  character(len=2), allocatable :: elem(:)
  character(len=240) :: fname = ' '
  logical :: u = .true. ! True/False for unrestricted/restricted SCF
 end type frag

 type(frag), allocatable :: frags(:)
end module frag_info

program main
 use print_id, only: iout
 implicit none
 integer :: i
 character(len=240) :: gau_path, gjfname
 character(len=64), parameter :: error_warn = &
  'ERROR in subroutine frag_comb_wfn: wrong command line arguments!'
 logical :: alive

 i = iargc()
 if(i /= 1) then
  write(iout,'(/,A)') error_warn
  write(iout,'(A,/)') "Example: frag_comb_wfn water_dimer.gjf &"
  stop
 end if

 gjfname = ' '
 call getarg(1, gjfname)
 if(index(gjfname,'.gjf') == 0) then
  write(iout,'(/,A)') error_warn
  write(iout,'(A)') "'.gjf' suffix not found in filename "//TRIM(gjfname)
  stop
 end if

 inquire(file=TRIM(gjfname), exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine frag_comb_wfn: file does not exist.'
  write(iout,'(A)') 'Filename = '//TRIM(gjfname)
  stop
 end if

 call get_gau_path(gau_path)
 call frag_guess_wfn(gau_path, gjfname)
 stop
end program main 

! perform fragment-guess wavefunction calculations, one by one fragment
subroutine frag_guess_wfn(gau_path, gjfname)
 use print_id, only: iout
 use frag_info, only: nfrag, frags
 implicit none
 integer :: i, j, k, m, fid, ifrag, iatom
 integer :: charge, mult, natom
 integer :: mem, nproc
 integer, allocatable :: cm(:)
 real(kind=8), allocatable :: coor(:,:), tmp_coor(:,:)
 character(len=2), allocatable :: elem(:)
 character(len=8) :: method
 character(len=17) :: basis ! 17 is for longest 6-311++G(3df,2pd)
 character(len=240) :: buf
 character(len=240), intent(in) :: gau_path, gjfname
 character(len=1200) :: longbuf

 buf = ' '; longbuf = ' '
 call read_mem_and_nproc_from_gjf(gjfname, mem, nproc)
 write(iout,'(A,I0,A)') '%mem=', mem, 'MB'
 write(iout,'(A,I0)') '%nprocshared=', nproc
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
 if(nfrag > 999) then
  write(iout,'(A)') 'ERROR in subroutine frag_guess_wfn: nfrag>999, too large.'
  stop
 end if

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
 allocate(elem(natom))
 elem = '  '

 do i = 1, natom, 1
  read(fid,'(A)',iostat=j) buf
  if(j/=0 .or. LEN_TRIM(buf)==0) exit

  j = index(buf,'(')
  k = index(buf,'=')
  m = index(buf,')')
  if(j*k*m == 0) then
   write(iout,'(A)') 'ERROR in subroutine frag_guess_wfn: wrong format in file '//TRIM(gjfname)
   stop
  end if
  read(buf(1:j-1),*) elem(i)
  read(buf(k+1:m-1),*) ifrag

  if(ifrag > nfrag) then
   write(iout,'(A)') 'ERROR in subroutine frag_guess_wfn: the number of fragments&
                    & detected in Cartesian coordinates is'
   write(iout,'(A)') 'larger than that found in charges and spin multiplicities.'
   write(iout,'(2(A,I0))') 'ifrag = ', ifrag, ', nfrag=', nfrag
   stop
  end if

  read(buf(m+1:),*) coor(1:3,i)
  frags(ifrag)%natom = frags(ifrag)%natom + 1
  iatom = frags(ifrag)%natom

  if(iatom == 1) then
   allocate(frags(ifrag)%coor(3,1))
   frags(ifrag)%coor(:,1) = coor(:,i)

   frags(ifrag)%elem = [elem(i)]
   frags(ifrag)%atm_map = [1]
  else
   tmp_coor = frags(ifrag)%coor
   deallocate(frags(ifrag)%coor)
   allocate(frags(ifrag)%coor(3,iatom))
   frags(ifrag)%coor(:,1:iatom-1) = tmp_coor
   deallocate(tmp_coor)
   frags(ifrag)%coor(:,iatom) = coor(:,i)

   frags(ifrag)%elem = [frags(ifrag)%elem, elem(i)]
   frags(ifrag)%atm_map =[frags(ifrag)%atm_map, [i]]
  end if

  j = LEN_TRIM(buf)
  if(buf(j:j)=='r' .or. buf(j:j)=='R') frags(ifrag)%u = .false.
 end do ! for i

 close(fid)

 do i = 1, nfrag, 1
  write(frags(i)%fname,'(I3.3,A)') i, 'frg-'//TRIM(gjfname)
  call gen_frag_gjf(frags(i), method, basis, mem, nproc)
 end do ! for i

 do i = 1, nfrag, 1
  call perform_scf_and_read_e(gau_path, frags(i)%fname, frags(i)%e, frags(i)%ssquare)
  write(iout,'(A,F18.9)') 'frags(i)%e = ', frags(i)%e
 end do ! for i

 deallocate(frags)
 return
end subroutine frag_guess_wfn

! read memory and nprocshared from a given .gjf file
subroutine read_mem_and_nproc_from_gjf(gjfname, mem, np)
 use print_id, only: iout
 implicit none
 integer :: i, j, k, fid
 integer, intent(out) :: mem, np
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 ! default settings
 mem = 1000 ! 1000 MB
 np = 1     ! 1 core

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call lower(buf)
  if(buf(1:1) == '#') exit

  j = LEN_TRIM(buf)
  k = index(buf,'=')
  if(buf(1:4) == '%mem') then
   read(buf(k+1:j-2),*) mem
   select case(buf(j-1:j))
   case('gb')
    mem = 1000*mem ! you like 1024? I prefer 1000
   case('mb')
   case('gw')
    mem = 1000*8*mem
   case('mw')
    mem = 8*mem
   case default
    write(iout,'(A)') 'ERROR in subroutine read_mem_and_nproc_from_gjf: memory&
                     & unit cannot be recognized.'
    write(iout,'(A)') 'unit = '//TRIM(buf(j-1:j))
    stop
   end select
  else if(buf(1:6) == '%nproc') then
   read(buf(k+1:),*) np
  end if
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_mem_and_nproc_from_gjf: incomplete&
                   & file '//TRIM(gjfname)
  stop
 end if

 close(fid)
 return
end subroutine read_mem_and_nproc_from_gjf

! read the method and basis set from a string
! Note: please call subroutine lower before calling this subroutine,
!       in order to transform all letters to lower case
subroutine read_method_and_basis_from_buf(buf, method, basis)
 use print_id, only: iout
 implicit none
 integer :: i, j
 character(len=1200), intent(in) :: buf
 character(len=8), intent(out) :: method
 character(len=21), intent(out) :: basis

 j = index(buf, '/')
 if(j == 0) then
  write(iout,'(A)') "ERROR in subroutine read_method_and_basis_from_buf: no&
                   & '/' symbol found in '#' line. Fail to identify the method."
  stop
 end if

 i = index(buf(1:j-1), ' ', back=.true.)
 method = ' '
 method = buf(i+1:j-1)

 i = index(buf(j+1:), ' ')
 basis = ' '
 basis = buf(j+1:j+i-1)
 return
end subroutine read_method_and_basis_from_buf

! read nfrag (the number of fragments) from a string
! Note: please call subroutine lower before calling this subroutine,
!       in order to transform all letters to lower case
subroutine read_nfrag_from_buf(buf, nfrag)
 use print_id, only: iout
 implicit none
 integer :: i, j, k
 integer, intent(out) :: nfrag
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

! generate a .gjf file from a type frag
subroutine gen_frag_gjf(frag0, method, basis, mem, nproc)
 use frag_info, only: frag
 use print_id, only: iout
 implicit none
 integer :: i, fid
 integer, intent(in) :: mem, nproc
 character(len=8) :: method0
 character(len=8), intent(in) :: method
 character(len=21), intent(in) :: basis
 type(frag), intent(in) :: frag0

 if(index(method,'ro3lyp')/=0 .or. index(method,'roo3lyp')/=0) then
  write(iout,'(A)') 'ERROR in subroutine gen_frag_gjf: R- and RO-type O3LYP&
                   & functional is not supported.'
  write(iout,'(A)') 'Because this functional name is confusing.'
  stop
 end if

 method0 = ' '
 if(method(1:1)=='u' .or. (method(1:1)=='r' .and. method(2:2)/='o')) then
  method0 = method(2:)
 else if(method(1:2) == 'ro') then
  method0 = method(3:)
 else
  method0 = method
 end if

 if(frag0%u) then
  method0 = 'u'//method0
 else
  if(frag0%mult == 1) then
   method0 = 'r'//method0
  else
   method0 = 'ro'//method0
  end if
 end if

 open(newunit=fid,file=TRIM(frag0%fname),status='replace')

 i = index(frag0%fname, '.gjf', back=.true.)
 write(fid,'(A)') '%chk='//frag0%fname(1:i-1)//'.chk'
 write(fid,'(A,I0,A)') '%mem=', mem, 'MB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(A)') '#p '//TRIM(method0)//'/'//TRIM(basis)//&
                               ' nosymm int=nobasistransform'

 write(fid,'(/,A,/)') 'auto-generated file by frag_guess_wfn of MOKIT'
 write(fid,'(I0,1X,I0)') frag0%charge, frag0%mult

 do i = 1, frag0%natom, 1
  write(fid,'(A2,3X,3F18.10)') frag0%elem(i), frag0%coor(1:3,i)
 end do ! for i

 write(fid,'(/)',advance='no')
 close(fid)
 return
end subroutine gen_frag_gjf

