! written by jxzou at 20231109: partition a large molecule into monomers
!  and generate initial MOs as well as amplitudes from fragement calculations

program frag_qc
 implicit none
 integer :: i
 character(len=24) :: data_string
 character(len=240) :: gjfname = ' '

 call getarg(1, gjfname)
 call require_file_exist(gjfname)

 call fdate(data_string)
 write(6,'(A)') TRIM(data_string)


 call fdate(data_string)
 write(6,'(/,A)') TRIM(data_string)
end program frag_qc
