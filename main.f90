program test
	!
	use read_qe_wfc,        only: print_header, read_single_wfc, read_upto_n_wfc, read_all_wfc, &
	                              overlap, ik, nbnd, ispin, npol, ngw, igwx, xk, gamma_only
	use iso_fortran_env,    only: dp=> real64
	!
	implicit none
	!
	integer                   :: ibnd=40
	complex(dp), allocatable  :: evc0(:), evc1(:,:)
	character(4)              :: kpoint
	character(100)            :: filename
	!
	integer :: t1
	write(kpoint,'(I4)') 0
	filename = "wfc"//trim(adjustl(kpoint))//".dat"
	!
	call print_header(filename)
	!call read_single_wfc(ibnd, filename, evc0)
	call read_upto_n_wfc(ibnd, filename, evc1)
	!
	do t1 = 1, ibnd
		write(*,*) "norm 1,",t1," =", overlap(evc1(:,1),evc1(:,t1))
	enddo
	!
end program test
