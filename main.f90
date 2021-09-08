program test
	!
	use read_wfc,        only: read_single_wfc
	use iso_fortran_env, only: dp=> real64
	!
	implicit none
	!
	integer                   :: ibnd=1
	complex(dp), allocatable  :: evc1(:), evc2(:)
	integer                   :: ik, nbnd, ispin, npol, ngw, igwx
	real(dp)                  :: xk(3)
	logical                   :: gamma_only
	character(4)              :: kpoint
	!
	write(kpoint,'(I4)') 1
	call read_single_wfc(ibnd, "wfc"//trim(adjustl(kpoint))//".dat", evc1, ik, nbnd, ispin, npol, ngw, igwx, xk, gamma_only)
	!
	write(*,*) "ik         = ", ik
	write(*,*) "nbnd       = ", nbnd
	write(*,*) "ispin      = ", ispin
	write(*,*) "npol       = ", npol
	write(*,*) "ngw        = ", ngw
	write(*,*) "igwx       = ", igwx
	write(*,*) "xk         = ", xk
	write(*,*) "gamma_only = ", gamma_only
	write(*,*) "norm 1     = ", 2.d0*dot_product(evc1(2:),evc1(2:)) + conjg(evc1(1)) * evc1(1)
	write(*,*) "norm 1     = ", dot_product(evc1,evc1)
	write(*,*) ""
	!
	ibnd = 2
	call read_single_wfc(ibnd, "wfc1.dat", evc2, ik, nbnd, ispin, npol, ngw, igwx, xk, gamma_only)
	write(*,*) "ik         = ", ik
	write(*,*) "nbnd       = ", nbnd
	write(*,*) "ispin      = ", ispin
	write(*,*) "npol       = ", npol
	write(*,*) "ngw        = ", ngw
	write(*,*) "igwx       = ", igwx
	write(*,*) "xk         = ", xk
	write(*,*) "gamma_only = ", gamma_only
	write(*,*) "norm 2     = ", 2.d0*dot_product(evc2(2:),evc2(2:)) + conjg(evc2(1)) * evc2(1)
	write(*,*) "norm 2     = ", dot_product(evc2,evc2)
	write(*,*) ""
	write(*,*) "overlap    = ", 2.d0*dble(dot_product(evc2(2:),evc1(2:))) + conjg(evc2(1)) * evc1(1)
	write(*,*) "overlap    = ", dot_product(evc2,evc1)
end program test
