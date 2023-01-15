! Modified from "How to read a selected wave function from a wfc#k.dat file 
! inside the `prefix.save` directrory"
!   https://gitlab.com/QEF/q-e/-/snippets/1869202
!   Authored by Pietro Delugas
!
! Modified by TQ
!
module read_qe_wfc
	!
	! this module read the wfc from QE (>= v6.3 ?) files, i.e., wfc**.dat
	! common usage:
	!    call print_header(filename)
	!    call read_single_wfc(ibnd, filename, evc)
	!    call read_upto_n_wfc(ibnd, filename, evc)
	!    call read_all_wfc(filename, evc)
	!    ...
	! see main.f90 for an example
	!
	use iso_fortran_env, only: dp=> real64
	!
	implicit none
	!
	integer   :: iunit = 97
	integer   :: ik, nbnd, ispin, npol, ngw, igwx
	real(dp)  :: xk(3)
	logical   :: gamma_only
	!
	!
	private   :: read_header
	!
	!
	contains
	!
	!
	subroutine read_header(filename)
		!
		implicit none
		!
		character(len=*), intent(in)   :: filename
		!
		! working variable
		integer   :: tmp_i, tmp_unit
		real(dp)  :: scalef
		real(dp)  :: b1(3), b2(3), b3(3)
		logical   :: file_exist
		!
		! check if file exists or not
		inquire(file = trim(filename), number = tmp_unit, exist = file_exist)
		if ( .not. file_exist ) then
			write(*,*) "Error: the file ",trim(filename), " does not exist!"
			stop 1
		endif
		if (tmp_unit .ne. -1) close(tmp_unit)
		!
		! open file
		open(unit = iunit, file = trim(filename), form = 'unformatted', status = 'old')
		!
		! read header
		read(iunit) ik, xk, ispin, gamma_only, scalef
		read(iunit) ngw, igwx, npol, nbnd
		read(iunit) b1, b2, b3
		!
		! avoid reading miller indices of G vectors below E_cut for this kpoint
		! if needed allocate  an integer array of dims (1:3,1:igwx)
		! NOTE: the above comment is copied from https://gitlab.com/QEF/q-e/-/snippets/1869202
		! not sure what it means
		!
		read (iunit) tmp_i
		!
	end subroutine read_header
	!
	!
	subroutine read_header_only(filename)
		!
		implicit none
		!
		character(len=*), intent(in)   :: filename
		!
		call read_header(filename)
		close(iunit)
		!
	end subroutine read_header_only
	!
	!
	subroutine print_header(filename)
		!
		implicit none
		!
		character(len=*), intent(in)   :: filename
		!
		call read_header_only(filename)
		!
		write(*,*) "filename   = ", trim(filename)
		write(*,*) "ik         = ", ik
		write(*,*) "nbnd       = ", nbnd
		write(*,*) "ispin      = ", ispin
		write(*,*) "npol       = ", npol
		write(*,*) "ngw        = ", ngw
		write(*,*) "igwx       = ", igwx
		write(*,*) "xk         = ", xk
		write(*,*) "gamma_only = ", gamma_only
		!
		close(iunit)
		!
	end subroutine print_header
	!
	!
	subroutine read_single_wfc(ibnd, filename, evc)
		!
		implicit none
		!
		integer                 , intent(in)   :: ibnd
		character(len=*)        , intent(in)   :: filename
		complex(dp), allocatable, intent(out)  :: evc(:)
		!
		! internal variable
		integer                  :: t1
		complex(dp), allocatable :: tmp_r(:)
		!
		call read_header(filename)
		!
		if (ibnd > nbnd) then
			write(*,'("Error, ibnd(",I5,") > nbnd(",I5,")")') ibnd, nbnd
			stop 2
		endif
		!
		if (allocated(evc)) deallocate(evc)
		allocate (evc(npol*igwx))
		allocate (tmp_r(npol*igwx))
		!
		do t1 = 1,ibnd-1
			! this line is originally read in a real number
			read(iunit) tmp_r(1:npol*igwx)
		enddo
		!
		read(iunit) evc(1:npol*igwx)
		!
		deallocate(tmp_r)
		close(iunit)
		!
	end subroutine read_single_wfc
	!
	!
	subroutine read_upto_n_wfc(ibnd, filename, evc)
		!
		implicit none
		!
		integer                 , intent(in)   :: ibnd
		character(len=*)        , intent(in)   :: filename
		complex(dp), allocatable, intent(out)  :: evc(:,:)
		!
		! internal variable
		integer   :: t1, end_bnd
		!
		call read_header(filename)
		!
		if (ibnd > nbnd) then
			write(*,'("Warning: ibnd(",I5,") > nbnd(",I5,"), only read up to nbnd wfcs!")') ibnd, nbnd
			end_bnd = nbnd
		else
			end_bnd = ibnd
		endif
		!
		if (allocated(evc)) deallocate(evc)
		allocate (evc(npol*igwx, end_bnd))
		!
		do t1 = 1, end_bnd
			read(iunit) evc(1:npol*igwx,t1)
		enddo
		!
		close(iunit)
		!
	end subroutine read_upto_n_wfc
	!
	!
	subroutine read_all_wfc(filename, evc)
		!
		implicit none
		!
		character(len=*)        , intent(in)   :: filename
		complex(dp), allocatable, intent(out)  :: evc(:,:)
		!
		!call read_header_only(filename)
		! because the address of nbnd is passed, no need to worry it has no valued yet
		call read_upto_n_wfc(nbnd, filename, evc)
		!
	end subroutine read_all_wfc
	!
	!
	function overlap(evc1,evc2)
		!
		implicit none
		!
		complex(dp), dimension(:) :: evc1, evc2
		complex(dp)               :: overlap
		!
		if (gamma_only) then
			overlap = 2.d0*dble(dot_product(evc1,evc2)) - evc1(1)*evc2(1)
		else
			overlap = dot_product(evc1,evc2)
		endif
		!
	end function overlap
	!
end module read_qe_wfc
