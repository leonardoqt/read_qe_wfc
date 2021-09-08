! Modified from "How to read a selected wave function from a wfc#k.dat file 
! inside the `prefix.save` directrory"
!   https://gitlab.com/QEF/q-e/-/snippets/1869202
!   Authored by Pietro Delugas
!
! Modified by TQ
!
module read_wfc
	!
	use iso_fortran_env, only: dp=> real64
	!
	implicit none
	!
	contains
	!
	subroutine read_single_wfc(ibnd, filename, evc, ik, nbnd, ispin, npol, ngw, igwx, xk, gamma_only)
		!
		implicit none
		!
		integer                 , intent(in)   :: ibnd
		character(len=*)        , intent(in)   :: filename
		complex(dp), allocatable, intent(out)  :: evc(:)
		integer                 , intent(out)  :: ik, nbnd, ispin, npol, ngw, igwx
		real(dp)                , intent(out)  :: xk(3)
		logical                                :: gamma_only
		!
		! internal variable
		integer   :: t1, t2, t3, tmp_i, iunit=97
		real(dp)  :: scalef
		real(dp)  :: b1(3), b2(3), b3(3)
		complex(dp), allocatable :: tmp_r(:)
		!
		! open file
		write(*,*) "reading file", filename
		open(unit = iunit, file = trim(filename), form = 'unformatted', status = 'old')
		!
		! read header
		read(iunit) ik, xk, ispin, gamma_only, scalef
		read(iunit) ngw, igwx, npol, nbnd
		read(iunit) b1, b2, b3
		!
		! avoid reading miller indices of G vectors below E_cut for this kpoint
		! if needed allocate  an integer array of dims (1:3,1:igwx)
		!
		read (iunit) tmp_i
		!
		! read wfc
		!
		if (ibnd > nbnd) then
			write(*,*) "Error, ibnd > nbnd, (",nbnd,")"
			stop
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
end module read_wfc
