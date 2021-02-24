subroutine definefilter (mm, nx, delta, k2, g)
  use mconstant
  implicit none

  integer, intent(in) :: mm, nx
  real(sp), intent(in) :: delta
  real(sp), dimension(nx/2+1,nx,nx), intent(in) :: k2
  real(sp), dimension(nx/2+1,nx,nx), intent(out) :: g


  if ( mm .eq. 1 ) then
    where (k2 .ge. (pi/delta) * (pi/delta))
            g = 0.
    elsewhere
            g = 1.
    endwhere
  else if ( mm .eq. 0 ) then 
    g=exp(-k2*delta**2/24.)
  else
    write(*,*) 'Unknown filter type! Stopped in definefilter'
    stop
  end if
  
end subroutine definefilter
