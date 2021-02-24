!!!! This subroutine is needed for the helical wave components, 
!!!! which are anti-hermitian in wave space.
subroutine antihermit(c,lx1,ly,lz)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly,lz
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: c

  integer :: iy,iz

  ! imposing the anti-symmetry constraint on the kx=0 plane, excluding the line kx=kz=0.
  do iy=1,ly
  do iz=lz/2+2,lz
  ! using MOD to handle iy=1; for iy>1, MOD(...)+1 = ly+2-iy
  c(1,iy,iz) = - CONJG(c(1,MOD(ly+1-iy,ly)+1,lz+2-iz))
  enddo
  enddo

  ! imposing the constraint on the line kx=kz=0.
  do iy=ly/2+2,ly
  c(1,iy,1) = - CONJG(c(1,ly+2-iy,1))
  enddo
 
  return
end subroutine antihermit


