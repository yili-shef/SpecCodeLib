! Difference between dealiasing and symmetrize is at the last statement.
subroutine dealiasing(c,k2,lx1,ly,lz)
  use mconstant
  implicit none

  integer, intent(in)  :: lx1,ly,lz
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: c
  real(sp),    dimension(lx1,ly,lz), intent(in)    :: k2


  integer :: iy,iz

! zero ky=-ly/2 modes
  c(:,ly/2+1,:) = (0.0_SP,0.0_SP)
! zero kz=-lz/2 modes
  c(:,:,lz/2+1) = (0.0_SP,0.0_SP)
! zero kx=lx modes
  c(lx1,:,:)=(0.0_SP,0.0_SP)
! zero kx=ky=kz=0 mode
  c(1,1,1) = (0.0_sp,0.0_sp)

  ! imposing the symmetry constraint on the kx=0 plane, excluding the line kx=kz=0.
  do iy=1,ly
  do iz=lz/2+2,lz
  ! using MOD to handle iy=1; for iy>1, MOD(...)+1 = ly+2-iy
  c(1,iy,iz) = CONJG(c(1,MOD(ly+1-iy,ly)+1,lz+2-iz))
  enddo
  enddo

  ! imposing the constraint on the line kx=kz=0.
  do iy=ly/2+2,ly
  c(1,iy,1) = CONJG(c(1,ly+2-iy,1))
  enddo

  where(k2 .ge. real(ly*ly,sp)/9) c = (0.0,0.0)

end subroutine dealiasing
