subroutine setzero (c, k2, lx1, ly, lz)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly,lz
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: c
  real(sp),    dimension(lx1,ly,lz), intent(in) :: k2

! zero ky=-ly/2 modes
  c(:,ly/2+1,:) = (0.0_SP,0.0_SP)
! zero kz=-lz/2 modes
  c(:,:,lz/2+1) = (0.0_SP,0.0_SP)
! zero kx=lx modes
  c(lx1,:,:)=(0.0_SP,0.0_SP)
! zero kx=ky=kz=0 mode
  c(1,1,1) = (0.0_sp,0.0_sp)

  where (k2 .ge. (lx1-1)**2) c=cmplx(0._sp,0._sp)

end subroutine setzero    
