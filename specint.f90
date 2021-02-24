subroutine specint (ui,upi,xp,yp,zp,kx,ky,kz,lx1,ly,lz)
  use mconstant
  implicit none

  integer, intent(in) :: lx1,ly,lz
  real(sp), intent(in) :: xp,yp,zp
  complex(sp), dimension(lx1,ly,lz), intent(in) :: ui
  real(sp), dimension(lx1,ly,lz), intent(in) :: kx,ky,kz
  real(sp), intent(out) :: upi

  real(sp), dimension(lx1,ly,lz) :: ut

  ut=2.*real( ui * exp( eye*(xp*kx+yp*ky+zp*kz) ) )
  ut(1,:,:)=ut(1,:,:)/2.
  upi=sum(ut)
  ! use the sum function directly
  ! the roundoff error could potentially be big, since we are summing a lot of
  ! terms together.

  ! this method makes use of the knowledge that the magnitude of the modes at
  ! the same wavenumber is about the same to reduce the possibility of adding a
  ! small number to a large one. test shows that this way gives the same results
  ! as directly using sum function.
  ! notice that this method requires k2 as a dummy variable.
!  upi=0.
!  do i=1,lx1-1
!    tr=sum(ut,mask=(abs(sqrt(k2)-i).lt.0.5))
!    upi=upi+tr
!  end do


  ! adding the modes one by one. the order are from large wavenumbers to small
  ! ones, but in terms of components. this gives a smaller value than previous two
  ! methods for ux at ifile=500 and xp=yp=zp=2. another way of adding, starting 
  ! from small wavenumbers, ending at large wavenumbers, gives an even smaller result.
!   upi=0.
!   do ix=lx1,1,-1
!   do iy=ly,1,-1
!   do iz=lz,1,-1
!   tt=eye*(kx(ix,iy,iz)*xp+ky(ix,iy,iz)*yp+kz(ix,iy,iz)*zp)
!   tt=ui(ix,iy,iz)*exp(tt)
!   upi=upi + tt
!   end do
!   end do
!   end do
!   do ix=lx1,2,-1
!   do iy=ly,1,-1
!   do iz=lz,1,-1
!   tt=-eye*(kx(ix,iy,iz)*xp+ky(ix,iy,iz)*yp+kz(ix,iy,iz)*zp)
!   tt=conjg(ui(ix,iy,iz))*exp(tt)
!   upi=upi + tt
!   end do
!   end do
!   end do

  
end subroutine specint      
