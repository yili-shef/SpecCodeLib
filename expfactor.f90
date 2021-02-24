program efcheck
! check the evaluation of the exponential integration facor in
! hyperviscosity simulation
  use mconstant

  implicit none
  integer, parameter :: nx=128, q=8,lx=nx/2,lx2=lx*lx, lx1=nx/2+1
  real(sp), dimension(lx1,nx,nx) :: k2
  real(sp), dimension(lx1) :: kx
  real(sp), dimension(nx) :: ky,kz

  real(sp) :: rnukq, ef1, ef2, ef3, dt
  integer  :: mx,my,mz

  rnukq=200.
  dt = 0.0025
  mx=58
  my=2
  mz=65

  call wavenumber(kx,ky,kz,k2,lx1,nx,nx)
  write(*,*) 'k2', k2(mx,my,mz)
  ef1=exp(-(k2(mx,my,mz)/lx2)**q*dt*rnukq)
  ef2=exp(-(k2(mx,my,mz)*0.00047344)**q*dt)
  ef3=exp(-(sqrt(k2(mx,my,mz))*0.021758)**(2*q)*dt)
  write(*,*) ef1,ef2,ef3
end program efcheck
