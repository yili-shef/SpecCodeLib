subroutine energy_budget (ux,uy,uz,kx,ky,kz,k2,nu,tek,piek,ek,dek, &
                          lx1,ly,lz,nx,ny,nz)
  use mconstant
  implicit none

  integer,  intent(in) :: lx1,ly,lz
  integer,  intent(in) :: nx,ny,nz
  real(sp), intent(in) :: nu
  complex(sp),dimension(lx1,ly,lz),intent(in)  :: ux,uy,uz
  real(sp),   dimension(lx1,ly,lz),intent(in)  :: kx,ky,kz,k2
  real(sp),   dimension(lx1-1),    intent(out) :: tek,ek,dek,piek

  complex(sp),dimension(lx1,ly,lz) :: wx,wy,wz
  real(sp),   dimension(lx1,ly,lz) :: atmp
  integer :: i

  wx=eye*(ky*uz-kz*uy)
  wy=eye*(kz*ux-kx*uz)
  wz=eye*(kx*uy-ky*ux)
 
  atmp=2._sp*(.5_sp*(ux*conjg(ux)+uy*conjg(uy)+uz*conjg(uz)))
  atmp(1,:,:)=.5_sp*atmp(1,:,:) 
  do i=1,lx1-1
     ek(i)=sum(atmp,mask=(abs(sqrt(k2)-i).lt.0.5_sp))
     dek(i)=2._sp*nu*real(i,sp)*real(i,sp)*ek(i)
  end do

  call convec_dns(ux,uy,uz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
  call projection(wx,wy,wz,kx,ky,kz,lx1,ly,lz)

  ! note the definition of t is the nonlinear transfer for the kinetic energy of
  ! the kth mode |u(vector k)|^2/2.
  atmp=2._sp*real(wx*conjg(ux)+wy*conjg(uy)+wz*conjg(uz))
  atmp(1,:,:)=.5_sp*atmp(1,:,:)

  do i=lx1-1,1,-1
     tek(i)=sum(atmp,mask=(abs(sqrt(k2)-i).lt.0.5_sp))
     piek(i)=sum(tek(i:lx1-1))
  end do

end subroutine energy_budget
