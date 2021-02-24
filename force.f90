! Combine the right hand side of the Navier-Stokes equation with the forcing
! term. When only the forcing term is needed, set wx wy wz =0 upon input.
subroutine force(wx,wy,wz,vx,vy,vz,k2,lx1,ly,lz,fkmax,eps)
  use mconstant
  implicit none

  integer,  intent(in) :: lx1,ly,lz
  real(sp), intent(in) :: fkmax,eps
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: wx,wy,wz
  complex(sp), dimension(lx1,ly,lz), intent(in)  :: vx,vy,vz
  real(sp),    dimension(lx1,ly,lz), intent(in)  :: k2
  
  real(sp) :: kinet, a
  real(sp), dimension(lx1,ly,lz) :: tmp

  where (k2 .lt. fkmax*fkmax)
    tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
  elsewhere
    tmp = 0.
  end where
  tmp(1,:,:) = .5_sp * tmp(1,:,:)
  kinet=sum(tmp)

  a=eps/(2._sp*kinet)
  where (k2 .lt. fkmax*fkmax)
    wx=wx+a*vx
    wy=wy+a*vy
    wz=wz+a*vz
  endwhere

end subroutine force
