subroutine force_heli(fx,fy,fz,vx,vy,vz,wx,wy,wz,k2,lx1,ly,lz,fkmax,eps,eta)
  use mconstant
  implicit none

  integer,  intent(in) :: lx1,ly,lz
  real(sp), intent(in) :: fkmax,eps,eta
  complex(sp), dimension(lx1,ly,lz), intent(out) :: fx,fy,fz
  complex(sp), dimension(lx1,ly,lz), intent(in)  :: vx,vy,vz,wx,wy,wz
  real(sp),    dimension(lx1,ly,lz), intent(in)  :: k2
  
  real(sp) :: den, kinet, enst, uw, a, b, fk2
  real(sp), dimension(lx1,ly,lz) :: tmp

  ! input both energy and helicity
  fk2=fkmax*fkmax
  where (k2 .lt. fk2)
  tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
  elsewhere
          tmp=0.
  end where
  tmp(1,:,:) = .5_sp * tmp(1,:,:)
  kinet=sum(tmp)

  where (k2 .lt. fk2)
  tmp = wx*conjg(wx) + wy*conjg(wy) + wz*conjg(wz)
  elsewhere
          tmp=0.
  end where
  tmp(1,:,:) = .5_sp * tmp(1,:,:)
  enst=sum(tmp)

  where (k2 .lt. fk2)
  tmp=vx*conjg(wx)+vy*conjg(wy)+vz*conjg(wz)+conjg(vx)*wx+conjg(vy)*wy+conjg(vz)*wz
  elsewhere
          tmp=0.
  end where
  tmp(1,:,:) = .5_sp * tmp(1,:,:)
  uw=sum(tmp)
  
  den=4._sp*kinet*enst-uw*uw
  a=(4._sp*eps*enst-eta*uw)/(2._sp*den)
  b=(eta*kinet-eps*uw)/den

  where (k2.lt.fk2)
    fx=a*vx+b*wx
    fy=a*vy+b*wy
    fz=a*vz+b*wz    
  elsewhere
            fx=0.
            fy=0.
            fz=0.
  endwhere
  fx(1,1,1)=0._sp
  fy(1,1,1)=0._sp
  fz(1,1,1)=0._sp

end subroutine force_heli
