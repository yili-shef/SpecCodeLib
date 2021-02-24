subroutine pressure_dns(vx,vy,vz,k2,kx,ky,kz,lx1,ly,lz,nx,ny,nz,p)
  use mconstant
  use mfftwplan3d
  implicit none

  integer, intent(in) :: lx1,ly,lz,nx,ny,nz
  complex(sp),dimension(lx1,ly,lz),intent(in)    :: vx,vy,vz
  real(sp),   dimension(lx1,ly,lz),intent(in)    :: k2
  real(sp),   dimension(lx1),      intent(in)    :: kx
  real(sp),   dimension(ly),       intent(in)    :: ky
  real(sp),   dimension(lz),       intent(in)    :: kz
  complex(sp),dimension(lx1,ly,lz),intent(out)   :: p

  complex(sp),dimension(lx1,ly,lz) :: axy,ayx,aa
  real(sp) :: const, ignore_me
  integer  :: ii,jj,kk

  const=1._sp/real(nx*ny*nz,sp)

  do kk=1,lz
  do jj=1,ly
    axy(:,ly,lz)=eye*kx*vx(:,ly,lz) ! a11
  end do
  end do
  call rfftwnd_f77_one_complex_to_real(c2r3d,axy,ignore_me)
  aa=cmplx(real(axy,sp)*real(axy,sp), aimag(axy)*aimag(axy))

  do kk=1,lz
  do ii=1,lx1
    axy(ii,:,kk)=eye*ky*vy(ii,:,kk) ! a22
  end do
  end do
  call rfftwnd_f77_one_complex_to_real(c2r3d,axy,ignore_me)
  aa=aa+cmplx(real(axy,sp)*real(axy,sp), aimag(axy)*aimag(axy))

  do jj=1,ly
  do ii=1,lx1
    axy(ii,jj,:)=eye*kz*vz(ii,jj,:) ! a22
  end do
  end do
  call rfftwnd_f77_one_complex_to_real(c2r3d,axy,ignore_me)
  aa=aa+cmplx(real(axy,sp)*real(axy,sp), aimag(axy)*aimag(axy))

  do kk=1,lz
  do jj=1,ly
    axy(:,jj,kk)=eye*kx*vy(:,jj,kk) !a12
  end do
  do ii=1,lx1
    ayx(ii,:,kk)=eye*ky*vx(ii,:,kk) !a21
  end do
  end do
  call rfftwnd_f77_one_complex_to_real(c2r3d,axy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,ayx,ignore_me)
  aa=aa+2.*cmplx(real(axy,sp)*real(ayx,sp), aimag(axy)*aimag(ayx))

  do jj=1,ly
  do kk=1,lz
    axy(:,jj,kk)=eye*kx*vz(:,jj,kk) !a13
  end do
  do ii=1,lx1
    ayx(ii,jj,:)=eye*kz*vx(ii,jj,:) !a31
  end do
  end do
  call rfftwnd_f77_one_complex_to_real(c2r3d,axy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,ayx,ignore_me)
  aa=aa+2.*cmplx(real(axy,sp)*real(ayx,sp), aimag(axy)*aimag(ayx))

  do ii=1,lx1
  do kk=1,lz
    axy(ii,:,kk)=eye*ky*vz(ii,:,kk) !a23
  end do
  do jj=1,ly
    ayx(ii,jj,:)=eye*kz*vy(ii,jj,:) !a32
  end do
  end do
  call rfftwnd_f77_one_complex_to_real(c2r3d,axy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,ayx,ignore_me)
  aa=aa+2.*cmplx(real(axy,sp)*real(ayx,sp), aimag(axy)*aimag(ayx))

  aa=aa*const
  call rfftwnd_f77_one_real_to_complex(r2c3d,aa,ignore_me)
  p=p/k2
  call symmetrize(p,k2,lx1,ly,lz)

end subroutine pressure_dns
