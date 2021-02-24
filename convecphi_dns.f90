subroutine convecphi_dns(vx,vy,vz,phi,wphi,kx,ky,kz,lx1,ly,lz,nx,ny,nz)
  use mconstant
  use mfftwplan3d
  implicit none

  integer, intent(in) :: lx1,ly,lz,nx,ny,nz
  complex(sp), dimension(lx1,ly,lz), intent(in)    :: vx,vy,vz,phi
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: wphi
  real(sp), dimension(lx1), intent(in) :: kx
  real(sp), dimension(ly),  intent(in) :: ky
  real(sp), dimension(lz),  intent(in) :: kz

  complex(sp), dimension(lx1,ly,lz) :: vt,phit
  real(sp) :: const, ignore_me
  integer :: ii, jj, kk

  const=1._sp/real(nx*ny*nz,sp)

  phit = phi
  call rfftwnd_f77_one_complex_to_real(c2r3d,phit,ignore_me)

  vt=vx
  call rfftwnd_f77_one_complex_to_real(c2r3d,vt,ignore_me)
  vt = cmplx( real(vt) * real(phit), aimag(vt) * aimag(phit) )
  call rfftwnd_f77_one_real_to_complex(r2c3d,vt,ignore_me)
  do ii = 1, lx1
    wphi(ii,:,:) = eye * kx(ii) * vt(ii,:,:)
  end do

  vt=vy
  call rfftwnd_f77_one_complex_to_real(c2r3d,vt,ignore_me)
  vt = cmplx( real(vt) * real(phit), aimag(vt) * aimag(phit) )
  call rfftwnd_f77_one_real_to_complex(r2c3d,vt,ignore_me)
  do jj = 1, ly
    wphi(:,jj,:) = wphi(:,jj,:) + eye * ky(jj) * vt(:,jj,:)
  end do


  vt=vz
  call rfftwnd_f77_one_complex_to_real(c2r3d,vt,ignore_me)
  vt = cmplx( real(vt) * real(phit), aimag(vt) * aimag(phit) )
  call rfftwnd_f77_one_real_to_complex(r2c3d,vt,ignore_me)
  do kk = 1, lz
    wphi(:,:,kk) = wphi(:,:,kk) + eye * kz(kk) * vt(:,:,kk)
  end do

  wphi=wphi*const

end subroutine convecphi_dns
