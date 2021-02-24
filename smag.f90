subroutine smag (vx,vy,vz,t11,t12,t13,t22,t23,t33,kx,ky,kz,lx1,ly,lz, &
                 nx,ny,nz,delta)
  use mconstant
  use mfftwplan3d
  implicit none

  real(sp), parameter  :: c1=0.026 !cut-off filter
  integer,  intent(in) :: lx1,ly,lz,nx,ny,nz
  real(sp), intent(in) :: delta
  complex(sp), dimension(lx1,ly,lz), intent(in)  :: vx,vy,vz
  complex(sp), dimension(lx1,ly,lz), intent(out) :: t11,t12,t13,t22,t23,t33
  real(sp),    dimension(lx1), intent(in)  :: kx
  real(sp),    dimension(ly),  intent(in)  :: ky
  real(sp),    dimension(lz),  intent(in)  :: kz

  real(sp) :: ignore_me, const, sr, si
  integer  :: ii,jj,kk

  const=1._sp/real(nx*ny*nz,sp)

  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    t11(ii,jj,kk) = eye * kx(ii) * vx(ii,jj,kk)
    t22(ii,jj,kk) = eye * ky(jj) * vy(ii,jj,kk)
    t33(ii,jj,kk) = eye * kz(kk) * vz(ii,jj,kk)
    t12(ii,jj,kk) = eye * ( kx(ii) * vy(ii,jj,kk) + ky(jj) * vx(ii,jj,kk) ) * .5
    t13(ii,jj,kk) = eye * ( kx(ii) * vz(ii,jj,kk) + kz(kk) * vx(ii,jj,kk) ) * .5
    t23(ii,jj,kk) = eye * ( ky(jj) * vz(ii,jj,kk) + kz(kk) * vy(ii,jj,kk) ) * .5
  end do
  end do
  end do

  call rfftwnd_f77_one_complex_to_real(c2r3d,t11,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,t12,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,t13,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,t22,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,t23,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,t33,ignore_me)

  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1 - 1
    sr = real( t11(ii,jj,kk) )**2 + real( t22(ii,jj,kk) )**2 + real( t33(ii,jj,kk) )**2       &
         + 2. * (real( t12(ii,jj,kk) )**2 + real( t13(ii,jj,kk) )**2 + real( t23(ii,jj,kk) )**2)   
    sr = sqrt(2. * sr)

    si = aimag( t11(ii,jj,kk) )**2 + aimag( t22(ii,jj,kk) )**2 + aimag( t33(ii,jj,kk) )**2       &
         + 2. * (aimag( t12(ii,jj,kk) )**2 + aimag( t13(ii,jj,kk) )**2 + aimag( t23(ii,jj,kk) )**2)   
    si = sqrt(2. * si)

    t11(ii,jj,kk) = cmplx( sr * real( t11(ii,jj,kk) ), si * aimag( t11(ii,jj,kk) ) )
    t22(ii,jj,kk) = cmplx( sr * real( t22(ii,jj,kk) ), si * aimag( t22(ii,jj,kk) ) )
    t33(ii,jj,kk) = cmplx( sr * real( t33(ii,jj,kk) ), si * aimag( t33(ii,jj,kk) ) )
    t12(ii,jj,kk) = cmplx( sr * real( t12(ii,jj,kk) ), si * aimag( t12(ii,jj,kk) ) )
    t13(ii,jj,kk) = cmplx( sr * real( t13(ii,jj,kk) ), si * aimag( t13(ii,jj,kk) ) )
    t23(ii,jj,kk) = cmplx( sr * real( t23(ii,jj,kk) ), si * aimag( t23(ii,jj,kk) ) )

    t11(ii,jj,kk) = -2. * c1 * delta * delta * t11(ii,jj,kk)*const
    t12(ii,jj,kk) = -2. * c1 * delta * delta * t12(ii,jj,kk)*const
    t13(ii,jj,kk) = -2. * c1 * delta * delta * t13(ii,jj,kk)*const
    t22(ii,jj,kk) = -2. * c1 * delta * delta * t22(ii,jj,kk)*const
    t23(ii,jj,kk) = -2. * c1 * delta * delta * t23(ii,jj,kk)*const
    t33(ii,jj,kk) = -2. * c1 * delta * delta * t33(ii,jj,kk)*const

  end do
  end do
  end do
  
  call rfftwnd_f77_one_real_to_complex(r2c3d,t11,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,t12,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,t13,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,t22,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,t23,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,t33,ignore_me)


end subroutine smag
