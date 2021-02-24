! use work spaces to reduce the use of memory.
subroutine convecdnsphi2c(va,vb,phi,wphi,kx,ky,kz,lx1,ly,lz,const,wspace1,wspace2)
  use mconstant
  use mfftwplan3d
  implicit none

  integer, intent(in) :: lx1,ly,lz
  real(sp), intent(in) :: const
  complex(sp), dimension(lx1,ly,lz), intent(in)    :: va,vb,phi
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: wphi
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: wspace1, wspace2
  real(sp), dimension(lx1), intent(in) :: kx
  real(sp), dimension(ly),  intent(in) :: ky
  real(sp), dimension(lz),  intent(in) :: kz

  real(sp) :: ignore_me
  integer :: ii, jj, kk

  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    call crayatoxyz (va(ii,jj,kk), vb(ii,jj,kk), kx(ii), ky(jj), kz(kk), &
                     wspace1(ii,jj,kk), ignore_me, ignore_me)
    wspace2(ii,jj,kk) = eye * kx(ii) * phi(ii,jj,kk)
  end do
  end do
  end do
  call rfftwnd_f77_one_complex_to_real(c2r3d,wspace1,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wspace2,ignore_me)
  wphi = cmplx( real(wspace1) * real(wspace2), aimag(wspace1) * aimag(wspace2) )

  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    call crayatoxyz (va(ii,jj,kk), vb(ii,jj,kk), kx(ii), ky(jj), kz(kk), &
                     ignore_me, wspace1(ii,jj,kk), ignore_me)
    wspace2(ii,jj,kk) = eye * ky(jj) * phi(ii,jj,kk)
  end do
  end do
  end do
  call rfftwnd_f77_one_complex_to_real(c2r3d,wspace1,  ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wspace2,ignore_me)
  wphi = wphi + cmplx( real(wspace1) * real(wspace2), aimag(wspace1) * aimag(wspace2) )

  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    call crayatoxyz (va(ii,jj,kk), vb(ii,jj,kk), kx(ii), ky(jj), kz(kk), &
                     ignore_me, ignore_me, wspace1(ii,jj,kk))
    wspace2(ii,jj,kk) = eye * kz(kk) * phi(ii,jj,kk)
  end do
  end do
  end do
  call rfftwnd_f77_one_complex_to_real(c2r3d,wspace1,  ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wspace2,ignore_me)
  wphi = wphi + cmplx( real(wspace1) * real(wspace2), aimag(wspace1) * aimag(wspace2) )

  wphi=wphi*const
  call rfftwnd_f77_one_real_to_complex(r2c3d,wphi,ignore_me)

end subroutine convecdnsphi2c

