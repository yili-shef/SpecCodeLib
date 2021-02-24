subroutine xyztocrayafull (ux, uy, uz, v1, v2, kx, ky, kz, lx1, ly, lz)
  use mconstant
  implicit none

  integer, intent(in) :: lx1, ly, lz
  complex(sp), dimension(lx1,ly,lz), intent(in)  :: ux, uy, uz
  complex(sp), dimension(lx1,ly,lz), intent(out) :: v1, v2 
  real(sp),    dimension(lx1),       intent(in)  :: kx
  real(sp),    dimension(ly),        intent(in)  :: ky
  real(sp),    dimension(lz),        intent(in)  :: kz

  real(sp) :: eax, eay, eaz, ebx, eby, ebz
  real(sp) :: kperp, k, kxii, kyjj, kzkk
  integer :: ii, jj, kk

  do kk = 1, lz
  kzkk = kz(kk)
  do jj = 1, ly
  kyjj = ky(jj)
  do ii = 1, lx1
    kxii = kx(ii)

    kperp = sqrt( kxii * kxii + kyjj * kyjj )
 
    if (kperp .lt. smallest) then
 
        eax = 1._sp; eay = 0._sp; eaz = 0._sp
        ebx = 0._sp; eby = 1._sp; ebz = 0._sp
 
    else
 
        k = sqrt( kperp * kperp + kzkk * kzkk )
 
        eax = kyjj / kperp; eay = - kxii / kperp; eaz = 0._sp
 
        ebx = kxii * kzkk / (kperp * k)
        eby = kyjj * kzkk / (kperp * k)
        ebz = - kperp / k
 
    end if
    v1(ii,jj,kk) = ux(ii,jj,kk) * eax + uy(ii,jj,kk) * eay + uz(ii,jj,kk) * eaz
    v2(ii,jj,kk) = ux(ii,jj,kk) * ebx + uy(ii,jj,kk) * eby + uz(ii,jj,kk) * ebz
  end do
  end do
  end do

end subroutine xyztocrayafull
