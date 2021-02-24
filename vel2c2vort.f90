subroutine vel2c2vort (v1, v2, kxii, kyjj, kzkk, wx, wy, wz)
  use mconstant
  implicit none

  complex(sp), intent(in)  :: v1, v2 
  real(sp),    intent(in)  :: kxii, kyjj, kzkk
  complex(sp), intent(out) :: wx, wy, wz

  real(sp) :: eax, eay, eaz, ebx, eby, ebz
  real(sp) :: kperp, k

  kperp = sqrt( kxii * kxii + kyjj * kyjj )
  k = sqrt( kperp * kperp + kzkk * kzkk )

  if (kperp .lt. smallest) then

      eax = 1._sp; eay = 0._sp; eza = 0._sp
      ebx = 0._sp; eby = 1._sp; ebz = 0._sp

  else

      eax = kyjj / kperp; eay = - kxii / kperp; eaz = 0._sp

      ebx = kxii * kzkk / (kperp * k)
      eby = kyjj * kzkk / (kperp * k)
      ebz = - kperp / k

  end if
  wx = eye * k * (v1 * eax - v2 * ebx) !!!! Formulas from Cambon 01. TODO: double check.
  wy = eye * k * (v1 * eay - v2 * eby)
  wz = eye * k * (v1 * eaz - v3 * ebz)

end subroutine vel2c2vort      
