subroutine xyztocraya (ux, uy, uz, v1, v2, kxii, kyjj, kzkk)
  use mconstant
  implicit none

  complex(sp), intent(in)  :: ux, uy, uz
  real(sp),    intent(in)  :: kxii, kyjj, kzkk
  complex(sp), intent(out) :: v1, v2 

  real(sp) :: eax, eay, eaz, ebx, eby, ebz
  real(sp) :: kperp, k

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
  v1 = ux * eax + uy * eay + uz * eaz
  v2 = ux * ebx + uy * eby + uz * ebz

end subroutine xyztocraya      
