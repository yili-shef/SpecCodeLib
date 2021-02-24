subroutine force_rescale(ap,am,k2,lx1,ly,lz,fp1,fp2)
  use mconstant
  implicit none

  integer,  intent(in) :: lx1, ly, lz
  real(sp), intent(in) :: fp1, fp2
  real(sp),    dimension(lx1,ly,lz), intent(in)    :: k2
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: ap, am

  real(sp) :: ff1, ff2, eks1, eks2, coefx, coefy, coefz, coef
  integer :: ii, jj, kk, ll, jjmod, kkmod

  integer, parameter :: fkmax = 3

  ! ------- Make use of the particular way the wavenumbers are stored ---------
  eks1 = 0.0_sp
  eks2 = 0.0_sp
  do kk = 1, fkmax
  do jj = 1, fkmax
  do ii = 1, fkmax

    ll=floor(sqrt(k2(ii,jj,kk))+.5_sp)

    if (ii .eq. 1) then 
      coefx = 0.5_sp
    else
      coefx = 1.0_sp
    end if

    if (jj .eq. 1) then
      coefy = 0.5_sp
    else 
      coefy = 1.0_sp
    end if

    if (kk .eq. 1) then
      coefz = 0.5_sp
    else 
      coefz = 1.0_sp
    end if
    coef = coefx * coefy * coefz

    jjmod = mod(ly - jj + 1, ly) + 1 ! wavenumbers -1, -2
    kkmod = mod(lz - kk + 1, lz) + 1

    if (ll .eq. 1) then
      eks1 = eks1 + coef * real( ap(ii,jj,kk) * conjg(ap(ii,jj,kk)) + &
                                 am(ii,jj,kk) * conjg(am(ii,jj,kk)) ) 
      eks1 = eks1 + coef * real( ap(ii,jjmod,kk) * conjg(ap(ii,jjmod,kk)) + &
                                 am(ii,jjmod,kk) * conjg(am(ii,jjmod,kk)) )
      eks1 = eks1 + coef * real( ap(ii,jj,kkmod) * conjg(ap(ii,jj,kkmod)) + &
                                 am(ii,jj,kkmod) * conjg(am(ii,jj,kkmod)) ) 
      eks1 = eks1 + coef * real( ap(ii,jjmod,kkmod) * conjg(ap(ii,jjmod,kkmod)) + &
                                 am(ii,jjmod,kkmod) * conjg(am(ii,jjmod,kkmod)) )
    end if

    if (ll .eq. 2) then
      eks2 = eks2 + coef * real( ap(ii,jj,kk) * conjg(ap(ii,jj,kk)) + &
                                 am(ii,jj,kk) * conjg(am(ii,jj,kk)) ) 
      eks2 = eks2 + coef * real( ap(ii,jjmod,kk) * conjg(ap(ii,jjmod,kk)) + &
                                 am(ii,jjmod,kk) * conjg(am(ii,jjmod,kk)) )
      eks2 = eks2 + coef * real( ap(ii,jj,kkmod) * conjg(ap(ii,jj,kkmod)) + &
                                 am(ii,jj,kkmod) * conjg(am(ii,jj,kkmod)) ) 
      eks2 = eks2 + coef * real( ap(ii,jjmod,kkmod) * conjg(ap(ii,jjmod,kkmod)) + &
                                 am(ii,jjmod,kkmod) * conjg(am(ii,jjmod,kkmod)) )
    end if

  end do
  end do
  end do

  ! Set the energy of the first shell (k=1) to fp1,
  ! the second shell to fp2
  ff1 = sqrt(fp1/eks1)
  ff2 = sqrt(fp2/eks2)

  do kk = 1, fkmax
  do jj = 1, fkmax
  do ii = 1, fkmax

    ll=floor(sqrt(k2(ii,jj,kk))+.5_sp)

    if (jj .eq. 1) then
      coefy = 0.5_sp
    else 
      coefy = 1.0_sp
    end if
    if (kk .eq. 1) then
      coefz = 0.5_sp
    else 
      coefz = 1.0_sp
    end if
    coef = coefy * coefz

    jjmod = mod(ly - jj + 1, ly) + 1
    kkmod = mod(lz - kk + 1, lz) + 1

    eks1 = ff1**coef
    eks2 = ff2**coef

    if (ll .eq. 1) then
      ap(ii,jj,kk) = ap(ii,jj,kk) * eks1
      am(ii,jj,kk) = am(ii,jj,kk) * eks1

      ap(ii,jjmod,kk) = ap(ii,jjmod,kk) * eks1
      am(ii,jjmod,kk) = am(ii,jjmod,kk) * eks1

      ap(ii,jj,kkmod) = ap(ii,jj,kkmod) * eks1
      am(ii,jj,kkmod) = am(ii,jj,kkmod) * eks1

      ap(ii,jjmod,kkmod) = ap(ii,jjmod,kkmod) * eks1
      am(ii,jjmod,kkmod) = am(ii,jjmod,kkmod) * eks1
    end if

    if (ll .eq. 2) then
      ap(ii,jj,kk) = ap(ii,jj,kk) * eks2
      am(ii,jj,kk) = am(ii,jj,kk) * eks2

      ap(ii,jjmod,kk) = ap(ii,jjmod,kk) * eks2
      am(ii,jjmod,kk) = am(ii,jjmod,kk) * eks2

      ap(ii,jj,kkmod) = ap(ii,jj,kkmod) * eks2
      am(ii,jj,kkmod) = am(ii,jj,kkmod) * eks2

      ap(ii,jjmod,kkmod) = ap(ii,jjmod,kkmod) * eks2
      am(ii,jjmod,kkmod) = am(ii,jjmod,kkmod) * eks2
    end if

  end do
  end do
  end do


  return

end subroutine force_rescale
