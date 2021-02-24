subroutine phiforce(phi,wphi,k2,lx1,ly,lz,fkmax,epsphi)
  use mconstant
  implicit none

  integer,  intent(in) :: lx1,ly,lz
  real(sp), intent(in) :: fkmax,epsphi
  complex(sp), dimension(lx1,ly,lz), intent(in) :: phi
  real(sp),    dimension(lx1,ly,lz), intent(in) :: k2
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: wphi

  real(sp) :: kinet, a, tmp
  integer  :: ii, jj, kk

  kinet = 0.
  a = fkmax*fkmax
  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    if ( k2(ii,jj,kk) .lt. a) then 
        tmp = phi(ii,jj,kk) * conjg(phi(ii,jj,kk)) 
        if ( ii .eq. 1) then
            kinet = kinet + tmp/2.
        else
            kinet = kinet + tmp
        end if 
    end if

  end do
  end do
  end do

  tmp=epsphi/(2._sp*kinet)
  where ( k2 .lt. a)
      wphi = wphi + tmp * phi
  end where

end subroutine phiforce      
