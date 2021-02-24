subroutine forcephi(phi,wphi,k2,lx1,ly,lz,fkmax,epsphi)
  use mconstant
  implicit none

  integer,  intent(in) :: lx1,ly,lz
  real(sp), intent(in) :: fkmax,epsphi
  complex(sp), dimension(lx1,ly,lz), intent(in) :: phi
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: wphi
  real(sp),    dimension(lx1,ly,lz), intent(in)  :: k2

  real(sp),    dimension(lx1,ly,lz) :: tmp
  real(sp) :: kinet, a

  where (k2 .lt. fkmax*fkmax)
    tmp = phi * conjg(phi)
  elsewhere
    tmp = 0.
  end where
  tmp(1,:,:) = .5_sp * tmp(1,:,:)
  kinet=sum(tmp)

  a=epsphi/(2._sp*kinet)
  where (k2 .lt. fkmax*fkmax)
    wphi=wphi+a*phi
  endwhere

end subroutine forcephi      
