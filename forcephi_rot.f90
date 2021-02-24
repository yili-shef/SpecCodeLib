subroutine forcephi_rot(phi,wphi,k2,lx1,ly,lz,fk1,fk2,epsphi)
  use mconstant
  implicit none

  integer,  intent(in) :: lx1,ly,lz
  real(sp), intent(in) :: fk1,fk2,epsphi
  complex(sp), dimension(lx1,ly,lz), intent(in) :: phi
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: wphi
  real(sp),    dimension(lx1,ly,lz), intent(in)  :: k2

  real(sp),    dimension(lx1,ly,lz) :: tmp
  real(sp) :: kinet, a

  where (k2 .lt. fk2*fk2 .and. k2 .ge. fk1*fk1)
    tmp = phi * conjg(phi)
  elsewhere
    tmp = 0.
  end where
  tmp(1,:,:) = .5_sp * tmp(1,:,:)
  kinet=sum(tmp)

  a=epsphi/(2._sp*kinet)
  where (k2 .lt. fk2*fk2 .and. k2 .ge. fk1*fk1)
    wphi=wphi+a*phi
  endwhere

end subroutine forcephi_rot
