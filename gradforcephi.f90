subroutine gradforcephi(uz,phirms,urms,epsphi,wphi)
  use mconstant
  implicit none

  real(sp), parameter :: cphiuz = -0.3 !TODO find the correct value.

  real(sp), intent(in) :: epsphi, urms, phirms
  complex(sp), dimension(:,:,:), intent(in) :: uz
  complex(sp), dimension(:,:,:), intent(inout) :: wphi

  real(sp) :: G

  G = -epsphi/(urms*phirms*cphiuz)

  wphi = wphi - G*uz

end subroutine gradforcephi  
