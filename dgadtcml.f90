subroutine dgadtcml(t,gaij,dgadt,ny)
  use mconstant
  implicit none

  integer, parameter :: three = 3
  real(dp), parameter :: gmma = 0.01

  integer,  intent(in) :: ny
  real(dp), intent(in) :: t
  real(dp), dimension(three,three,2), intent(in) :: gaij
  real(dp), dimension(three,three,2), intent(out) :: dgadt

  real(dp), dimension(three,three) :: cinvij, dbij, id
  real(dp) :: trcinvij, traa

  id = 0.; id(1,1) = 1.; id(2,2) = 1.; id(3,3) = 1.

  ! Nonlinea term
  dbij = matmul( gaij(:,:,2), gaij(:,:,2) )
  traa = dbij(1,1) + dbij(2,2) + dbij(3,3)
 
  ! Cauchy-Green tensor
  cinvij = matmul( transpose(gaij(:,:,1)), gaij(:,:,1) )
  trcinvij = cinvij(1,1) + cinvij(2,2) + cinvij(3,3)

  ! rate of changes of gij
  dgadt(:,:,1) = - matmul( gaij(:,:,1), gaij(:,:,2) ) ! -  (trcinvij / 3) * (gaij(:,:,1) - id)
 
  ! rate of changes of aij
  dgadt(:,:,2) = - dbij + traa / trcinvij * cinvij - (trcinvij / 3) * gaij(:,:,2)

end subroutine dgadtcml
