subroutine dcadtcml(t,caij,dcadt,ny)
  use mconstant
  implicit none

  integer, parameter :: three = 3
  real(dp), parameter :: alpha = 0.1

  integer,  intent(in) :: ny
  real(dp), intent(in) :: t
  real(dp), dimension(three,three,2), intent(in) :: caij
  real(dp), dimension(three,three,2), intent(out) :: dcadt

  real(dp), dimension(three,three) :: aaij, id
  real(dp) :: trcinv, traa

  id = 0.; id(1,1) = 1.; id(2,2) = 1.; id(3,3) = 1.

  ! Nonlinea term
  aaij = matmul( caij(:,:,2), caij(:,:,2) )
  traa = aaij(1,1) + aaij(2,2) + aaij(3,3)
 
  ! Trace of cij^-1
  trcinv = caij(1,1,1) + caij(2,2,1) + caij(3,3,1)

  ! Rate of changes of cij^-1
  dcadt(:,:,1) = - matmul( transpose(caij(:,:,2)), caij(:,:,1) ) &
                 - matmul( caij(:,:,1), caij(:,:,2) ) &
                 + alpha * ( (trcinv / 3) * caij(:,:,1)  &
                 - matmul(caij(:,:,1), caij(:,:,1)) ) ! Damping term
 
  ! Rate of changes of aij
  dcadt(:,:,2) = - aaij + traa / trcinv * caij(:,:,1) - (trcinv / 3) * caij(:,:,2)

end subroutine dcadtcml
