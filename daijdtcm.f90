subroutine daijdtcm(t,aij,daijdt,ny)
  use mconstant
  implicit none

  integer, parameter :: three = 3
  real(dp), parameter :: gmma = 0.01

  integer,  intent(in) :: ny
  real(dp), intent(in) :: t
  real(dp), dimension(three,three), intent(in) :: aij
  real(dp), dimension(three,three), intent(out) :: daijdt

  real(dp), dimension(three,three) :: cinvij, dbij, id
  real(dp) :: anorm, trcinvij, traa
  ! real(dp) :: f1, f2

  ! f1 = 2*y(1) - y(1) * y(1)
  ! f2 = -2*y(1)
  ! f1 = -f1
  ! f1 = f2
  ! dydx(1) = (4 * exp(f1) - exp(-2*f1)) / (2 * exp(f1) + exp(-2*f1)) * y(1) * y(1) ! &
  !           ! - (1./3.) * (2*exp(f2) + exp(-2*f2))*y(1)

  id = 0.; id(1,1) = 1.; id(2,2) = 1.; id(3,3) = 1.

  ! Nonlinea term
  dbij = matmul(aij, aij)
  traa = dbij(1,1) + dbij(2,2) + dbij(3,3)
 
  ! Matrix exponential for pressure
  cinvij = gmma * aij - .5 * gmma * gmma * (dbij - traa * id / 3)
  cinvij = transpose(cinvij)

  daijdt = matmul(cinvij, transpose(cinvij)) 
  anorm = sqrt( daijdt(1,1) + daijdt(2,2) + daijdt(3,3) )

  call taylorch(cinvij, anorm, three)

  cinvij = matmul(cinvij, transpose(cinvij))

  trcinvij = cinvij(1,1) + cinvij(2,2) + cinvij(3,3)
 
  ! rate of changes of aij
  daijdt = -dbij + traa / trcinvij * cinvij - (trcinvij / 3) * aij 

end subroutine daijdtcm
