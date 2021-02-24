SUBROUTINE sgsm_p (p,t11,t12,t13,t22,t23,t33,kx,ky,kz,lx1,ly,lz,delta)
  USE mconstant
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lx1,ly,lz
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN) :: p
  REAL(SP), DIMENSION(lx1,ly,lz), INTENT(IN) :: kx,ky,kz
  REAL(SP), INTENT(IN) :: delta

  REAL(SP), DIMENSION(lx1,ly,lz), INTENT(OUT) :: t11,t12,t13
  REAL(SP), DIMENSION(lx1,ly,lz), INTENT(OUT) :: t22,t23,t33

  REAL(SP) :: coef

  coef=-delta*delta/12

  t11=-coef*kx*kx*p
  t12=-coef*kx*ky*p
  t13=-coef*kx*kz*p
  t22=-coef*ky*ky*p
  t23=-coef*ky*kz*p
  t33=-coef*kz*kz*p

END SUBROUTINE sgsm_p      
