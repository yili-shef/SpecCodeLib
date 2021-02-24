! Adopted from Numerical Recipe. Using the fifth
! order Cash-Karp Runge-Kutta method to advance 
! the solution.

SUBROUTINE rkck(y,dydx,rwkspc,nrwk,x,ny,h,yout,yerr)
  USE mconstant
  USE miderivs
  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: ny,nrwk
  REAL(DP), DIMENSION(nrwk), INTENT(IN) :: rwkspc
  REAL(DP), DIMENSION(ny), INTENT(IN) :: y, dydx
  REAL(DP), INTENT(IN) :: x, h
  REAL(DP), DIMENSION(ny), INTENT(OUT) :: yout, yerr

  REAL(DP), DIMENSION(ny) :: ak2,ak3,ak4,ak5,ak6,ytemp
  REAL(DP), PARAMETER :: A2=0.2_DP, A3=0.3_DP, A4=0.6_DP, A5=1.0_DP, &
    A6=0.875_DP, B21=0.2_DP, B31=3.0_DP/40.0_DP, B32=9.0_DP/40.0_DP, &
    B41=0.3_DP, B42=-0.9_DP, B43=1.2_DP, B51=-11.0_DP/54.0_DP, &
    B52=2.5_DP, B53=-70.0_DP/27.0_DP, B54=35.0_DP/27.0_DP, &
    B61=1631.0_DP/55296.0_DP, B62=175.0_DP/512.0_DP, &
    B63=575.0_DP/13824.0_DP, B64=44275.0_DP/110592.0_DP, &
    B65=253.0_DP/4096.0_DP, C1=37.0_DP/378.0_DP, &
    C3=250.0_DP/621.0_DP, C4=125.0_DP/594.0_DP,  &
    C6=512.0_DP/1771.0_DP, DC1=C1-2825.0_DP/27648.0_DP, &
    DC3=C3-18575.0_DP/48384.0_DP, DC4=C4-13525.0_DP/55296.0_DP, &
    DC5=-277.0_DP/14336.0_DP, DC6=C6-0.25_DP

  ytemp=y+B21*h*dydx                                !First step
  CALL derivs(x+A2*h,ytemp,ak2,ny,rwkspc,nrwk)                   !Second step
  ytemp=y+h*(B31*dydx+B32*ak2)
  CALL derivs(x+A3*h,ytemp,ak3,ny,rwkspc,nrwk)                     !Third step
  ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
  CALL derivs(x+A4*h,ytemp,ak4,ny,rwkspc,nrwk)                     !Fourth step
  ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
  CALL derivs(x+A5*h,ytemp,ak5,ny,rwkspc,nrwk)                     !Fifth step
  ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
  CALL derivs(x+A6*h,ytemp,ak6,ny,rwkspc,nrwk)                     !Sixth step
  yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
  yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
END SUBROUTINE rkck
