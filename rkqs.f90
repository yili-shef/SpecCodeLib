SUBROUTINE rkqs(y,dydx,rwkspc,nrwk,x,ny,htry,eps,yscal,hdid,hnext)
  USE mconstant
  USE miderivs
  IMPLICIT NONE
  
  INTEGER,  INTENT(IN) :: ny,nrwk
  REAL(DP), DIMENSION(nrwk), INTENT(IN) :: rwkspc
  REAL(DP), DIMENSION(ny), INTENT(INOUT) :: y
  REAL(DP), DIMENSION(ny), INTENT(IN) :: dydx,yscal
  REAL(DP), INTENT(INOUT) :: x
  REAL(DP), INTENT(IN) :: htry,eps
  REAL(DP), INTENT(OUT) :: hdid,hnext

  REAL(DP) :: errmax, h, htemp, xnew
  REAL(DP), DIMENSION(ny) :: yerr,ytemp
  REAL(DP), PARAMETER :: SAFETY=0.9_DP, PGROW=-0.2_DP, PSHRNK=-0.25_DP, &
    ERRCON=1.89E-4_DP

  h=htry
  DO
    CALL rkck(y,dydx,rwkspc,nrwk,x,ny,h,ytemp,yerr)
    errmax=MAXVAL(ABS(yerr(:)/yscal(:)))/eps
    IF (errmax<=1.0_DP) EXIT
    htemp=SAFETY*h*(errmax**PSHRNK)
    h=SIGN(MAX(ABS(htemp),0.1_DP*ABS(h)),h)
    xnew=x+h
    IF (xnew==x) THEN
            WRITE(*,*) 'stepsize underflow in rkqs. stop'
            STOP
    END IF
  END DO
  IF (errmax > ERRCON) THEN
          hnext=SAFETY*h*(errmax**PGROW)
  ELSE
          hnext=5.0_DP*h
  END IF
  hdid=h
  x=x+h
  y=ytemp
END SUBROUTINE rkqs
