SUBROUTINE rkcaller(ystart,rwkspc,nrwk,ny,x1,x2,eps,h1,hmin,hnext,nok,nbad)
  USE mconstant
  USE miderivs
  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: ny,nrwk
  REAL(DP), DIMENSION(ny), INTENT(INOUT) :: ystart
  REAL(DP), DIMENSION(nrwk), INTENT(IN) :: rwkspc
  REAL(DP), INTENT(IN) :: x1,x2,eps,h1,hmin
  REAL(DP), INTENT(OUT) :: hnext
  INTEGER, INTENT(INOUT) :: nok,nbad

  INTEGER, PARAMETER :: MAXSTP=10000
  INTEGER :: nstp
  REAL(DP) :: h, hdid, x, xsav
  REAL(DP), DIMENSION(ny) :: dydx, y, yscal

  x=x1
  h=SIGN(h1,x2-x1)
  y=ystart
  DO nstp=1,MAXSTP
    CALL derivs(x,y,dydx,ny,rwkspc,nrwk)
    yscal(:)=ABS(y(:))+ABS(h*dydx(:))+TINY
    IF ((x+h-x2)*(x+h-x1)>0.0_DP) h=x2-x
    CALL rkqs(y,dydx,rwkspc,nrwk,x,ny,h,eps,yscal,hdid,hnext)
    IF (hdid==h) THEN
            nok=nok+1
    ELSE
            nbad=nbad+1
    END IF
    IF ((x-x2)*(x2-x1)>=0.0_DP) THEN
            ystart=y
            RETURN
    END IF
    IF (ABS(hnext)<hmin) THEN
            WRITE(*,*) 'Stepsize smaller than hmin in rkcaller. Stop.'
            STOP
    END IF
    h=hnext
  END DO
  WRITE(*,*) 'Too many steps in rkcaller. Stop.'
  STOP
END SUBROUTINE rkcaller
