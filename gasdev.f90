! Adopted from Numerical Recipe Page 280, with slight change.
FUNCTION gasdev(idum)
  USE mconstant
  INTEGER :: idum
  REAL(sp) :: gasdev
  INTEGER :: iset
  REAL(sp) :: fac,gset,rsq,v1,v2
  REAL :: ran1
  SAVE :: iset,gset
  DATA iset/0/

  IF (idum .lt. 0) iset=0
  IF (iset .eq. 0) THEN
          rsq =10_sp
          DO WHILE (rsq .GE. 1. .OR. rsq .EQ. 0.)
            v1=2._sp*ran1(idum)-1._sp
            v2=2._sp*ran1(idum)-1._sp
            rsq=v1*v1+v2*v2
          END DO
          fac=SQRT(-2._sp*LOG(rsq)/rsq)
          gset=v1*fac
          gasdev=v2*fac
          iset=1
  ELSE
          gasdev=gset
          iset=0
  END IF
  RETURN
END FUNCTION gasdev
