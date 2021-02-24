PROGRAM postdeform
  USE mconstant
  IMPLICIT NONE

  INTEGER,  PARAMETER :: if1=1,if2=20
  REAL(DP), PARAMETER :: dt=0.017

  
  REAL(DP)     :: Qre,Rre,Qres,Rres,Q03
  REAL(DP),    DIMENSION(3,3):: Apre,Bre1,Bre2,Bre,Bren,Brens,Apres
  REAL(DP),    DIMENSION(3)  :: evle,evlen

  REAL(DP), DIMENSION(3,3) :: evtr,Cre,Apt
  REAL(DP), DIMENSION(3)   :: fv1,fv2

  INTEGER  :: i,ifile,id1,id2,id3,id4,ip,idum,matz,ierr,idnty(3,3),mm
  REAL(DP) :: Qt,Rt,tmp1,tmp2
  REAL(DP) :: fp2p,fp1p,fp2,fp1,fm2p,fm1p,fm2,fm1,f02p,f01p,f02,f01
  REAL(DP) :: time

  INTEGER :: VIsDouble,Debug
  CHARACTER(1) :: NULLCHR

  INTEGER :: TecIni, TecZne, TecDat, TecEnd

  VIsDouble=0
  Debug=1
  NULLCHR=CHAR(0)

  mm=40000

! This routine has the same function as postdeform_c.f90. Namely, it also
! compares the analytic and numerical solutions to the Restricted Euler dynamcis
! and deformation equations. The difference is that, first, here the analytic solutions
! used are not in dimensionless form, as oppose to in postdeform_c.f90, second,
! here the numericaly scheme is a second order one while that in
! postdeform_c.f90 is a fourth order one.
  
  idnty=0
  idnty(1,1)=1; idnty(2,2)=1; idnty(3,3)=1
 
  Apre=0.
  Apre(1,1)=1.; Apre(2,2)=1.; Apre(3,3)=-2.
!  Apre(1,2)=1.; Apre(2,1)=1.
  Bren=idnty
  Qre=-.5*SUM(MATMUL(Apre,Apre),MASK=idnty .EQ. 1)
  Rre=-(1./3.)*SUM(MATMUL(Apre,MATMUL(Apre,Apre)),MASK=idnty .EQ. 1)
  Q03=27.*Rre*Rre/4. + Qre*Qre*Qre
  WRITE(*,*) 'Qre,Rre,Q03', Qre, Rre, Q03
  
  Qt=Qre; Rt=Rre; Apt=Apre
  IF(Q03 .GT. 0.) THEN
    tmp2=(2./3.)*fp2p(Rt,Q03)*Qt**2
    tmp1=(2./3.)*fp1p(Rt,Q03)*Qt**2
    Bre1=idnty*tmp2-Apt*fp2(Rt,Q03)
    Bre2=Apt*fp1(Rt,Q03)-idnty*tmp1
    Bre1(:,:)=Bre1(:,:)/(fp1(Rt,Q03)*tmp2-fp2(Rt,Q03)*tmp1)
    Bre2(:,:)=Bre2(:,:)/(fp1(Rt,Q03)*tmp2-fp2(Rt,Q03)*tmp1)
  ELSE IF(Q03 .LT. 0.) THEN
    tmp2=(2./3.)*fm2p(Rt,Q03)*Qt**2
    tmp1=(2./3.)*fm1p(Rt,Q03)*Qt**2
    Bre1=idnty*tmp2-Apt*fm2(Rt,Q03)
    Bre2=Apt*fm1(Rt,Q03)-idnty*tmp1
    Bre1(:,:)=Bre1(:,:)/(fm1(Rt,Q03)*tmp2-fm2(Rt,Q03)*tmp1)
    Bre2(:,:)=Bre2(:,:)/(fm1(Rt,Q03)*tmp2-fm2(Rt,Q03)*tmp1)
  ELSE
    tmp2=(2./3.)*f02p(Rt,Q03)*Qt**2
    tmp1=(2./3.)*f01p(Rt,Q03)*Qt**2
    Bre1=idnty*tmp2-Apt*f02(Rt,Q03)
    Bre2=Apt*f01(Rt,Q03)-idnty*tmp1
    Bre1(:,:)=Bre1(:,:)/(f01(Rt,Q03)*tmp2-f02(Rt,Q03)*tmp1)
    Bre2(:,:)=Bre2(:,:)/(f01(Rt,Q03)*tmp2-f02(Rt,Q03)*tmp1)
  END IF

  time=0. 
  DO ifile=if1,if2

    Rt=Rre
    IF(Q03 .GT. 0.) THEN
      Bre=Bre1*fp1(Rt,Q03)+Bre2*fp2(Rt,Q03)
    ELSE IF(Q03 .LT. 0.) THEN
      Bre=Bre1*fm1(Rt,Q03)+Bre2*fm2(Rt,Q03)
    ELSE
      Bre=Bre1*f01(Rt)+Bre2*f02(Rt)
    END IF
    Cre=MATMUL(Bre,TRANSPOSE(Bre))
    matz=5
    CALL rs(3,3,Cre,evle,matz,evtr,fv1,fv2,ierr)

    Cre=MATMUL(Bren,TRANSPOSE(Bren))
    matz=5
    CALL rs(3,3,Cre,evlen,matz,evtr,fv1,fv2,ierr)

    WRITe(16,*) time, evle(1),evlen(1)

    Qres=Qre+dt*(-3.*Rre)
    Rres=Rre+dt*(2.*Qre*Qre/3.)
    Brens=Bren+dt*MATMUL(Bren,Apre)
    Apres=Apre+dt*(-MATMUL(Apre,Apre)-(2./3.)*idnty*Qre)
    
    Qre=.5*(Qre+Qres)+.5*dt*(-3.*Rres)
    Rre=.5*(Rre+Rres)+.5*dt*(2.*Qres*Qres/3.)
    Bren=.5*(Bren+Brens)+.5*dt*MATMUL(Brens,Apres)
    Apre=.5*(Apre+Apres)+.5*dt*(-MATMUL(Apres,Apres)-(2./3.)*idnty*Qres)

    time=time+dt
  END DO

  WRITE(*,*) 'finished'

END PROGRAM postdeform      


FUNCTION fp2p(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fp2p
  REAL(DP), INTENT(IN)  :: r,q03

  fp2p=.5*(((SQRT(q03)+3.*SQRT(3.)/2.*r)**2)**(-1./3.)+((SQRT(q03)-3.*SQRT(3.)/2.*r)**2)**(-1./3.))
END FUNCTION fp2p

FUNCTION fp1p(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fp1p
  REAL(DP), INTENT(IN)  :: r,q03

  fp1p=.25*SQRT(3.)*(((SQRT(q03)+3.*SQRT(3.)/2.*r)**2)**(-1./3.)-((SQRT(q03)-3.*SQRT(3.)/2.*r)**2)**(-1./3.))
END FUNCTION fp1p

FUNCTION fp2(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fp2
  REAL(DP), INTENT(IN)  :: r,q03
  REAL(DP) :: tt1,tt2

  tt1=SQRT(q03)+3.*SQRT(3.)/2.*r
  tt2=SQRT(q03)-3.*SQRT(3.)/2.*r
  fp2=(1./SQRT(3.))*(SIGN((ABS(tt1))**(1./3.),tt1)-SIGN((ABS(tt2))**(1./3.),tt2))
END FUNCTION fp2

FUNCTION fp1(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fp1
  REAL(DP), INTENT(IN)  :: r,q03
  REAL(DP) :: tt1,tt2

  tt1=SQRT(q03)+3.*SQRT(3.)/2.*r
  tt2=SQRT(q03)-3.*SQRT(3.)/2.*r
  fp1=.5*(SIGN((ABS(tt1))**(1./3.),tt1)+SIGN((ABS(tt2))**(1./3.),tt2))
END FUNCTION fp1

FUNCTION fm2p(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fm2p
  REAL(DP), INTENT(IN)  :: r,q03
  
  REAL(DP) :: tt

  tt=ABS(q03)+27.*r*r/4.
  fm2p=1.5*SQRT(3.)*r*tt**(-5./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r/SQRT(ABS(q03)))) &
       +.5*SQRT(3.)*SQRT(ABS(q03))*tt**(-5./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r/SQRT(ABS(q03))))
  
END FUNCTION fm2p

FUNCTION fm1p(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fm1p
  REAL(DP), INTENT(IN) :: r,q03

  REAL(DP) :: tt

  tt=ABS(q03)+27.*r*r/4.
  fm1p=(9./4.)*r*tt**(-5./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r/SQRT(ABS(q03)))) &
       -(SQRT(3.)/2.)*SQRT(ABS(q03))*tt**(-5./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r/SQRT(ABS(q03))))
  
END FUNCTION fm1p

FUNCTION fm2(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fm2
  REAL(DP), INTENT(IN) :: r,q03

  REAL(DP) :: tt

  tt=ABS(q03)+27.*r*r/4.
  fm2=(2./SQRT(3.))*tt**(1./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r/SQRT(ABS(q03))))
  
END FUNCTION fm2


FUNCTION fm1(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fm1
  REAL(DP), INTENT(IN) :: r,q03

  REAL(DP) :: tt

  tt=ABS(q03)+27.*r*r/4.
  fm1=tt**(1./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r/SQRT(ABS(q03))))
  
END FUNCTION fm1

FUNCTION f02p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: f02p
  REAL(DP), INTENT(IN) :: r

  f02p=2.**(1./3.)/9. * (r**2)**(-1./3.)
END FUNCTION f02p


FUNCTION f01p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: f01p
  REAL(DP), INTENT(IN) :: r

  f01p=-4./9. * SIGN((ABS(r))**(-5./3.),r)
END FUNCTION f01p

FUNCTION f02(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: f02
  REAL(DP), INTENT(IN) :: r

  f02=2.**(1./3.)/3*SIGN((ABS(r))**(1./3.),r)
END FUNCTION f02

FUNCTION f01(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: f01
  REAL(DP), INTENT(IN) :: r

  f01=2.**(1./3.)*((1.5*SQRT(3.)*r)**2)**(-1./3.)

END FUNCTION f01
