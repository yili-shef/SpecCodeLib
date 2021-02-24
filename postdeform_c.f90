PROGRAM postdeform
  USE mconstant
  IMPLICIT NONE

  INTEGER,  PARAMETER :: nprtcl=1024*1024
  INTEGER,  PARAMETER :: if1=1,if2=20
  REAL(DPP), PARAMETER :: dt=0.017

  
  REAL(DPP) :: Qre1,Rre1,Qre2,Rre2,Qre3,Rre3,Qres,Rres
  REAL(DPP),    DIMENSION(nprtcl)    :: Qre,Rre,Q03,Qi
  REAL(DPP),    DIMENSION(3,3,nprtcl):: Apre,Bre,Bren,Bre1,Bre2
  REAL(DPP),    DIMENSION(3,3):: Bren1,Apre1,Bren2,Apre2,Bren3,Apre3,Brens,Apres
  REAL(DPP),    DIMENSION(3,nprtcl)  :: evle,evlen

  REAL(DPP), DIMENSION(3,3) :: evtr,Cre,Apt
  REAL(DPP), DIMENSION(3)   :: fv1,fv2

  INTEGER  :: i,ifile,id1,id2,id3,id4,ip,idum,matz,ierr,idnty(3,3),mm
  REAL(DPP) :: Qt,Rt,tmp1,tmp2
  REAL(DPP) :: fp2p,fp1p,fp2,fp1,fm2p,fm1p,fm2,fm1,f02p,f01p,f02,f01
  REAL(DPP) :: time

  REAL(DP) :: Qref(nprtcl), Rref(nprtcl), Apref(3,3,nprtcl)
  INTEGER :: VIsDouble,Debug
  CHARACTER(1) :: NULLCHR

  INTEGER :: TecIni, TecZne, TecDat, TecEnd

  VIsDouble=0
  Debug=1
  NULLCHR=CHAR(0)

  mm=40000

  idnty=0
  idnty(1,1)=1; idnty(2,2)=1; idnty(3,3)=1
  OPEN(20,FILE='./xBdata/Aini.dat',FORM='unformatted')
  READ(20) Qre,Rre,Apre
  CLOSE(20)
  
  DO ip=1,nprtcl
  Bren(:,:,ip)=idnty
  END DO



! This routine compares the results of the eigenvalues of the Cauchy-Green
! tensor obtained from Restricted Euler dynamics, either by solving the
! equations numerically, or by evaluating the analytic solutions.  
  Q03=27.*Rre*Rre/4. + Qre*Qre*Qre
  Qi=Qre
  DO ip=1,nprtcl
    Qt=Qre(ip); Rt=Rre(ip); Apt=Apre(:,:,ip)
    IF(Q03(ip) .GT. 0.) THEN
      Qt=Qt/(ABS(Q03(ip)))**(1./3.); Rt=Rt/(ABS(Q03(ip)))**.5; Apt=Apt/(ABS(Q03(ip)))**(1./6.)
!      tmp2=(2./3.)*fp2p(Rt)*Qt**2
!      tmp1=(2./3.)*fp1p(Rt)*Qt**2
      tmp2=(2./3.)*fp2p(Rt)
      tmp1=(2./3.)*fp1p(Rt)
      Bre1(:,:,ip)=idnty*tmp2-Apt*fp2(Rt)
      Bre2(:,:,ip)=Apt*fp1(Rt)-idnty*tmp1
      Bre1(:,:,ip)=Bre1(:,:,ip)/(fp1(Rt)*tmp2-fp2(Rt)*tmp1)
      Bre2(:,:,ip)=Bre2(:,:,ip)/(fp1(Rt)*tmp2-fp2(Rt)*tmp1)
    ELSE IF(Q03(ip) .LT. 0.) THEN
      Qt=Qt/(ABS(Q03(ip)))**(1./3.); Rt=Rt/(ABS(Q03(ip)))**.5; Apt=Apt/(ABS(Q03(ip)))**(1./6.)
      tmp2=(2./3.)*fm2p(Rt)*Qt**2
      tmp1=(2./3.)*fm1p(Rt)*Qt**2
      Bre1(:,:,ip)=idnty*tmp2-Apt*fm2(Rt)
      Bre2(:,:,ip)=Apt*fm1(Rt)-idnty*tmp1
      Bre1(:,:,ip)=Bre1(:,:,ip)/(fm1(Rt)*tmp2-fm2(Rt)*tmp1)
      Bre2(:,:,ip)=Bre2(:,:,ip)/(fm1(Rt)*tmp2-fm2(Rt)*tmp1)
    ELSE
      Qt=Qt/(ABS(Qi(ip))); Rt=Rt/(ABS(Qi(ip)))**1.5; Apt=Apt/(ABS(Qi(ip)))**.5
      tmp2=(2./3.)*f02p(Rt)*Qt**2
      tmp1=(2./3.)*f01p(Rt)*Qt**2
      Bre1(:,:,ip)=idnty*tmp2-Apt*f02(Rt)
      Bre2(:,:,ip)=Apt*f01(Rt)-idnty*tmp1
      Bre1(:,:,ip)=Bre1(:,:,ip)/(f01(Rt)*tmp2-f02(Rt)*tmp1)
      Bre2(:,:,ip)=Bre2(:,:,ip)/(f01(Rt)*tmp2-f02(Rt)*tmp1)
    END IF
  END DO

  time=0. 
  DO ifile=if1,if2


    DO ip=1,nprtcl
      Rt=Rre(ip)
      IF(Q03(ip) .GT. 0.) THEN
        Rt=Rt/(ABS(Q03(ip)))**.5
        Bre(:,:,ip)=Bre1(:,:,ip)*fp1(Rt)+Bre2(:,:,ip)*fp2(Rt)
      ELSE IF(Q03(ip) .LT. 0.) THEN
        Rt=Rt/(ABS(Q03(ip)))**.5
        Bre(:,:,ip)=Bre1(:,:,ip)*fm1(Rt)+Bre2(:,:,ip)*fm2(Rt)
      ELSE
        Rt=Rt/(ABS(Qi(ip)))**1.5
        Bre(:,:,ip)=Bre1(:,:,ip)*f01(Rt)+Bre2(:,:,ip)*f02(Rt)
      END IF
      Cre=MATMUL(Bre(:,:,ip),TRANSPOSE(Bre(:,:,ip)))
      matz=5
      CALL rs(3,3,Cre,evle(:,ip),matz,evtr,fv1,fv2,ierr)
    END DO 

    DO ip=1,nprtcl
      Cre=MATMUL(Bren(:,:,ip),TRANSPOSE(Bren(:,:,ip)))
      matz=5
      CALL rs(3,3,Cre,evlen(:,ip),matz,evtr,fv1,fv2,ierr)
    END DO

    IF (ifile .EQ. if2) THEN
            DO ip=1,mm
              IF (Q03(ip) .GT. 0.) THEN
                      WRITE(15,*) LOG(evle(1,ip)),LOG(evlen(1,ip))
              ELSE IF (Q03(ip) .LT. 0.) THEN
                      WRITE(16,*) LOG(evle(1,ip)),LOG(evlen(1,ip))
              ELSE
                      WRITE(17,*) LOG(evle(1,ip)),LOG(evlen(1,ip))
              END IF
            END DO
    END IF


    DO ip=1,nprtcl
    Qre1=Qre(ip)+.5*dt*(-3.*Rre(ip))
    Rre1=Rre(ip)+.5*dt*(2.*Qre(ip)*Qre(ip)/3.)
      Bren1=Bren(:,:,ip)+.5*dt*MATMUL(Bren(:,:,ip),Apre(:,:,ip))
      Apre1=Apre(:,:,ip)+.5*dt*(-MATMUL(Apre(:,:,ip),Apre(:,:,ip))-(2./3.)*idnty*Qre(ip))

    Qre2=Qre(ip)+.5*dt*(-3.*Rre1)
    Rre2=Rre(ip)+.5*dt*(2.*Qre1*Qre1/3.)
      Bren2=Bren(:,:,ip)+.5*dt*MATMUL(Bren1,Apre1)
      Apre2=Apre(:,:,ip)+.5*dt*(-MATMUL(Apre1,Apre1)-(2./3.)*idnty*Qre1)
    
    Qre3=Qre(ip)+dt*(-3.*Rre2)
    Rre3=Rre(ip)+dt*(2.*Qre2*Qre2/3.)
      Bren3=Bren(:,:,ip)+dt*MATMUL(Bren2,Apre2)
      Apre3=Apre(:,:,ip)+dt*(-MATMUL(Apre2,Apre2)-(2./3.)*idnty*Qre2)
    
    Qres=Qre(ip)+(dt/6.)*(-3.*Rre(ip))+(2.*dt/3.)*(-3.*(Rre1+Rre2)/2.)+(dt/6.)*(-3.*Rre3)
    Rres=Rre(ip)+(dt/6.)*(2.*Qre(ip)*Qre(ip)/3.)+(2.*dt/3.)*(2./3.*(Qre1+Qre2)**2/4.)+(dt/6.)*(2.*Qre3*Qre3/3.)
      Brens=Bren(:,:,ip)+(dt/6.)*MATMUL(Bren(:,:,ip),Apre(:,:,ip)) &
                   +(2.*dt/3.)*MATMUL((Bren1+Bren2)/2.,(Apre1+Apre2)/2.) &
                   +(dt/6.)*MATMUL(Bren3,Apre3)
      Apres=Apre(:,:,ip)+(dt/6.)*(-MATMUL(Apre(:,:,ip),Apre(:,:,ip))-(2./3.)*idnty*Qre(ip)) &
                   +(2.*dt/3.)*(-MATMUL(.5*(Apre1+Apre2),.5*(Apre1+Apre2)) & 
                                -(2./3.)*idnty*(Qre1+Qre2)/2.)   &
                   +(dt/6.)*(-MATMUL(Apre3,Apre3)-(2./3.)*idnty*Qre3)
    Qre(ip)=Qres; Rre(ip)=Rres; Bren(:,:,ip)=Brens; Apre(:,:,ip)=Apres
    END DO

    time=time+dt
  END DO

  WRITE(*,*) 'finished'

END PROGRAM postdeform      


FUNCTION fp2p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fp2p
  REAL(DPP), INTENT(IN)  :: r

  REAL(DPP) :: tt1,tt2,tt3,tt4

  tt1=1.+3.*SQRT(3.)/2.*r
  tt2=1.-3.*SQRT(3.)/2.*r
!  tt3=SIGN((ABS(tt1))**(1./3.),tt1)
!  tt4=SIGN((ABS(tt2))**(1./3.),tt2)
  IF(tt1 .GE. 0.) THEN
          tt3=tt1**(1./3.)
  ELSE
          tt3=-(-tt1)**(1./3.)
  END IF
  IF(tt2 .GE. 0.) THEN
          tt4=tt2**(1./3.)
  ELSE
          tt4=-(-tt2)**(1./3.)
  END IF
            
  fp2p=.5*(tt3**2+tt4**2)
END FUNCTION fp2p

FUNCTION fp1p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fp1p
  REAL(DPP), INTENT(IN)  :: r

  REAL(DPP) :: tt1,tt2,tt3,tt4

  tt1=1.+3.*SQRT(3.)/2.*r
  tt2=1.-3.*SQRT(3.)/2.*r
!  tt3=SIGN((ABS(tt1))**(1./3.),tt1)
!  tt4=SIGN((ABS(tt2))**(1./3.),tt2)
  IF(tt1 .GE. 0.) THEN
          tt3=tt1**(1./3.)
  ELSE
          tt3=-(-tt1)**(1./3.)
  END IF
  IF(tt2 .GE. 0.) THEN
          tt4=tt2**(1./3.)
  ELSE
          tt4=-(-tt2)**(1./3.)
  END IF
            

  fp1p=.25*SQRT(3.)*(tt4**2-tt3**2)
END FUNCTION fp1p

FUNCTION fp2(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fp2
  REAL(DPP), INTENT(IN)  :: r
  REAL(DPP) :: tt1,tt2,tt3,tt4

  tt1=1.+3.*SQRT(3.)/2.*r
  tt2=1.-3.*SQRT(3.)/2.*r
  IF(tt1 .GE. 0.) THEN
          tt3=tt1**(1./3.)
  ELSE
          tt3=-(-tt1)**(1./3.)
  END IF
  IF(tt2 .GE. 0.) THEN
          tt4=tt2**(1./3.)
  ELSE
          tt4=-(-tt2)**(1./3.)
  END IF
            
!  fp2=(1./SQRT(3.))*(SIGN((ABS(tt1))**(1./3.),tt1)-SIGN((ABS(tt2))**(1./3.),tt2))
  fp2=(1./SQRT(3.))*(tt3-tt4)
END FUNCTION fp2

FUNCTION fp1(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fp1
  REAL(DPP), INTENT(IN)  :: r
  REAL(DPP) :: tt1,tt2,tt3,tt4

  tt1=1.+3.*SQRT(3.)/2.*r
  tt2=1.-3.*SQRT(3.)/2.*r
  IF(tt1 .GE. 0.) THEN
          tt3=tt1**(1./3.)
  ELSE
          tt3=-(-tt1)**(1./3.)
  END IF
  IF(tt2 .GE. 0.) THEN
          tt4=tt2**(1./3.)
  ELSE
          tt4=-(-tt2)**(1./3.)
  END IF
            
  fp1=.5*(tt3+tt4)
END FUNCTION fp1

FUNCTION fm2p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fm2p
  REAL(DPP), INTENT(IN)  :: r
  
  REAL(DPP) :: tt

  tt=1.+27.*r*r/4.
  fm2p=1.5*SQRT(3.)*r*tt**(-5./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r)) &
       +tt**(-5./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r))
  
END FUNCTION fm2p

FUNCTION fm1p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fm1p
  REAL(DPP), INTENT(IN) :: r

  REAL(DPP) :: tt

  tt=1.+27.*r*r/4.
  fm1p=(9./4.)*r*tt**(-5./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r)) &
       -(SQRT(3.)/2.)*tt**(-5./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r))
  
END FUNCTION fm1p

FUNCTION fm2(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fm2
  REAL(DPP), INTENT(IN) :: r

  REAL(DPP) :: tt

  tt=1.+27.*r*r/4.
  fm2=(2./SQRT(3.))*tt**(1./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r))
  
END FUNCTION fm2


FUNCTION fm1(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fm1
  REAL(DPP), INTENT(IN) :: r

  REAL(DPP) :: tt

  tt=1.+27.*r*r/4.
  fm1=tt**(1./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r))
  
END FUNCTION fm1

FUNCTION f02p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: f02p
  REAL(DPP), INTENT(IN) :: r

  f02p=2.**(1./3.)/9. * (r**2)**(-1./3.)
END FUNCTION f02p


FUNCTION f01p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: f01p
  REAL(DPP), INTENT(IN) :: r

  f01p=-4./9. * SIGN((ABS(r))**(-5./3.),r)
END FUNCTION f01p

FUNCTION f02(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: f02
  REAL(DPP), INTENT(IN) :: r

  f02=2.**(1./3.)/3*SIGN((ABS(r))**(1./3.),r)
END FUNCTION f02

FUNCTION f01(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: f01
  REAL(DPP), INTENT(IN) :: r

  f01=2.**(1./3.)*((1.5*SQRT(3.)*r)**2)**(-1./3.)

END FUNCTION f01
