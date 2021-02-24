PROGRAM eigen
  USE mconstant
  IMPLICIT NONE

  INTEGER,  PARAMETER :: nprtcl=1024*1024
  INTEGER,  PARAMETER :: if1=1,if2=20 !40 for filtered case, 20 for DNS
  REAL(DP), PARAMETER :: dt=0.017

  REAL(DP), DIMENSION(3,3,nprtcl) :: B,Bre,Bre1,Bre2,Ap,evtr,evtrre
  REAL(DP), DIMENSION(3,nprtcl)   :: evle,evlere
  REAL(DP), DIMENSION(nprtcl)     :: Q03,Qre,Rre,Q,R,Qi
  INTEGER,  DIMENSION(nprtcl)     :: marker
  REAL(SP), DIMENSION(3,3,nprtcl) :: Apf
  REAL(SP), DIMENSION(nprtcl)     :: trash

  REAL(DP), DIMENSION(3,3) :: C,Cre,Apt
  REAL(DP), DIMENSION(3)   :: fv1,fv2

  INTEGER  :: i,ii,kk,ifile,ip,idum,matz,ierr,mm,nn,idnty(3,3)
  REAL(DP) :: time,Qt,Rt
  REAL(DP) :: fp2p,fp1p,fp2,fp1,fm2p,fm1p,fm2,fm1,f02p,f01p,f02,f01
  REAL(DP) :: tmp1, tmp2
  
  idnty=0
  idnty(1,1)=1; idnty(2,2)=1; idnty(3,3)=1

  OPEN(20,FILE='./xBdata/Are.dat',FORM='unformatted')
  READ(20) Qre,Rre
  CLOSE(20)
  OPEN(20,FILE='./xBdata/A.dat',FORM='unformatted')
  READ(20) trash, trash, Apf
  CLOSE(20)
  Ap=REAL(Apf,DP)

  Q03=27.*Rre*Rre/4. + Qre*Qre*Qre
  Qi=Qre
  DO ip=1,nprtcl
    Qt=Qre(ip); Rt=Rre(ip); Apt=Ap(:,:,ip)
    IF(Q03(ip) .GT. 0.) THEN
      Qt=Qt/(ABS(Q03(ip)))**(1./3.); Rt=Rt/(ABS(Q03(ip)))**.5; Apt=Apt/(ABS(Q03(ip)))**(1./6.)
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


  
  OPEN(23,FILE='./xBdata/evle_re.dat',FORM='unformatted')
  OPEN(22,FILE='./xBdata/evle_dns.dat',FORM='unformatted')
  OPEN(21,FILE='./xBdata/evtr_re.dat',FORM='unformatted')
  OPEN(20,FILE='./xBdata/evtr_dns.dat',FORM='unformatted')
  OPEN(19,FILE='./xBdata/B.dat',FORM='unformatted')
  OPEN(18,FILE='./xBdata/Are.dat',FORM='unformatted')


  time=0.
  DO ifile=if1,if2
    READ(19) B
    READ(18) Qre,Rre


    DO ip=1,nprtcl


      C=MATMUL(TRANSPOSE(B(:,:,ip)),B(:,:,ip))
      matz=5
      CALL rs(3,3,C,evle(:,ip),matz,evtr(:,:,ip),fv1,fv2,ierr)

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
      Cre=MATMUL(TRANSPOSE(Bre(:,:,ip)),Bre(:,:,ip))
      matz=5
      CALL rs(3,3,Cre,evlere(:,ip),matz,evtrre(:,:,ip),fv1,fv2,ierr)

    END DO

    time=time+dt
    
    WRITE(20) evtr
    WRITE(21) evtrre
    WRITE(22) evle
    WRITE(23) evlere
  END DO
  CLOSE(18)
  CLOSE(19)
  CLOSE(20)
  CLOSE(21)
  CLOSE(22)
  CLOSE(23)


  WRITE(*,*) 'finished'

END PROGRAM eigen      


FUNCTION fp2p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fp2p
  REAL(DP), INTENT(IN)  :: r

  REAL(DP) :: tt1,tt2,tt3,tt4

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

  REAL(DP) :: fp1p
  REAL(DP), INTENT(IN)  :: r

  REAL(DP) :: tt1,tt2,tt3,tt4

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

  REAL(DP) :: fp2
  REAL(DP), INTENT(IN)  :: r
  REAL(DP) :: tt1,tt2,tt3,tt4

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

  REAL(DP) :: fp1
  REAL(DP), INTENT(IN)  :: r
  REAL(DP) :: tt1,tt2,tt3,tt4

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

  REAL(DP) :: fm2p
  REAL(DP), INTENT(IN)  :: r
  
  REAL(DP) :: tt

  tt=1.+27.*r*r/4.
  fm2p=1.5*SQRT(3.)*r*tt**(-5./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r)) &
       +tt**(-5./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r))
  
END FUNCTION fm2p

FUNCTION fm1p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fm1p
  REAL(DP), INTENT(IN) :: r

  REAL(DP) :: tt

  tt=1.+27.*r*r/4.
  fm1p=(9./4.)*r*tt**(-5./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r)) &
       -(SQRT(3.)/2.)*tt**(-5./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r))
  
END FUNCTION fm1p

FUNCTION fm2(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fm2
  REAL(DP), INTENT(IN) :: r

  REAL(DP) :: tt

  tt=1.+27.*r*r/4.
  fm2=(2./SQRT(3.))*tt**(1./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r))
  
END FUNCTION fm2


FUNCTION fm1(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fm1
  REAL(DP), INTENT(IN) :: r

  REAL(DP) :: tt

  tt=1.+27.*r*r/4.
  fm1=tt**(1./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r))
  
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
