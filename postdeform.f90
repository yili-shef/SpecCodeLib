PROGRAM postdeform
  USE mconstant
  IMPLICIT NONE

  INTEGER,  PARAMETER :: nprtcl=1024*1024
  INTEGER,  PARAMETER :: if1=1,if2=20
  REAL(DP), PARAMETER :: dt=0.017

  REAL(DP),    DIMENSION(nprtcl)    :: Q,R,Qre,Rre,Qres,Rres,Q0,Qi
  REAL(DP),    DIMENSION(3,3,nprtcl):: B,Ap,Bre,Bre1,Bre2

  REAL(DP), DIMENSION(3,3) :: C,evtr,Cre
  REAL(DP), DIMENSION(3)   :: evle,fv1,fv2

  INTEGER  :: i,ifile,id1,id2,id3,id4,ip,idum,matz,ierr,idnty(3,3)
  REAL(DP) :: detB,detC,time,sdetB,sdetC,slogev1,slogev2,slogev3,sQ,sR,Qt,Rt
  REAL(DP) :: fp2p,fp1p,fp2,fp1,fm2p,fm1p,fm2,fm1,f02p,f01p,f02,f01
  
  idnty=0
  idnty(1,1)=1; idnty(2,2)=1; idnty(3,3)=1
  OPEN(18,FILE='./xBdata/x.dat',FORM='unformatted')
  OPEN(19,FILE='./xBdata/B.dat',FORM='unformatted')
  OPEN(20,FILE='./xBdata/A.dat',FORM='unformatted')
  OPEN(16,FILE='time.dat')

  time=0.
  DO ifile=if1,if2
    READ(19) B

    slogev1=0.; slogev2=0.; slogev3=0.
    DO ip=1,nprtcl

    !    CALL dtmnnt(B(:,:,ip),3,detB) !Check if the determinant of B equals 1.
    C=MATMUL(B(:,:,ip),TRANSPOSE(B(:,:,ip)))
    CALL rs(3,3,C,evle,matz,evtr,fv1,fv2,ierr)
    slogev1=slogev1+LOG(evle(1))
    slogev2=slogev2+LOG(evle(2))
    slogev3=slogev3+LOG(evle(3))

    END DO
    slogev1=slogev1/REAL(nprtcl)
    slogev2=slogev2/REAL(nprtcl)
    slogev3=slogev3/REAL(nprtcl)

    WRITE(16,'(20E12.3)') time,slogev1,slogev2,slogev3
    time=time+dt
    
  END DO
  CLOSE(16)
  CLOSE(18)
  CLOSE(19)
  CLOSE(20)

  OPEN(20,FILE='./xBdata/A.dat',FORM='unformatted')
  READ(20) Q,R,Ap
  CLOSE(20)
  
  Qre=Q; Rre=R
  OPEN(15,FILE='./xBdata/QRre.dat',FORM='unformatted')
  DO ifile=if1, if2 
    WRITE(15) Qre,Rre
    Qres=Qre+dt*(-3.*Rre)
    Rres=Rre+dt*(2.*Qre*Qre/3.)
    Qre=.5*(Qre+Qres)+.5*dt*(-3.*Rres)
    Rre=.5*(Rre+Rres)+.5*dt*(2.*Qres*Qres/3.)
  END DO
  CLOSE(15)

  OPEN(15,FILE='./xBdata/QRre.dat',FORM='unformatted')
    READ(15) Qre,Rre
  CLOSE(15)
  Qi=Qre
  Q0=(27.*Rre*Rre/4. + Qre*Qre*Qre)**(1./3.)
  time=0. 
  DO ip=1,nprtcl
    IF(Q0(ip) .GT. 0.) THEN
      Qt=Qre(ip)/Q0(ip); Rt=Rre(ip)/(Q0(ip))**1.5
      Bre1(:,:,ip)=idnty*Qt**2*fp2p(Rt)-1.5*Ap(:,:,ip)/SQRT(Q0(ip))*fp2(Rt)
      Bre2(:,:,ip)=1.5*Ap(:,:,ip)/SQRT(Q0(ip))*fp1(Rt)-idnty*Qt**2*fp1p(Rt)
    ELSE IF(Q0(ip) .LT. 0.) THEN
      Qt=-Qre(ip)/Q0(ip); Rt=Rre(ip)/(ABS(Q0(ip)))**1.5
      Bre1(:,:,ip)=idnty*Qt**2*fm2p(Rt)-1.5*Ap(:,:,ip)/SQRT(ABS(Q0(ip)))*fm2(Rt)
      Bre2(:,:,ip)=1.5*Ap(:,:,ip)/SQRT(ABS(Q0(ip)))*fm1(Rt)-idnty*Qt**2*fm1p(Rt)
    ELSE
      Qt=Qre(ip)/(ABS(Qi(ip))+TINY); Rt=Rre(ip)/((ABS(Qi(ip)))**1.5+TINY)
      Bre1(:,:,ip)=idnty*Qt**2*f02p(Rt)-1.5*Ap(:,:,ip)/SQRT(ABS(Qi(ip)))*f02(Rt)
      Bre2(:,:,ip)=1.5*Ap(:,:,ip)/SQRT(ABS(Qi(ip))+TINY)*f01(Rt)-idnty*Qt**2*f01p(Rt)
    END IF
  END DO

  OPEN(16,FILE='./xBdata/QRre.dat',FORM='unformatted')
  OPEN(15,FILE='time_re.dat')
  DO ifile=if1,if2
    READ(16) Qre,Rre

    sdetB=0.; sdetC=0.
    slogev1=0.; slogev2=0.; slogev3=0.
    DO ip=1,nprtcl
      IF(Q0(ip) .GT. 0.) THEN
        Rt=Rre(ip)/(Q0(ip))**1.5
        Bre(:,:,ip)=Bre1(:,:,ip)*fp1(Rt)+Bre2(:,:,ip)*fp2(Rt)
      ELSE IF(Q0(ip) .LT. 0.) THEN
        Rt=Rre(ip)/(ABS(Q0(ip)))**1.5
        Bre(:,:,ip)=Bre1(:,:,ip)*fm1(Rt)+Bre2(:,:,ip)*fm2(Rt)
      ELSE
        Rt=Rre(ip)/(ABS(Qi(ip)))**1.5
        Bre(:,:,ip)=Bre1(:,:,ip)*f01(Rt)+Bre2(:,:,ip)*f02(Rt)
      END IF
      CALL dtmnnt(Bre(:,:,ip),3,detB)
      Cre=MATMUL(Bre(:,:,ip),TRANSPOSE(Bre(:,:,ip)))
      CALL dtmnnt(Cre,3,detC)
      CALL rs(3,3,Cre,evle,matz,evtr,fv1,fv2,ierr)
      sdetB=sdetB+detB
      sdetC=sdetC+detC
      slogev1=slogev1+LOG(evle(1))
      slogev2=slogev2+LOG(evle(2))
      slogev3=slogev3+LOG(evle(3))
    END DO 
    sdetB=sdetB/REAL(nprtcl)
    sdetC=sdetC/REAL(nprtcl)
    slogev1=slogev1/REAL(nprtcl)
    slogev2=slogev2/REAL(nprtcl)
    slogev3=slogev3/REAL(nprtcl)
    sQ=SUM(Qre)/REAL(nprtcl)
    sR=SUM(Rre)/REAL(nprtcl)
    WRITE(15,'(20E12.3)') time,sdetB,sdetC,slogev1,slogev2,slogev3,sQ,sR

    time=time+dt
  END DO
  CLOSE(15)
  CLOSE(16)

  WRITE(*,*) 'finished'

END PROGRAM postdeform      


FUNCTION fp2p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fp2p
  REAL(DP), INTENT(IN)  :: r

  fp2p=.5*((1.+3.*SQRT(3.)/2.*r)**(-2./3.)+(1.-3.*SQRT(3.)/2.*r)**(-2./3.))
END FUNCTION fp2p

FUNCTION fp1p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fp1p
  REAL(DP), INTENT(IN)  :: r

  fp1p=.25*SQRT(3.)*((1.+3.*SQRT(3.)/2.*r)**(-2./3.)-(1.-3.*SQRT(3.)/2.*r)**(-2./3.))
END FUNCTION fp1p

FUNCTION fp2(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fp2
  REAL(DP), INTENT(IN)  :: r

  fp2=(1./SQRT(3.))*((1.+3.*SQRT(3.)/2.*r)**(1./3.)-(1.-3.*SQRT(3.)/2.*r)**(1./3.))
END FUNCTION fp2

FUNCTION fp1(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fp1
  REAL(DP), INTENT(IN)  :: r

  fp1=.5*((1.+3.*SQRT(3.)/2.*r)**(1./3.)+(1.-3.*SQRT(3.)/2.*r)**(1./3.))
END FUNCTION fp1

FUNCTION fm2p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fm2p
  REAL(DP), INTENT(IN)  :: r

  fm2p=1.5*SQRT(3.)*r*(1.+27.*r*r/4.)**(-5./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r)) &
       +(1.+27.*r*r/4.)**(-5./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r))
  
END FUNCTION fm2p

FUNCTION fm1p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fm1p
  REAL(DP), INTENT(IN) :: r


  fm1p=(9./4.)*r*(1.+27.*r*r/4.)**(-5./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r)) &
       -(SQRT(3.)/2.)*(1.+27.*r*r/4.)**(-5./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r))
  
END FUNCTION fm1p

FUNCTION fm2(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fm2
  REAL(DP), INTENT(IN) :: r


  fm2=(2./SQRT(3.))*(1.+27.*r*r/4.)**(1./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r))
  
END FUNCTION fm2


FUNCTION fm1(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: fm1
  REAL(DP), INTENT(IN) :: r


  fm1=(1.+27.*r*r/4.)**(1./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r))
  
END FUNCTION fm1

FUNCTION f02p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: f02p
  REAL(DP), INTENT(IN) :: r

  f02p=2.**(1./3.)/9. * r**(-2./3.)
END FUNCTION f02p


FUNCTION f01p(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: f01p
  REAL(DP), INTENT(IN) :: r

  f01p=-4./9. * r**(-5./3.)
END FUNCTION f01p

FUNCTION f02(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: f02
  REAL(DP), INTENT(IN) :: r

  f02=2.**(1./3.)/3*r**(1./3.)
END FUNCTION f02

FUNCTION f01(r)
  USE mconstant
  IMPLICIT NONE

  REAL(DP) :: f01
  REAL(DP), INTENT(IN) :: r

  f01=2.**(1./3.)*(1.5*SQRT(3.)*r)**(-2./3.)

END FUNCTION f01
