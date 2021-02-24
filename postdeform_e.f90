PROGRAM postdeform
  USE mconstant
  IMPLICIT NONE

  INTEGER,  PARAMETER :: nprtcl=1024*1024
  INTEGER,  PARAMETER :: if1=1,if2=20
  REAL(DPP), PARAMETER :: dt=0.017

  INTEGER,  PARAMETER :: nev=70
  INTEGER,  DIMENSION(nev,nev) :: pev1,pev2,pev3
  REAL(DPP), DIMENSION(nev,nev) :: rpev1,rpev2,rpev3
  INTEGER,  DIMENSION(nev) :: mpev1,mpev2,mpev3,mpev1re,mpev2re,mpev3re
  REAL(DPP), DIMENSION(nev) :: rmpev1,rmpev2,rmpev3,rmpev1re,rmpev2re,rmpev3re
  REAL(DPP), DIMENSION(nev) :: ev1,ev2,ev3
  
  REAL(DP), DIMENSION(3,3,nprtcl) :: Bf,Apf
  REAL(DP), DIMENSION(nprtcl) :: Qref,Rref
  REAL(DPP), DIMENSION(3,3,nprtcl) :: B,Bre,Bre1,Bre2,Ap
  REAL(DPP), DIMENSION(nprtcl)     :: Qre,Rre,Qres,Rres
  INTEGER,  DIMENSION(nprtcl)     :: marker
  REAL(DPP), DIMENSION(nprtcl)    :: Q03

  REAL(DPP), DIMENSION(3,3) :: C,Cre,evtr,evtrre,Apt
  REAL(DPP), DIMENSION(3)   :: evle,evlere,fv1,fv2

  INTEGER  :: i,ii,kk,ifile,id1,id2,id3,id4,ip,idum,matz,ierr,mm,nn,idnty(3,3)
  REAL(DPP) :: lev1,lev2,lev3,levre1,levre2,levre3
  REAL(DPP) :: time,ev1b,ev2b,ev3b,ev,evre
  REAL(DPP) :: mev1s,mev1e,mev2s,mev2e,mev3s,mev3e
  REAL(DPP) :: xev1s,xev1e,xev2s,xev2e,xev3s,xev3e
  REAL(DPP) :: mev1,xev1,mev2,xev2,mev3,xev3
  REAL(DPP) :: meanev1,meanev2,meanev3,rmsev1,rmsev2,rmsev3
  REAL(DPP) :: meanev1re,meanev2re,meanev3re,rmsev1re,rmsev2re,rmsev3re
  REAL(DPP) :: rhoev1,rhoev2,rhoev3,Qt,Rt
  REAL(DPP) :: fp2p,fp1p,fp2,fp1,fm2p,fm1p,fm2,fm1,f02p,f01p,f02,f01
  REAL(DPP) :: tmp1, tmp2
  
  mev1e=-14.;xev1e=.2;mev2e=-7.1;xev2e=7.1;mev3e=-.2;xev3e=14.
  mev1s=-.5;xev1s=.2;mev2s=-.35;xev2s=.35;mev3s=-.2;xev3s=.5

  idnty=0
  idnty(1,1)=1; idnty(2,2)=1; idnty(3,3)=1
  OPEN(20,FILE='./xBdata/A.dat',FORM='unformatted')
  READ(20) Qref,Rref,Apf
  CLOSE(20)
  Qre=DBLE(Qref); Rre=DBLE(Rref); Ap=DBLE(Apf)

  Q03=27.*Rre*Rre/4. + Qre*Qre*Qre
  DO ip=1,nprtcl
    Qt=Qre(ip); Rt=Rre(ip); Apt=Ap(:,:,ip)
    IF(Q03(ip) .GT. 0.) THEN
      tmp2=(2./3.)*fp2p(Rt,Q03(ip))*Qt**2
      tmp1=(2./3.)*fp1p(Rt,Q03(ip))*Qt**2
      Bre1(:,:,ip)=idnty*tmp2-Apt*fp2(Rt,Q03(ip))
      Bre2(:,:,ip)=Apt*fp1(Rt,Q03(ip))-idnty*tmp1
      Bre1(:,:,ip)=Bre1(:,:,ip)/(fp1(Rt,Q03(ip))*tmp2-fp2(Rt,Q03(ip))*tmp1)
      Bre2(:,:,ip)=Bre2(:,:,ip)/(fp1(Rt,Q03(ip))*tmp2-fp2(Rt,Q03(ip))*tmp1)
    ELSE IF(Q03(ip) .LT. 0.) THEN
      tmp2=(2./3.)*fm2p(Rt,Q03(ip))*Qt**2
      tmp1=(2./3.)*fm1p(Rt,Q03(ip))*Qt**2
      Bre1(:,:,ip)=idnty*tmp2-Apt*fm2(Rt,Q03(ip))
      Bre2(:,:,ip)=Apt*fm1(Rt,Q03(ip))-idnty*tmp1
      Bre1(:,:,ip)=Bre1(:,:,ip)/(fm1(Rt,Q03(ip))*tmp2-fm2(Rt,Q03(ip))*tmp1)
      Bre2(:,:,ip)=Bre2(:,:,ip)/(fm1(Rt,Q03(ip))*tmp2-fm2(Rt,Q03(ip))*tmp1)
    ELSE
      tmp2=(2./3.)*f02p(Rt)*Qt**2
      tmp1=(2./3.)*f01p(Rt)*Qt**2
      Bre1(:,:,ip)=idnty*tmp2-Apt*f02(Rt)
      Bre2(:,:,ip)=Apt*f01(Rt)-idnty*tmp1
      Bre1(:,:,ip)=Bre1(:,:,ip)/(f01(Rt)*tmp2-f02(Rt)*tmp1)
      Bre2(:,:,ip)=Bre2(:,:,ip)/(f01(Rt)*tmp2-f02(Rt)*tmp1)
    END IF
  END DO


  OPEN(18,FILE='./xBdata/QRre.dat',FORM='unformatted')
  DO ifile=if1,if2
    READ(18) Qref,Rref
  END DO
  Qre=DBLE(Qref); Rre=DBLE(Rref)
  CLOSE(18)

  ev1b=(xev1e-mev1e)/nev; ev2b=(xev2e-mev2e)/nev; ev3b=(xev3e-mev3e)/nev
  DO ip=1,nprtcl
    Rt=Rre(ip)
    IF(Q03(ip) .GT. 0.) THEN
      Bre(:,:,ip)=Bre1(:,:,ip)*fp1(Rt,Q03(ip))+Bre2(:,:,ip)*fp2(Rt,Q03(ip))
    ELSE IF(Q03(ip) .LT. 0.) THEN
      Bre(:,:,ip)=Bre1(:,:,ip)*fm1(Rt,Q03(ip))+Bre2(:,:,ip)*fm2(Rt,Q03(ip))
    ELSE
      Bre(:,:,ip)=Bre1(:,:,ip)*f01(Rt)+Bre2(:,:,ip)*f02(Rt)
    END IF
    Cre=MATMUL(Bre(:,:,ip),TRANSPOSE(Bre(:,:,ip)))
    matz=5
    CALL rs(3,3,Cre,evlere,matz,evtr,fv1,fv2,ierr)

    levre1=(LOG(evlere(1))-mev1e)/ev1b
    levre2=(LOG(evlere(2))-mev2e)/ev2b
    levre3=(LOG(evlere(3))-mev3e)/ev3b
    IF ((levre1 .GE. 0. .AND. levre1 .LT. nev) .AND. &
        (levre2 .GE. 0. .AND. levre2 .LT. nev) .AND. &
        (levre3 .GE. 0. .AND. levre3 .LT. nev)) THEN
      marker(ip)=1
    ELSE
      marker(ip)=0
    END IF

  END DO
  
  OPEN(19,FILE='./xBdata/B.dat',FORM='unformatted')
  OPEN(17,FILE='deform_rho.dat')
  OPEN(16,FILE='deform_jpdf1.dat')
  OPEN(15,FILE='deform_jpdf2.dat')
  OPEN(14,FILE='deform_jpdf3.dat')
  OPEN(13,FILE='pdf_time_re.dat')
  OPEN(12,FILE='pdf_time_dns.dat')

  WRITE(17,*)'VARIABLES = "time" "mev1" "mev2" "mev3" "sev1"' 
  WRITE(17,*)'"sev2" "sev3" "mev1re" "mev2re" "mev3re"'
  WRITE(17,*)'"sev1re" "sev2re" "sev3re"'
  WRITE(17,*)'"rhoev1" "rhoev2" "rhoev3"'
  WRITE(16,'(''VARIABLES = "ev1", "evre1", "jpdf1"'')')
  WRITE(15,'(''VARIABLES = "ev2", "evre2", "jpdf2"'')')
  WRITE(14,'(''VARIABLES = "ev3", "evre3", "jpdf3"'')')
  WRITE(13,'(''Variables = "ev1" "pev1" "ev2" "pev2" "ev3" "pev3"'')')
  WRITE(12,'(''Variables = "ev1" "pev1" "ev2" "pev2" "ev3" "pev3"'')')

  OPEN(18,FILE='./xBdata/QRre.dat',FORM='unformatted')
    READ(18) Qref,Rref
  CLOSE(18)
  Qre=DBLE(Qref); Rre=DBLE(Rref)

  time=0.
  DO ifile=if1,if2
    READ(19) Bf
    B=DBLE(Bf)

    xev1=xev1s+(xev1e-xev1s)/(if2-if1)*(ifile-if1)
    xev2=xev2s+(xev2e-xev2s)/(if2-if1)*(ifile-if1)
    xev3=xev3s+(xev3e-xev3s)/(if2-if1)*(ifile-if1)
    mev1=mev1s+(mev1e-mev1s)/(if2-if1)*(ifile-if1)
    mev2=mev2s+(mev2e-mev2s)/(if2-if1)*(ifile-if1)
    mev3=mev3s+(mev3e-mev3s)/(if2-if1)*(ifile-if1)

    ev1b=(xev1-mev1)/nev; ev2b=(xev2-mev2)/nev; ev3b=(xev3-mev3)/nev
    ev1=(/(mev1+(i-.5)*ev1b,i=1,nev)/)
    ev2=(/(mev2+(i-.5)*ev2b,i=1,nev)/)
    ev3=(/(mev3+(i-.5)*ev3b,i=1,nev)/)

    pev1=0; pev2=0; pev3=0
    mpev1=0; mpev2=0; mpev3=0
    mpev1re=0; mpev2re=0; mpev3re=0
    meanev1=0.; meanev1re=0.; rmsev1=0.; rmsev1re=0.
    meanev2=0.; meanev2re=0.; rmsev2=0.; rmsev2re=0.
    meanev3=0.; meanev3re=0.; rmsev3=0.; rmsev3re=0.
    rhoev1=0.; rhoev2=0.; rhoev3=0.
    ii=0
    DO ip=1,nprtcl

      IF (marker(ip) .EQ. 0) CYCLE

      C=MATMUL(B(:,:,ip),TRANSPOSE(B(:,:,ip)))
      matz=5
      CALL rs(3,3,C,evle,matz,evtr,fv1,fv2,ierr)

      Rt=Rre(ip)
      IF(Q03(ip) .GT. 0.) THEN
        Bre(:,:,ip)=Bre1(:,:,ip)*fp1(Rt,Q03(ip))+Bre2(:,:,ip)*fp2(Rt,Q03(ip))
      ELSE IF(Q03(ip) .LT. 0.) THEN
        Bre(:,:,ip)=Bre1(:,:,ip)*fm1(Rt,Q03(ip))+Bre2(:,:,ip)*fm2(Rt,Q03(ip))
      ELSE
        Bre(:,:,ip)=Bre1(:,:,ip)*f01(Rt)+Bre2(:,:,ip)*f02(Rt)
      END IF
      Cre=MATMUL(Bre(:,:,ip),TRANSPOSE(Bre(:,:,ip)))
      matz=5
      CALL rs(3,3,Cre,evlere,matz,evtr,fv1,fv2,ierr)

      lev1=(LOG(evle(1))-mev1)/ev1b
      levre1=(LOG(evlere(1))-mev1)/ev1b
      lev2=(LOG(evle(2))-mev2)/ev2b
      levre2=(LOG(evlere(2))-mev2)/ev2b
      lev3=(LOG(evle(3))-mev3)/ev3b
      levre3=(LOG(evlere(3))-mev3)/ev3b
      IF ((lev1 .GE. 0. .AND. lev1 .LT. nev) .AND. (levre1 .GE. 0. .AND. levre1 .LT. nev) .AND. &
          (lev2 .GE. 0. .AND. lev2 .LT. nev) .AND. (levre2 .GE. 0. .AND. levre2 .LT. nev) .AND. &
          (lev3 .GE. 0. .AND. lev3 .LT. nev) .AND. (levre3 .GE. 0. .AND. levre3 .LT. nev)) THEN
              mm=FLOOR(lev1)+1
              nn=FLOOR(levre1)+1
              mpev1(mm)=mpev1(mm)+1
              mpev1re(nn)=mpev1re(nn)+1
              pev1(mm,nn)=pev1(mm,nn)+1
              
              mm=FLOOR(lev2)+1
              nn=FLOOR(levre2)+1
              mpev2(mm)=mpev2(mm)+1
              mpev2re(nn)=mpev2re(nn)+1
              pev2(mm,nn)=pev2(mm,nn)+1
              
              mm=FLOOR(lev3)+1
              nn=FLOOR(levre3)+1
              mpev3(mm)=mpev3(mm)+1
              mpev3re(nn)=mpev3re(nn)+1
              pev3(mm,nn)=pev3(mm,nn)+1

              meanev1=meanev1+LOG(evle(1))
              meanev2=meanev2+LOG(evle(2))
              meanev3=meanev3+LOG(evle(3))
              meanev1re=meanev1re+LOG(evlere(1))
              meanev2re=meanev2re+LOG(evlere(2))
              meanev3re=meanev3re+LOG(evlere(3))

              rmsev1=rmsev1+(LOG(evle(1)))**2
              rmsev2=rmsev2+(LOG(evle(2)))**2
              rmsev3=rmsev3+(LOG(evle(3)))**2
              rmsev1re=rmsev1re+(LOG(evlere(1)))**2
              rmsev2re=rmsev2re+(LOG(evlere(2)))**2
              rmsev3re=rmsev3re+(LOG(evlere(3)))**2

              rhoev1=rhoev1+LOG(evle(1))*LOG(evlere(1))
              rhoev2=rhoev2+LOG(evle(2))*LOG(evlere(2))
              rhoev3=rhoev3+LOG(evle(3))*LOG(evlere(3))

              ii=ii+1
      END IF

    END DO
    WRITE(*,*) 'check', REAL(ii)/nprtcl
    rpev1=pev1/REAL(ii)
    rpev2=pev2/REAL(ii)
    rpev3=pev3/REAL(ii)
    rmpev1=mpev1/REAL(ii)
    rmpev2=mpev2/REAL(ii)
    rmpev3=mpev3/REAL(ii)
    rmpev1re=mpev1re/REAL(ii)
    rmpev2re=mpev2re/REAL(ii)
    rmpev3re=mpev3re/REAL(ii)
    
    rpev1=rpev1/ev1b/ev1b 
    rpev2=rpev2/ev2b/ev2b
    rpev3=rpev3/ev3b/ev3b

    rmpev1=rmpev1/ev1b
    rmpev2=rmpev2/ev2b
    rmpev3=rmpev3/ev3b

    rmpev1re=rmpev1re/ev1b
    rmpev2re=rmpev2re/ev2b
    rmpev3re=rmpev3re/ev3b

    meanev1=meanev1/REAL(ii)
    meanev2=meanev2/REAL(ii)
    meanev3=meanev3/REAL(ii)
    meanev1re=meanev1re/REAL(ii)
    meanev2re=meanev2re/REAL(ii)
    meanev3re=meanev3re/REAL(ii)
    
    rmsev1=rmsev1/REAL(ii)
    rmsev2=rmsev2/REAL(ii)
    rmsev3=rmsev3/REAL(ii)
    rmsev1re=rmsev1re/REAL(ii)
    rmsev2re=rmsev2re/REAL(ii)
    rmsev3re=rmsev3re/REAL(ii)

    rmsev1=rmsev1-meanev1*meanev1
    rmsev2=rmsev2-meanev2*meanev2
    rmsev3=rmsev3-meanev3*meanev3
    rmsev1re=rmsev1re-meanev1re*meanev1re
    rmsev2re=rmsev2re-meanev2re*meanev2re
    rmsev3re=rmsev3re-meanev3re*meanev3re

    rmsev1=SQRT(rmsev1)
    rmsev2=SQRT(rmsev2)
    rmsev3=SQRT(rmsev3)
    rmsev1re=SQRT(rmsev1re)
    rmsev2re=SQRT(rmsev2re)
    rmsev3re=SQRT(rmsev3re)

    rhoev1=rhoev1/REAL(ii)-meanev1*meanev1re
    rhoev2=rhoev2/REAL(ii)-meanev2*meanev2re
    rhoev3=rhoev3/REAL(ii)-meanev3*meanev3re

    rhoev1=rhoev1/(rmsev1*rmsev1re+TINY)
    rhoev2=rhoev2/(rmsev2*rmsev2re+TINY)
    rhoev3=rhoev3/(rmsev3*rmsev3re+TINY)

    WRITE(17,'(20E12.3)') time,meanev1,meanev2,meanev3,rmsev1,rmsev2,rmsev3, &
                          meanev1re,meanev2re,meanev3re,rmsev1re,rmsev2re,         &
                          rmsev3re,rhoev1,rhoev2,rhoev3

    WRITE(16,'(''ZONE T="t='',F8.4,''", I='',I6,'' J='',I6,'',F=POINT'')') time,nev,nev
    WRITE(15,'(''ZONE T="t='',F8.4,''", I='',I6,'' J='',I6,'',F=POINT'')') time,nev,nev
    WRITE(14,'(''ZONE T="t='',F8.4,''", I='',I6,'' J='',I6,'',F=POINT'')') time,nev,nev
    DO ii=1,nev
    DO i=1,nev
      WRITE(16,'(20E12.3)') ev1(i),ev1(ii),rpev1(i,ii)
      WRITE(15,'(20E12.3)') ev2(i),ev2(ii),rpev2(i,ii)
      WRITE(14,'(20E12.3)') ev3(i),ev3(ii),rpev3(i,ii)
    END DO
    END DO
    WRITE(13,*) 'ZONE T="t=',time,'"'
    WRITE(12,*) 'ZONE T="t=',time,'"'
    DO i=1,nev
      WRITE(13,'(20E12.3)') ev1(i),rmpev1re(i),ev2(i),rmpev2re(i),ev3(i),rmpev3re(i) 
      WRITE(12,'(20E12.3)') ev1(i),rmpev1(i),  ev2(i),rmpev2(i),  ev3(i),rmpev3(i) 
    END DO


    
    DO i=1,nprtcl
    Qres(i)=Qre(i)+dt*(-3.*Rre(i))
    Rres(i)=Rre(i)+dt*(2.*Qre(i)*Qre(i)/3.)
    Qre(i)=.5*(Qre(i)+Qres(i))+.5*dt*(-3.*Rres(i))
    Rre(i)=.5*(Rre(i)+Rres(i))+.5*dt*(2.*Qres(i)*Qres(i)/3.)
    END DO



    time=time+dt
    
  END DO
  CLOSE(12)
  CLOSE(13)
  CLOSE(14)
  CLOSE(15)
  CLOSE(16)
  CLOSE(17)
  CLOSE(19)


  WRITE(*,*) 'finished'

END PROGRAM postdeform      


FUNCTION fp2p(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fp2p
  REAL(DPP), INTENT(IN)  :: r,q03

  REAL(DPP) :: tt1,tt2,tt3,tt4

  tt1=SQRT(q03)+3.*SQRT(3.)/2.*r
  tt2=SQRT(q03)-3.*SQRT(3.)/2.*r
  tt3=SIGN((ABS(tt1))**(1./3.),tt1)
  tt4=SIGN((ABS(tt2))**(1./3.),tt2)
  fp2p=.5*(tt3**(-2)+tt4**(-2))
END FUNCTION fp2p

FUNCTION fp1p(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fp1p
  REAL(DPP), INTENT(IN)  :: r,q03

  REAL(DPP) :: tt1,tt2,tt3,tt4

  tt1=SQRT(q03)+3.*SQRT(3.)/2.*r
  tt2=SQRT(q03)-3.*SQRT(3.)/2.*r
  tt3=SIGN((ABS(tt1))**(1./3.),tt1)
  tt4=SIGN((ABS(tt2))**(1./3.),tt2)

  fp1p=.25*SQRT(3.)*(tt3**(-2)-tt4**(-2))
END FUNCTION fp1p

FUNCTION fp2(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fp2
  REAL(DPP), INTENT(IN)  :: r,q03
  REAL(DPP) :: tt1,tt2

  tt1=SQRT(q03)+3.*SQRT(3.)/2.*r
  tt2=SQRT(q03)-3.*SQRT(3.)/2.*r
  fp2=(1./SQRT(3.))*(SIGN((ABS(tt1))**(1./3.),tt1)-SIGN((ABS(tt2))**(1./3.),tt2))
END FUNCTION fp2

FUNCTION fp1(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fp1
  REAL(DPP), INTENT(IN)  :: r,q03
  REAL(DPP) :: tt1,tt2

  tt1=SQRT(q03)+3.*SQRT(3.)/2.*r
  tt2=SQRT(q03)-3.*SQRT(3.)/2.*r
  fp1=.5*(SIGN((ABS(tt1))**(1./3.),tt1)+SIGN((ABS(tt2))**(1./3.),tt2))
END FUNCTION fp1

FUNCTION fm2p(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fm2p
  REAL(DPP), INTENT(IN)  :: r,q03
  
  REAL(DPP) :: tt

  tt=ABS(q03)+27.*r*r/4.
  fm2p=1.5*SQRT(3.)*r*tt**(-5./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r/SQRT(ABS(q03)))) &
       +.5*SQRT(3.)*SQRT(ABS(q03))*tt**(-5./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r/SQRT(ABS(q03))))
  
END FUNCTION fm2p

FUNCTION fm1p(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fm1p
  REAL(DPP), INTENT(IN) :: r,q03

  REAL(DPP) :: tt

  tt=ABS(q03)+27.*r*r/4.
  fm1p=(9./4.)*r*tt**(-5./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r/SQRT(ABS(q03)))) &
       -(SQRT(3.)/2.)*SQRT(ABS(q03))*tt**(-5./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r/SQRT(ABS(q03))))
  
END FUNCTION fm1p

FUNCTION fm2(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fm2
  REAL(DPP), INTENT(IN) :: r,q03

  REAL(DPP) :: tt

  tt=ABS(q03)+27.*r*r/4.
  fm2=(2./SQRT(3.))*tt**(1./6.)*SIN((1./3.)*ATAN(1.5*SQRT(3.)*r/SQRT(ABS(q03))))
  
END FUNCTION fm2


FUNCTION fm1(r,q03)
  USE mconstant
  IMPLICIT NONE

  REAL(DPP) :: fm1
  REAL(DPP), INTENT(IN) :: r,q03

  REAL(DPP) :: tt

  tt=ABS(q03)+27.*r*r/4.
  fm1=tt**(1./6.)*COS((1./3.)*ATAN(1.5*SQRT(3.)*r/SQRT(ABS(q03))))
  
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
