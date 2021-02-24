PROGRAM jpdf
  USE mconstant
  IMPLICIT NONE

  INTEGER,  PARAMETER :: nprtcl=1024*1024
  INTEGER,  PARAMETER :: if1=1,if2=68 !40 for filtered case, 20 for DNS
  REAL(DP), PARAMETER :: dt=0.005

  INTEGER,  PARAMETER :: nev=100
  INTEGER,  DIMENSION(nev,nev) :: pev1,pev2,pev3
  INTEGER,  DIMENSION(nev) :: mpev1,mpev2,mpev3,mpev1re,mpev2re,mpev3re
  REAL(DP), DIMENSION(nev,nev) :: rpev1,rpev2,rpev3
  REAL(DP), DIMENSION(nev) :: rmpev1,rmpev2,rmpev3,rmpev1re,rmpev2re,rmpev3re
  REAL(DP), DIMENSION(nev) :: ev1,ev2,ev3

  REAL(DP), DIMENSION(3,3,nprtcl) :: evtr,evtrre
  REAL(DP), DIMENSION(3,nprtcl)  :: evle,evlere
  INTEGER,  DIMENSION(nprtcl)     :: marker

  INTEGER  :: i,ii,kk,ifile,ip,mm,nn
  REAL(DP) :: lev1,lev2,lev3,levre1,levre2,levre3
  REAL(DP) :: time,ev1b,ev2b,ev3b,ev,evre
  REAL(DP) :: mev1s,mev1e,mev2s,mev2e,mev3s,mev3e
  REAL(DP) :: xev1s,xev1e,xev2s,xev2e,xev3s,xev3e
  REAL(DP) :: mev1,xev1,mev2,xev2,mev3,xev3
  REAL(DP) :: meanev1,meanev2,meanev3,rmsev1,rmsev2,rmsev3
  REAL(DP) :: meanev1re,meanev2re,meanev3re,rmsev1re,rmsev2re,rmsev3re
  REAL(DP) :: rhoev1,rhoev2,rhoev3
  REAL(DP) :: tmp1, tmp2
  
!  mev1e=-40.; xev1e=.2; mev2e=-20.1; xev2e=20.1; mev3e=-.2; xev3e=40.  !For filtered
  mev1e=-20.; xev1e=.2; mev2e=-10.1; xev2e=10.1; mev3e=-.2; xev3e=20. !For DNS
  mev1s=-4;   xev1s=.2; mev2s=-2.1;  xev2s=2.1;  mev3s=-.2; xev3s=4

  OPEN(18,FILE='./xBdata/evle_re.dat',FORM='unformatted')
  DO ifile=if1,if2
    READ(18) evlere
  END DO
  CLOSE(18)

!  ev1b=(xev1e-mev1e)/nev; ev2b=(xev2e-mev2e)/nev; ev3b=(xev3e-mev3e)/nev
!  DO ip=1,nprtcl

!    levre1=(LOG(evlere(1,ip))-mev1e)/ev1b
!    levre2=(LOG(evlere(2,ip))-mev2e)/ev2b
!    levre3=(LOG(evlere(3,ip))-mev3e)/ev3b
!    IF ((levre1 .GE. 0. .AND. levre1 .LT. nev) .AND. &
!        (levre2 .GE. 0. .AND. levre2 .LT. nev) .AND. &
!        (levre3 .GE. 0. .AND. levre3 .LT. nev)) THEN
!      marker(ip)=1
!    ELSE
!      marker(ip)=0
!    END IF

!  END DO

!  OPEN(15,FILE='marker.dat')
!  WRITE(15,*) marker
!  CLOSE(15)
!  STOP
  OPEN(15,FILE='marker.dat')
  READ(15,*) marker
  CLOSE(15)
  
  OPEN(21,FILE='./xBdata/evle_re.dat',FORM='unformatted')
  OPEN(20,FILE='./xBdata/evle_dns.dat',FORM='unformatted')
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


  time=0.
  DO ifile=if1,if2

    READ(21) evlere
    READ(20) evle

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


      lev1  =(LOG(evle(1,ip))  -mev1)/ev1b
      levre1=(LOG(evlere(1,ip))-mev1)/ev1b
      lev2  =(LOG(evle(2,ip))  -mev2)/ev2b
      levre2=(LOG(evlere(2,ip))-mev2)/ev2b
      lev3  =(LOG(evle(3,ip))  -mev3)/ev3b
      levre3=(LOG(evlere(3,ip))-mev3)/ev3b
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

              meanev1=meanev1+LOG(evle(1,ip))
              meanev2=meanev2+LOG(evle(2,ip))
              meanev3=meanev3+LOG(evle(3,ip))
              meanev1re=meanev1re+LOG(evlere(1,ip))
              meanev2re=meanev2re+LOG(evlere(2,ip))
              meanev3re=meanev3re+LOG(evlere(3,ip))

              rmsev1=rmsev1+(LOG(evle(1,ip)))**2
              rmsev2=rmsev2+(LOG(evle(2,ip)))**2
              rmsev3=rmsev3+(LOG(evle(3,ip)))**2
              rmsev1re=rmsev1re+(LOG(evlere(1,ip)))**2
              rmsev2re=rmsev2re+(LOG(evlere(2,ip)))**2
              rmsev3re=rmsev3re+(LOG(evlere(3,ip)))**2

              rhoev1=rhoev1+LOG(evle(1,ip))*LOG(evlere(1,ip))
              rhoev2=rhoev2+LOG(evle(2,ip))*LOG(evlere(2,ip))
              rhoev3=rhoev3+LOG(evle(3,ip))*LOG(evlere(3,ip))

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


    time=time+dt
    
  END DO
  CLOSE(12)
  CLOSE(13)
  CLOSE(14)
  CLOSE(15)
  CLOSE(16)
  CLOSE(17)
  CLOSE(18)
  CLOSE(19)
  CLOSE(20)
  CLOSE(21)


  WRITE(*,*) 'finished'

END PROGRAM jpdf      

