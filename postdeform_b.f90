PROGRAM postdeform
  USE mconstant
  IMPLICIT NONE

  INTEGER,  PARAMETER :: nprtcl=1024*1024
  INTEGER,  PARAMETER :: if1=1,if2=2
  REAL(DPP), PARAMETER :: dt=0.017

  REAL(DP), DIMENSION(3,3,nprtcl) :: Apref,Bf
  REAL(DP), DIMENSION(nprtcl) :: Qref,Rref
  
  REAL(DPP), DIMENSION(3,3,nprtcl):: B
  REAL(DPP),    DIMENSION(nprtcl)    :: Qre,Rre,Qres,Rres
  REAL(DPP),    DIMENSION(3,3,nprtcl):: Apre,Bren,Brens,Apres
  REAL(DPP),    DIMENSION(3,nprtcl)  :: evle,evlen

  REAL(DPP), DIMENSION(3,3) :: evtr,Cre,C
  REAL(DPP), DIMENSION(3)   :: fv1,fv2

  INTEGER  :: i,ifile,ip,idum,matz,ierr,mm,idnty(3,3)
  REAL(DPP) :: time

  INTEGER :: VIsDouble,Debug
  CHARACTER(1) :: NULLCHR

  INTEGER :: TecIni, TecZne, TecDat, TecEnd

  VIsDouble=1
  Debug=1
  NULLCHR=CHAR(0)

  mm=40000

  idnty=0
  idnty(1,1)=1; idnty(2,2)=1; idnty(3,3)=1
  OPEN(20,FILE='./xBdata/A.dat',FORM='unformatted')
  READ(20) Qref,Rref,Apref
  CLOSE(20)
  Qre=DBLE(Qref); Rre=DBLE(Rref); Apre=DBLE(Apref)
  
  DO ip=1,nprtcl
  Bren(:,:,ip)=idnty
  END DO

! The purpose of this routine is to compare the eigenvalues of the Cauchy-Green
! tensor numerically calculated from the DNS data and Restricted Euler dynamics.
! The Restricted Euler equations for velocity gradient, coupled with the
! equations for the deformation tensor, are solved numerically, other than using
! the analytic solution for the deformation tensor.
  OPEN(19,FILE='./xBdata/B.dat',FORM='unformatted')

  time=0.
  DO ifile=if1,if2
    READ(19) Bf
    B=DBLE(Bf)

    DO ip=1,nprtcl

      C=MATMUL(B(:,:,ip),TRANSPOSE(B(:,:,ip)))
      matz=5
      CALL rs(3,3,C,evle(:,ip),matz,evtr,fv1,fv2,ierr)
    END DO
    DO ip=1,nprtcl
      Cre=MATMUL(Bren(:,:,ip),TRANSPOSE(Bren(:,:,ip)))
      matz=5
      CALL rs(3,3,Cre,evlen(:,ip),matz,evtr,fv1,fv2,ierr)
    END DO

    IF (ifile .EQ. if2) THEN
    DO ip=1,mm
    WRITE(15,*) LOG(evle(1,ip)),LOG(evlen(1,ip))
    END DO
      WRITE(*,*) 'data output'
      i = TecIni('evle'//NULLCHR, 'evle evlen'//NULLCHR, 'evl1.plt'//NULLCHR, &
                 '.'//NULLCHR, Debug, VIsDouble)
      i = TecZne('Zone'//NULLCHR,mm,1,1,'BLOCK'//NULLCHR,NULLCHR)
  
      i = TecDat(mm,LOG(evle(1,1:mm)),0)
      i = TecDat(mm,LOG(evlen(1,1:mm)),0)
  
      WRITE(*,*) 'ending data output'
      i = TecEnd()
    END IF


    Qres=Qre+dt*(-3.*Rre)
    Rres=Rre+dt*(2.*Qre*Qre/3.)
    DO ip=1,nprtcl
    Brens(:,:,ip)=Bren(:,:,ip)+dt*MATMUL(Bren(:,:,ip),Apre(:,:,ip))
    Apres(:,:,ip)=Apre(:,:,ip)+dt*(-MATMUL(Apre(:,:,ip),Apre(:,:,ip))-(2./3.)*idnty*Qre(ip))
    END DO
    
    Qre=.5*(Qre+Qres)+.5*dt*(-3.*Rres)
    Rre=.5*(Rre+Rres)+.5*dt*(2.*Qres*Qres/3.)
    DO ip=1,nprtcl
    Bren(:,:,ip)=.5*(Bren(:,:,ip)+Brens(:,:,ip))+.5*dt*MATMUL(Brens(:,:,ip),Apres(:,:,ip))
    Apre(:,:,ip)=.5*(Apre(:,:,ip)+Apres(:,:,ip))+.5*dt*(-MATMUL(Apres(:,:,ip),Apres(:,:,ip))-(2./3.)*idnty*Qres(ip))
    END DO

    time=time+dt
    
  END DO
  CLOSE(19)


  WRITE(*,*) 'finished'

END PROGRAM postdeform      

