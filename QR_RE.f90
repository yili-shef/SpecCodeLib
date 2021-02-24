PROGRAM materialdeform
  USE mconstant
  IMPLICIT NONE

  INTEGER,  PARAMETER :: nprtcl=1024*1024
  INTEGER,  PARAMETER :: if1=1,if2=20
  REAL(DP), PARAMETER :: dt=0.017

  REAL(DP), DIMENSION(nprtcl)     :: Qre,Rre
  REAL(SP), DIMENSION(nprtcl)     :: Qresp,Rresp

  REAL(DP) :: Qre1,Qre2,Qre3,Rre1,Rre2,Rre3,Qres,Rres

  INTEGER :: i,ifile,ip

  OPEN(20,FILE='./xBdata/A.dat',FORM='unformatted')
  READ(20) Qresp,Rresp
  CLOSE(20)
  Qre=REAL(Qresp,DP)
  Rre=REAL(Rresp,DP)

  OPEN(21,FILE='./xBdata/Are.dat',FORM='unformatted')
  
  WRITE(21) Qre,Rre
  
  DO ifile = if1+1,if2
    WRITE(*,*) ifile

    DO ip=1,nprtcl
      Qre1=Qre(ip)+.5*dt*(-3.*Rre(ip))
      Rre1=Rre(ip)+.5*dt*(2.*Qre(ip)*Qre(ip)/3.)
      
      Qre2=Qre(ip)+.5*dt*(-3.*Rre1)
      Rre2=Rre(ip)+.5*dt*(2.*Qre1*Qre1/3.)
      
      Qre3=Qre(ip)+dt*(-3.*Rre2)
      Rre3=Rre(ip)+dt*(2.*Qre2*Qre2/3.)

      Qres=Qre(ip)+(dt/6.)*(-3.*Rre(ip))+(2.*dt/3.)*(-3.*(Rre1+Rre2)/2.)+(dt/6.)*(-3.*Rre3)
      Rres=Rre(ip)+(dt/6.)*(2.*Qre(ip)*Qre(ip)/3.)+(2.*dt/3.)*(2./3.*(Qre1+Qre2)**2/4.)+(dt/6.)*(2.*Qre3*Qre3/3.)
      
      Qre(ip)=Qres; Rre(ip)=Rres
    END DO

    WRITE(21) Qre,Rre
  END DO

  CLOSE(21)


  WRITE(*,*) 'finished'
  STOP

END PROGRAM materialdeform      
