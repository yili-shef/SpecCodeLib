PROGRAM materialdeform
  USE mconstant
  IMPLICIT NONE

  INTEGER,   PARAMETER :: nx=128,ny=128,nz=128, nprtcl=1024*1024
  INTEGER,   PARAMETER :: if1=1,if2=211
  REAL(DP), PARAMETER :: dt=0.017

  REAL(DP), DIMENSION(3,3,nprtcl) :: B,B1,B2,B3,Ap,App
  REAL(SP),  DIMENSION(3,3,nprtcl) :: Apf
  REAL(SP),  DIMENSION(nprtcl)     :: trash

  INTEGER :: idnty(3,3),i,ifile,ip

  idnty=0
  idnty(1,1)=1; idnty(2,2)=1; idnty(3,3)=1

  OPEN(19,FILE='./xBdata/B.dat',FORM='unformatted')
  OPEN(20,FILE='./xBdata/A.dat',FORM='unformatted')
  
  DO i=1,nprtcl
  B(:,:,i)=idnty
  END DO
  WRITE(19) B
  
  READ(20) trash,trash,Apf
  Ap=DBLE(Apf)

  DO ifile = if1+1,if2
    WRITE(*,*) ifile
    READ(20) trash,trash,Apf
    App=DBLE(Apf)

    DO ip=1, nprtcl

      B1(:,:,ip)=B(:,:,ip)+ .5*dt*MATMUL(B(:,:,ip),Ap(:,:,ip))
      B2(:,:,ip)=B(:,:,ip)+.25*dt*MATMUL(B1(:,:,ip),Ap(:,:,ip)+App(:,:,ip))
      B3(:,:,ip)=B(:,:,ip)+ .5*dt*MATMUL(B2(:,:,ip),Ap(:,:,ip)+App(:,:,ip))
      B(:,:,ip) =B(:,:,ip)+(dt/6.)*(MATMUL(B(:,:,ip),Ap(:,:,ip))+ &
                 MATMUL(B1(:,:,ip)+B2(:,:,ip),Ap(:,:,ip)+App(:,:,ip))+ &
                 MATMUL(B3(:,:,ip),App(:,:,ip)))
    END DO
    WRITE(19) B
    Ap=App
  END DO

  CLOSE(20)
  CLOSE(19)


  WRITE(*,*) 'finished'
  STOP

END PROGRAM materialdeform      
