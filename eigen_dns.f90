PROGRAM eigen
  USE mconstant
  IMPLICIT NONE

  INTEGER,  PARAMETER :: nprtcl=1024*1024
  INTEGER,  PARAMETER :: if1=1,if2=211 !40 for filtered case, 20 for DNS
  REAL(DP), PARAMETER :: dt=0.017

  REAL(DP), DIMENSION(3,3,nprtcl) :: B
  REAL(DP), DIMENSION(3,nprtcl)   :: evle

  REAL(DP), DIMENSION(3,3) :: C,Apt
  REAL(DP), DIMENSION(3)   :: fv1,fv2

  INTEGER  :: i,ii,kk,ifile,ip,matz,ierr
  REAL(DP) :: time
  
  
  OPEN(22,FILE='./xBdata/evle_dns.dat',FORM='unformatted')
  OPEN(19,FILE='./xBdata/B.dat',FORM='unformatted')


  time=0.
  DO ifile=if1,if2
    READ(19) B

    DO ip=1,nprtcl


      C=MATMUL(TRANSPOSE(B(:,:,ip)),B(:,:,ip))
      matz=5
      CALL rs(3,3,C,evle(:,ip),matz,Apt,fv1,fv2,ierr)


    END DO

    time=time+dt
    
    WRITE(22) evle
  END DO
  CLOSE(19)
  CLOSE(22)


  WRITE(*,*) 'finished'

END PROGRAM eigen      

