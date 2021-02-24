SUBROUTINE sgshelidiss_spec(uxi,uyi,uzi,tau11,tau22,tau33,tau12,tau13,tau23, &
                            kx,ky,kz,k2,SGSDHk,G,lx1,ly,lz,lx)
  USE mconstant
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lx1,ly,lz,lx

  REAL(SP),    DIMENSION(lx),        INTENT(OUT) :: SGSDHk
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: uxi,uyi,uzi
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: tau11,tau22,tau33,tau12,tau13,tau23
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz,k2,G
  
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: wx,wy,wz,R11,R12,R13,R22,R23,R33
  REAL(SP),    DIMENSION(lx1,ly,lz) :: atmp

  INTEGER :: i

  R11=uxi*G
  R12=uyi*G
  R13=uzi*G
  wx=eye*(ky*R13-kz*R12)
  wy=eye*(kz*R11-kx*R13)
  wz=eye*(kx*R12-ky*R11)
  CALL symmpart(wx,wy,wz,R11,R12,R13,R22,R23,R33,lx1,ly,lz,kx,ky,kz)
  atmp=-2._SP*REAL(tau11*CONJG(R11)+tau22*CONJG(R22)+tau33*CONJG(R33)+ &
         2._SP*(tau12*CONJG(R12)+tau13*CONJG(R13)+tau23*CONJG(R23)))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  atmp=2.*atmp
  DO i = 1, lx 
    SGSDHk(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i).LT.0.5_SP))
  END DO
  
END SUBROUTINE sgshelidiss_spec    
