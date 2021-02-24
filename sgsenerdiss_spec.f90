SUBROUTINE sgsenerdiss_spec(uxi,uyi,uzi,tau11,tau22,tau33,tau12,tau13,tau23, &
                  kx,ky,kz,k2,SGSDEk,G,lx1,ly,lz,lx)
  USE mconstant 
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lx1,ly,lz,lx
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: uxi,uyi,uzi
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: tau11,tau22,tau33,tau12,tau13,tau23
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz,k2,G
  REAL(SP),    DIMENSION(lx),        INTENT(OUT) :: SGSDEk

  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: ux,uy,uz,S11,S12,S13,S22,S23,S33
  REAL(SP),    DIMENSION(lx1,ly,lz) :: atmp

  INTEGER :: i
  

  ux=uxi*G
  uy=uyi*G
  uz=uzi*G
  CALL symmpart(ux,uy,uz,S11,S12,S13,S22,S23,S33,lx1,ly,lz,kx,ky,kz)
  atmp=-2._SP*REAL(tau11*CONJG(S11)+tau22*CONJG(S22)+tau33*CONJG(S33)+ &
         2._SP*(tau12*CONJG(S12)+tau13*CONJG(S13)+tau23*CONJG(S23)))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=1,lx
    SGSDEk(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i).LT.0.5_SP))
  END DO
END SUBROUTINE sgsenerdiss_spec
