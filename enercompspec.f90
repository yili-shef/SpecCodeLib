SUBROUTINE enercompspec (ux,uy,uz,kx,ky,kz,k2,E11,E22,E33,V11,V22,V33,lx1,ly,lz,nx,ny,nz)
  USE mconstant
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lx1,ly,lz
  INTEGER, INTENT(IN) :: nx,ny,nz
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: ux,uy,uz
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz,k2
  REAL(SP),    DIMENSION(lx1-1),     INTENT(OUT) :: E11,E22,E33,V11,V22,V33

  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: wx,wy,wz
  REAL(SP),    DIMENSION(lx1,ly,lz) :: atmp
  INTEGER :: i

  atmp=2._SP*(.5_SP*ux*CONJG(ux))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=1,lx1-1
    E11(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i).LT.0.5_SP))
  END DO
  
  atmp=2._SP*(.5_SP*uy*CONJG(uy))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=1,lx1-1
    E22(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i).LT.0.5_SP))
  END DO

  atmp=2._SP*(.5_SP*uz*CONJG(uz))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=1,lx1-1
    E33(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i).LT.0.5_SP))
  END DO

  wx=eye*(ky*uz-kz*uy)
  wy=eye*(kz*ux-kx*uz)
  wz=eye*(kx*uy-ky*ux)

  atmp=2._SP*(.5_SP*wx*CONJG(wx))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=1,lx1-1
    V11(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i).LT.0.5_SP))
  END DO
  
  atmp=2._SP*(.5_SP*wy*CONJG(wy))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=1,lx1-1
    V22(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i).LT.0.5_SP))
  END DO

  atmp=2._SP*(.5_SP*wz*CONJG(wz))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=1,lx1-1
    V33(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i).LT.0.5_SP))
  END DO

END SUBROUTINE enercompspec      
