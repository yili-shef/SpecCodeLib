SUBROUTINE helicity_budget(ux,uy,uz,kx,ky,kz,k2,nu,THk,PiHk,Hk,DHk, &
                           lx1,ly,lz,nx,ny,nz)
  USE mconstant
  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: lx1,ly,lz
  INTEGER,  INTENT(IN) :: nx,ny,nz
  REAL(SP), INTENT(IN) :: nu
  COMPLEX(SP),DIMENSION(lx1,ly,lz),INTENT(IN)  :: ux,uy,uz
  REAL(SP),   DIMENSION(lx1,ly,lz),INTENT(IN)  :: kx,ky,kz,k2
  REAL(SP),   DIMENSION(lx1-1),    INTENT(OUT) :: THk,Hk,DHk,PiHk

  COMPLEX(SP),DIMENSION(lx1,ly,lz) :: wx,wy,wz,wxo,wyo,wzo
  REAL(SP),   DIMENSION(lx1,ly,lz) :: atmp
  INTEGER :: i
  
  wx=eye*(ky*uz-kz*uy)
  wy=eye*(kz*ux-kx*uz)
  wz=eye*(kx*uy-ky*ux)
 
  atmp=2._SP*(wx*CONJG(ux)+wy*CONJG(uy)+wz*CONJG(uz))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=1,lx1-1
     Hk(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i).LT..5_SP))
     DHk(i)=2._SP*nu*REAL(i,SP)*REAL(i,SP)*Hk(i)
  END DO

  wxo=wx
  wyo=wy
  wzo=wz
  CALL convec_dns(ux,uy,uz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)

  atmp=2._SP*2._SP*REAL(wxo*CONJG(wx)+wyo*CONJG(wy)+wzo*CONJG(wz))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=lx1-1,1,-1
    THk(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i).LT..5_SP))
    PiHk(i)=SUM(THk(i:lx1-1))
  END DO

END SUBROUTINE helicity_budget
