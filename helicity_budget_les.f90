SUBROUTINE helicity_budget_les(ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz,k2,nu,THk,PiHk,Hk,DHk,SGSDHk, &
                           lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb)
  USE mconstant
  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb
  REAL(SP), INTENT(IN) :: nu
  COMPLEX(SP),DIMENSION(lx1,ly,lz),INTENT(IN)  :: ux,uy,uz
  COMPLEX(SP),DIMENSION(lx1,ly,lz),INTENT(IN)  :: tau11,tau12,tau13,tau22,tau23,tau33
  REAL(SP),   DIMENSION(lx1,ly,lz),INTENT(IN)  :: kx,ky,kz,k2
  REAL(SP),   DIMENSION(lx1-1),    INTENT(OUT) :: THk,Hk,DHk,PiHk,SGSDHk

  COMPLEX(SP),DIMENSION(lx1,ly,lz) :: wx,wy,wz,wxo,wyo,wzo
  COMPLEX(SP),DIMENSION(lx1,ly,lz) :: r11,r12,r13,r22,r23,r33
  REAL(SP),   DIMENSION(lx1,ly,lz) :: atmp
  INTEGER :: i
  
  wx=eye*(ky*uz-kz*uy)
  wy=eye*(kz*ux-kx*uz)
  wz=eye*(kx*uy-ky*ux)
  r11=eye*kx*wx
  r22=eye*ky*wy
  r33=eye*kz*wz
  r12=.5_SP*eye*(ky*wx+kx*wy)
  r13=.5_SP*eye*(kz*wx+kx*wz)
  r23=.5_SP*eye*(ky*wz+kz*wy)
 
  atmp=2._SP*(wx*CONJG(ux)+wy*CONJG(uy)+wz*CONJG(uz))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=1,lx1-1
     Hk(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).lt..5_SP))
     DHk(i)=2._SP*nu*REAL(i,SP)*REAL(i,SP)*Hk(i)
  END DO

  atmp=-2._SP*2._SP*REAL(tau11*CONJG(r11)+tau22*CONJG(r22)+tau33*CONJG(r33)+ &
                         2._SP*(tau12*CONJG(r12)+tau13*CONJG(r13)+tau23*CONJG(r23)))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=1,lx1-1
    SGSDHk(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).LT..5_SP))
  END DO
  
  wxo=wx
  wyo=wy
  wzo=wz
  CALL convec(ux,uy,uz,wx,wy,wz,lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb)

  atmp=2._SP*2._SP*REAL(wxo*CONJG(wx)+wyo*CONJG(wy)+wzo*CONJG(wz))
  atmp(1,:,:)=.5_SP*atmp(1,:,:)
  DO i=lx1-1,1,-1
    THk(i)=SUM(atmp,mask=(ABS(SQRT(k2)-i-0.5_SP*oneless).LT..5_SP))
    PiHk(i)=SUM(THk(i:lx1-1))
  END DO

END SUBROUTINE helicity_budget_les
