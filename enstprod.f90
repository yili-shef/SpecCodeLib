SUBROUTINE enstroprod (ux,uy,uz,kx,ky,kz,enspro,lx1,ly,lz,lx,nx,ny,nz)
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lx1,ly,lz,lx,nx,ny,nz

  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: ux,uy,uz
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz
  REAL(SP),    DIMENSION(nx,ny,nz),  INTENT(OUT) :: enspro

  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: wx,wy,wz,s11,s12,s13,s22,s23,s33
  REAL(SP) :: ignore_me

  wx=eye*(ky*uz-kz*uy)
  wy=eye*(kz*ux-kx*uz)
  wz=eye*(kx*uy-ky*ux)
  s11=eye*kx*ux
  s22=eye*ky*uy
  s33=eye*kz*uz
  s12=.5_SP*eye*(ky*ux+kx*uy)
  s13=.5_SP*eye*(kz*ux+kx*uz)
  s23=.5_SP*eye*(ky*uz+kz*uy)
  
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,wx,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,wy,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,wz,ignore_me)

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,s11,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,s12,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,s13,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,s22,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,s23,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,s33,ignore_me)

  enspro(1:nx:2,:,:) = REAL(s11(1:lx,:,:))*REAL(wx(1:lx,:,:))*REAL(wx(1:lx,:,:)) &
                      +REAL(s22(1:lx,:,:))*REAL(wy(1:lx,:,:))*REAL(wy(1:lx,:,:)) &
                      +REAL(s33(1:lx,:,:))*REAL(wz(1:lx,:,:))*REAL(wz(1:lx,:,:)) &
                      +2._SP*(REAL(s12(1:lx,:,:))*REAL(wx(1:lx,:,:))*REAL(wy(1:lx,:,:)) &
                             +REAL(s13(1:lx,:,:))*REAL(wx(1:lx,:,:))*REAL(wz(1:lx,:,:)) &
                             +REAL(s23(1:lx,:,:))*REAL(wy(1:lx,:,:))*REAL(wz(1:lx,:,:)) &
                             )

  enspro(2:nx:2,:,:) = AIMAG(s11(1:lx,:,:))*AIMAG(wx(1:lx,:,:))*AIMAG(wx(1:lx,:,:)) &
                      +AIMAG(s22(1:lx,:,:))*AIMAG(wy(1:lx,:,:))*AIMAG(wy(1:lx,:,:)) &
                      +AIMAG(s33(1:lx,:,:))*AIMAG(wz(1:lx,:,:))*AIMAG(wz(1:lx,:,:)) &
                      +2._SP*(AIMAG(s12(1:lx,:,:))*AIMAG(wx(1:lx,:,:))*AIMAG(wy(1:lx,:,:)) &
                             +AIMAG(s13(1:lx,:,:))*AIMAG(wx(1:lx,:,:))*AIMAG(wz(1:lx,:,:)) &
                             +AIMAG(s23(1:lx,:,:))*AIMAG(wy(1:lx,:,:))*AIMAG(wz(1:lx,:,:)) &
                             )

END SUBROUTINE enstroprod 
