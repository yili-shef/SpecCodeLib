SUBROUTINE enstrophy (ux,uy,uz,kx,ky,kz,enstr,lx1,ly,lz,lx,nx,ny,nz)
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lx1,lx,ly,lz,nx,ny,nz

  COMPLEX(DP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: ux,uy,uz
  REAL(DP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz
  REAL(DP),    DIMENSION(nx,ny,nz),  INTENT(OUT) :: enstr

  COMPLEX(DP), DIMENSION(lx1,ly,lz) :: wx, wy, wz
  REAL(DP) :: ignore_me


  wx=eye*(ky*uz-kz*uy)
  wy=eye*(kz*ux-kx*uz)
  wz=eye*(kx*uy-ky*ux)

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,wx,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,wy,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,wz,ignore_me)

  enstr(1:nx:2,:,:)=.5_DP*(REAL(wx(1:lx,:,:))*REAL(wx(1:lx,:,:))   &
                          +REAL(wy(1:lx,:,:))*REAL(wy(1:lx,:,:))   &
                          +REAL(wz(1:lx,:,:))*REAL(wz(1:lx,:,:))   &
                          )

  enstr(2:nx:2,:,:)=.5_DP*(AIMAG(wx(1:lx,:,:))*AIMAG(wx(1:lx,:,:))   &
                          +AIMAG(wy(1:lx,:,:))*AIMAG(wy(1:lx,:,:))   &
                          +AIMAG(wz(1:lx,:,:))*AIMAG(wz(1:lx,:,:))   &
                          )

END SUBROUTINE enstrophy      
