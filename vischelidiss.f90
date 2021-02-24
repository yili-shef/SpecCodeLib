SUBROUTINE vischelidiss (ux,uy,uz,wx,wy,wz,kx,ky,kz,nu,vischd,lx1,lx,ly,lz,nx,ny,nz)
  USE mconstant
  USE mfftwplan3d

  INTEGER,  INTENT(IN) :: lx1,lx,ly,lz,nx,ny,nz
  REAL(SP), INTENT(IN) :: nu

  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: ux,uy,uz,wx,wy,wz
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz
  REAL(SP),    DIMENSION(nx,ny,nz),  INTENT(OUT) :: vischd

  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: s11,s12,s13,s22,s23,s33
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: r11,r12,r13,r22,r23,r33
  REAL(SP) :: ignore_me

  s11=eye*kx*ux
  s22=eye*ky*uy
  s33=eye*kz*uz
  s12=.5_SP*eye*(ky*ux+kx*uy)
  s13=.5_SP*eye*(kz*ux+kx*uz)
  s23=.5_SP*eye*(ky*uz+kz*uy)

  r11=eye*kx*wx
  r22=eye*ky*wy
  r33=eye*kz*wz
  r12=.5_SP*eye*(ky*wx+kx*wy)
  r13=.5_SP*eye*(kz*wx+kx*wz)
  r23=.5_SP*eye*(ky*wz+kz*wy)

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,s11,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,s12,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,s13,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,s22,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,s23,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,s33,ignore_me)

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,r11,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,r12,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,r13,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,r22,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,r23,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,r33,ignore_me)

  vischd(1:nx:2,:,:)= REAL(s11(1:lx,:,:))*REAL(r11(1:lx,:,:)) &
                     +REAL(s22(1:lx,:,:))*REAL(r22(1:lx,:,:)) &
                     +REAL(s33(1:lx,:,:))*REAL(r33(1:lx,:,:)) &
                     +2._SP*(REAL(s12(1:lx,:,:))*REAL(r12(1:lx,:,:)) &
                            +REAL(s13(1:lx,:,:))*REAL(r13(1:lx,:,:)) &
                            +REAL(s23(1:lx,:,:))*REAL(r23(1:lx,:,:)) &
                            )
                    

  vischd(2:nx:2,:,:)= AIMAG(s11(1:lx,:,:))*AIMAG(r11(1:lx,:,:)) &
                     +AIMAG(s22(1:lx,:,:))*AIMAG(r22(1:lx,:,:)) &
                     +AIMAG(s33(1:lx,:,:))*AIMAG(r33(1:lx,:,:)) &
                     +2._SP*(AIMAG(s12(1:lx,:,:))*AIMAG(r12(1:lx,:,:)) &
                            +AIMAG(s13(1:lx,:,:))*AIMAG(r13(1:lx,:,:)) &
                            +AIMAG(s23(1:lx,:,:))*AIMAG(r23(1:lx,:,:)) &
                            ) 

  vischd=4._SP*nu*vischd
END SUBROUTINE vischelidiss      

