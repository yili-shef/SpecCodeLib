!!TODO: to be worked: tranlate from helidiss to enerdiss
SUBROUTINE sgsenerdiss (ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33, &
                        kx,ky,kz,sgshd,lx1,lx,ly,lz,nx,ny,nz)

  USE mconstant
  USE mfftwplan3d  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lx1,lx,ly,lz,nx,ny,nz

  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: ux,uy,uz
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: tau11,tau12,tau13,tau22,tau23,tau33
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz
  REAL(SP),    DIMENSION(nx,ny,nz),  INTENT(OUT) :: sgshd

  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: t11,t12,t13,t22,t23,t33
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: r11,r12,r13,r22,r23,r33
  REAL(SP) :: ignore_me


  t11=eye*(ky*uz-kz*uy)
  t12=eye*(kz*ux-kx*uz)
  t13=eye*(kx*uy-ky*ux)
  r11=eye*kx*t11
  r22=eye*ky*t12
  r33=eye*kz*t13
  r12=.5_SP*eye*(ky*t11+kx*t12)
  r13=.5_SP*eye*(kz*t11+kx*t13)
  r23=.5_SP*eye*(ky*t13+kz*t12)
  t11=tau11
  t12=tau12
  t13=tau13
  t22=tau22
  t23=tau23
  t33=tau33

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,r11,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,r12,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,r13,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,r22,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,r23,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,r33,ignore_me)

  
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,t11,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,t12,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,t13,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,t22,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,t23,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,t33,ignore_me)

  sgshd(1:nx:2,:,:)=-2._SP*(REAL(t11(1:lx,:,:))*REAL(r11(1:lx,:,:)) &
                           +REAL(t22(1:lx,:,:))*REAL(r22(1:lx,:,:)) &
                           +REAL(t33(1:lx,:,:))*REAL(r33(1:lx,:,:)) &
                           +2._SP*(REAL(t12(1:lx,:,:))*REAL(r12(1:lx,:,:)) &
                                  +REAL(t13(1:lx,:,:))*REAL(r13(1:lx,:,:)) &
                                  +REAL(t23(1:lx,:,:))*REAL(r23(1:lx,:,:)) &
                                  )&
                           )

  sgshd(2:nx:2,:,:)=-2._SP*(AIMAG(t11(1:lx,:,:))*AIMAG(r11(1:lx,:,:)) &
                           +AIMAG(t22(1:lx,:,:))*AIMAG(r22(1:lx,:,:)) &
                           +AIMAG(t33(1:lx,:,:))*AIMAG(r33(1:lx,:,:)) &
                           +2._SP*(AIMAG(t12(1:lx,:,:))*AIMAG(r12(1:lx,:,:)) &
                                  +AIMAG(t13(1:lx,:,:))*AIMAG(r13(1:lx,:,:)) &
                                  +AIMAG(t23(1:lx,:,:))*AIMAG(r23(1:lx,:,:)) &
                                  ) &
                           )

END SUBROUTINE sgsenerdiss 
