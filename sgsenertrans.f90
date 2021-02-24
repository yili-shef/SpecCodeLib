SUBROUTINE sgsenertrans (ux,uy,uz,tau11,tau12,tau13,tau22,tau23,tau33, &
                         kx,ky,kz,sgset,lx1,ly,lz,lx,nx,ny,nz)
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: lx1,ly,lz,lx,nx,ny,nz

  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: ux,uy,uz
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: tau11,tau12,tau13
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: tau22,tau23,tau33
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz
  REAL(SP),    DIMENSION(nx,ny,nz),  INTENT(OUT) :: sgset

  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: ut, dtaut, dtauts
  REAL(SP) :: ignore_me


  dtauts=0._SP
  ut=ux
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,ut,ignore_me)

  dtaut=eye*kx*tau11
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,dtaut,ignore_me)
  dtauts=dtauts+dtaut
  
  dtaut=eye*ky*tau12
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,dtaut,ignore_me)
  dtauts=dtauts+dtaut

  dtaut=eye*kz*tau13
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,dtaut,ignore_me)
  dtauts=dtauts+dtaut
 
  sgset(1:lx:2,:,:)=REAL(ut(1:lx,:,:))*REAL(dtauts(1:lx,:,:))
  sgset(2:lx:2,:,:)=AIMAG(ut(1:lx,:,:))*AIMAG(dtauts(1:lx,:,:))
  
  dtauts=0._SP
  ut=uy
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,ut,ignore_me)

  dtaut=eye*kx*tau12
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,dtaut,ignore_me)
  dtauts=dtauts+dtaut
  
  dtaut=eye*ky*tau22
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,dtaut,ignore_me)
  dtauts=dtauts+dtaut

  dtaut=eye*kz*tau23
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,dtaut,ignore_me)
  dtauts=dtauts+dtaut

  sgset(1:lx:2,:,:)=sgset(1:lx:2,:,:)+REAL(ut(1:lx,:,:))*REAL(dtauts(1:lx,:,:))
  sgset(2:lx:2,:,:)=sgset(2:lx:2,:,:)+AIMAG(ut(1:lx,:,:))*AIMAG(dtauts(1:lx,:,:))

  dtauts=0._SP
  ut=uz
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,ut,ignore_me)

  dtaut=eye*kx*tau13
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,dtaut,ignore_me)
  dtauts=dtauts+dtaut
  
  dtaut=eye*ky*tau23
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,dtaut,ignore_me)
  dtauts=dtauts+dtaut

  dtaut=eye*kz*tau33
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,dtaut,ignore_me)
  dtauts=dtauts+dtaut

  sgset(1:lx:2,:,:)=sgset(1:lx:2,:,:)+REAL(ut(1:lx,:,:))*REAL(dtauts(1:lx,:,:))
  sgset(2:lx:2,:,:)=sgset(2:lx:2,:,:)+AIMAG(ut(1:lx,:,:))*AIMAG(dtauts(1:lx,:,:))

END SUBROUTINE sgsenertrans      
