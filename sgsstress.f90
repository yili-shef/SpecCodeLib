! Calculate the Fourier transform of the SGS stresses a priori from DNS data. The
! returned variables tauij are the Fourier transform of the (real) SGS stresses.
! This is different from the MPI version, which returns the real SGS stresses.
SUBROUTINE sgsstress (vxi,vyi,vzi,tau11,tau22,tau33,tau12,tau13,tau23,G,lx1,ly,lz,lx,nx,ny,nz) 
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE
 
  INTEGER, INTENT(IN) :: lx1,ly,lz,lx,nx,ny,nz
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: vxi, vyi, vzi
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(OUT) :: tau11,tau12,tau13,tau22,tau23,tau33
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: G

  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: vx, vy, vz, arrt
  REAL(SP) :: const, ignore_me
  INTEGER  :: i

  const=1._SP/REAL(nx*ny*nz,SP)
  vx=vxi 
  vy=vyi
  vz=vzi

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,vx,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,vy,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,vz,ignore_me)

  DO i = 1, lx
  arrt(i,:,:) = CMPLX(REAL(vx(i,:,:))*REAL(vx(i,:,:)),AIMAG(vx(i,:,:))*AIMAG(vx(i,:,:)))
  END DO
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,arrt,ignore_me)
  arrt=G*const*arrt
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,arrt,ignore_me)
  tau11(1:lx,:,:)=arrt(1:lx,:,:)

  DO i = 1, lx
  arrt(i,:,:) = CMPLX(REAL(vy(i,:,:))*REAL(vy(i,:,:)),AIMAG(vy(i,:,:))*AIMAG(vy(i,:,:)))
  END DO
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,arrt,ignore_me)
  arrt=G*const*arrt
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,arrt,ignore_me)
  tau22(1:lx,:,:)=arrt(1:lx,:,:)

  DO i = 1, lx
  arrt(i,:,:) = CMPLX(REAL(vz(i,:,:))*REAL(vz(i,:,:)),AIMAG(vz(i,:,:))*AIMAG(vz(i,:,:)))
  END DO
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,arrt,ignore_me)
  arrt=G*const*arrt
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,arrt,ignore_me)
  tau33(1:lx,:,:)=arrt(1:lx,:,:)

  DO i = 1, lx
  arrt(i,:,:) = CMPLX(REAL(vx(i,:,:))*REAL(vy(i,:,:)),AIMAG(vx(i,:,:))*AIMAG(vy(i,:,:)))
  END DO
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,arrt,ignore_me)
  arrt=G*const*arrt
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,arrt,ignore_me)
  tau12(1:lx,:,:)=arrt(1:lx,:,:)

  DO i = 1, lx
  arrt(i,:,:) = CMPLX(REAL(vx(i,:,:))*REAL(vz(i,:,:)),AIMAG(vx(i,:,:))*AIMAG(vz(i,:,:)))
  END DO
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,arrt,ignore_me)
  arrt=G*const*arrt
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,arrt,ignore_me)
  tau13(1:lx,:,:)=arrt(1:lx,:,:)

  DO i = 1, lx
  arrt(i,:,:) = CMPLX(REAL(vy(i,:,:))*REAL(vz(i,:,:)),AIMAG(vy(i,:,:))*AIMAG(vz(i,:,:)))
  END DO
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,arrt,ignore_me)
  arrt=G*const*arrt
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,arrt,ignore_me)
  tau23(1:lx,:,:)=arrt(1:lx,:,:)

  vx=G*vxi
  vy=G*vyi
  vz=G*vzi

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,vx,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,vy,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,vz,ignore_me)

  tau11=tau11-CMPLX(REAL(vx)*REAL(vx),AIMAG(vx)*AIMAG(vx))
  tau22=tau22-CMPLX(REAL(vy)*REAL(vy),AIMAG(vy)*AIMAG(vy))
  tau33=tau33-CMPLX(REAL(vz)*REAL(vz),AIMAG(vz)*AIMAG(vz))
  tau12=tau12-CMPLX(REAL(vx)*REAL(vy),AIMAG(vx)*AIMAG(vy))
  tau13=tau13-CMPLX(REAL(vx)*REAL(vz),AIMAG(vx)*AIMAG(vz))
  tau23=tau23-CMPLX(REAL(vy)*REAL(vz),AIMAG(vy)*AIMAG(vz))

  CALL rfftwnd_f77_one_real_to_complex(R2C3D,tau11,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,tau12,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,tau13,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,tau22,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,tau23,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,tau33,ignore_me)
  tau11=const*tau11
  tau12=const*tau12
  tau13=const*tau13
  tau22=const*tau22
  tau23=const*tau23
  tau33=const*tau33

END SUBROUTINE sgsstress

