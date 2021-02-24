SUBROUTINE pressure(vx,vy,vz,k2,kx,ky,kz,lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,p)
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb
  COMPLEX(SP),DIMENSION(lx1,ly,lz),INTENT(IN)    :: vx,vy,vz
  REAL(SP),   DIMENSION(lx1,ly,lz),INTENT(IN)    :: k2,kx,ky,kz
  COMPLEX(SP),DIMENSION(lx1,ly,lz),INTENT(OUT)   :: p

  COMPLEX(SP),DIMENSION(lx1,ly,lz) :: Axy,Ayx
  COMPLEX(SP),DIMENSION(lxb1,lyb,lzb) :: Axyb, Ayxb, AAb
  REAL(SP) :: const, ignore_me

  const=1._SP/REAL(nxb*nyb*nzb,SP)

  Axy=eye*kx*vx ! A11
  CALL padd(Axy,Axyb,lx1,ly,lz,lxb1,lyb,lzb)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,Axyb,ignore_me)
  AAb=CMPLX(REAL(Axyb,SP)*REAL(Axyb,SP), AIMAG(Axyb)*AIMAG(Axyb))

  Axy=eye*ky*vy ! A22
  CALL padd(Axy,Axyb,lx1,ly,lz,lxb1,lyb,lzb)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,Axyb,ignore_me)
  AAb=AAb+CMPLX(REAL(Axyb,SP)*REAL(Axyb,SP), AIMAG(Axyb)*AIMAG(Axyb))

  Axy=eye*kz*vz ! A22
  CALL padd(Axy,Axyb,lx1,ly,lz,lxb1,lyb,lzb)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,Axyb,ignore_me)
  AAb=AAb+CMPLX(REAL(Axyb,SP)*REAL(Axyb,SP), AIMAG(Axyb)*AIMAG(Axyb))

  Axy=eye*kx*vy !A12
  Ayx=eye*ky*vx !A21
  CALL padd(Axy,Axyb,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(Ayx,Ayxb,lx1,ly,lz,lxb1,lyb,lzb)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,Axyb,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,Ayxb,ignore_me)
  AAb=AAb+2.*CMPLX(REAL(Axyb,SP)*REAL(Ayxb,SP), AIMAG(Axyb)*AIMAG(Ayxb))

  Axy=eye*kx*vz !A13
  Ayx=eye*kz*vx !A31
  CALL padd(Axy,Axyb,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(Ayx,Ayxb,lx1,ly,lz,lxb1,lyb,lzb)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,Axyb,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,Ayxb,ignore_me)
  AAb=AAb+2.*CMPLX(REAL(Axyb,SP)*REAL(Ayxb,SP), AIMAG(Axyb)*AIMAG(Ayxb))

  Axy=eye*ky*vz !A23
  Ayx=eye*kz*vy !A32
  CALL padd(Axy,Axyb,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(Ayx,Ayxb,lx1,ly,lz,lxb1,lyb,lzb)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,Axyb,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,Ayxb,ignore_me)
  AAb=AAb+2.*CMPLX(REAL(Axyb,SP)*REAL(Ayxb,SP), AIMAG(Axyb)*AIMAG(Ayxb))

  AAb=AAb*const
  CALL rfftwnd_f77_one_real_to_complex(R2C3DB,AAb,ignore_me)
  CALL unpadd(AAb,p,lx1,ly,lz,lxb1,lyb,lzb)
  p=p/k2
  CALL symmetrize(p,k2,lx1,ly,lz)

END SUBROUTINE pressure
