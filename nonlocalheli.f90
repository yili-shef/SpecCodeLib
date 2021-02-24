SUBROUTINE nonlocalheli(vx,vy,vz,wx,wy,wz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz, &
                        lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,delta,eps,eta)
  USE mconstant
  USE mprmtr
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb
  REAL(SP), INTENT(IN) :: delta,eps,eta
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: vx,vy,vz,wx,wy,wz
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(OUT) :: tau11,tau12,tau13,tau22,tau23,tau33
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz

  COMPLEX(SP), DIMENSION(lxb1,lyb,lzb) :: Sb,S11b,S12b,S13b,S22b,S23b,S33b
  COMPLEX(SP), DIMENSION(lxb1,lyb,lzb) :: R11b,R12b,R13b,R22b,R23b,R33b
  REAL(SP) :: ignore_me,const, c2

  const=1._SP/REAL(nxb*nyb*nzb,SP)
  c2=beta_nonlocal*eta*delta/eps

  CALL symmpart(vx,vy,vz,tau11,tau12,tau13,tau22,tau23,tau33,lx1,ly,lz,kx,ky,kz)
  CALL padd(tau11,S11b,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(tau12,S12b,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(tau13,S13b,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(tau22,S22b,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(tau23,S23b,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(tau33,S33b,lx1,ly,lz,lxb1,lyb,lzb)
  CALL symmpart(wx,wy,wz,tau11,tau12,tau13,tau22,tau23,tau33,lx1,ly,lz,kx,ky,kz)
  CALL padd(tau11,R11b,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(tau12,R12b,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(tau13,R13b,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(tau22,R22b,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(tau23,R23b,lx1,ly,lz,lxb1,lyb,lzb)
  CALL padd(tau33,R33b,lx1,ly,lz,lxb1,lyb,lzb)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,S11b,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,S12b,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,S13b,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,S22b,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,S23b,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,S33b,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,R11b,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,R12b,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,R13b,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,R22b,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,R23b,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3DB,R33b,ignore_me)

  Sb = CMPLX(                                                                     &
             SQRT(                                                                &
                  2._SP*(REAL(S11b,SP)**2+REAL(S22b,SP)**2+REAL(S33b,SP)**2       &
                  +2._SP*(REAL(S12b,SP)**2+REAL(S13b,SP)**2+REAL(S23b,SP)**2))    &
                 )                                                       &
            ,                                                            &
             SQRT(                                                       &
                  2._SP*(AIMAG(S11b)**2+AIMAG(S22b)**2+AIMAG(S33b)**2    &
                  +2._SP*(AIMAG(S12b)**2+AIMAG(S13b)**2+AIMAG(S23b)**2)) &
                 )                                                       &
            )
  S11b=-2.*c1*delta**2*CMPLX(REAL(Sb,SP)*REAL(S11b,SP),AIMAG(Sb)*AIMAG(S11b))
  S12b=-2.*c1*delta**2*CMPLX(REAL(Sb,SP)*REAL(S12b,SP),AIMAG(Sb)*AIMAG(S12b))
  S13b=-2.*c1*delta**2*CMPLX(REAL(Sb,SP)*REAL(S13b,SP),AIMAG(Sb)*AIMAG(S13b))
  S22b=-2.*c1*delta**2*CMPLX(REAL(Sb,SP)*REAL(S22b,SP),AIMAG(Sb)*AIMAG(S22b))
  S23b=-2.*c1*delta**2*CMPLX(REAL(Sb,SP)*REAL(S23b,SP),AIMAG(Sb)*AIMAG(S23b))
  S33b=-2.*c1*delta**2*CMPLX(REAL(Sb,SP)*REAL(S33b,SP),AIMAG(Sb)*AIMAG(S33b))

  S11b=S11b-c2*delta**3*CMPLX(REAL(Sb,SP)*REAL(R11b,SP),AIMAG(Sb)*AIMAG(R11b))
  S12b=S12b-c2*delta**3*CMPLX(REAL(Sb,SP)*REAL(R12b,SP),AIMAG(Sb)*AIMAG(R12b))
  S13b=S13b-c2*delta**3*CMPLX(REAL(Sb,SP)*REAL(R13b,SP),AIMAG(Sb)*AIMAG(R13b))
  S22b=S22b-c2*delta**3*CMPLX(REAL(Sb,SP)*REAL(R22b,SP),AIMAG(Sb)*AIMAG(R22b))
  S23b=S23b-c2*delta**3*CMPLX(REAL(Sb,SP)*REAL(R23b,SP),AIMAG(Sb)*AIMAG(R23b))
  S33b=S33b-c2*delta**3*CMPLX(REAL(Sb,SP)*REAL(R33b,SP),AIMAG(Sb)*AIMAG(R33b))

  S11b=S11b*const
  S12b=S12b*const
  S13b=S13b*const
  S22b=S22b*const
  S23b=S23b*const
  S33b=S33b*const
  CALL rfftwnd_f77_one_real_to_complex(R2C3DB,S11b,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3DB,S12b,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3DB,S13b,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3DB,S22b,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3DB,S23b,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3DB,S33b,ignore_me)

  CALL unpadd(S11b,tau11,lx1,ly,lz,lxb1,lyb,lzb)
  CALL unpadd(S12b,tau12,lx1,ly,lz,lxb1,lyb,lzb)
  CALL unpadd(S13b,tau13,lx1,ly,lz,lxb1,lyb,lzb)
  CALL unpadd(S22b,tau22,lx1,ly,lz,lxb1,lyb,lzb)
  CALL unpadd(S23b,tau23,lx1,ly,lz,lxb1,lyb,lzb)
  CALL unpadd(S33b,tau33,lx1,ly,lz,lxb1,lyb,lzb)

END SUBROUTINE nonlocalheli
