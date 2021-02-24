SUBROUTINE localheli(vx,vy,vz,wx,wy,wz,tau11,tau12,tau13,tau22,tau23,tau33,kx,ky,kz, &
                     lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb,delta)
  USE mconstant
  USE mprmtr
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER,  INTENT(IN) :: lx1,ly,lz,lxb1,lyb,lzb,nxb,nyb,nzb
  REAL(SP), INTENT(IN) :: delta
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: vx,vy,vz,wx,wy,wz
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(OUT) :: tau11,tau12,tau13,tau22,tau23,tau33
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz

  COMPLEX(SP), DIMENSION(lxb1,lyb,lzb) :: Sb,S11b,S12b,S13b,S22b,S23b,S33b
  COMPLEX(SP), DIMENSION(lxb1,lyb,lzb) :: SRdRb,R11b,R12b,R13b,R22b,R23b,R33b
  REAL(SP) :: ignore_me,const

  const=1._SP/REAL(nxb*nyb*nzb,SP)
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
 
  SRdRb = cmplx(                                                         &
                sqrt(                                                    &
                     2*(REAL(R11b,SP)**2+REAL(R22b,SP)**2+REAL(R33b,SP)**2        &
                     +2.*(REAL(R12b,SP)**2+REAL(R13b,SP)**2+REAL(R23b,SP)**2))    &
                    )                                                    &
               ,                                                         &
                sqrt(                                                    &
                     2*(aimag(R11b)**2+aimag(R22b)**2+aimag(R33b)**2     &
                     +2.*(aimag(R12b)**2+aimag(R13b)**2+aimag(R23b)**2)) &
                    )                                                    &
               )
  SRdRb = cmplx(                                                                                &
                 ( REAL(S11b,SP)*REAL(R11b,SP)+REAL(S22b,SP)*REAL(R22b,SP)+REAL(S33b,SP)*REAL(R33b,SP)            &
                   +2*(REAL(S12b,SP)*REAL(R12b,SP)+REAL(S13b,SP)*REAL(R13b,SP)+REAL(S23b,SP)*REAL(R23b,SP))       &
                 )                                                                              &
                 /REAL(SRdRb,SP)                                                                   &
               ,                                                                                &
                 ( aimag(S11b)*aimag(R11b)+aimag(S22b)*aimag(R22b)+aimag(S33b)*aimag(R33b)      &
                   +2*(aimag(S12b)*aimag(R12b)+aimag(S13b)*aimag(R13b)+aimag(S23b)*aimag(R23b)) &
                 )                                                                              &
                 /aimag(SRdRb)                                                                  &
               )

  S11b=-2.*c1*delta**2*CMPLX(REAL(Sb,SP)*REAL(S11b,SP),aimag(Sb)*aimag(S11b))
  S12b=-2.*c1*delta**2*CMPLX(REAL(Sb,SP)*REAL(S12b,SP),aimag(Sb)*aimag(S12b))
  S13b=-2.*c1*delta**2*CMPLX(REAL(Sb,SP)*REAL(S13b,SP),aimag(Sb)*aimag(S13b))
  S22b=-2.*c1*delta**2*CMPLX(REAL(Sb,SP)*REAL(S22b,SP),aimag(Sb)*aimag(S22b))
  S23b=-2.*c1*delta**2*CMPLX(REAL(Sb,SP)*REAL(S23b,SP),aimag(Sb)*aimag(S23b))
  S33b=-2.*c1*delta**2*CMPLX(REAL(Sb,SP)*REAL(S33b,SP),aimag(Sb)*aimag(S33b))

  S11b=S11b-beta*delta**3*CMPLX(REAL(SRdRb,SP)*REAL(R11b,SP),aimag(SRdRb)*aimag(R11b))
  S12b=S12b-beta*delta**3*CMPLX(REAL(SRdRb,SP)*REAL(R12b,SP),aimag(SRdRb)*aimag(R12b))
  S13b=S13b-beta*delta**3*CMPLX(REAL(SRdRb,SP)*REAL(R13b,SP),aimag(SRdRb)*aimag(R13b))
  S22b=S22b-beta*delta**3*CMPLX(REAL(SRdRb,SP)*REAL(R22b,SP),aimag(SRdRb)*aimag(R22b))
  S23b=S23b-beta*delta**3*CMPLX(REAL(SRdRb,SP)*REAL(R23b,SP),aimag(SRdRb)*aimag(R23b))
  S33b=S33b-beta*delta**3*CMPLX(REAL(SRdRb,SP)*REAL(R33b,SP),aimag(SRdRb)*aimag(R33b))

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

END SUBROUTINE localheli
