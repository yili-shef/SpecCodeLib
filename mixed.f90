subroutine mixed (vx,vy,vz,t11,t12,t13,t22,t23,t33,kx,ky,kz,lx1,ly,lz,nx,ny,nz,delta)
  use mconstant
  use mfftwplan3d

  integer, intent(in) :: lx1,ly,lz,nx,ny,nz
  real(sp), intent(in) :: delta
  complex(sp), dimension(lx1,ly,lz), intent(in) :: vx,vy,vz
  complex(sp), dimension(lx1,ly,lz), intent(out) :: t11,t12,t13,t22,t23,t33
  real(sp), dimension(lx1,ly,lz), intent(in) :: kx,ky,kz

  complex(sp), dimension(lx1,ly,lz) :: S
  complex(sp), dimension(lx1,ly,lz) :: a11, a12, a13
  complex(sp), dimension(lx1,ly,lz) :: a21, a22, a23
  complex(sp), dimension(lx1,ly,lz) :: a31, a32, a33
  real(sp) :: ignore_me,const

  real(sp), parameter :: cs2=0.026, cnl=1./12.

  const=1._sp/real(nx*ny*nz,sp)

  a11=eye*kx*vx; a12=eye*kx*vy; a13=eye*kx*vz
  a21=eye*ky*vx; a22=eye*ky*vy; a23=eye*ky*vz
  a31=eye*kz*vx; a32=eye*kz*vy; a33=eye*kz*vz

  CALL rfftwnd_f77_one_complex_to_real(C2R3D,a11,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,a12,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,a13,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,a21,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,a22,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,a23,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,a31,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,a32,ignore_me)
  CALL rfftwnd_f77_one_complex_to_real(C2R3D,a33,ignore_me)

  S = cmplx( real(a11,sp)**2+real(a22,sp)**2+real(a33,sp)**2, &
            aimag(a11)**2+aimag(a22)**2+aimag(a33)**2)
  S = S + cmplx( (real(a12,sp)+real(a21,sp))**2/2., &
                 (aimag(a12)+aimag(a21))**2/2.)
  S = S + cmplx( (real(a23,sp)+real(a32,sp))**2/2., &
                 (aimag(a23)+aimag(a32))**2/2.)
  S = S + cmplx( (real(a13,sp)+real(a31,sp))**2/2., &
                 (aimag(a13)+aimag(a31))**2/2.)
  S = 2.*S
  S = cmplx( sqrt(real(S)), sqrt(aimag(S)))

  t11=-cs2*delta**2*cmplx(real(S)*2.*real(a11),aimag(S)*2.*aimag(a11))
  t22=-cs2*delta**2*cmplx(real(S)*2.*real(a22),aimag(S)*2.*aimag(a22))
  t33=-cs2*delta**2*cmplx(real(S)*2.*real(a33),aimag(S)*2.*aimag(a33))
  t12=-cs2*delta**2*cmplx(real(S)*(real(a12)+real(a21)),aimag(S)*(aimag(a12)+aimag(a21)))
  t13=-cs2*delta**2*cmplx(real(S)*(real(a13)+real(a31)),aimag(S)*(aimag(a13)+aimag(a31)))
  t23=-cs2*delta**2*cmplx(real(S)*(real(a23)+real(a32)),aimag(S)*(aimag(a23)+aimag(a32)))

  t11=t11+cnl*delta**2*cmplx(real(a11)*real(a11),aimag(a11)*aimag(a11))
  t11=t11+cnl*delta**2*cmplx(real(a21)*real(a21),aimag(a21)*aimag(a21))
  t11=t11+cnl*delta**2*cmplx(real(a31)*real(a31),aimag(a31)*aimag(a31))

  t22=t22+cnl*delta**2*cmplx(real(a12)*real(a12),aimag(a12)*aimag(a12))
  t22=t22+cnl*delta**2*cmplx(real(a22)*real(a22),aimag(a22)*aimag(a22))
  t22=t22+cnl*delta**2*cmplx(real(a32)*real(a32),aimag(a32)*aimag(a32))

  t33=t33+cnl*delta**2*cmplx(real(a13)*real(a13),aimag(a13)*aimag(a13))
  t33=t33+cnl*delta**2*cmplx(real(a23)*real(a23),aimag(a23)*aimag(a23))
  t33=t33+cnl*delta**2*cmplx(real(a33)*real(a33),aimag(a33)*aimag(a33))


  t12=t12+cnl*delta**2*cmplx(real(a11)*real(a12),aimag(a11)*aimag(a12))
  t12=t12+cnl*delta**2*cmplx(real(a21)*real(a22),aimag(a21)*aimag(a22))
  t12=t12+cnl*delta**2*cmplx(real(a31)*real(a32),aimag(a31)*aimag(a32))

  t13=t13+cnl*delta**2*cmplx(real(a11)*real(a13),aimag(a11)*aimag(a13))
  t13=t13+cnl*delta**2*cmplx(real(a21)*real(a23),aimag(a21)*aimag(a23))
  t13=t13+cnl*delta**2*cmplx(real(a31)*real(a33),aimag(a31)*aimag(a33))

  t23=t23+cnl*delta**2*cmplx(real(a12)*real(a13),aimag(a12)*aimag(a13))
  t23=t23+cnl*delta**2*cmplx(real(a22)*real(a23),aimag(a22)*aimag(a23))
  t23=t23+cnl*delta**2*cmplx(real(a32)*real(a33),aimag(a32)*aimag(a33))

  t11=t11*const; t12=t12*const; t13=t13*const
  t22=t22*const; t23=t23*const; t33=t33*const

  CALL rfftwnd_f77_one_real_to_complex(R2C3D,t11,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,t12,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,t13,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,t22,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,t23,ignore_me)
  CALL rfftwnd_f77_one_real_to_complex(R2C3D,t33,ignore_me)

  a11=-(t11+t22+t33)/3.
  t11=t11+a11
  t22=t22+a11
  t33=t33+a11
end subroutine mixed      
