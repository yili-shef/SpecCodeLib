SUBROUTINE dynsmag_frc (ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,G_test,  &
                        lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,LM,cs2,  &
                        update)
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER,  INTENT(IN)  :: lx1,ly,lz,nx,ny,nz,lx
  REAL(SP), INTENT(IN)  :: delta,delta_test
  LOGICAL,  INTENT(IN)  :: update
  REAL(SP), INTENT(OUT) :: cs2,MM,LM

  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: ux,uy,uz
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz,G_test
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(OUT) :: t11,t12,t13,t22,t23,t33

  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: L11,L12,L13,L22,L23,L33
  COMPLEX(SP), DIMENSIOn(lx1,ly,lz) :: L1,L2,L3,M1,M2,M3
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: uxt,uyt,uzt,S
  REAL(SP) :: const,ignore_me

  const=1._SP/REAL(nx*ny*nz,SP)
  
  IF (update) THEN
    ! L_ij
    t11=ux
    t22=uy
    t33=uz
    t12=G_test*ux
    t13=G_test*uy
    t23=G_test*uz
 
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t33,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t23,ignore_me)
 
    L11=CMPLX(REAL(t11)*REAL(t11),AIMAG(t11)*AIMAG(t11))
    L12=CMPLX(REAL(t11)*REAL(t22),AIMAG(t11)*AIMAG(t22))
    L13=CMPLX(REAL(t11)*REAL(t33),AIMAG(t11)*AIMAG(t33))
    L22=CMPLX(REAL(t22)*REAL(t22),AIMAG(t22)*AIMAG(t22))
    L23=CMPLX(REAL(t22)*REAL(t33),AIMAG(t22)*AIMAG(t33))
    L33=CMPLX(REAL(t33)*REAL(t33),AIMAG(t33)*AIMAG(t33))
 
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L11,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L12,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L13,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L22,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L23,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L33,ignore_me)
    
    L11=G_test*L11*const
    L12=G_test*L12*const
    L13=G_test*L13*const
    L22=G_test*L22*const
    L23=G_test*L23*const
    L33=G_test*L33*const
    
    L1=eye*(kx*L11+ky*L12+kz*L13)
    L2=eye*(kx*L12+ky*L22+kz*L23)
    L3=eye*(kx*L13+ky*L23+kz*L33)
 
    L11=-CMPLX(REAL(t12)*REAL(t12),AIMAG(t12)*AIMAG(t12))
    L12=-CMPLX(REAL(t12)*REAL(t13),AIMAG(t12)*AIMAG(t13))
    L13=-CMPLX(REAL(t12)*REAL(t23),AIMAG(t12)*AIMAG(t23))
    L22=-CMPLX(REAL(t13)*REAL(t13),AIMAG(t13)*AIMAG(t13))
    L23=-CMPLX(REAL(t13)*REAL(t23),AIMAG(t13)*AIMAG(t23))
    L33=-CMPLX(REAL(t23)*REAL(t23),AIMAG(t23)*AIMAG(t23))
    
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L11,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L12,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L13,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L22,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L23,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L33,ignore_me)

    L1=L1+eye*(kx*L11+ky*L12+kz*L13)*const
    L2=L2+eye*(kx*L12+ky*L22+kz*L23)*const
    L3=L3+eye*(kx*L13+ky*L23+kz*L33)*const

    CALL rfftwnd_f77_one_complex_to_real(C2R3D,L1,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,L2,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,L3,ignore_me)

    CALL symmpart(ux,uy,uz,t11,t12,t13,t22,t23,t33,lx1,ly,lz,kx,ky,kz)
 
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t33,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t23,ignore_me)
 
    S=CMPLX(SQRT(2._SP*(REAL(t11)*REAL(t11)+REAL(t22)*REAL(t22)+REAL(t33)*REAL(t33)                     &
                        +2._SP*(REAL(t12)*REAL(t12)+REAL(t13)*REAL(t13)+REAL(t23)*REAL(t23))))          &
            , &
            SQRT(2._SP*(AIMAG(t11)*AIMAG(t11)+AIMAG(t22)*AIMAG(t22)+AIMAG(t33)*AIMAG(t33)               &
                       +2._SP*(AIMAG(t12)*AIMAG(t12)+AIMAG(t13)*AIMAG(t13)+AIMAG(t23)*AIMAG(t23))))     &
           )

    L11=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t11),AIMAG(S)*AIMAG(t11))
    L12=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t12),AIMAG(S)*AIMAG(t12))
    L13=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t13),AIMAG(S)*AIMAG(t13))
    L22=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t22),AIMAG(S)*AIMAG(t22))
    L23=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t23),AIMAG(S)*AIMAG(t23))
    L33=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t33),AIMAG(S)*AIMAG(t33))
    
    t11=L11
    t12=L12
    t13=L13
    t22=L22
    t23=L23
    t33=L33
    
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t11,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t12,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t13,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t22,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t23,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t33,ignore_me)
    
    t11=G_test*t11*const
    t12=G_test*t12*const
    t13=G_test*t13*const
    t22=G_test*t22*const
    t23=G_test*t23*const
    t33=G_test*t33*const
    
    M1=eye*(kx*t11+ky*t12+kz*t13)
    M2=eye*(kx*t12+ky*t22+kz*t23)
    M3=eye*(kx*t13+ky*t23+kz*t33)
    
    uxt=ux*G_test
    uyt=uy*G_test
    uzt=uz*G_test
 
    CALL symmpart(uxt,uyt,uzt,t11,t12,t13,t22,t23,t33,lx1,ly,lz,kx,ky,kz)
 
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t33,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t23,ignore_me)
                       
    S=CMPLX(SQRT(2._SP*(REAL(t11)*REAL(t11)+REAL(t22)*REAL(t22)+REAL(t33)*REAL(t33)                     &
                        +2._SP*(REAL(t12)*REAL(t12)+REAL(t13)*REAL(t13)+REAL(t23)*REAL(t23))))          &
           , &
            SQRT(2._SP*(AIMAG(t11)*AIMAG(t11)+AIMAG(t22)*AIMAG(t22)+AIMAG(t33)*AIMAG(t33)               &
                        +2._SP*(AIMAG(t12)*AIMAG(t12)+AIMAG(t13)*AIMAG(t13)+AIMAG(t23)*AIMAG(t23))))    &
           )
 
    t11=-2._SP*delta_test*delta_test*CMPLX(REAL(S)*REAL(t11),AIMAG(S)*AIMAG(t11))
    t12=-2._SP*delta_test*delta_test*CMPLX(REAL(S)*REAL(t12),AIMAG(S)*AIMAG(t12))
    t13=-2._SP*delta_test*delta_test*CMPLX(REAL(S)*REAL(t13),AIMAG(S)*AIMAG(t13))
    t22=-2._SP*delta_test*delta_test*CMPLX(REAL(S)*REAL(t22),AIMAG(S)*AIMAG(t22))
    t23=-2._SP*delta_test*delta_test*CMPLX(REAL(S)*REAL(t23),AIMAG(S)*AIMAG(t23))
    t33=-2._SP*delta_test*delta_test*CMPLX(REAL(S)*REAL(t33),AIMAG(S)*AIMAG(t33))
    
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t11,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t12,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t13,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t22,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t23,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t33,ignore_me)
    
    M1=M1+eye*(kx*t11+ky*t12+kz*t13)*const
    M2=M2+eye*(kx*t12+ky*t22+kz*t23)*const
    M3=M3+eye*(kx*t13+ky*t23+kz*t33)*const
 
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,M1,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,M2,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,M3,ignore_me)

END SUBROUTINE dynsmag_frc      
