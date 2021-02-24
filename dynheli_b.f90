SUBROUTINE dynheli_b (ux,uy,uz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz, &
                      G_test,lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,PM,&
                      PP,LM,LP,c1,c2,update)

  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER,  INTENT(IN)  :: lx1,ly,lz,nx,ny,nz,lx
  REAL(SP), INTENT(IN)  :: delta,delta_test
  LOGICAL,  INTENT(IN)  :: update
  REAL(SP), INTENT(OUT) :: c1,c2,MM,PM,PP,LM,LP
  
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: ux,uy,uz,wx,wy,wz
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz,G_test
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(OUT) :: t11,t12,t13,t22,t23,t33

  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: L11,L12,L13,L22,L23,L33
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: S11,S12,S13,S22,S23,S33
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: R11,R12,R13,R22,R23,R33
  COMPLEX(SP), DIMENSIOn(lx1,ly,lz) :: L1,L2,L3,M1,M2,M3,P1,P2,P3
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: uxt,uyt,uzt,S,S_test
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

    ! M_ij and P_ij
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

    S11=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t11),AIMAG(S)*AIMAG(t11))
    S12=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t12),AIMAG(S)*AIMAG(t12))
    S13=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t13),AIMAG(S)*AIMAG(t13))
    S22=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t22),AIMAG(S)*AIMAG(t22))
    S23=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t23),AIMAG(S)*AIMAG(t23))
    S33=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t33),AIMAG(S)*AIMAG(t33))
    
    t11=S11
    t12=S12
    t13=S13
    t22=S22
    t23=S23
    t33=S33
    
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
    
    CALL symmpart(wx,wy,wz,t11,t12,t13,t22,t23,t33,lx1,ly,lz,kx,ky,kz)
 
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t33,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t23,ignore_me)
 
    R11=delta*delta*delta*CMPLX(REAL(S)*REAL(t11),AIMAG(S)*AIMAG(t11))
    R12=delta*delta*delta*CMPLX(REAL(S)*REAL(t12),AIMAG(S)*AIMAG(t12))
    R13=delta*delta*delta*CMPLX(REAL(S)*REAL(t13),AIMAG(S)*AIMAG(t13))
    R22=delta*delta*delta*CMPLX(REAL(S)*REAL(t22),AIMAG(S)*AIMAG(t22))
    R23=delta*delta*delta*CMPLX(REAL(S)*REAL(t23),AIMAG(S)*AIMAG(t23))
    R33=delta*delta*delta*CMPLX(REAL(S)*REAL(t33),AIMAG(S)*AIMAG(t33))

    t11=R11
    t12=R12
    t13=R13
    t22=R22
    t23=R23
    t33=R33
    
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

    P1=eye*(kx*t11+ky*t12+kz*t13)
    P2=eye*(kx*t12+ky*t22+kz*t23)
    P3=eye*(kx*t13+ky*t23+kz*t33)

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
                       
    S_test=CMPLX(SQRT(2._SP*(REAL(t11)*REAL(t11)+REAL(t22)*REAL(t22)+REAL(t33)*REAL(t33)                     &
                             +2._SP*(REAL(t12)*REAL(t12)+REAL(t13)*REAL(t13)+REAL(t23)*REAL(t23))))          &
                , &
                 SQRT(2._SP*(AIMAG(t11)*AIMAG(t11)+AIMAG(t22)*AIMAG(t22)+AIMAG(t33)*AIMAG(t33)               &
                             +2._SP*(AIMAG(t12)*AIMAG(t12)+AIMAG(t13)*AIMAG(t13)+AIMAG(t23)*AIMAG(t23))))    &
                )
 
    L11=-2._SP*delta_test*delta_test*CMPLX(REAL(S_test)*REAL(t11),AIMAG(S_test)*AIMAG(t11))
    L12=-2._SP*delta_test*delta_test*CMPLX(REAL(S_test)*REAL(t12),AIMAG(S_test)*AIMAG(t12))
    L13=-2._SP*delta_test*delta_test*CMPLX(REAL(S_test)*REAL(t13),AIMAG(S_test)*AIMAG(t13))
    L22=-2._SP*delta_test*delta_test*CMPLX(REAL(S_test)*REAL(t22),AIMAG(S_test)*AIMAG(t22))
    L23=-2._SP*delta_test*delta_test*CMPLX(REAL(S_test)*REAL(t23),AIMAG(S_test)*AIMAG(t23))
    L33=-2._SP*delta_test*delta_test*CMPLX(REAL(S_test)*REAL(t33),AIMAG(S_test)*AIMAG(t33))
    
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L11,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L12,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L13,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L22,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L23,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L33,ignore_me)
    
    M1=M1+eye*(kx*L11+ky*L12+kz*L13)*const
    M2=M2+eye*(kx*L12+ky*L22+kz*L23)*const
    M3=M3+eye*(kx*L13+ky*L23+kz*L33)*const
 
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,M1,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,M2,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,M3,ignore_me)

    uxt=wx*G_test
    uyt=wy*G_test
    uzt=wz*G_test
 
    CALL symmpart(uxt,uyt,uzt,t11,t12,t13,t22,t23,t33,lx1,ly,lz,kx,ky,kz)
 
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t33,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,t23,ignore_me)
    
    L11=-delta_test**3*CMPLX(REAL(S_test)*REAL(t11),AIMAG(S_test)*AIMAG(t11))
    L12=-delta_test**3*CMPLX(REAL(S_test)*REAL(t12),AIMAG(S_test)*AIMAG(t12))
    L13=-delta_test**3*CMPLX(REAL(S_test)*REAL(t13),AIMAG(S_test)*AIMAG(t13))
    L22=-delta_test**3*CMPLX(REAL(S_test)*REAL(t22),AIMAG(S_test)*AIMAG(t22))
    L23=-delta_test**3*CMPLX(REAL(S_test)*REAL(t23),AIMAG(S_test)*AIMAG(t23))
    L33=-delta_test**3*CMPLX(REAL(S_test)*REAL(t33),AIMAG(S_test)*AIMAG(t33))

    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L11,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L12,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L13,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L22,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L23,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,L33,ignore_me)
    
    P1=P1+eye*(kx*L11+ky*L12+kz*L13)*const
    P2=P2+eye*(kx*L12+ky*L22+kz*L23)*const
    P3=P3+eye*(kx*L13+ky*L23+kz*L33)*const
 
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,P1,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,P2,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,P3,ignore_me)

    t11=CMPLX(REAL(L1)*REAL(M1)+REAL(L2)*REAL(M2)+REAL(L3)*REAL(M3),      &
              AIMAG(L1)*AIMAG(M1)+AIMAG(L2)*AIMAG(M2)+AIMAG(L3)*AIMAG(M3) &
             )
    LM=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    t11=CMPLX(REAL(L1)*REAL(P1)+REAL(L2)*REAL(P2)+REAL(L3)*REAL(P3),      &
              AIMAG(L1)*AIMAG(P1)+AIMAG(L2)*AIMAG(P2)+AIMAG(L3)*AIMAG(P3) &
             )
    LP=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    t11=CMPLX(REAL(P1)*REAL(P1)+REAL(P2)*REAL(P2)+REAL(P3)*REAL(P3),      &
              AIMAG(P1)*AIMAG(P1)+AIMAG(P2)*AIMAG(P2)+AIMAG(P3)*AIMAG(P3) &
             )
    PP=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    t11=CMPLX(REAL(M1)*REAL(M1)+REAL(M2)*REAL(M2)+REAL(M3)*REAL(M3),      &
              AIMAG(M1)*AIMAG(M1)+AIMAG(M2)*AIMAG(M2)+AIMAG(M3)*AIMAG(M3) &
             )
    MM=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    t11=CMPLX(REAL(P1)*REAL(M1)+REAL(P2)*REAL(M2)+REAL(P3)*REAL(M3),      &
              AIMAG(P1)*AIMAG(M1)+AIMAG(P2)*AIMAG(M2)+AIMAG(P3)*AIMAG(M3) &
             )
    PM=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const


    c1=(LM*PP-LP*PM)/(MM*PP-PM*PM)
    c2=(MM*LP-PM*LM)/(MM*PP-PM*PM)
    
    t11=(-c1*S11-c2*R11)*const
    t12=(-c1*S12-c2*R12)*const
    t13=(-c1*S13-c2*R13)*const
    t22=(-c1*S22-c2*R22)*const
    t23=(-c1*S23-c2*R23)*const
    t33=(-c1*S33-c2*R33)*const
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t11,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t12,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t13,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t22,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t23,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t33,ignore_me)
             
  ELSE
    CALL symmpart(ux,uy,uz,S11,S12,S13,S22,S23,S33,lx1,ly,lz,kx,ky,kz)
    CALL symmpart(wx,wy,wz,R11,R12,R13,R22,R23,R33,lx1,ly,lz,kx,ky,kz)

    CALL rfftwnd_f77_one_complex_to_real(C2R3D,S11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,S22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,S33,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,S12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,S13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,S23,ignore_me)
    
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,R11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,R22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,R33,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,R12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,R13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,R23,ignore_me)

    S=CMPLX(SQRT(2._SP*(REAL(S11)*REAL(S11)+REAL(S22)*REAL(S22)+REAL(S33)*REAL(S33)                     &
                        +2._SP*(REAL(S12)*REAL(S12)+REAL(S13)*REAL(S13)+REAL(S23)*REAL(S23))))          &
            , &
            SQRT(2._SP*(AIMAG(S11)*AIMAG(S11)+AIMAG(S22)*AIMAG(S22)+AIMAG(S33)*AIMAG(S33)               &
                       +2._SP*(AIMAG(S12)*AIMAG(S12)+AIMAG(S13)*AIMAG(S13)+AIMAG(S23)*AIMAG(S23))))     &
           )

    t11=-c1*2._SP*delta*delta*CMPLX(REAL(S)*REAL(S11),AIMAG(S)*AIMAG(S11))
    t12=-c1*2._SP*delta*delta*CMPLX(REAL(S)*REAL(S12),AIMAG(S)*AIMAG(S12))
    t13=-c1*2._SP*delta*delta*CMPLX(REAL(S)*REAL(S13),AIMAG(S)*AIMAG(S13))
    t22=-c1*2._SP*delta*delta*CMPLX(REAL(S)*REAL(S22),AIMAG(S)*AIMAG(S22))
    t23=-c1*2._SP*delta*delta*CMPLX(REAL(S)*REAL(S23),AIMAG(S)*AIMAG(S23))
    t33=-c1*2._SP*delta*delta*CMPLX(REAL(S)*REAL(S33),AIMAG(S)*AIMAG(S33))
    
    t11=t11-c2*delta*delta*delta*CMPLX(REAL(S)*REAL(R11),AIMAG(S)*AIMAG(R11))
    t12=t12-c2*delta*delta*delta*CMPLX(REAL(S)*REAL(R12),AIMAG(S)*AIMAG(R12))
    t13=t13-c2*delta*delta*delta*CMPLX(REAL(S)*REAL(R13),AIMAG(S)*AIMAG(R13))
    t22=t22-c2*delta*delta*delta*CMPLX(REAL(S)*REAL(R22),AIMAG(S)*AIMAG(R22))
    t23=t23-c2*delta*delta*delta*CMPLX(REAL(S)*REAL(R23),AIMAG(S)*AIMAG(R23))
    t33=t33-c2*delta*delta*delta*CMPLX(REAL(S)*REAL(R33),AIMAG(S)*AIMAG(R33))
    
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t11,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t12,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t13,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t22,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t23,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t33,ignore_me)
    t11=t11*const
    t12=t12*const
    t13=t13*const
    t22=t22*const
    t23=t23*const
    t33=t33*const
  END IF
END SUBROUTINE dynheli_b
