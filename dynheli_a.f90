SUBROUTINE dynheli_a (ux,uy,uz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz, &
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
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: M11,M12,M13,M22,M23,M33
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: P11,P12,P13,P22,P23,P33
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: S11,S12,S13,S22,S23,S33
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: R11,R12,R13,R22,R23,R33
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
 
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,L11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,L22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,L33,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,L12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,L13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,L23,ignore_me)
 
    L11=L11-CMPLX(REAL(t12)*REAL(t12),AIMAG(t12)*AIMAG(t12))
    L12=L12-CMPLX(REAL(t12)*REAL(t13),AIMAG(t12)*AIMAG(t13))
    L13=L13-CMPLX(REAL(t12)*REAL(t23),AIMAG(t12)*AIMAG(t23))
    L22=L22-CMPLX(REAL(t13)*REAL(t13),AIMAG(t13)*AIMAG(t13))
    L23=L23-CMPLX(REAL(t13)*REAL(t23),AIMAG(t13)*AIMAG(t23))
    L33=L33-CMPLX(REAL(t23)*REAL(t23),AIMAG(t23)*AIMAG(t23))

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
    
    M11=S11
    M12=S12
    M13=S13
    M22=S22
    M23=S23
    M33=S33
    
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,M11,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,M12,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,M13,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,M22,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,M23,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,M33,ignore_me)
    
    M11=G_test*M11*const
    M12=G_test*M12*const
    M13=G_test*M13*const
    M22=G_test*M22*const
    M23=G_test*M23*const
    M33=G_test*M33*const
    
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,M11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,M22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,M33,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,M12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,M13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,M23,ignore_me)
 
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

    P11=R11
    P12=R12
    P13=R13
    P22=R22
    P23=R23
    P33=R33
    
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,P11,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,P12,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,P13,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,P22,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,P23,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,P33,ignore_me)
    
    P11=G_test*P11*const
    P12=G_test*P12*const
    P13=G_test*P13*const
    P22=G_test*P22*const
    P23=G_test*P23*const
    P33=G_test*P33*const

    CALL rfftwnd_f77_one_complex_to_real(C2R3D,P11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,P22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,P33,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,P12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,P13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,P23,ignore_me)
 
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
 
    M11=M11-2._SP*delta_test*delta_test*CMPLX(REAL(S_test)*REAL(t11),AIMAG(S_test)*AIMAG(t11))
    M12=M12-2._SP*delta_test*delta_test*CMPLX(REAL(S_test)*REAL(t12),AIMAG(S_test)*AIMAG(t12))
    M13=M13-2._SP*delta_test*delta_test*CMPLX(REAL(S_test)*REAL(t13),AIMAG(S_test)*AIMAG(t13))
    M22=M22-2._SP*delta_test*delta_test*CMPLX(REAL(S_test)*REAL(t22),AIMAG(S_test)*AIMAG(t22))
    M23=M23-2._SP*delta_test*delta_test*CMPLX(REAL(S_test)*REAL(t23),AIMAG(S_test)*AIMAG(t23))
    M33=M33-2._SP*delta_test*delta_test*CMPLX(REAL(S_test)*REAL(t33),AIMAG(S_test)*AIMAG(t33))
 
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
    
    P11=P11-delta_test**3*CMPLX(REAL(S_test)*REAL(t11),AIMAG(S_test)*AIMAG(t11))
    P12=P12-delta_test**3*CMPLX(REAL(S_test)*REAL(t12),AIMAG(S_test)*AIMAG(t12))
    P13=P13-delta_test**3*CMPLX(REAL(S_test)*REAL(t13),AIMAG(S_test)*AIMAG(t13))
    P22=P22-delta_test**3*CMPLX(REAL(S_test)*REAL(t22),AIMAG(S_test)*AIMAG(t22))
    P23=P23-delta_test**3*CMPLX(REAL(S_test)*REAL(t23),AIMAG(S_test)*AIMAG(t23))
    P33=P33-delta_test**3*CMPLX(REAL(S_test)*REAL(t33),AIMAG(S_test)*AIMAG(t33))

    t11=CMPLX(REAL(M11)*REAL(M11)+REAL(M22)*REAL(M22)+REAL(M33)*REAL(M33)+        &
              2._SP*(REAL(M12)*REAL(M12)+REAL(M13)*REAL(M13)+REAL(M23)*REAL(M23)) &
             , &
              AIMAG(M11)*AIMAG(M11)+AIMAG(M22)*AIMAG(M22)+AIMAG(M33)*AIMAG(M33)+        &
              2._SP*(AIMAG(M12)*AIMAG(M12)+AIMAG(M13)*AIMAG(M13)+AIMAG(M23)*AIMAG(M23)) &
             )
    MM=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    t11=CMPLX(REAL(P11)*REAL(P11)+REAL(P22)*REAL(P22)+REAL(P33)*REAL(P33)+        &
              2._SP*(REAL(P12)*REAL(P12)+REAL(P13)*REAL(P13)+REAL(P23)*REAL(P23)) &
             , &
              AIMAG(P11)*AIMAG(P11)+AIMAG(P22)*AIMAG(P22)+AIMAG(P33)*AIMAG(P33)+        &
              2._SP*(AIMAG(P12)*AIMAG(P12)+AIMAG(P13)*AIMAG(P13)+AIMAG(P23)*AIMAG(P23)) &
             )
    PP=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    t11=CMPLX(REAL(P11)*REAL(M11)+REAL(P22)*REAL(M22)+REAL(P33)*REAL(M33)+        &
              2._SP*(REAL(P12)*REAL(M12)+REAL(P13)*REAL(M13)+REAL(P23)*REAL(M23)) &
             , &
              AIMAG(P11)*AIMAG(M11)+AIMAG(P22)*AIMAG(M22)+AIMAG(P33)*AIMAG(M33)+        &
              2._SP*(AIMAG(P12)*AIMAG(M12)+AIMAG(P13)*AIMAG(M13)+AIMAG(P23)*AIMAG(M23)) &
             )
    PM=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    t11=CMPLX(REAL(L11)*REAL(M11)+REAL(L22)*REAL(M22)+REAL(L33)*REAL(M33)+        &
              2._SP*(REAL(L12)*REAL(M12)+REAL(L13)*REAL(M13)+REAL(L23)*REAL(M23)) &
             , &
              AIMAG(L11)*AIMAG(M11)+AIMAG(L22)*AIMAG(M22)+AIMAG(L33)*AIMAG(M33)+        &
              2._SP*(AIMAG(L12)*AIMAG(M12)+AIMAG(L13)*AIMAG(M13)+AIMAG(L23)*AIMAG(M23)) &
             )
    LM=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    t11=CMPLX(REAL(L11)*REAL(P11)+REAL(L22)*REAL(P22)+REAL(L33)*REAL(P33)+        &
              2._SP*(REAL(L12)*REAL(P12)+REAL(L13)*REAL(P13)+REAL(L23)*REAL(P23)) &
             , &
              AIMAG(L11)*AIMAG(P11)+AIMAG(L22)*AIMAG(P22)+AIMAG(L33)*AIMAG(P33)+        &
              2._SP*(AIMAG(L12)*AIMAG(P12)+AIMAG(L13)*AIMAG(P13)+AIMAG(L23)*AIMAG(P23)) &
             )
    LP=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

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



END SUBROUTINE dynheli_a
