SUBROUTINE nonlinear (ux,uy,uz,wx,wy,wz,t11,t12,t13,t22,t23,t33,kx,ky,kz,G_test,    &
                      lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,MM,NM,NN,LM,LN,c1,c2,  &
                      update)
  USE mconstant
  USE mfftwplan3d
  IMPLICIT NONE

  INTEGER,  INTENT(IN)  :: lx1,ly,lz,nx,ny,nz,lx
  REAL(SP), INTENT(IN)  :: delta,delta_test
  LOGICAL,  INTENT(IN)  :: update
  REAL(SP), INTENT(OUT) :: c1,c2,MM,NM,NN,LM,LN
  
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(IN)  :: ux,uy,uz,wx,wy,wz
  REAL(SP),    DIMENSION(lx1,ly,lz), INTENT(IN)  :: kx,ky,kz,G_test
  COMPLEX(SP), DIMENSION(lx1,ly,lz), INTENT(OUT) :: t11,t12,t13,t22,t23,t33

  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: L11,L12,L13,L22,L23,L33
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: M11,M12,M13,M22,M23,M33
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: N11,N12,N13,N22,N23,N33
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: S11,S12,S13,S22,S23,S33
  COMPLEX(SP), DIMENSION(lx1,ly,lz) :: AA11,AA12,AA13,AA22,AA23,AA33
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

    ! M_ij and N_ij
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

    M11=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t11),AIMAG(S)*AIMAG(t11))
    M12=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t12),AIMAG(S)*AIMAG(t12))
    M13=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t13),AIMAG(S)*AIMAG(t13))
    M22=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t22),AIMAG(S)*AIMAG(t22))
    M23=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t23),AIMAG(S)*AIMAG(t23))
    M33=2._SP*delta*delta*CMPLX(REAL(S)*REAL(t33),AIMAG(S)*AIMAG(t33))
    
    S11=t11
    S12=t12
    S13=t13
    S22=t22
    S23=t23
    S33=t33
    
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

    uxt=.5_SP*wx
    uyt=.5_SP*wy
    uzt=.5_SP*wz
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,uxt,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,uyt,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,uzt,ignore_me)

    AA11=CMPLX( REAL(t11)*REAL(t11)+REAL(t12-uzt)*REAL(t12-uzt)+REAL(t13+uyt)*REAL(t13+uyt),      &
                AIMAG(t11)*AIMAG(t11)+AIMAG(t12-uzt)*AIMAG(t12-uzt)+AIMAG(t13+uyt)*AIMAG(t13+uyt) &
              )
    AA22=CMPLX( REAL(t12+uzt)*REAL(t12+uzt)+REAL(t22)*REAL(t22)+REAL(t23-uxt)*REAL(t23-uxt),      &
                AIMAG(t12+uzt)*AIMAG(t12+uzt)+AIMAG(t22)*AIMAG(t22)+AIMAG(t23-uxt)*AIMAG(t23-uxt) &
              )
    AA33=CMPLX( REAL(t13-uyt)*REAL(t13-uyt)+REAL(t12+uxt)*REAL(t12+uxt)+REAL(t33)*REAL(t33),      &
                AIMAG(t13-uyt)*AIMAG(t13-uyt)+AIMAG(t12+uxt)*AIMAG(t12+uxt)+AIMAG(t33)*AIMAG(t33) &
              )
    AA12=CMPLX( REAL(t11)*REAL(t12+uzt)+REAL(t12-uzt)*REAL(t22)+REAL(t13+uyt)*REAL(t23-uxt),      &
                AIMAG(t11)*AIMAG(t12+uzt)+AIMAG(t12-uzt)*AIMAG(t22)+AIMAG(t13+uyt)*AIMAG(t23-uxt) &
              )
    AA13=CMPLX( REAL(t11)*REAL(t13-uyt)+REAL(t12-uzt)*REAL(t23+uxt)+REAL(t13+uyt)*REAL(t33),      &
                AIMAG(t11)*AIMAG(t13-uyt)+AIMAG(t12-uzt)*AIMAG(t23+uxt)+AIMAG(t13+uyt)*AIMAG(t33) &
              )
    AA23=CMPLX( REAL(t12+uzt)*REAL(t13-uyt)+REAL(t22)*REAL(t23+uxt)+REAL(t23-uxt)*REAL(t33),      &
                AIMAG(t12+uzt)*AIMAG(t13-uyt)+AIMAG(t22)*AIMAG(t23+uxt)+AIMAG(t23-uxt)*AIMAG(t33) &
              )

    N11=delta*delta*AA11
    N12=delta*delta*AA12
    N13=delta*delta*AA13
    N22=delta*delta*AA22
    N23=delta*delta*AA23
    N33=delta*delta*AA33

    CALL rfftwnd_f77_one_real_to_complex(R2C3D,N11,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,N12,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,N13,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,N22,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,N23,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,N33,ignore_me)

    N11=G_test*N11*const
    N12=G_test*N12*const
    N13=G_test*N13*const
    N22=G_test*N22*const
    N23=G_test*N23*const
    N33=G_test*N33*const
    
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,N11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,N22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,N33,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,N12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,N13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,N23,ignore_me)
 
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
 
    uxt=.5_SP*wx*G_test
    uyt=.5_SP*wy*G_test
    uzt=.5_SP*wz*G_test
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,uxt,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,uyt,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,uzt,ignore_me)

    N11=delta_test*delta_test*                                                                   &
        CMPLX( REAL(t11)*REAL(t11)+REAL(t12-uzt)*REAL(t12-uzt)+REAL(t13+uyt)*REAL(t13+uyt),      &
               AIMAG(t11)*AIMAG(t11)+AIMAG(t12-uzt)*AIMAG(t12-uzt)+AIMAG(t13+uyt)*AIMAG(t13+uyt) &
             ) - N11
    N22=delta_test*delta_test*                                                                   &
        CMPLX( REAL(t12+uzt)*REAL(t12+uzt)+REAL(t22)*REAL(t22)+REAL(t23-uxt)*REAL(t23-uxt),      &
               AIMAG(t12+uzt)*AIMAG(t12+uzt)+AIMAG(t22)*AIMAG(t22)+AIMAG(t23-uxt)*AIMAG(t23-uxt) &
             ) - N22
    N33=delta_test*delta_test*                                                                   &
        CMPLX( REAL(t13-uyt)*REAL(t13-uyt)+REAL(t12+uxt)*REAL(t12+uxt)+REAL(t33)*REAL(t33),      &
               AIMAG(t13-uyt)*AIMAG(t13-uyt)+AIMAG(t12+uxt)*AIMAG(t12+uxt)+AIMAG(t33)*AIMAG(t33) &
             ) - N33
    N12=delta_test*delta_test*                                                                   &
        CMPLX( REAL(t11)*REAL(t12+uzt)+REAL(t12-uzt)*REAL(t22)+REAL(t13+uyt)*REAL(t23-uxt),      &
               AIMAG(t11)*AIMAG(t12+uzt)+AIMAG(t12-uzt)*AIMAG(t22)+AIMAG(t13+uyt)*AIMAG(t23-uxt) &
             ) - N12
    N13=delta_test*delta_test*                                                                   &
        CMPLX( REAL(t11)*REAL(t13-uyt)+REAL(t12-uzt)*REAL(t23+uxt)+REAL(t13+uyt)*REAL(t33),      &
               AIMAG(t11)*AIMAG(t13-uyt)+AIMAG(t12-uzt)*AIMAG(t23+uxt)+AIMAG(t13+uyt)*AIMAG(t33) &
             ) - N13
    N23=delta_test*delta_test*                                                                   &
        CMPLX( REAL(t12+uzt)*REAL(t13-uyt)+REAL(t22)*REAL(t23+uxt)+REAL(t23-uxt)*REAL(t33),      &
               AIMAG(t12+uzt)*AIMAG(t13-uyt)+AIMAG(t22)*AIMAG(t23+uxt)+AIMAG(t23-uxt)*AIMAG(t33) &
             ) - N23


    S_test=(N11+N22+N33)/3._SP
    N11=N11-S_test
    N22=N22-S_test
    N33=N33-S_test

    !coefficients
    t11=CMPLX(REAL(M11)*REAL(M11)+REAL(M22)*REAL(M22)+REAL(M33)*REAL(M33)+        &
              2._SP*(REAL(M12)*REAL(M12)+REAL(M13)*REAL(M13)+REAL(M23)*REAL(M23)) &
             , &
              AIMAG(M11)*AIMAG(M11)+AIMAG(M22)*AIMAG(M22)+AIMAG(M33)*AIMAG(M33)+        &
              2._SP*(AIMAG(M12)*AIMAG(M12)+AIMAG(M13)*AIMAG(M13)+AIMAG(M23)*AIMAG(M23)) &
             )
    MM=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    t11=CMPLX(REAL(N11)*REAL(N11)+REAL(N22)*REAL(N22)+REAL(N33)*REAL(N33)+        &
              2._SP*(REAL(N12)*REAL(N12)+REAL(N13)*REAL(N13)+REAL(N23)*REAL(N23)) &
             , &
              AIMAG(N11)*AIMAG(N11)+AIMAG(N22)*AIMAG(N22)+AIMAG(N33)*AIMAG(N33)+        &
              2._SP*(AIMAG(N12)*AIMAG(N12)+AIMAG(N13)*AIMAG(N13)+AIMAG(N23)*AIMAG(N23)) &
             )
    NN=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    t11=CMPLX(REAL(N11)*REAL(M11)+REAL(N22)*REAL(M22)+REAL(N33)*REAL(M33)+        &
              2._SP*(REAL(N12)*REAL(M12)+REAL(N13)*REAL(M13)+REAL(N23)*REAL(M23)) &
             , &
              AIMAG(N11)*AIMAG(M11)+AIMAG(N22)*AIMAG(M22)+AIMAG(N33)*AIMAG(M33)+        &
              2._SP*(AIMAG(N12)*AIMAG(M12)+AIMAG(N13)*AIMAG(M13)+AIMAG(N23)*AIMAG(M23)) &
             )
    NM=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    t11=CMPLX(REAL(L11)*REAL(M11)+REAL(L22)*REAL(M22)+REAL(L33)*REAL(M33)+        &
              2._SP*(REAL(L12)*REAL(M12)+REAL(L13)*REAL(M13)+REAL(L23)*REAL(M23)) &
             , &
              AIMAG(L11)*AIMAG(M11)+AIMAG(L22)*AIMAG(M22)+AIMAG(L33)*AIMAG(M33)+        &
              2._SP*(AIMAG(L12)*AIMAG(M12)+AIMAG(L13)*AIMAG(M13)+AIMAG(L23)*AIMAG(M23)) &
             )
    LM=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    t11=CMPLX(REAL(L11)*REAL(N11)+REAL(L22)*REAL(N22)+REAL(L33)*REAL(N33)+        &
              2._SP*(REAL(L12)*REAL(N12)+REAL(L13)*REAL(N13)+REAL(L23)*REAL(N23)) &
             , &
              AIMAG(L11)*AIMAG(N11)+AIMAG(L22)*AIMAG(N22)+AIMAG(L33)*AIMAG(N33)+        &
              2._SP*(AIMAG(L12)*AIMAG(N12)+AIMAG(L13)*AIMAG(N13)+AIMAG(L23)*AIMAG(N23)) &
             )
    LN=(SUM(REAL(t11(1:lx,:,:)))+SUM(AIMAG(t11(1:lx,:,:))))*const

    c1=(LM*NN-LN*NM)/(MM*NN-NM*NM)
    c2=(MM*LN-NM*LM)/(MM*NN-NM*NM)

    !stresses
    S_test=(AA11+AA22+AA33)/3._SP
    AA11=AA11-S_test
    AA22=AA22-S_test
    AA33=AA33-S_test

    t11=(-2._SP*c1*delta*delta*S*S11+c2*delta*delta*AA11)*const
    t12=(-2._SP*c1*delta*delta*S*S12+c2*delta*delta*AA12)*const
    t13=(-2._SP*c1*delta*delta*S*S13+c2*delta*delta*AA13)*const
    t22=(-2._SP*c1*delta*delta*S*S22+c2*delta*delta*AA22)*const
    t23=(-2._SP*c1*delta*delta*S*S23+c2*delta*delta*AA23)*const
    t33=(-2._SP*c1*delta*delta*S*S33+c2*delta*delta*AA33)*const
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t11,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t12,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t13,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t22,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t23,ignore_me)
    CALL rfftwnd_f77_one_real_to_complex(R2C3D,t33,ignore_me)

  ELSE
    uxt=.5_SP*wx
    uyt=.5_SP*wy
    uzt=.5_SP*wz
    CALL symmpart(ux,uy,uz,S11,S12,S13,S22,S23,S33,lx1,ly,lz,kx,ky,kz)

    CALL rfftwnd_f77_one_complex_to_real(C2R3D,S11,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,S22,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,S33,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,S12,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,S13,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,S23,ignore_me)
    
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,uxt,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,uyt,ignore_me)
    CALL rfftwnd_f77_one_complex_to_real(C2R3D,uzt,ignore_me)

    AA11=CMPLX( REAL(S11)*REAL(S11)+REAL(S12-uzt)*REAL(S12-uzt)+REAL(S13+uyt)*REAL(S13+uyt),      &
                AIMAG(S11)*AIMAG(S11)+AIMAG(S12-uzt)*AIMAG(S12-uzt)+AIMAG(S13+uyt)*AIMAG(S13+uyt) &
              )
    AA22=CMPLX( REAL(S12+uzt)*REAL(S12+uzt)+REAL(S22)*REAL(S22)+REAL(S23-uxt)*REAL(S23-uxt),      &
                AIMAG(S12+uzt)*AIMAG(S12+uzt)+AIMAG(S22)*AIMAG(S22)+AIMAG(S23-uxt)*AIMAG(S23-uxt) &
              )
    AA33=CMPLX( REAL(S13-uyt)*REAL(S13-uyt)+REAL(S12+uxt)*REAL(S12+uxt)+REAL(S33)*REAL(S33),      &
                AIMAG(S13-uyt)*AIMAG(S13-uyt)+AIMAG(S12+uxt)*AIMAG(S12+uxt)+AIMAG(S33)*AIMAG(S33) &
              )
    AA12=CMPLX( REAL(S11)*REAL(S12+uzt)+REAL(S12-uzt)*REAL(S22)+REAL(S13+uyt)*REAL(S23-uxt),      &
                AIMAG(S11)*AIMAG(S12+uzt)+AIMAG(S12-uzt)*AIMAG(S22)+AIMAG(S13+uyt)*AIMAG(S23-uxt) &
              )
    AA13=CMPLX( REAL(S11)*REAL(S13-uyt)+REAL(S12-uzt)*REAL(S23+uxt)+REAL(S13+uyt)*REAL(S33),      &
                AIMAG(S11)*AIMAG(S13-uyt)+AIMAG(S12-uzt)*AIMAG(S23+uxt)+AIMAG(S13+uyt)*AIMAG(S33) &
              )
    AA23=CMPLX( REAL(S12+uzt)*REAL(S13-uyt)+REAL(S22)*REAL(S23+uxt)+REAL(S23-uxt)*REAL(S33),      &
                AIMAG(S12+uzt)*AIMAG(S13-uyt)+AIMAG(S22)*AIMAG(S23+uxt)+AIMAG(S23-uxt)*AIMAG(S33) &
              )
              
    S_test=(AA11+AA22+AA33)/3._SP
    AA11=AA11-S_test
    AA22=AA22-S_test
    AA33=AA33-S_test
    
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
    
    t11=t11+c2*delta*delta*AA11
    t12=t12+c2*delta*delta*AA12
    t13=t13+c2*delta*delta*AA13
    t22=t22+c2*delta*delta*AA22
    t23=t23+c2*delta*delta*AA23
    t33=t33+c2*delta*delta*AA33
    
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

END SUBROUTINE nonlinear      
