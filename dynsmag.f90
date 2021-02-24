subroutine dynsmag (ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,g_test,  &
                    lx1,ly,lz,nx,ny,nz,lx,delta,delta_test,mm,lm,cs2,  &
                    update)
  use mconstant
  use mfftwplan3d
  implicit none

  integer,  intent(in)  :: lx1,ly,lz,nx,ny,nz,lx
  real(sp), intent(in)  :: delta,delta_test
  logical,  intent(in)  :: update
  real(sp), intent(out) :: cs2,mm,lm
  
  complex(sp), dimension(lx1,ly,lz), intent(in)  :: ux,uy,uz
  real(sp),    dimension(lx1,ly,lz), intent(in)  :: kx,ky,kz,g_test
  complex(sp), dimension(lx1,ly,lz), intent(out) :: t11,t12,t13,t22,t23,t33

  complex(sp), dimension(lx1,ly,lz) :: L11,L12,L13,L22,L23,L33
  complex(sp), dimension(lx1,ly,lz) :: m11,m12,m13,m22,m23,m33
  complex(sp), dimension(lx1,ly,lz) :: s11,s12,s13,s22,s23,s33
  complex(sp), dimension(lx1,ly,lz) :: s
  real(sp) :: const,ignore_me, delta2, delta_test2

  const=1._sp/real(nx*ny*nz,sp)
  delta2 = delta * delta 
  delta_test2 = delta_test * delta_test
  
  if (update) then
    ! L_ij
    t11=ux
    t22=uy
    t33=uz
    t12=g_test*ux
    t13=g_test*uy
    t23=g_test*uz
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,t11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t23,ignore_me)
 
    L11=cmplx(real(t11)*real(t11),aimag(t11)*aimag(t11))
    L12=cmplx(real(t11)*real(t22),aimag(t11)*aimag(t22))
    L13=cmplx(real(t11)*real(t33),aimag(t11)*aimag(t33))
    L22=cmplx(real(t22)*real(t22),aimag(t22)*aimag(t22))
    L23=cmplx(real(t22)*real(t33),aimag(t22)*aimag(t33))
    L33=cmplx(real(t33)*real(t33),aimag(t33)*aimag(t33))
 
    call rfftwnd_f77_one_real_to_complex(r2c3d,L11,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,L12,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,L13,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,L22,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,L23,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,L33,ignore_me)
    
    L11=g_test*L11*const
    L12=g_test*L12*const
    L13=g_test*L13*const
    L22=g_test*L22*const
    L23=g_test*L23*const
    L33=g_test*L33*const
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,L11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,L22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,L33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,L12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,L13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,L23,ignore_me)
 
    L11=L11-cmplx(real(t12)*real(t12),aimag(t12)*aimag(t12))
    L12=L12-cmplx(real(t12)*real(t13),aimag(t12)*aimag(t13))
    L13=L13-cmplx(real(t12)*real(t23),aimag(t12)*aimag(t23))
    L22=L22-cmplx(real(t13)*real(t13),aimag(t13)*aimag(t13))
    L23=L23-cmplx(real(t13)*real(t23),aimag(t13)*aimag(t23))
    L33=L33-cmplx(real(t23)*real(t23),aimag(t23)*aimag(t23))
    ! Lij is found

    ! m_ij
    t11=eye*kx*ux
    t22=eye*ky*uy
    t33=eye*kz*uz
    t12=.5_SP*eye*(ky*ux + kx*uy)
    t13=.5_SP*eye*(kz*ux + kx*uz)
    t23=.5_SP*eye*(ky*uz + kz*uy)
    ! sij in complex

    m11=g_test*t11
    m12=g_test*t12
    m13=g_test*t13
    m22=g_test*t22
    m23=g_test*t23
    m33=g_test*t33
    ! test filtered sij in complex

    call rfftwnd_f77_one_complex_to_real(c2r3d,t11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,t23,ignore_me)
 
    s=cmplx(sqrt(2._sp*(real(t11)*real(t11)+real(t22)*real(t22)+real(t33)*real(t33)                     &
                        +2._sp*(real(t12)*real(t12)+real(t13)*real(t13)+real(t23)*real(t23))))          &
            , &
            sqrt(2._sp*(aimag(t11)*aimag(t11)+aimag(t22)*aimag(t22)+aimag(t33)*aimag(t33)               &
                       +2._sp*(aimag(t12)*aimag(t12)+aimag(t13)*aimag(t13)+aimag(t23)*aimag(t23))))     &
           )

    s11=2._sp*delta2*cmplx(real(s)*real(t11),aimag(s)*aimag(t11))*const
    s12=2._sp*delta2*cmplx(real(s)*real(t12),aimag(s)*aimag(t12))*const
    s13=2._sp*delta2*cmplx(real(s)*real(t13),aimag(s)*aimag(t13))*const
    s22=2._sp*delta2*cmplx(real(s)*real(t22),aimag(s)*aimag(t22))*const
    s23=2._sp*delta2*cmplx(real(s)*real(t23),aimag(s)*aimag(t23))*const
    s33=2._sp*delta2*cmplx(real(s)*real(t33),aimag(s)*aimag(t33))*const
    
    call rfftwnd_f77_one_real_to_complex(r2c3d,s11,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,s12,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,s13,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,s22,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,s23,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,s33,ignore_me)
    ! 2Delta^2|s|sij in complex, will be used to calculate tij in complex and the first part of Mij
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,m11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,m22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,m33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,m12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,m13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,m23,ignore_me)
    ! test-filtered sij in real
                       
    s=cmplx(sqrt(2._sp*(real(m11)*real(m11)+real(m22)*real(m22)+real(m33)*real(m33)                     &
                        +2._sp*(real(m12)*real(m12)+real(m13)*real(m13)+real(m23)*real(m23))))          &
           , &
            sqrt(2._sp*(aimag(m11)*aimag(m11)+aimag(m22)*aimag(m22)+aimag(m33)*aimag(m33)               &
                        +2._sp*(aimag(m12)*aimag(m12)+aimag(m13)*aimag(m13)+aimag(m23)*aimag(m23))))    &
           )

    t11=-2._sp*delta_test2*cmplx(real(s)*real(m11),aimag(s)*aimag(m11))
    t12=-2._sp*delta_test2*cmplx(real(s)*real(m12),aimag(s)*aimag(m12))
    t13=-2._sp*delta_test2*cmplx(real(s)*real(m13),aimag(s)*aimag(m13))
    t22=-2._sp*delta_test2*cmplx(real(s)*real(m22),aimag(s)*aimag(m22))
    t23=-2._sp*delta_test2*cmplx(real(s)*real(m23),aimag(s)*aimag(m23))
    t33=-2._sp*delta_test2*cmplx(real(s)*real(m33),aimag(s)*aimag(m33))
    ! 2nd part of Mij
 
    m11=g_test*s11
    m12=g_test*s12
    m13=g_test*s13
    m22=g_test*s22
    m23=g_test*s23
    m33=g_test*s33
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,m11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,m22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,m33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,m12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,m13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,m23,ignore_me)
    ! first part of Mij in real

    m11=m11+t11
    m12=m12+t12
    m13=m13+t13
    m22=m22+t22
    m23=m23+t23
    m33=m33+t33
    
    t11=cmplx(real(m11)*real(m11)+real(m22)*real(m22)+real(m33)*real(m33)+        &
              2._sp*(real(m12)*real(m12)+real(m13)*real(m13)+real(m23)*real(m23)) &
             , &
              aimag(m11)*aimag(m11)+aimag(m22)*aimag(m22)+aimag(m33)*aimag(m33)+        &
              2._sp*(aimag(m12)*aimag(m12)+aimag(m13)*aimag(m13)+aimag(m23)*aimag(m23)) &
             )
    mm=sum(real(t11(1:lx,:,:)))+sum(aimag(t11(1:lx,:,:)))

    t11=cmplx(real(L11)*real(m11)+real(L22)*real(m22)+real(L33)*real(m33)+        &
              2._sp*(real(L12)*real(m12)+real(L13)*real(m13)+real(L23)*real(m23)) &
             , &
              aimag(L11)*aimag(m11)+aimag(L22)*aimag(m22)+aimag(L33)*aimag(m33)+        &
              2._sp*(aimag(L12)*aimag(m12)+aimag(L13)*aimag(m13)+aimag(L23)*aimag(m23)) &
             )
    lm=sum(real(t11(1:lx,:,:)))+sum(aimag(t11(1:lx,:,:)))

    cs2=lm/mm

    t11=-cs2*s11
    t12=-cs2*s12
    t13=-cs2*s13
    t22=-cs2*s22
    t23=-cs2*s23
    t33=-cs2*s33
    ! tij in complex

  else

    s11=eye*kx*ux
    s22=eye*ky*uy
    s33=eye*kz*uz
    s12=.5_sp*eye*(ky*ux + kx*uy)
    s13=.5_sp*eye*(kz*ux + kx*uz)
    s23=.5_sp*eye*(ky*uz + kz*uy)


    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    
    s=cmplx(sqrt(2._sp*(real(s11)*real(s11)+real(s22)*real(s22)+real(s33)*real(s33)                     &
                        +2._sp*(real(s12)*real(s12)+real(s13)*real(s13)+real(s23)*real(s23))))          &
            , &
            sqrt(2._sp*(aimag(s11)*aimag(s11)+aimag(s22)*aimag(s22)+aimag(s33)*aimag(s33)               &
                       +2._sp*(aimag(s12)*aimag(s12)+aimag(s13)*aimag(s13)+aimag(s23)*aimag(s23))))     &
           )

    t11=-cs2*2._sp*delta2*cmplx(real(s)*real(s11),aimag(s)*aimag(s11))*const
    t12=-cs2*2._sp*delta2*cmplx(real(s)*real(s12),aimag(s)*aimag(s12))*const
    t13=-cs2*2._sp*delta2*cmplx(real(s)*real(s13),aimag(s)*aimag(s13))*const
    t22=-cs2*2._sp*delta2*cmplx(real(s)*real(s22),aimag(s)*aimag(s22))*const
    t23=-cs2*2._sp*delta2*cmplx(real(s)*real(s23),aimag(s)*aimag(s23))*const
    t33=-cs2*2._sp*delta2*cmplx(real(s)*real(s33),aimag(s)*aimag(s33))*const

    call rfftwnd_f77_one_real_to_complex(r2c3d,t11,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t12,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t13,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t22,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t23,ignore_me)
    call rfftwnd_f77_one_real_to_complex(r2c3d,t33,ignore_me)

  end if
    
end subroutine dynsmag      
