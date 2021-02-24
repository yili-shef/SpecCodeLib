subroutine csgsexpa (vx,vy,vz,t11,t12,t13,t22,t23,t33,kx,ky,kz,lx1,ly,lz,nx,ny,nz,delta)
  use mconstant
  use mfftwplan3d
  implicit none

  integer,  intent(in) :: lx1,ly,lz,nx,ny,nz
  real(sp), intent(in) :: delta
  real(sp),    dimension(lx1,ly,lz), intent(in)  :: kx,ky,kz
  complex(sp), dimension(lx1,ly,lz), intent(in)  :: vx,vy,vz
  complex(sp), dimension(lx1,ly,lz), intent(out) :: t11,t12,t13,t22,t23,t33

  complex(sp), dimension(lx1,ly,lz) :: t21,t31,t32

  real(sp) :: ignore_me, const

  real(sp), parameter :: cs2=0.01
  real(sp), parameter :: ctau=1.

  integer :: ii,jj,kk,ll,mm

  integer, parameter :: three=3
  ! Stuff for sgpadm, matrix exponential routine
  ! integer,  parameter :: ideg=6
  ! integer :: ns, nsmax


  real(sp), dimension(three,three) :: A, RA
  real(sp) :: ANorm, S12, S13,S23, SS


  const = 1./real(nx*ny*nz,sp)

  ! Using tij to store velocity gradient temporarily.
  t11=eye*kx*vx; t12=eye*kx*vy; t13=eye*kx*vz
  t21=eye*ky*vx; t22=eye*ky*vy; t23=eye*ky*vz
  t31=eye*kz*vx; t32=eye*kz*vy; t33=eye*kz*vz

  call rfftwnd_f77_one_complex_to_real(C2R3D,t11,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t12,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t13,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t21,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t22,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t23,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t31,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t32,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t33,ignore_me)

  do ii=1,lx1-1
    do jj=1,ly
      do kk=1,lz
        A(1,1)=real(t11(ii,jj,kk)); A(1,2)=real(t12(ii,jj,kk)); A(1,3)=real(t13(ii,jj,kk))
        A(2,1)=real(t21(ii,jj,kk)); A(2,2)=real(t22(ii,jj,kk)); A(2,3)=real(t23(ii,jj,kk))
        A(3,1)=real(t31(ii,jj,kk)); A(3,2)=real(t32(ii,jj,kk)); A(3,3)=real(t33(ii,jj,kk))
        
        S12=.5*(A(1,2)+A(2,1)); S13=.5*(A(1,3)+A(3,1)); S23=.5*(A(2,3)+A(3,2))
        SS=A(1,1)**2+A(2,2)**2+A(3,3)**2

        ANorm = SS+A(1,2)**2+A(1,3)**2+A(2,1)**2+A(2,3)**2 +A(3,1)**2+A(3,2)**2

        SS=2*(SS+2*(S12*S12+S13*S13+S23*S23))
        SS=cs2*delta*delta*SS

        ANorm = ctau/sqrt(ANorm)
        A=-A*Anorm

        call taylorch(A,ctau,three)
        A=matmul(A,transpose(A))

        RA=SS*A


        A(1,1)=aimag(t11(ii,jj,kk)); A(1,2)=aimag(t12(ii,jj,kk)); A(1,3)=aimag(t13(ii,jj,kk))
        A(2,1)=aimag(t21(ii,jj,kk)); A(2,2)=aimag(t22(ii,jj,kk)); A(2,3)=aimag(t23(ii,jj,kk))
        A(3,1)=aimag(t31(ii,jj,kk)); A(3,2)=aimag(t32(ii,jj,kk)); A(3,3)=aimag(t33(ii,jj,kk))

        S12=.5*(A(1,2)+A(2,1)); S13=.5*(A(1,3)+A(3,1)); S23=.5*(A(2,3)+A(3,2))
        SS=A(1,1)**2+A(2,2)**2+A(3,3)**2

        ANorm = SS+A(1,2)**2+A(1,3)**2+A(2,1)**2+A(2,3)**2 +A(3,1)**2+A(3,2)**2

        SS=2*(SS+2*(S12*S12+S13*S13+S23*S23))
        SS=cs2*delta*delta*SS

        ANorm = ctau/sqrt(ANorm)
        A=-A*Anorm

        call taylorch(A,ctau,three)
        A=matmul(A,transpose(A))
        A=SS*A

        t11(ii,jj,kk)=cmplx(RA(1,1),A(1,1)) * const
        t12(ii,jj,kk)=cmplx(RA(1,2),A(1,2)) * const
        t13(ii,jj,kk)=cmplx(RA(1,3),A(1,3)) * const
        t22(ii,jj,kk)=cmplx(RA(2,2),A(2,2)) * const
        t23(ii,jj,kk)=cmplx(RA(2,3),A(2,3)) * const
        t33(ii,jj,kk)=cmplx(RA(3,3),A(3,3)) * const

        ss = real( -( t11(ii,jj,kk) + t22(ii,jj,kk) + t33(ii,jj,kk) ) / 3. )
        anorm = aimag( -( t11(ii,jj,kk) + t22(ii,jj,kk) + t33(ii,jj,kk) ) / 3. ) 
        t11(ii,jj,kk) = t11(ii,jj,kk) + cmplx(ss, anorm)
        t22(ii,jj,kk) = t22(ii,jj,kk) + cmplx(ss, anorm)
        t33(ii,jj,kk) = t33(ii,jj,kk) + cmplx(ss, anorm)

      end do
    end do
  end do


  call rfftwnd_f77_one_real_to_complex(R2C3D,t11,ignore_me)
  call rfftwnd_f77_one_real_to_complex(R2C3D,t12,ignore_me)
  call rfftwnd_f77_one_real_to_complex(R2C3D,t13,ignore_me)
  call rfftwnd_f77_one_real_to_complex(R2C3D,t22,ignore_me)
  call rfftwnd_f77_one_real_to_complex(R2C3D,t23,ignore_me)
  call rfftwnd_f77_one_real_to_complex(R2C3D,t33,ignore_me)

!  t21=-(t11+t22+t33)/3._sp

!  t11=t11+t21
!  t22=t22+t21
!  t33=t33+t21
end subroutine csgsexpa
