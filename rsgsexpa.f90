subroutine rsgsexpa (vx,vy,vz,t11,t12,t13,t22,t23,t33,kx,ky,kz,lx1,ly,lz,nx,ny,nz,delta)
  use mconstant
  use mfftwplan3d
  implicit none

  integer,  intent(in) :: lx1,ly,lz,nx,ny,nz
  real(sp), intent(in) :: delta
  real(sp),    dimension(lx1), intent(in)  :: kx
  real(sp),    dimension(ly ), intent(in)  :: ky
  real(sp),    dimension(lz ), intent(in)  :: kz
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


  ! Using tij to store velocity gradient temporarily.
  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    t11(ii,jj,kk)=eye*kx(ii)*vx(ii,jj,kk); t21(ii,jj,kk)=eye*kx(ii)*vy(ii,jj,kk); t31(ii,jj,kk)=eye*kx(ii)*vz(ii,jj,kk)
    t12(ii,jj,kk)=eye*ky(jj)*vx(ii,jj,kk); t22(ii,jj,kk)=eye*ky(jj)*vy(ii,jj,kk); t32(ii,jj,kk)=eye*ky(jj)*vz(ii,jj,kk)
    t13(ii,jj,kk)=eye*kz(kk)*vx(ii,jj,kk); t23(ii,jj,kk)=eye*kz(kk)*vy(ii,jj,kk); t33(ii,jj,kk)=eye*kz(kk)*vz(ii,jj,kk)
  end do
  end do
  end do

  call rfftwnd_f77_one_complex_to_real(C2R3D,t11,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t12,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t13,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t21,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t22,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t23,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t31,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t32,ignore_me)
  call rfftwnd_f77_one_complex_to_real(C2R3D,t33,ignore_me)

  do kk=1,nz
    do jj=1,ny
      do ii=1,nx

        if ( mod(ii, 2) .eq. 1) then
                ll = (ii+1) / 2
          A(1,1)=real(t11(ll,jj,kk)); A(1,2)=real(t12(ll,jj,kk)); A(1,3)=real(t13(ll,jj,kk))
          A(2,1)=real(t21(ll,jj,kk)); A(2,2)=real(t22(ll,jj,kk)); A(2,3)=real(t23(ll,jj,kk))
          A(3,1)=real(t31(ll,jj,kk)); A(3,2)=real(t32(ll,jj,kk)); A(3,3)=real(t33(ll,jj,kk))
        else
                ll = ii / 2
          A(1,1)=aimag(t11(ll,jj,kk)); A(1,2)=aimag(t12(ll,jj,kk)); A(1,3)=aimag(t13(ll,jj,kk))
          A(2,1)=aimag(t21(ll,jj,kk)); A(2,2)=aimag(t22(ll,jj,kk)); A(2,3)=aimag(t23(ll,jj,kk))
          A(3,1)=aimag(t31(ll,jj,kk)); A(3,2)=aimag(t32(ll,jj,kk)); A(3,3)=aimag(t33(ll,jj,kk))
        end if

        
        S12=.5*(A(1,2)+A(2,1)); S13=.5*(A(1,3)+A(3,1)); S23=.5*(A(2,3)+A(3,2))
        SS=A(1,1)**2+A(2,2)**2+A(3,3)**2

        ANorm = SS+A(1,2)**2+A(1,3)**2+A(2,1)**2+A(2,3)**2 +A(3,1)**2+A(3,2)**2

        SS=2*(SS+2*(S12*S12+S13*S13+S23*S23))
        SS=cs2*delta*delta*SS

        ANorm = ctau/sqrt(ANorm)

        A = -A*Anorm

        call taylorch(A,ctau,three)
        A=matmul(A,transpose(A))

        if ( mod(ii,2) .eq. 1) then
          RA = SS * A
        else
          A = SS * A

          ll = ii / 2

          t11(ll,jj,kk)=cmplx(RA(1,1),A(1,1))
          t12(ll,jj,kk)=cmplx(RA(1,2),A(1,2))
          t13(ll,jj,kk)=cmplx(RA(1,3),A(1,3))
          t22(ll,jj,kk)=cmplx(RA(2,2),A(2,2))
          t23(ll,jj,kk)=cmplx(RA(2,3),A(2,3))
          t33(ll,jj,kk)=cmplx(RA(3,3),A(3,3))
         
          ss = real( -( t11(ll,jj,kk) + t22(ll,jj,kk) + t33(ll,jj,kk) ) / 3. )
          anorm = aimag( -( t11(ll,jj,kk) + t22(ll,jj,kk) + t33(ll,jj,kk) ) / 3. ) 

          t11(ll,jj,kk) = t11(ll,jj,kk) + cmplx(ss, anorm)
          t22(ll,jj,kk) = t22(ll,jj,kk) + cmplx(ss, anorm)
          t33(ll,jj,kk) = t33(ll,jj,kk) + cmplx(ss, anorm)

        end if

      end do
    end do
  end do

end subroutine rsgsexpa
