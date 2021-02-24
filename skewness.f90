subroutine skewness (ux,kx,s,nx,ny,nz)
  use mconstant
  use mfftwplan3d
  implicit none

  integer, intent(in) :: nx,ny,nz
  complex(sp), dimension(nx/2+1,ny,nz), intent(in) :: ux
  real(sp),    dimension(nx/2+1),     intent(in) :: kx
  real(sp), intent(out) :: s

  complex(sp), dimension(nx/2+1,ny,nz) :: duxdx
  real(sp) :: ignore_me,const,tmp

  integer :: jj,kk

  const=1._sp/real(nx*ny*nz,sp)
  do kk=1,nz
  do jj=1,ny
    duxdx(:,jj,kk)=eye*kx*ux(:,jj,kk)
  end do
  end do

  call rfftwnd_f77_one_complex_to_real(c2r3d,duxdx,ignore_me)
  tmp=(   sum(real(duxdx(1:nx/2,:,:)))  +  sum(aimag(duxdx(1:nx/2,:,:)))   )*const
  
  s=(                            &
     sum(                        &
         (real(duxdx(1:nx/2,:,:)))**3         &
        )                        &
     +                           &
     sum(                        &
         (aimag(duxdx(1:nx/2,:,:)))**3        &
        )                        &
    )*const                      &
    /                            &
    (                            &
      (                          &
        sum(                     &
            (real(duxdx(1:nx/2,:,:)))**2      &
           )                     &
        +                        &
        sum(                     &
            (aimag(duxdx(1:nx/2,:,:)))**2     &
           )                     &
      )*const                    &
    )**1.5_sp

end subroutine skewness      
