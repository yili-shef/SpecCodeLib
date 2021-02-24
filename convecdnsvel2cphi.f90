subroutine convecdnsvel2cphi(va,vb,phi,kx,ky,kz,wa,wb,wphi,lx1,ly,lz,const)
  use mconstant
  use mfftwplan3d
  implicit none

  integer,  intent(in) :: lx1,ly,lz
  real(sp), intent(in) :: const
  complex(sp), dimension(lx1,ly,lz), intent(in)    :: va, vb, phi
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: wa, wb, wphi
  real(sp),    dimension(lx1), intent(in) :: kx
  real(sp),    dimension(ly),  intent(in) :: ky
  real(sp),    dimension(lz),  intent(in) :: kz

  complex(sp), dimension(lx1,ly,lz) :: vxt,vyt,vzt,wz
  real(sp) :: ignore_me
  integer :: ii,jj,kk

  ! convection term in NS equation
  call crayatoxyzfull(va, vb, kx, ky, kz, vxt, vyt, vzt, lx1, ly, lz)
  do kk = 1, lz
  do jj = 1, ly
  do ii = 1, lx1
    wa(ii,jj,kk) = eye * ( ky(jj) * vzt(ii,jj,kk) - kz(kk) * vyt(ii,jj,kk) )
    wb(ii,jj,kk) = eye * ( kz(kk) * vxt(ii,jj,kk) - kx(ii) * vzt(ii,jj,kk) )
    wz(ii,jj,kk) = eye * ( kx(ii) * vyt(ii,jj,kk) - ky(jj) * vxt(ii,jj,kk) )
  end do
  end do
  end do

  call rfftwnd_f77_one_complex_to_real(c2r3d,vxt,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,vyt,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,vzt,ignore_me)

  call rfftwnd_f77_one_complex_to_real(c2r3d,wa,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wb,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

  ! Lamb vector
  wphi=cmplx(  real (vyt)*real (wz)-real (vzt)*real (wb) &
           ,  aimag(vyt)*aimag(wz)-aimag(vzt)*aimag(wb) &
           )
  wz=cmplx(   real (vzt)*real (wa)-real (vxt)*real (wz) &
           ,  aimag(vzt)*aimag(wa)-aimag(vxt)*aimag(wz) &
           )
  wa=cmplx(   real (vxt)*real (wb)-real (vyt)*real (wa) &
           ,  aimag(vxt)*aimag(wb)-aimag(vyt)*aimag(wa) &
           )

  ! convection of passive scalar
  do ii = 1, lx1
    wb(ii,:,:) = eye * kx(ii) * phi(ii,:,:)
  end do
  call rfftwnd_f77_one_complex_to_real(c2r3d,wb,ignore_me)
  wb = cmplx( real(vxt) * real(wb), aimag(vxt) * aimag(wb) )  

  vxt = wphi * const ! wphi is released from temporary storage. vxt now store scaled Lamb vector x-comp

  do jj = 1, ly
    wphi(:,jj,:) = eye * ky(jj) * phi(:,jj,:)
  end do
  call rfftwnd_f77_one_complex_to_real(c2r3d,wphi,ignore_me)
  wphi = wb + cmplx( real(vyt) * real(wphi), aimag(vyt) * aimag(wphi) )  

  do kk = 1, lz
    wb(:,:,kk) = eye * kz(kk) * phi(:,:,kk)
  end do
  call rfftwnd_f77_one_complex_to_real(c2r3d,wb,ignore_me)
  wphi = wphi + cmplx( real(vzt) * real(wb), aimag(vzt) * aimag(wb) )  

  wphi=wphi*const
  call rfftwnd_f77_one_real_to_complex(r2c3d,wphi,ignore_me)


  vyt = wz*const
  vzt = wa*const
  call rfftwnd_f77_one_real_to_complex(r2c3d,vxt,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,vyt,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,vzt,ignore_me)
  call xyztocrayafull (vxt, vyt, vzt, wa, wb, kx, ky, kz, lx1, ly, lz)

end subroutine convecdnsvel2cphi
