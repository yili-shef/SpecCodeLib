subroutine convecdns2c(va,vb,kx,ky,kz,wa,wb,lx1,ly,lz,const)
  use mconstant
  use mfftwplan3d
  implicit none

  integer,  intent(in) :: lx1,ly,lz
  real(sp), intent(in) :: const
  complex(sp), dimension(lx1,ly,lz), intent(in)    :: va, vb
  complex(sp), dimension(lx1,ly,lz), intent(inout) :: wa, wb
  real(sp),    dimension(lx1), intent(in) :: kx
  real(sp),    dimension(ly),  intent(in) :: ky
  real(sp),    dimension(lz),  intent(in) :: kz

  complex(sp), dimension(lx1,ly,lz) :: vxt,vyt,vzt,wz,tmp
  real(sp) :: ignore_me
  integer :: ii,jj,kk

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

  tmp=cmplx(  real (vyt)*real (wz)-real (vzt)*real (wb) &
           ,  aimag(vyt)*aimag(wz)-aimag(vzt)*aimag(wb) &
           )
  wz=cmplx(   real (vzt)*real (wa)-real (vxt)*real (wz) &
           ,  aimag(vzt)*aimag(wa)-aimag(vxt)*aimag(wz) &
           )
  wa=cmplx(   real (vxt)*real (wb)-real (vyt)*real (wa) &
           ,  aimag(vxt)*aimag(wb)-aimag(vyt)*aimag(wa) &
           )

  vxt = tmp*const
  vyt = wz*const
  vzt = wa*const
  call rfftwnd_f77_one_real_to_complex(r2c3d,vxt,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,vyt,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,vzt,ignore_me)
  call xyztocrayafull (vxt, vyt, vzt, wa, wb, kx, ky, kz, lx1, ly, lz)

end subroutine convecdns2c
