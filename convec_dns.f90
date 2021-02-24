subroutine convec_dns(vx,vy,vz,wx,wy,wz,lx1,ly,lz,nx,ny,nz)
  use mconstant
  use mfftwplan3d
  implicit none

  integer, intent(in) :: lx1,ly,lz,nx,ny,nz
  complex(sp),dimension(lx1,ly,lz),intent(in)    :: vx,vy,vz
  complex(sp),dimension(lx1,ly,lz),intent(inout) :: wx,wy,wz

  complex(sp), dimension(lx1,ly,lz) :: vxt,vyt,vzt,tmp
  real(sp) :: const, ignore_me

  integer, parameter :: numth = 4

  const=1._sp/real(nx*ny*nz,sp)
  !$omp parallel workshare num_threads(numth)
  vxt=vx
  vyt=vy
  vzt=vz
  !$omp end parallel workshare

  call rfftwnd_f77_one_complex_to_real(c2r3d,vxt,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,vyt,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,vzt,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)

  !$omp parallel workshare num_threads(numth)
  tmp=cmplx(  real (vyt)*real (wz)-real (vzt)*real (wy)               &
           ,  aimag(vyt)*aimag(wz)-aimag(vzt)*aimag(wy)                       &
           )
  wz=cmplx(   real (vzt)*real (wx)-real (vxt)*real (wz)               &
           ,  aimag(vzt)*aimag(wx)-aimag(vxt)*aimag(wz)                       &
           )
  wx=cmplx(   real (vxt)*real (wy)-real (vyt)*real (wx)               &
           ,  aimag(vxt)*aimag(wy)-aimag(vyt)*aimag(wx)                       &
           )

  vxt=tmp*const
  vyt=wz*const
  vzt=wx*const
  !$omp end parallel workshare 

  call rfftwnd_f77_one_real_to_complex(r2c3d,vxt,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,vyt,ignore_me)
  call rfftwnd_f77_one_real_to_complex(r2c3d,vzt,ignore_me)

  !$omp parallel workshare num_threads(numth)
  wx=vxt
  wy=vyt
  wz=vzt
  !$omp end parallel workshare 

end subroutine convec_dns
