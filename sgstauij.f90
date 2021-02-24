! Calculate the Fourier transform of the SGS stresses a priori from DNS data. The
! returned variables tauij are the Fourier transform of the (real) SGS stresses.
! This is different from the MPI version, which returns the real SGS stresses.
subroutine sgstauij (vi,vj,tij,g,nx,ny,nz) 
  use mconstant
  use mfftwplan3d
  implicit none
 
  integer, intent(in) :: nx,ny,nz
  complex(sp), dimension(nx/2+1,ny,nz), intent(in)  :: vi, vj
  complex(sp), dimension(nx/2+1,ny,nz), intent(out) :: tij
  real(sp),    dimension(nx/2+1,ny,nz), intent(in)  :: g

  complex(sp), dimension(nx/2+1,ny,nz) :: ui, uj
  real(sp) :: const, ignore_me

  const=1._sp/real(nx*ny*nz,sp)
  ui=vi 
  uj=vj

  call rfftwnd_f77_one_complex_to_real(c2r3d,ui,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,uj,ignore_me)

  tij=cmplx(real(ui)*real(uj),aimag(ui)*aimag(uj))
  call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
  tij=g*const*tij
  call rfftwnd_f77_one_complex_to_real(c2r3d,tij,ignore_me)


  ui=g*vi
  uj=g*vj

  call rfftwnd_f77_one_complex_to_real(c2r3d,ui,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,uj,ignore_me)

  tij=tij-cmplx(real(ui)*real(uj),aimag(ui)*aimag(uj))

  call rfftwnd_f77_one_real_to_complex(r2c3d,tij,ignore_me)
  tij=const*tij

end subroutine sgstauij

