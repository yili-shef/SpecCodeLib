! Calculate the Fourier transform of the SGS stresses a priori from DNS data. The
! returned variables tauij are the real SGS stresses: diagonal elements
subroutine rsgstauii (vi,vir2c,vic2r,tij,g,nx,ny,nz) 
  use mconstant
  implicit none
 
  integer,    intent(in) :: nx,ny,nz
  integer(8), intent(in) :: vic2r,vir2c
  complex(sp), dimension(nx/2+1,ny,nz)  :: vi
  complex(sp), dimension(nx/2+1,ny,nz), intent(out) :: tij
  real(sp),    dimension(nx/2+1,ny,nz), intent(in)  :: g

  complex(sp), dimension(nx/2+1,ny,nz) :: ui
  real(sp) :: const

  const=1._sp/real(nx*ny*nz,sp)

  ui = vi

  call dfftw_execute(vic2r)

  vi=cmplx(real(vi, kind = sp)*real(vi, kind = sp), aimag(vi)*aimag(vi))
  call dfftw_execute(vir2c)
  vi=g*const*vi
  call dfftw_execute(vic2r)
  tij = vi


  vi = g*ui

  call dfftw_execute(vic2r)

  tij=tij-cmplx(real(vi, kind = sp)*real(vi, kind = sp),aimag(vi)*aimag(vi))

end subroutine rsgstauii
