! Calculate the Fourier transform of the SGS stresses a priori from DNS data. The
! returned variables tauij are the real SGS stresses.

! Fourier in real out
subroutine rsgstauij (vi,vj,vir2c,vic2r,vjr2c,vjc2r,tij,g,nx,ny,nz) 
  use mconstant
  implicit none
 
  integer,    intent(in) :: nx,ny,nz
  integer(8), intent(in) :: vic2r,vir2c,vjc2r,vjr2c
  complex(sp), dimension(nx/2+1,ny,nz)  :: vi, vj
  complex(sp), dimension(nx/2+1,ny,nz), intent(out) :: tij
  real(sp),    dimension(nx/2+1,ny,nz), intent(in)  :: g

  complex(sp), dimension(nx/2+1,ny,nz) :: ui, uj
  real(sp) :: const

  const=1._sp/real(nx*ny*nz,sp)

  ui = vi; uj = vj




  call dfftw_execute(vic2r)

!write(*,*) 'check point 6'

  call dfftw_execute(vjc2r)


!write(*,*) 'check point 7'

  vi=cmplx(real(vi, kind = sp)*real(vj, kind = sp), aimag(vi)*aimag(vj))
  call dfftw_execute(vir2c)
  
!write(*,*) 'check point 8'

  vi=g*const*vi
  call dfftw_execute(vic2r)
  tij = vi


!write(*,*) 'check point 9'

  vi = g*ui
  vj = g*uj

  call dfftw_execute(vic2r)
  
!write(*,*) 'check point 10'

  call dfftw_execute(vjc2r)

  tij=tij-cmplx(real(vi, kind=sp)*real(vj, kind=sp),aimag(vi)*aimag(vj))

end subroutine rsgstauij
