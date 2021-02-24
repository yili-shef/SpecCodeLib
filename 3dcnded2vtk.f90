program cnded2vtk
  use mconstant
  implicit none

  integer, parameter :: npnt = 40
  real(sp), dimension(npnt) :: xco, yco, zco
  real(sp), dimension(npnt, npnt, npnt) :: pdf, cnded3d, weightedcnded

  integer :: ii, jj, kk
  character(80) :: fpath
  integer :: vtkopenfile, vtk_rg_header, vtkscalar, vtkclosefile

  open(17, file='cndsgsed-align-st-3d-16dx-mtlm128heli.dat')
  read(17,*)
  do kk = 1, npnt
  do jj = 1, npnt
  do ii = 1, npnt
    read(17,'(15E15.5)') xco(ii), yco(jj), zco(kk), pdf(ii,jj,kk), cnded3d(ii,jj,kk), &
                         weightedcnded(ii,jj,kk)
  end do
  end do
  end do
  close(17)

  ! Change suffix and appending 0 at the end for C output.
  fpath = "cndsgsed-align-st-3d-16dx-mtlm128heli.vtk0"
  kk = vtkopenfile (fpath( 1:len_trim(fpath) )) 
  kk = vtk_rg_header (xco, yco, zco, sp, npnt, npnt, npnt)
  fpath = "cndsgsed0"
  kk = vtkscalar (cnded3d, sp, npnt, npnt, npnt, fpath(1 : len_trim(fpath)) )
  fpath = "alignpdf0"
  kk = vtkscalar (pdf, sp, npnt, npnt, npnt, fpath(1 : len_trim(fpath)) )
  fpath = "weightedcnded0"
  kk = vtkscalar (weightedcnded, sp, npnt, npnt, npnt, fpath(1 : len_trim(fpath)) )
  kk = vtkclosefile()

end program cnded2vtk      
