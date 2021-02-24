program xtrctu
  implicit none

  integer, parameter :: nx = 128, ny = 128, nz = 33
  integer, parameter :: lh = nx/2 +1, ld = lh * 2

  real(8), parameter :: L_x = 32.d0
  real(8), parameter :: dx = L_x/nx              
  real(8), parameter :: dy = dx, dz = dx
  real(8), parameter :: L_y = dy * ny, L_z = dz * nz

  real(8), dimension(ld,ny,nz) :: u, trsh
  real(8), dimension(nx,ny,nz) :: ur
  real(8), dimension(nx) :: xgrid
  real(8), dimension(ny) :: ygrid
  real(8), dimension(nz) :: zgrid


  integer :: ii, jj, kk
  character(80) :: fpath
  integer :: vtkopenfile, vtk_rg_header, vtkscalar, vtkclosefile

  open(21, file = 'vel.out', form = 'unformatted')
  read(21) u, trsh, trsh, trsh, trsh, trsh, zgrid
  close(21)
  ur = u(1:nx,:,:)
  xgrid = dx * (/(ii, ii = 1, nx)/)
  ygrid = dy * (/(jj, jj = 1, ny)/)
  zgrid = dz * (/(kk, kk = 1, nz)/)

  ! Change suffix and appending 0 at the end for C output.
  fpath = "u-ABLLESJHU.vtk0"
  kk = vtkopenfile (fpath( 1:len_trim(fpath) )) 
  kk = vtk_rg_header (xgrid, ygrid, zgrid, 8, nx, ny, nz)
  fpath = "u0"
  kk = vtkscalar (ur, 8, nx, ny, nz, fpath(1 : len_trim(fpath)) )
  kk = vtkclosefile()

end program xtrctu      
