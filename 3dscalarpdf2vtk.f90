program staligntovtk
  use mconstant
  implicit none

  integer :: npntx, npnty, npntz
  real(sp), allocatable, dimension(:) :: xco, yco, zco
  real(sp), allocatable, dimension(:, :, :) :: pdf

  integer :: ii, jj, kk
  character(80) :: fnm, scalarname
  integer :: vtkopenfile, vtk_rg_header, vtkscalar, vtkclosefile

  kk = iargc()
  if ( kk .ne. 5 ) then
      write(*,*) 
      write(*,*) ' >>>> Wrong number of arguments. Usage: '
      write(*,*)
      write(*,*) ' ./3dscalarpdf2vtk.x npntx npnty npntz fnm scalarname'
      write(*,*) 
      write(*,*) ' Stopped'
      stop
  end if

  call getarg(1, fnm)
  read(fnm, '(I6)') npntx
  call getarg(2, fnm)
  read(fnm, '(I6)') npnty
  call getarg(3, fnm)
  read(fnm, '(I6)') npntz

  call getarg(4, fnm)
  fnm = adjustl(fnm)

  call getarg(5, scalarname)
  scalarname = adjustl(scalarname)

  allocate( xco(npntx), yco(npnty), zco(npntz) )
  allocate( pdf(npntx, npnty, npntz) )

  open(17, file=fnm(1:len_trim(fnm))//'.dat')
  read(17,*)
  do kk = 1, npntz
  do jj = 1, npnty
  do ii = 1, npntx
    read(17,*) xco(ii), yco(jj), zco(kk), pdf(ii,jj,kk)
  end do
  end do
  end do
  close(17)

  ! Change suffix and appending 0 at the end for C output.
  fnm = fnm(1:len_trim(fnm))//".vtk0"
  scalarname = scalarname( 1:len_trim(scalarname) )//"0"
  write(*,*) fnm(1:len_trim(fnm))
  write(*,*) scalarname(1:len_trim(scalarname))

  kk = vtkopenfile( fnm( 1:len_trim(fnm) ) )
  kk = vtk_rg_header(xco, yco, zco, sp, npntx, npnty, npntz)
  !kk = vtkscalar(pdf, sp, npntx, npnty, npntz, len_trim(scalarname), &
  !               scalarname(1 : len_trim(scalarname)) )
  kk = vtkscalar(pdf, sp, npntx, npnty, npntz, scalarname(1 : len_trim(scalarname)) )
  kk = vtkclosefile()

  deallocate(xco, yco, zco, pdf)

  write(*,*) 'Finished'

end program staligntovtk      
