program gridxyz
    use mconstant

    integer :: nx, ny, nz, ii, jj, kk, ll

    real(sp), allocatable, dimension(:) :: xp, yp, zp

    character(80) :: prefix, str, str1
    real(sp) :: dx, dy, dz

    if ( iargc() .ne. 3) stop 'Usage: ./gridxyz.x nx nfilestart prefix'

    call getarg(1, str)
    read(str, '(I20)') nx

    call getarg(2, str)
    str = adjustl(str)

    call getarg(3, prefix)
    prefix = adjustl(prefix)


    ny = nx; nz = nx
    allocate( xp(nx*ny*nz), yp(nx*ny*nz), zp(nx*ny*nz) )


    dx = 2 * pi / nx
    dy = 2 * pi / ny
    dz = 2 * pi / nz

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      ll = (kk-1) * nx * ny + (jj-1) * nx + ii

      xp(ll) = (ii-1)*dx
      yp(ll) = (jj-1)*dy
      zp(ll) = (kk-1)*dz

    end do
    end do
    end do

    str1 = './out/'//prefix(1:len_trim(prefix))//str(1:len_trim(str))//'.dat'
    open(16, file = str1(1:len_trim(str1)), form = 'unformatted')
      write(16) xp, yp, zp
    close(16)

    deallocate( xp, yp, zp )

end program gridxyz
