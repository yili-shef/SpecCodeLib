program rossbyno
    use mconstant
    use mwavenumber
    implicit none

    complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz
    real(sp),    allocatable, dimension(:)     :: kx, ky, kz, ek

    integer :: ii, jj, kk, ll, nx, ny, nz, lx, lx1, ly, lz
    real(sp) :: kxii, kyjj, kzkk, etot, intlen, uprime, omprime, omintlen, eee, enstot
    real(sp) :: kkk
    character(80) :: str


    nx = iargc()
    if (nx .ne. 2) stop 'Usage: ./rossbyno.x nx nfile'

    call getarg(1, str) 
    read(str, '(I20)') nx
    call getarg(2, str) 
    str = adjustl(str)

    ny = nx; nz = nx; lx = nx / 2; lx1 = lx + 1; ly = ny; lz = nz

    allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
    allocate( kx(lx1), ky(ly), kz(lz), ek(lx) )

    call wavenumber(kx,ky,kz,lx1,ly,lz) 

    open(20, file = './out/ux'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(20) ux
    close(20)
    open(20, file = './out/uy'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(20) uy
    close(20)
    open(20, file = './out/uz'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(20) uz
    close(20)

    ux = ux * conjg(ux) + uy * conjg(uy) + uz * conjg(uz)
    ux(1,:,:) = .5 * ux(1,:,:)

    etot = 0.
    eee = 0.
    enstot = 0.
    do kk = 1, lz
    kzkk = kz(kk)
    do jj = 1, ly
    kyjj = ky(jj)
    do ii = 1, lx1
    kxii = kx(ii)

      kkk = sqrt(kzkk * kzkk + kyjj * kyjj + kxii * kxii)

      ll = floor( kkk + .5 )
      if ( ll .ge. 1 .and. ll .le. lx ) then
          ! ek(ll) = ek(ll) + real( ux(ii,jj,kk) )
          etot = etot + real( ux(ii,jj,kk) )
          eee = eee + real( ux(ii,jj,kk) ) / kkk
          enstot = enstot + 2 * kkk * kkk * real( ux(ii,jj,kk) )
      end if

    end do
    end do
    end do

    !etot = sum(ek) ! < ui ui > / 2

    !enstot = 0.
    !eee = 0.
    !do ll = 1, lx
    !    enstot = enstot + 2 * ll * ll * ek(ll)  ! < omi omi > 
    !    eee = eee + ek(ll) / ll
    !end do

    intlen = eee / etot
    omprime = sqrt( enstot / 3 )

    uprime = sqrt( 2 * etot / 3 )
    omintlen = uprime / intlen

    write(*,*) 'uprime   intlen   omintlen   omprime'
    write(*,*)  uprime,  intlen,  omintlen,  omprime


    deallocate( ux, uy, uz, kx, ky, kz )

end program rossbyno      
