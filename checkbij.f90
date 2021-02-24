program checkbij
  use mconstant

  integer, parameter :: nx = 128
  integer, parameter :: ny = nx, nz = nx, lx = nx/2, lx1 = lx+1, ly = ny, lz = nz

  complex(sp), dimension(lx1,ly,lz) :: b11, b12, b13, b21, b22, b23, b31, b32, b33

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  ! ----------------------------------------------------------------------------------

  real(sp) :: cc(evsize,evsize) 
  integer :: ii, jj, kk, ll
  character(80) :: fnsuffix

  if ( iarg() .ne. 1) stop 'Usage: checkbij.x fnsuffix'

  call getarg(1,fnsuffix)
  fnsuffix = adjustl(fnsuffix)

  call readbij

  call caleig
  !call caldet


contains 

  subroutine readbij
    
    open(30, file = './out/b11'//trim(fnsuffix)//'.dat', form = 'unformatted')
        read(30) b11
    close(30)
    open(30, file = './out/b12'//trim(fnsuffix)//'.dat', form = 'unformatted')
        read(30) b12
    close(30)
    open(30, file = './out/b13'//trim(fnsuffix)//'.dat', form = 'unformatted')
        read(30) b13
    close(30)
    open(30, file = './out/b21'//trim(fnsuffix)//'.dat', form = 'unformatted')
        read(30) b21
    close(30)
    open(30, file = './out/b22'//trim(fnsuffix)//'.dat', form = 'unformatted')
        read(30) b22
    close(30)
    open(30, file = './out/b23'//trim(fnsuffix)//'.dat', form = 'unformatted')
        read(30) b23
    close(30)
    open(30, file = './out/b31'//trim(fnsuffix)//'.dat', form = 'unformatted')
        read(30) b31
    close(30)
    open(30, file = './out/b32'//trim(fnsuffix)//'.dat', form = 'unformatted')
        read(30) b32
    close(30)
    open(30, file = './out/b33'//trim(fnsuffix)//'.dat', form = 'unformatted')
        read(30) b33
    close(30)
  end subroutine readbij

  subroutine caldet

    real(sp) :: maxdet, det

    maxdet = 0._sp

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

    if ( mod(ii, 2) .eq. 1 ) then
        ll = (ii + 1)/2

        cc(1,1)=real(b11(ll,jj,kk))
        cc(1,2)=real(b12(ll,jj,kk))
        cc(1,3)=real(b13(ll,jj,kk))
        cc(2,1)=real(b21(ll,jj,kk))
        cc(2,2)=real(b22(ll,jj,kk))
        cc(2,3)=real(b23(ll,jj,kk))
        cc(3,1)=real(b31(ll,jj,kk))
        cc(3,2)=real(b32(ll,jj,kk))
        cc(3,3)=real(b33(ll,jj,kk))
    else
        ll = ii / 2

        cc(1,1)=aimag(b11(ll,jj,kk))
        cc(1,2)=aimag(b12(ll,jj,kk))
        cc(1,3)=aimag(b13(ll,jj,kk))
        cc(2,1)=aimag(b21(ll,jj,kk))
        cc(2,2)=aimag(b22(ll,jj,kk))
        cc(2,3)=aimag(b23(ll,jj,kk))
        cc(3,1)=aimag(b31(ll,jj,kk))
        cc(3,2)=aimag(b32(ll,jj,kk))
        cc(3,3)=aimag(b33(ll,jj,kk))
    end if

    det = cc(1,1) * ( cc(2,2) * cc(3,3) - cc(3,2) * cc(2,3) ) + &
          cc(1,2) * ( cc(3,1) * cc(2,3) - cc(2,1) * cc(3,3) ) + &
          cc(1,3) * ( cc(2,1) * cc(3,2) - cc(3,1) * cc(2,2) )

    maxdet = max(abs(det), maxdet)

    end do
    end do
    end do

    write(*,*) 'maxdet = ', maxdet

  end subroutine caldet

  subroutine caleig

    real(dp), dimension(evsize) :: evcc 
    real(dp), dimension(evsize,evsize) :: evtrcc 

    real(dp) :: meanal, meanbe, meanga, varal, varbe, varga

    integer, parameter :: npnt = nx * ny * nz

    meanal = 0._dp; meanbe = 0._dp; meanga = 0._dp
    varal  = 0._dp; varbe  = 0._dp; varga  = 0._dp

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

    if ( mod(ii, 2) .eq. 1 ) then
        ll = (ii + 1)/2

        cc(1,1)=real(b11(ll,jj,kk))
        cc(1,2)=real(b12(ll,jj,kk))
        cc(1,3)=real(b13(ll,jj,kk))
        cc(2,1)=real(b21(ll,jj,kk))
        cc(2,2)=real(b22(ll,jj,kk))
        cc(2,3)=real(b23(ll,jj,kk))
        cc(3,1)=real(b31(ll,jj,kk))
        cc(3,2)=real(b32(ll,jj,kk))
        cc(3,3)=real(b33(ll,jj,kk))
    else
        ll = ii / 2

        cc(1,1)=aimag(b11(ll,jj,kk))
        cc(1,2)=aimag(b12(ll,jj,kk))
        cc(1,3)=aimag(b13(ll,jj,kk))
        cc(2,1)=aimag(b21(ll,jj,kk))
        cc(2,2)=aimag(b22(ll,jj,kk))
        cc(2,3)=aimag(b23(ll,jj,kk))
        cc(3,1)=aimag(b31(ll,jj,kk))
        cc(3,2)=aimag(b32(ll,jj,kk))
        cc(3,3)=aimag(b33(ll,jj,kk))
    end if

    call dsyevr("V", "A", "U", evsize, cc, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                nfound, evcc, evtrcc, evsize, isuppz, work, lwork, iwork, liwork, info)

    meanal = meanal + evcc(3)
    meanbe = meanbe + evcc(2)
    meanga = meanga + evcc(1)
    
    varal = varal + evcc(3) * evcc(3) 
    varbe = varbe + evcc(2) * evcc(2) 
    varga = varga + evcc(1) * evcc(1)
    
    end do
    end do
    end do

    meanal = meanal / npnt
    meanbe = meanbe / npnt
    meanga = meanga / npnt
    varal = varal / npnt - meanal * meanal 
    varbe = varbe / npnt - meanbe * meanbe
    varga = varga / npnt - meanga * meanga

    write(*,*) 'meanal', meanal, 'meanbe', meanbe, 'meanga', meanga
    write(*,*) 'varal', varal, 'varbe', varbe, 'varga', varga

  end subroutine caleig

end program checkbij      
