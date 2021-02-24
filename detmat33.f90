program detmat
    use mconstant
    use mfftwplan3d
    implicit none

    complex(sp), allocatable, dimension(:,:,:) :: b11, b12, b13
    complex(sp), allocatable, dimension(:,:,:) :: b21, b22, b23
    complex(sp), allocatable, dimension(:,:,:) :: b31, b32, b33

    real(sp), dimension(3,3) :: bb
    real(sp) :: det, ignore_me, detmax, detmin
    integer :: ll, nx, ny, nz, lx, ly, lz, lx1, ifile, ii, jj, kk
    character(80) :: str


    ll = iargc()
    if (ll .ne. 2) stop 'usage: detmat.x nx ifile'

    call getarg(1, str)
    read(str, '(I20)') nx
    call getarg(2, str)
    read(str, '(I20)') ifile

    lx = nx/2; lx1 = lx + 1; ly = nx; lz = nx
    ny = nx; nz = nx

    call fftwplan3de(nx,ny,nz)
    write(*,*) 'after fftwplan3d'

    allocate(b11(lx1,ly,lz), b12(lx1,ly,lz), b13(lx1,ly,lz))
    allocate(b21(lx1,ly,lz), b22(lx1,ly,lz), b23(lx1,ly,lz))
    allocate(b31(lx1,ly,lz), b32(lx1,ly,lz), b33(lx1,ly,lz))

    open(20, file = './out/b11'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(20) b11
    close(20)
    open(20, file = './out/b12'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(20) b12
    close(20)
    open(20, file = './out/b13'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(20) b13
    close(20)
    open(20, file = './out/b21'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(20) b21
    close(20)
    open(20, file = './out/b22'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(20) b22
    close(20)
    open(20, file = './out/b23'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(20) b23
    close(20)
    open(20, file = './out/b31'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(20) b31
    close(20)
    open(20, file = './out/b32'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(20) b32
    close(20)
    open(20, file = './out/b33'//str(1:len_trim(str))//'.dat', form = 'unformatted')
      read(20) b33
    close(20)

!    call rfftwnd_f77_one_complex_to_real(c2r3d,b11,ignore_me)
!    call rfftwnd_f77_one_complex_to_real(c2r3d,b12,ignore_me)
!    call rfftwnd_f77_one_complex_to_real(c2r3d,b13,ignore_me)
!    call rfftwnd_f77_one_complex_to_real(c2r3d,b21,ignore_me)
!    call rfftwnd_f77_one_complex_to_real(c2r3d,b22,ignore_me)
!    call rfftwnd_f77_one_complex_to_real(c2r3d,b23,ignore_me)
!    call rfftwnd_f77_one_complex_to_real(c2r3d,b31,ignore_me)
!    call rfftwnd_f77_one_complex_to_real(c2r3d,b32,ignore_me)
!    call rfftwnd_f77_one_complex_to_real(c2r3d,b33,ignore_me)


    detmax = -20.
    detmin = 20.
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

        if ( mod(ii,2) .eq. 1 )  then
            ll = (ii + 1) / 2
  
            bb(1,1) = real(b11(ll,jj,kk))
            bb(1,2) = real(b12(ll,jj,kk))
            bb(1,3) = real(b13(ll,jj,kk))
            bb(2,1) = real(b21(ll,jj,kk))
            bb(2,2) = real(b22(ll,jj,kk))
            bb(2,3) = real(b23(ll,jj,kk))
            bb(3,1) = real(b31(ll,jj,kk))
            bb(3,2) = real(b32(ll,jj,kk))
            bb(3,3) = real(b33(ll,jj,kk))
        else
            ll = ii / 2

            bb(1,1) = aimag(b11(ll,jj,kk))
            bb(1,2) = aimag(b12(ll,jj,kk))
            bb(1,3) = aimag(b13(ll,jj,kk))
            bb(2,1) = aimag(b21(ll,jj,kk))
            bb(2,2) = aimag(b22(ll,jj,kk))
            bb(2,3) = aimag(b23(ll,jj,kk))
            bb(3,1) = aimag(b31(ll,jj,kk))
            bb(3,2) = aimag(b32(ll,jj,kk))
            bb(3,3) = aimag(b33(ll,jj,kk))
        end if

        det = bb(1,1) * ( bb(2,2) * bb(3,3) - bb(2,3) * bb(3,2) ) &
             -bb(1,2) * ( bb(2,1) * bb(3,3) - bb(2,3) * bb(3,1) ) &
             +bb(1,3) * ( bb(2,1) * bb(3,2) - bb(2,2) * bb(3,1) )

        detmax = max(detmax, det)
        detmin = min(detmin, det)
 
    end do
    end do
    end do
    write(*,*) 'max det = ', detmax
    write(*,*) 'min det = ', detmin

    call destroyplan3d

    deallocate(b11,b12,b13,b21,b22,b23,b31,b32,b33)

end program detmat
