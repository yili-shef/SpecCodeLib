program findhwpoints
    use mconstant
    use mfftwplan3d
    implicit none
    
    real(sp), parameter :: rmsw = 6.3079 ! calculated from hr256.list
    integer, parameter :: nx = 256, ndel = 0, nrms = 3
 
    integer :: ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll, ip, nini
    real(sp) :: ignore_me, delta_c
 
    complex(sp) :: a12, a13, a21, a23, a31, a32
 
    complex(sp), allocatable, dimension(:,:,:) :: wx, wy, wz
    real(sp),    allocatable, dimension(:,:,:) :: g
    real(sp),    allocatable, dimension(:)     :: kx, ky, kz, xp, yp, zp

    character(80) :: str, strfn
    
    ny=nx; nz=nx
    lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz
    delta_c=ndel*2*pi/nx
 
    allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
    allocate(g (lx1,ly,lz),kx(lx1), ky(ly), kz(lz))
    allocate(xp(nx*ny*nz), yp(nx*ny*nz), zp(nx*ny*nz))
 
    call fftwplan3de(nx,ny,nz)
    write(*,*) 'after fftwplan'
 
    call wavenumber(kx,ky,kz,g,lx1,ly,lz)
    write(*,*) 'after wavenumber'
 
    ! Gaussian filter
    g=exp(-g*delta_c**2/24.)
 
    nini = 12
    open(35, file = 'findho.list')
 
        do while ( .not. eof(35) )
            read(35,*) str
 
            open(15,file='./out/ux'//str(1:len_trim(str)),form='unformatted')
              read(15) wx
            close(15)
            open(15,file='./out/uy'//str(1:len_trim(str)),form='unformatted')
              read(15) wy
            close(15)
            open(15,file='./out/uz'//str(1:len_trim(str)),form='unformatted')
              read(15) wz
            close(15)
            write(*,*) 'finishing reading data'
   
            wx = wx * g
            wy = wy * g
            wz = wz * g
            do kk = 1, lz
            do jj = 1, ly
            do ii = 1, lx1
              a12 = eye * ky(jj) * wx(ii,jj,kk)
              a13 = eye * kz(kk) * wx(ii,jj,kk)
              a21 = eye * kx(ii) * wy(ii,jj,kk)
              a23 = eye * kz(kk) * wy(ii,jj,kk)
              a31 = eye * kx(ii) * wz(ii,jj,kk) 
              a32 = eye * ky(jj) * wz(ii,jj,kk)
   
              wx(ii,jj,kk) = a32 - a23
              wy(ii,jj,kk) = a13 - a31
              wz(ii,jj,kk) = a21 - a12
            end do
            end do
            end do
            
            call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
            
            ip = 0
            do kk=1,nz
            do jj=1,ny
            do ii=1,nx
   
              if ( mod(ii, 2) .eq. 1 ) then
                ll = (ii + 1)/2
                ignore_me = real( wx(ll,jj,kk) ) ** 2 + real( wy(ll,jj,kk) ) ** 2 &
                          + real( wz(ll,jj,kk) ) ** 2
              else
                ll = ii / 2 
                ignore_me = aimag( wx(ll,jj,kk) ) ** 2 + aimag( wy(ll,jj,kk) ) ** 2 &
                          + aimag( wz(ll,jj,kk) ) ** 2
              endif
   
              ignore_me = sqrt(ignore_me)
   
              if (ignore_me .ge. nrms * rmsw)  then
                ip = ip + 1
                xp(ip) = (ii-1) * 2 * pi / nx
                yp(ip) = (jj-1) * 2 * pi / ny
                zp(ip) = (kk-1) * 2 * pi / nz
              end if
   
            end do
            end do
            end do
   
            write(*,*) 'nini =', nini
            write(*,*) 'threshold omega = ', nrms, '* ', rmsw
            write(*,*) 'number of points: ', ip
            write(*,*) 'percentage: ', real(ip) / (nx * ny * nz)
   
            write(strfn, '(I20)') nini
            strfn = adjustl(strfn)

            str = './out/'//strfn(1:len_trim(strfn))//'higho-'//str(1:len_trim(str))
            open(10, file = str(1:len_trim(str)), form = 'unformatted')
              write(10) xp(1:ip), yp(1:ip), zp(1:ip)
            close(10)

            nini = nini + 1
        end do
    close(35)
   
    deallocate(wx,wy,wz)
    deallocate(kx,ky,kz,g, xp, yp, zp)
            
    call destroyplan3d
   
    write(*,*) 'done.'
end program findhwpoints
