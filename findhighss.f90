program findhwpoints
    use mconstant
    use mfftwplan3d
    implicit none
    
    real(sp), parameter :: rmsw = 6.3079 ! calculated from hr256.list
    integer, parameter :: nx = 256, ndel = 0, nrms = 2.5
 
    integer :: ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll, ip, nini
    real(sp) :: ignore_me, delta_c
 
    complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
    real(sp),    allocatable, dimension(:,:,:) :: g
    real(sp),    allocatable, dimension(:)     :: kx, ky, kz, xp, yp, zp

    character(80) :: str, strfn
    
    ny=nx; nz=nx
    lx=nx/2;lx1=nx/2+1;ly=ny;lz=nz
    delta_c=ndel*2*pi/nx
 
    allocate(s11(lx1,ly,lz),s22(lx1,ly,lz),s33(lx1,ly,lz))
    allocate(s12(lx1,ly,lz),s13(lx1,ly,lz),s23(lx1,ly,lz))
    allocate(g (lx1,ly,lz),kx(lx1), ky(ly), kz(lz))
    allocate(xp(nx*ny*nz), yp(nx*ny*nz), zp(nx*ny*nz))
 
    call fftwplan3de(nx,ny,nz)
    write(*,*) 'after fftwplan'
 
    call wavenumber(kx,ky,kz,g,lx1,ly,lz)
    write(*,*) 'after wavenumber'
 
    ! Gaussian filter
    g=exp(-g*delta_c**2/24.)
 
    nini = 1
    open(35, file = 'findhs.list')
 
        do while ( .not. eof(35) )
            read(35,*) str
 
            open(15,file='./out/ux'//str(1:len_trim(str)),form='unformatted')
              read(15) s11
            close(15)
            open(15,file='./out/uy'//str(1:len_trim(str)),form='unformatted')
              read(15) s22
            close(15)
            open(15,file='./out/uz'//str(1:len_trim(str)),form='unformatted')
              read(15) s33
            close(15)
   
            s11 = s11 * g
            s22 = s22 * g
            s33 = s33 * g
            do kk = 1, lz
            do jj = 1, ly
            do ii = 1, lx1
              s12(ii,jj,kk) = .5 * eye * ( kx(ii) * s22(ii,jj,kk) + ky(jj) * s11(ii,jj,kk) )
              s13(ii,jj,kk) = .5 * eye * ( kx(ii) * s33(ii,jj,kk) + kz(kk) * s11(ii,jj,kk) )
              s23(ii,jj,kk) = .5 * eye * ( ky(jj) * s33(ii,jj,kk) + kz(kk) * s22(ii,jj,kk) )

              s11(ii,jj,kk) = eye * kx(ii) * s11(ii,jj,kk)
              s22(ii,jj,kk) = eye * ky(jj) * s22(ii,jj,kk)
              s33(ii,jj,kk) = eye * kz(kk) * s33(ii,jj,kk)
            end do
            end do
            end do
            
            call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
            
            ip = 0
            do kk=1,nz
            do jj=1,ny
            do ii=1,nx
   
              if ( mod(ii, 2) .eq. 1 ) then
                ll = (ii + 1)/2
                ignore_me = real( s11(ll,jj,kk) ) ** 2 + real( s22(ll,jj,kk) ) ** 2 &
                          + real( s33(ll,jj,kk) ) ** 2 + 2. * real( s12(ll,jj,kk) ) ** 2 &
                          + 2. * real( s13(ll,jj,kk) ) ** 2 + 2. * real( s23(ll,jj,kk) ) ** 2 
              else
                ll = ii / 2 
                ignore_me = aimag( s11(ll,jj,kk) ) ** 2 + aimag( s22(ll,jj,kk) ) ** 2 &
                          + aimag( s33(ll,jj,kk) ) ** 2 + 2. * aimag( s12(ll,jj,kk) ) ** 2 &
                          + 2. * aimag( s13(ll,jj,kk) ) ** 2 + 2. * aimag( s23(ll,jj,kk) ) ** 2 
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
   
            write(*,*) ip, real(ip) / (nx * ny * nz)
   
            write(strfn, '(I20)') nini
            strfn = adjustl(strfn)

            str = './out/'//strfn(1:len_trim(strfn))//'highs-'//str(1:len_trim(str))
            open(10, file = str(1:len_trim(str)), form = 'unformatted')
              write(10) xp(1:ip), yp(1:ip), zp(1:ip)
            close(10)

            nini = nini + 1
        end do
    close(35)
   
    deallocate(s11,s12,s13,s22,s23,s33)
    deallocate(kx,ky,kz,g, xp, yp, zp)
            
    call destroyplan3d
   
    write(*,*) 'done.'
end program findhwpoints
