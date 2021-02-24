program bijongrid
    use mconstant
    use mwavenumber
    use mfftwplan3d
    implicit none

    integer, parameter :: norder = 6
 
    integer :: nx, ny, nz, nprtcle, ii, jj, kk, ll, ip, nfile1, nfile2
    integer :: lx, lx1, ly, lz
 
    real(sp), allocatable, dimension(:,:,:) :: b, ap
    real(sp), allocatable, dimension(:,:,:) :: a11r, a12r, a13r
    real(sp), allocatable, dimension(:,:,:) :: a21r, a22r, a23r
    real(sp), allocatable, dimension(:,:,:) :: a31r, a32r, a33r
    complex(sp), allocatable, dimension(:,:,:) :: a11, a12, a13
    complex(sp), allocatable, dimension(:,:,:) :: a21, a22, a23
    complex(sp), allocatable, dimension(:,:,:) :: a31, a32, a33

    real(sp), allocatable, dimension(:) :: xp, yp, zp, kx, ky, kz

    integer, parameter :: idnty(3,3) = (/1, 0, 0, 0, 1, 0, 0, 0, 1/)
    real(sp), dimension(3,3) :: app, b1, b2, b3

    real(sp) :: y(3), dxyz(3), bg(norder,3)
    integer  :: ix(norder), iy(norder), iz(norder), lhnode(3)
    integer  :: ix1, iy1, iz1

    real(sp) :: dt, ignore_me
    character(80) :: prefix, str, fstr1, fstr2
 
    write(*,*)
    write(*,'(''>>> Deformation gradient on the grid points <<<'')' )
    write(*,*)
    write(*,'(''>>> Need velocity field data and fluid particle location data <<<'')')
    write(*,*)

    if ( iargc() .ne. 5) stop 'Usage: ./bijongrid.x nx dt firstfile-no lastfile-no prefix-for-xyz'

    call getarg(1, str)
    read(str, '(I20)') nx

    call getarg(2, str)
    read(str, '(F20.10)') dt

    call getarg(3, fstr1)
    fstr1 = adjustl(fstr1)
    read(fstr1, '(I20)') nfile1

    call getarg(4, fstr2)
    fstr2 = adjustl(fstr2)
    read(fstr2, '(I20)') nfile2
 
    call getarg(5, prefix)
    prefix = adjustl(prefix)

    ny = nx; nz = nx; lx = nx / 2; lx1 = lx + 1; ly = ny; lz = nz
    nprtcle = nx * ny * nz

    dxyz(1) = 2 * pi / nx; dxyz(2) = 2 * pi / ny; dxyz(3) = 2 * pi / nz

    call fftwplan3de(nx, ny, nz)

    allocate( b(3,3,nprtcle), xp(nprtcle), yp(nprtcle), zp(nprtcle) )
    allocate( a11(lx1,ly,lz), a12(lx1,ly,lz), a13(lx1,ly,lz) )
    allocate( a21(lx1,ly,lz), a22(lx1,ly,lz), a23(lx1,ly,lz) )
    allocate( a31(lx1,ly,lz), a32(lx1,ly,lz), a33(lx1,ly,lz) )
    allocate( a11r(nx,ny,nz), a12r(nx,ny,nz), a13r(nx,ny,nz) )
    allocate( a21r(nx,ny,nz), a22r(nx,ny,nz), a23r(nx,ny,nz) )
    allocate( a31r(nx,ny,nz), a32r(nx,ny,nz), a33r(nx,ny,nz) )
    allocate( ap(3,3,nprtcle) )
    allocate( kx(lx1), ky(ly), kz(lz) )
 
    call wavenumber(kx, ky, kz, lx1, ly, lz)

    do ii=1,nprtcle
        b(:,:,ii)=idnty
    end do
    
    ! reading the first file

    str = adjustl(fstr1) ! str to be used to in readdata
    str = str(1:len_trim(str))//'.dat'
    write(*,*) 'reading data file: ', str(1:len_trim(str))

    call readdata
    call uitoaijr
 
    ! interpolate on the first gradient field
    do ip=1, nprtcle
        y(1)=xp(ip)
        y(2)=yp(ip)
        y(3)=zp(ip)
       
        call pre_interp(y,dxyz,bg,lhnode)
        call ixiyiz(lhnode,nx,ny,nz,ix,iy,iz)

        app = 0._sp
        do  ii=1,norder
          ix1 = ix(ii)
          do  jj=1,norder
            iy1 = iy(jj)
            do  kk=1,norder
              iz1 = iz(kk)

              app(1,1)=app(1,1)+a11r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
              app(1,2)=app(1,2)+a12r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
              app(1,3)=app(1,3)+a13r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)

              app(2,1)=app(2,1)+a21r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
              app(2,2)=app(2,2)+a22r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
              app(2,3)=app(2,3)+a23r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)

              app(3,1)=app(3,1)+a31r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
              app(3,2)=app(3,2)+a32r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
              app(3,3)=app(3,3)+a33r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
            end do
          end do
        end do
        ap(:,:,ip) = app
    end do
 
    do ll = nfile1 + 1, nfile2
 
        write(str, '(I20)') ll
        str = adjustl(str)
        str = str(1:len_trim(str))//'.dat'
        write(*,*) 'reading data file: ', str(1:len_trim(str))
     
        call readdata
        call uitoaijr
 
        do ip = 1, nprtcle
            y(1) = xp(ip)
            y(2) = yp(ip)
            y(3) = zp(ip)
       
            call pre_interp(y,dxyz,bg,lhnode)
            call ixiyiz(lhnode,nx,ny,nz,ix,iy,iz)

            app = 0._sp
            do  ii=1,norder
              ix1 = ix(ii)
              do  jj=1,norder
                iy1 = iy(jj)
                do  kk=1,norder
                  iz1 = iz(kk)
         
                  app(1,1)=app(1,1)+a11r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
                  app(1,2)=app(1,2)+a12r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
                  app(1,3)=app(1,3)+a13r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
         
                  app(2,1)=app(2,1)+a21r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
                  app(2,2)=app(2,2)+a22r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
                  app(2,3)=app(2,3)+a23r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
         
                  app(3,1)=app(3,1)+a31r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
                  app(3,2)=app(3,2)+a32r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
                  app(3,3)=app(3,3)+a33r(ix1,iy1,iz1)*bg(ii,1)* bg(jj,2)*bg(kk,3)
                end do
              end do
            end do

            b1 = b(:,:,ip) +  .5_sp * dt * matmul( ap(:,:,ip), b(:,:,ip) )
            b2 = b(:,:,ip) + .25_sp * dt * matmul( ap(:,:,ip) + app, b1 )
            b3 = b(:,:,ip) +  .5_sp * dt * matmul( ap(:,:,ip) + app, b2 )
            b(:,:,ip) = b(:,:,ip) + (dt/6._sp) * ( matmul( ap(:,:,ip), b(:,:,ip) ) + &
                       matmul( ap(:,:,ip) + app, b1 + b2 ) + matmul( app, b3 ) )
           
            ap(:,:,ip) = app ! save for next step
        end do

    end do

    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

        ll = (kk-1) * nx * ny + (jj-1) * nx + ii
       
        a11r(ii,jj,kk) = b(1,1,ll)
        a12r(ii,jj,kk) = b(1,2,ll)
        a13r(ii,jj,kk) = b(1,3,ll)
        a21r(ii,jj,kk) = b(2,1,ll)
        a22r(ii,jj,kk) = b(2,2,ll)
        a23r(ii,jj,kk) = b(2,3,ll)
        a31r(ii,jj,kk) = b(3,1,ll)
        a32r(ii,jj,kk) = b(3,2,ll)
        a33r(ii,jj,kk) = b(3,3,ll)

    end do
    end do
    end do

    a11(1:lx,:,:) = cmplx( a11r(1:nx:2,:,:), a11r(2:nx:2,:,:) )
    a12(1:lx,:,:) = cmplx( a12r(1:nx:2,:,:), a12r(2:nx:2,:,:) )
    a13(1:lx,:,:) = cmplx( a13r(1:nx:2,:,:), a13r(2:nx:2,:,:) )

    a21(1:lx,:,:) = cmplx( a21r(1:nx:2,:,:), a21r(2:nx:2,:,:) )
    a22(1:lx,:,:) = cmplx( a22r(1:nx:2,:,:), a22r(2:nx:2,:,:) )
    a23(1:lx,:,:) = cmplx( a23r(1:nx:2,:,:), a23r(2:nx:2,:,:) )

    a31(1:lx,:,:) = cmplx( a31r(1:nx:2,:,:), a31r(2:nx:2,:,:) )
    a32(1:lx,:,:) = cmplx( a32r(1:nx:2,:,:), a32r(2:nx:2,:,:) )
    a33(1:lx,:,:) = cmplx( a33r(1:nx:2,:,:), a33r(2:nx:2,:,:) )

    fstr2 = 'b11from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
    open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
        write(25) a11
    close(25)
    fstr2 = 'b12from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
    open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
        write(25) a12
    close(25)
    fstr2 = 'b13from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
    open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
        write(25) a13
    close(25)
    fstr2 = 'b21from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
    open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
        write(25) a21
    close(25)
    fstr2 = 'b22from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
    open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
        write(25) a22
    close(25)
    fstr2 = 'b23from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
    open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
        write(25) a23
    close(25)
    fstr2 = 'b31from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
    open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
        write(25) a31
    close(25)
    fstr2 = 'b32from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
    open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
        write(25) a32
    close(25)
    fstr2 = 'b33from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
    open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
        write(25) a33
    close(25)


    call destroyplan3d

    deallocate(b, ap, kx, ky, kz, a11, a12, a13, a21, a22, a23, a31, a32, a33)
    deallocate(a11r, a12r, a13r, a21r, a22r, a23r, a31r, a32r, a33r)

    write(*,*) 'Finished.'
    stop

contains 

    subroutine readdata

        open(35, file = './out/'//prefix(1:len_trim(prefix))//str(1:len_trim(str)), form = 'unformatted')
          read(35) xp, yp, zp
        close(35)

        open(35, file = './out/ux'//str(1:len_trim(str)), form = 'unformatted')
          read(35) a11
        close(35)
        open(35, file = './out/uy'//str(1:len_trim(str)), form = 'unformatted')
          read(35) a22
        close(35)
        open(35, file = './out/uz'//str(1:len_trim(str)), form = 'unformatted')
          read(35) a33
        close(35)

    end subroutine readdata

    subroutine uitoaijr

        do kk = 1, lz
        do jj = 1, ly
        do ii = 1, lx1
            a12(ii,jj,kk) = eye * ky(jj) * a11(ii,jj,kk)
            a13(ii,jj,kk) = eye * kz(kk) * a11(ii,jj,kk)
            a11(ii,jj,kk) = eye * kx(ii) * a11(ii,jj,kk)

            a21(ii,jj,kk) = eye * kx(ii) * a22(ii,jj,kk)
            a23(ii,jj,kk) = eye * kz(kk) * a22(ii,jj,kk)
            a22(ii,jj,kk) = eye * ky(jj) * a22(ii,jj,kk)

            a31(ii,jj,kk) = eye * kx(ii) * a33(ii,jj,kk)
            a32(ii,jj,kk) = eye * ky(jj) * a33(ii,jj,kk)
            a33(ii,jj,kk) = eye * kz(kk) * a33(ii,jj,kk)
        end do
        end do
        end do

        call rfftwnd_f77_one_complex_to_real(c2r3d,a11,ignore_me)
        call rfftwnd_f77_one_complex_to_real(c2r3d,a12,ignore_me)
        call rfftwnd_f77_one_complex_to_real(c2r3d,a13,ignore_me)
        call rfftwnd_f77_one_complex_to_real(c2r3d,a21,ignore_me)
        call rfftwnd_f77_one_complex_to_real(c2r3d,a22,ignore_me)
        call rfftwnd_f77_one_complex_to_real(c2r3d,a23,ignore_me)
        call rfftwnd_f77_one_complex_to_real(c2r3d,a31,ignore_me)
        call rfftwnd_f77_one_complex_to_real(c2r3d,a32,ignore_me)
        call rfftwnd_f77_one_complex_to_real(c2r3d,a33,ignore_me)

        a11r(1:nx:2,:,:)=real(a11(1:lx,:,:),sp); a11r(2:nx:2,:,:)=aimag(a11(1:lx,:,:))
        a12r(1:nx:2,:,:)=real(a12(1:lx,:,:),sp); a12r(2:nx:2,:,:)=aimag(a12(1:lx,:,:))
        a13r(1:nx:2,:,:)=real(a13(1:lx,:,:),sp); a13r(2:nx:2,:,:)=aimag(a13(1:lx,:,:))

        a21r(1:nx:2,:,:)=real(a21(1:lx,:,:),sp); a21r(2:nx:2,:,:)=aimag(a21(1:lx,:,:))
        a22r(1:nx:2,:,:)=real(a22(1:lx,:,:),sp); a22r(2:nx:2,:,:)=aimag(a22(1:lx,:,:))
        a23r(1:nx:2,:,:)=real(a23(1:lx,:,:),sp); a23r(2:nx:2,:,:)=aimag(a23(1:lx,:,:))

        a33r(1:nx:2,:,:)=real(a33(1:lx,:,:),sp); a33r(2:nx:2,:,:)=aimag(a33(1:lx,:,:))
        a32r(1:nx:2,:,:)=real(a32(1:lx,:,:),sp); a32r(2:nx:2,:,:)=aimag(a32(1:lx,:,:))
        a31r(1:nx:2,:,:)=real(a31(1:lx,:,:),sp); a31r(2:nx:2,:,:)=aimag(a31(1:lx,:,:))

    end subroutine uitoaijr

end program bijongrid      
