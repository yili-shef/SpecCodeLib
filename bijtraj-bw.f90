program bijongrid
    use mconstant
    use mwavenumber
    use mfftwplan3d
    implicit none

    integer, parameter :: norder = 6
 
    integer :: nx, ny, nz, nprtcle, ii, jj, kk, ll, mm, ip, nfile1, nfile2
    integer :: lx, lx1, ly, lz, inc
 
    real(sp), allocatable, dimension(:,:,:) :: b, ap
    real(sp), allocatable, dimension(:,:,:) :: a11r, a12r, a13r
    real(sp), allocatable, dimension(:,:,:) :: a21r, a22r, a23r
    real(sp), allocatable, dimension(:,:,:) :: a31r, a32r, a33r
    real(sp), allocatable, dimension(:,:,:) :: kx, ky, kz
    complex(sp), allocatable, dimension(:,:,:) :: a11, a12, a13
    complex(sp), allocatable, dimension(:,:,:) :: a21, a22, a23
    complex(sp), allocatable, dimension(:,:,:) :: a31, a32, a33

    real(sp), allocatable, dimension(:) :: xp, yp, zp

    integer, parameter :: idnty(3,3) = (/1, 0, 0, 0, 1, 0, 0, 0, 1/)
    real(sp), dimension(3,3) :: app, b1, b2, b3

    real(sp) :: y(3), dxyz(3), bg(norder,3)
    integer  :: ix(norder), iy(norder), iz(norder), lhnode(3)
    integer  :: ix1, iy1, iz1, info

    real(sp) :: dt, ignore_me
    character(80) :: prefix, str, fstr1, fstr2, mark, inifname
 
    write(*,*)
    write(*,'(''>>> Deformation gradient on the grid points <<<'')' )
    write(*,*)
    write(*,'(''>>> Need velocity field data and fluid particle location data <<<'')')
    write(*,*)

    if ( iargc() .ne. 6) stop 'Usage: ./bijongrid.x nx dt 1stfile-no lastfile-no pref-for-xyz initfile'

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

    if ( nfile1 .gt. nfile2 ) then
        dt = -dt ! backward tracking
        mark = 'bw'
        inc = -1
    else
        mark = 'fw'
        inc = 1
    end if
    mark = adjustl(mark)
 
    call getarg(5, prefix)
    prefix = adjustl(prefix)

    call getarg(6,inifname)
    inifname = adjustl(inifname)

    ny = nx; nz = nx; lx = nx / 2; lx1 = lx + 1; ly = ny; lz = nz
    nprtcle = nx * ny * nz

    dxyz(1) = 2 * pi / nx; dxyz(2) = 2 * pi / ny; dxyz(3) = 2 * pi / nz

    call fftwplan3d(nx, ny, nz)

    allocate( b(3,3,nprtcle), xp(nprtcle), yp(nprtcle), zp(nprtcle) )
    allocate( a11(lx1,ly,lz), a12(lx1,ly,lz), a13(lx1,ly,lz) )
    allocate( a21(lx1,ly,lz), a22(lx1,ly,lz), a23(lx1,ly,lz) )
    allocate( a31(lx1,ly,lz), a32(lx1,ly,lz), a33(lx1,ly,lz) )
    allocate( a11r(nx,ny,nz), a12r(nx,ny,nz), a13r(nx,ny,nz) )
    allocate( a21r(nx,ny,nz), a22r(nx,ny,nz), a23r(nx,ny,nz) )
    allocate( a31r(nx,ny,nz), a32r(nx,ny,nz), a33r(nx,ny,nz) )
    allocate( ap(3,3,nprtcle) )
    allocate( kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz) )
 
    call wavenumber(kx, ky, kz, lx1, ly, lz)

    inquire( file = './out/b11'//inifname(1:len_trim(inifname))//'.dat', exist = info)

    if ( info .eq. .false. ) then

      ! start from spheres if initial field not exist
      do ii=1,nprtcle 
          b(:,:,ii)=idnty
      end do

      write(*,*) '!!!!!!!!!! WARNING: No initial field. Start from spheres !!!!!!!!!!!'

    else

      open(30, file = './out/b11'//inifname(1:len_trim(inifname))//'.dat', form = 'unformatted')
        read(30) a11
      close(30)
      open(30, file = './out/b12'//inifname(1:len_trim(inifname))//'.dat', form = 'unformatted')
        read(30) a12
      close(30)
      open(30, file = './out/b13'//inifname(1:len_trim(inifname))//'.dat', form = 'unformatted')
        read(30) a13
      close(30)
      open(30, file = './out/b21'//inifname(1:len_trim(inifname))//'.dat', form = 'unformatted')
        read(30) a21
      close(30)
      open(30, file = './out/b22'//inifname(1:len_trim(inifname))//'.dat', form = 'unformatted')
        read(30) a22
      close(30)
      open(30, file = './out/b23'//inifname(1:len_trim(inifname))//'.dat', form = 'unformatted')
        read(30) a23
      close(30)
      open(30, file = './out/b31'//inifname(1:len_trim(inifname))//'.dat', form = 'unformatted')
        read(30) a31
      close(30)
      open(30, file = './out/b32'//inifname(1:len_trim(inifname))//'.dat', form = 'unformatted')
        read(30) a32
      close(30)
      open(30, file = './out/b33'//inifname(1:len_trim(inifname))//'.dat', form = 'unformatted')
        read(30) a33
      close(30)

      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx
     
          ll = (kk-1) * nx * ny + (jj-1) * nx + ii

          if ( mod(ii, 2) .eq. 1 ) then
              mm = ii/2 + 1
         
              b(1,1,ll) = real( a11(mm,jj,kk),sp )
              b(1,2,ll) = real( a12(mm,jj,kk),sp )
              b(1,3,ll) = real( a13(mm,jj,kk),sp )
              b(2,1,ll) = real( a21(mm,jj,kk),sp )
              b(2,2,ll) = real( a22(mm,jj,kk),sp )
              b(2,3,ll) = real( a23(mm,jj,kk),sp )
              b(3,1,ll) = real( a31(mm,jj,kk),sp )
              b(3,2,ll) = real( a32(mm,jj,kk),sp )
              b(3,3,ll) = real( a33(mm,jj,kk),sp )
          else
              mm = ii/2
         
              b(1,1,ll) = aimag( a11(mm,jj,kk) )
              b(1,2,ll) = aimag( a12(mm,jj,kk) )
              b(1,3,ll) = aimag( a13(mm,jj,kk) )
              b(2,1,ll) = aimag( a21(mm,jj,kk) )
              b(2,2,ll) = aimag( a22(mm,jj,kk) )
              b(2,3,ll) = aimag( a23(mm,jj,kk) )
              b(3,1,ll) = aimag( a31(mm,jj,kk) )
              b(3,2,ll) = aimag( a32(mm,jj,kk) )
              b(3,3,ll) = aimag( a33(mm,jj,kk) )
          end if
     
      end do
      end do
      end do

    end if
    
    ! reading the first file

    str = adjustl(fstr1) ! NOTE: str to be used to in readdata 
    str = str(1:len_trim(str))//'.dat'
    write(*,*) 'bijtraj-bw.x: reading data file: ', str(1:len_trim(str))

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
 
    do ll = nfile1 + inc, nfile2, inc ! inc = 1 for fw, = -1 for bw
 
        ! output  (internal subprocedure)
        call writebij ! Need to use str to specify output file names

        write(str, '(I20)') ll
        str = adjustl(str)
        str = str(1:len_trim(str))//'.dat'
        write(*,*) 'bijtraj-bw.x: reading data file: ', str(1:len_trim(str))
     
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

    ! output last files
    call writebij

    call destroyplan3d

    deallocate(b, ap, kx, ky, kz, a11, a12, a13, a21, a22, a23, a31, a32, a33)
    deallocate(a11r, a12r, a13r, a21r, a22r, a23r, a31r, a32r, a33r)

    write(*,*) 'bijtraj-bw.x Finished.'
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

        a12 = eye * ky * a11
        a13 = eye * kz * a11
        a11 = eye * kx * a11

        a21 = eye * kx * a22
        a23 = eye * kz * a22
        a22 = eye * ky * a22

        a31 = eye * kx * a33
        a32 = eye * ky * a33
        a33 = eye * kz * a33

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

    subroutine writebij

        integer :: llll

        do kk = 1, nz
        do jj = 1, ny
        do ii = 1, nx
  
            llll = (kk-1) * nx * ny + (jj-1) * nx + ii
           
            a11r(ii,jj,kk) = b(1,1,llll)
            a12r(ii,jj,kk) = b(1,2,llll)
            a13r(ii,jj,kk) = b(1,3,llll)
            a21r(ii,jj,kk) = b(2,1,llll)
            a22r(ii,jj,kk) = b(2,2,llll)
            a23r(ii,jj,kk) = b(2,3,llll)
            a31r(ii,jj,kk) = b(3,1,llll)
            a32r(ii,jj,kk) = b(3,2,llll)
            a33r(ii,jj,kk) = b(3,3,llll)
  
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
  
        fstr2 = 'b11'//mark(1:len_trim(mark))//'from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
        open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
            write(25) a11
        close(25)
        fstr2 = 'b12'//mark(1:len_trim(mark))//'from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
        open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
            write(25) a12
        close(25)
        fstr2 = 'b13'//mark(1:len_trim(mark))//'from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
        open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
            write(25) a13
        close(25)
        fstr2 = 'b21'//mark(1:len_trim(mark))//'from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
        open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
            write(25) a21
        close(25)
        fstr2 = 'b22'//mark(1:len_trim(mark))//'from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
        open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
            write(25) a22
        close(25)
        fstr2 = 'b23'//mark(1:len_trim(mark))//'from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
        open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
            write(25) a23
        close(25)
        fstr2 = 'b31'//mark(1:len_trim(mark))//'from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
        open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
            write(25) a31
        close(25)
        fstr2 = 'b32'//mark(1:len_trim(mark))//'from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
        open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
            write(25) a32
        close(25)
        fstr2 = 'b33'//mark(1:len_trim(mark))//'from'//fstr1(1:len_trim(fstr1))//'-'//str(1:len_trim(str))
        open(25, file = './out/'//fstr2(1:len_trim(fstr2)), form = 'unformatted')
            write(25) a33
        close(25)

    end subroutine writebij

end program bijongrid      
