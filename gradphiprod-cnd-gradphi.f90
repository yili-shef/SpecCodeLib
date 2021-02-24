program gradphiprod
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer, parameter :: npnt=100
  real(sp), parameter ::  bnd = 10._sp, binw = 2*bnd / npnt
  real(dp), dimension(npnt) :: pdfx, pdfy, pdfz, cndprodx, cndprody, cndprodz
  
  integer :: nx,ny,nz,ii,jj,kk,ll,lx1,lx,ly,lz,ndel,nfile
  real(sp) :: ignore_me, const, delta_c
  real(sp) :: gg(3), aij(3,3), gprod(3)
  real(dp) :: rmsx, rmsy, rmsz, rmsx0

  complex(sp), allocatable, dimension(:,:,:) :: gpx, gpy, gpz
  complex(sp), allocatable, dimension(:,:,:) :: a11, a12, a13
  complex(sp), allocatable, dimension(:,:,:) :: a21, a22, a23
  complex(sp), allocatable, dimension(:,:,:) :: a31, a32, a33
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz
  character(80) :: str, flnm, path, str1
  
  nx=iargc()
  if (nx .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./gradphiprod-cnd-gradphi.x nx ndel filelist'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     ndel: filter scale = ndel * dx'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz

  const = 1._sp/(nx*ny*nz)

  allocate(   g(lx1,ly,lz) )
  allocate( gpx(lx1,ly,lz), gpy(lx1,ly,lz), gpz(lx1,ly,lz) )
  allocate(  kx(lx1,ly,lz),  ky(lx1,ly,lz),  kz(lx1,ly,lz) )
  allocate( a11(lx1,ly,lz), a12(lx1,ly,lz), a13(lx1,ly,lz) )
  allocate( a21(lx1,ly,lz), a22(lx1,ly,lz), a23(lx1,ly,lz) )
  allocate( a31(lx1,ly,lz), a32(lx1,ly,lz), a33(lx1,ly,lz) )

  ! length of displacement
  call getarg(2,str)
  read(str, '(I20)') ndel
  str=adjustl(str)

  delta_c=ndel*2*pi/nx

  ! file list 
  call getarg(3,flnm)
  flnm = adjustl(flnm)


  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  g=exp(-g*delta_c**2/24.)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  pdfx = 0._dp; pdfy = 0._dp; pdfz = 0._dp
  rmsx = 0._dp; rmsy = 0._dp; rmsz = 0._dp
  cndprodx = 0._dp; cndprody = 0._dp; cndprodz = 0._dp
  do while ( .true. )
    read(30,*,end=999) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/phi'//str1(1:len_trim(str1)),form='unformatted')
      read(15) gpz
    close(15)
    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) a11
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) a22
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) a33
    close(15)

    gpx = eye * kx * gpz
    gpy = eye * ky * gpz
    gpz = eye * kz * gpz
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,gpx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,gpy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,gpz,ignore_me)

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

    
    if ( nfile .eq. 0 ) rmsx0 = sqrt( sum( gpx(1:lx,:,:) * conjg(gpx(1:lx,:,:))) * const )
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx
      
        if ( mod(ii,2) .eq. 1) then
            ll = (ii + 1)/2
            gg(1) = real(gpx(ll,jj,kk),sp)
            gg(2) = real(gpy(ll,jj,kk),sp)
            gg(3) = real(gpz(ll,jj,kk),sp)
            aij(1,1) = real(a11(ll,jj,kk), sp)
            aij(1,2) = real(a12(ll,jj,kk), sp)
            aij(1,3) = real(a13(ll,jj,kk), sp)
            aij(2,1) = real(a21(ll,jj,kk), sp)
            aij(2,2) = real(a22(ll,jj,kk), sp)
            aij(2,3) = real(a23(ll,jj,kk), sp)
            aij(3,1) = real(a31(ll,jj,kk), sp)
            aij(3,2) = real(a32(ll,jj,kk), sp)
            aij(3,3) = real(a33(ll,jj,kk), sp)
        else
            ll = ii / 2
            gg(1) = aimag(gpx(ll,jj,kk))
            gg(2) = aimag(gpy(ll,jj,kk))
            gg(3) = aimag(gpz(ll,jj,kk))
            aij(1,1) = aimag(a11(ll,jj,kk))
            aij(1,2) = aimag(a12(ll,jj,kk))
            aij(1,3) = aimag(a13(ll,jj,kk))
            aij(2,1) = aimag(a21(ll,jj,kk))
            aij(2,2) = aimag(a22(ll,jj,kk))
            aij(2,3) = aimag(a23(ll,jj,kk))
            aij(3,1) = aimag(a31(ll,jj,kk))
            aij(3,2) = aimag(a32(ll,jj,kk))
            aij(3,3) = aimag(a33(ll,jj,kk))
        end if

        gprod = matmul(transpose(aij), gg)

    
        ignore_me = gg(1) / rmsx0
        rmsx = rmsx + gg(1)*gg(1)

        ll = 1 + floor((ignore_me + bnd) / binw)
        if (ll .ge. 1 .and. ll .le. npnt) then
          pdfx(ll) = pdfx(ll) + 1
          cndprodx(ll) = cndprodx(ll) + gprod(1)
        end if
 
        ignore_me = gg(2) / rmsx0
        rmsy = rmsy + gg(2)*gg(2)

        ll = 1 + floor((ignore_me + bnd) / binw)
        if (ll .ge. 1 .and. ll .le. npnt) then
          pdfy(ll) = pdfy(ll) + 1
          cndprody(ll) = cndprody(ll) + gprod(2)
        end if
 
        ignore_me = gg(3) / rmsx0
        rmsz = rmsz + gg(3)*gg(3)

        ll = 1 + floor((ignore_me + bnd) / binw)
        if (ll .ge. 1 .and. ll .le. npnt) then
          pdfz(ll) = pdfz(ll) + 1
          cndprodz(ll) = cndprodz(ll) + gprod(3)
        end if
 
    end do
    end do
    end do


    nfile = nfile + 1
  end do
  999 continue
  close(30)

  ignore_me = 1._sp / (nx*ny*nz) / nfile

  rmsx = sqrt( rmsx * ignore_me )
  rmsy = sqrt( rmsy * ignore_me )
  rmsz = sqrt( rmsz * ignore_me )

  cndprodx = cndprodx / (pdfx + mytiny)
  cndprody = cndprody / (pdfy + mytiny)
  cndprodz = cndprodz / (pdfz + mytiny)

  pdfx = pdfx * ignore_me
  pdfy = pdfy * ignore_me
  pdfz = pdfz * ignore_me

  write(*,*) 'normalization pdfx: ', sum(pdfx)
  write(*,*) 'normalization pdfy: ', sum(pdfy)
  write(*,*) 'normalization pdfz: ', sum(pdfz)

  pdfx = pdfx / binw
  pdfy = pdfy / binw
  pdfz = pdfz / binw

  
  path='gradphiprod-cnd-grad-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)

    write(15,'(''# Title = "'',4E12.4, ''"'')') rmsx, rmsy, rmsz
    write(15,'(''# variables = "xx","cndprodx","pdfgx","yy","cndprody","pdfgy","zz","cndprodz","pdfgz"'')')

    do lz=1,npnt

      ignore_me = -bnd + (lz - .5) * binw

      write(15,'(10E12.4)') ignore_me * rmsx0 / rmsx, cndprodx(lz),  pdfx(lz) * rmsx / rmsx0, &
      ignore_me * rmsx0 / rmsy, cndprody(lz), rmsy / rmsx0 * pdfy(lz), ignore_me * rmsx0 / rmsz, &
      cndprodz(lz), rmsz / rmsx0 * pdfz(lz)

    end do

  close(15)
  
  
  deallocate(gpx, gpy, gpz, a11, a12, a13, a21, a22, a23, a31, a32, a33, g, kx, ky, kz)
  
  call destroyplan3d

  write(*,*) 'gradphiprod-cnd-gradphi.x done.'

end program gradphiprod
