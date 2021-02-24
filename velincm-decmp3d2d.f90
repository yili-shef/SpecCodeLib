program velincmdecmp
  use mconstant
  use mfftwplan3d
  use mfftwplan2d
  implicit none
  
  integer :: nx,ny,nz,lx1,lx,ly,lz,nfile,ndel
  real(sp) :: ignore_me

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz
  real(sp),    allocatable, dimension(:,:,:) :: ur, dur, tmp

  complex(sp), allocatable, dimension(:,:) :: u2dx, u2dy, u2dz
  real(sp),    allocatable, dimension(:,:) :: u2dr, du2dr, tmp2d

  character(80) :: str, flnm, path

  ! == increments of 3D velocity over displacement in perpendicular direction ==
  ! longitudinal 
  real(dp) :: mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10  
  ! transverse of in plane component
  real(dp) :: mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, mv10  
  ! transverse of uz
  real(dp) :: mduzdt2, mduzdt3, mduzdt4, mduzdt5, mduzdt6, mduzdt7, mduzdt8, mduzdt9, mduzdt10

  ! == increments of 3D velocity over displacement in parallel direction ==
  ! longitudinal (i.e. uz)
  real(dp) :: mdudz2, mdudz3, mdudz4, mdudz5, mdudz6, mdudz7, mdudz8, mdudz9, mdudz10
  ! longitudinal (including contributions from ux and uy)
  real(dp) :: mdvdz2, mdvdz3, mdvdz4, mdvdz5, mdvdz6, mdvdz7, mdvdz8, mdvdz9, mdvdz10
  
  ! == increments of 2D velocity over displacement in perpendicular direction ==
  ! longtiduinal (from ux and uy)
  real(dp) :: mu2d2, mu2d3, mu2d4, mu2d5, mu2d6, mu2d7, mu2d8, mu2d9, mu2d10
  ! transverse from ux and uy
  real(dp) :: mv2d2, mv2d3, mv2d4, mv2d5, mv2d6, mv2d7, mv2d8, mv2d9, mv2d10
  ! transverse for uz
  real(dp) :: mu2dz2, mu2dz3, mu2dz4, mu2dz5, mu2dz6, mu2dz7, mu2dz8, mu2dz9, mu2dz10
  
  nx=iargc()
  if (nx .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./velincm-decmp3d2d.x nx filelist ndel'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     ndel: displacement = ndel * dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! file number 
  call getarg(2,flnm)
  flnm = adjustl(flnm)

  ! filter scale
  call getarg(3,str)
  read(str,'(I20)') ndel
  str=adjustl(str)

  ny = nx; nz = nx
  lx = nx/2; lx1 = nx/2 + 1; ly = ny; lz = nz

  allocate( ux(lx1,ly,lz), uy(lx1,ly,lz), uz(lx1,ly,lz) )
  allocate( ur(nx,ny,nz),  dur(nx,ny,nz), tmp(nx,ny,nz) )
  allocate( u2dx(lx1,ly), u2dy(lx1,ly), u2dz(lx1,ly) )
  allocate( u2dr(nx,ny),  du2dr(nx,ny), tmp2d(nx,ny) )

  write(*,*) 'fftwplan'
  call fftwplan3de(nx, ny, nz)
  call fftwplan2de(nx, ny)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  mu2  = 0._dp; mv2  = 0._dp;  mduzdt2  = 0._dp; mdudz2  = 0._dp; mdvdz2  = 0._dp
  mu3  = 0._dp; mv3  = 0._dp;  mduzdt3  = 0._dp; mdudz3  = 0._dp; mdvdz3  = 0._dp
  mu4  = 0._dp; mv4  = 0._dp;  mduzdt4  = 0._dp; mdudz4  = 0._dp; mdvdz4  = 0._dp
  mu5  = 0._dp; mv5  = 0._dp;  mduzdt5  = 0._dp; mdudz5  = 0._dp; mdvdz5  = 0._dp
  mu6  = 0._dp; mv6  = 0._dp;  mduzdt6  = 0._dp; mdudz6  = 0._dp; mdvdz6  = 0._dp
  mu7  = 0._dp; mv7  = 0._dp;  mduzdt7  = 0._dp; mdudz7  = 0._dp; mdvdz7  = 0._dp
  mu8  = 0._dp; mv8  = 0._dp;  mduzdt8  = 0._dp; mdudz8  = 0._dp; mdvdz8  = 0._dp
  mu9  = 0._dp; mv9  = 0._dp;  mduzdt9  = 0._dp; mdudz9  = 0._dp; mdvdz9  = 0._dp
  mu10 = 0._dp; mv10 = 0._dp;  mduzdt10 = 0._dp; mdudz10 = 0._dp; mdvdz10 = 0._dp

  mu2d2  = 0._dp; mv2d2  = 0._dp; mu2dz2  = 0._dp
  mu2d3  = 0._dp; mv2d3  = 0._dp; mu2dz3  = 0._dp
  mu2d4  = 0._dp; mv2d4  = 0._dp; mu2dz4  = 0._dp
  mu2d5  = 0._dp; mv2d5  = 0._dp; mu2dz5  = 0._dp
  mu2d6  = 0._dp; mv2d6  = 0._dp; mu2dz6  = 0._dp
  mu2d7  = 0._dp; mv2d7  = 0._dp; mu2dz7  = 0._dp
  mu2d8  = 0._dp; mv2d8  = 0._dp; mu2dz8  = 0._dp
  mu2d9  = 0._dp; mv2d9  = 0._dp; mu2dz9  = 0._dp
  mu2d10 = 0._dp; mv2d10 = 0._dp; mu2dz10 = 0._dp

  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) path
    write(*,*) 'data file ', path(1:len_trim(path))

    open(15,file='./out/ux'//path(1:len_trim(path)),form='unformatted')
      read(15) ux
    close(15)
    open(15,file='./out/uy'//path(1:len_trim(path)),form='unformatted')
      read(15) uy
    close(15)
    open(15,file='./out/uz'//path(1:len_trim(path)),form='unformatted')
      read(15) uz
    close(15)

    u2dx = ux(:,:,1) ! components with kz=0
    u2dy = uy(:,:,1)
    u2dz = uz(:,:,1)

    ux(:,:,1) = 0._sp ! removing kz=0 components
    uy(:,:,1) = 0._sp
    uz(:,:,1) = 0._sp
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)

    call rfftwnd_f77_one_complex_to_real(c2r2d,u2dx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r2d,u2dy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r2d,u2dz,ignore_me)
    
    
    ! =============== 3D components : perpendicular =================
    ! ux
    ur(1:nx:2,:,:) =  real( ux(1:lx,:,:) )
    ur(2:nx:2,:,:) = aimag( ux(1:lx,:,:) )

    ! longitudinal - perpendicular
    dur = abs(cshift(ur, ndel, 1) - ur)
    call calsum(mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10)

    ! transverse - perpendicular
    dur = abs(cshift(ur, ndel, 2) - ur)
    call calsum(mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, mv10)

    ! uy
    ur(1:nx:2,:,:) =  real( uy(1:lx,:,:) )
    ur(2:nx:2,:,:) = aimag( uy(1:lx,:,:) )

    ! longitudinal - perpendicular
    dur = abs(cshift(ur, ndel, 2) - ur)
    call calsum(mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10)

    ! Transverse - perpendicular
    dur = abs(cshift(ur, ndel, 1) - ur)
    call calsum(mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, mv10)

    ! transverse for uz 
    ur(1:nx:2,:,:) =  real( uz(1:lx,:,:) )
    ur(2:nx:2,:,:) = aimag( uz(1:lx,:,:) )

    dur = abs(cshift(ur, ndel, 2) - ur)
    call calsum(mduzdt2, mduzdt3, mduzdt4, mduzdt5, mduzdt6, mduzdt7, mduzdt8, mduzdt9, mduzdt10)
    dur = abs(cshift(ur, ndel, 1) - ur)
    call calsum(mduzdt2, mduzdt3, mduzdt4, mduzdt5, mduzdt6, mduzdt7, mduzdt8, mduzdt9, mduzdt10)

    ! =================== 2D component : perpendicular ===============

    ! u2dx 
    u2dr(1:nx:2,:) =  real( u2dx(1:lx,:) )
    u2dr(2:nx:2,:) = aimag( u2dx(1:lx,:) )

    ! longitudinal - perpendicu2dlar
    du2dr = abs(cshift(u2dr, ndel, 1) - u2dr)
    call calsum2d(mu2d2, mu2d3, mu2d4, mu2d5, mu2d6, mu2d7, mu2d8, mu2d9, mu2d10)

    ! transverse - perpendicular
    du2dr = abs(cshift(u2dr, ndel, 2) - u2dr)
    call calsum2d(mv2d2, mv2d3, mv2d4, mv2d5, mv2d6, mv2d7, mv2d8, mv2d9, mv2d10)

    ! u2dy 
    u2dr(1:nx:2,:) =  real( u2dy(1:lx,:) )
    u2dr(2:nx:2,:) = aimag( u2dy(1:lx,:) )

    ! longitudinal - perpendicular
    du2dr = abs(cshift(u2dr, ndel, 2) - u2dr)
    call calsum2d(mu2d2, mu2d3, mu2d4, mu2d5, mu2d6, mu2d7, mu2d8, mu2d9, mu2d10)

    ! transverse - perpendicular
    du2dr = abs(cshift(u2dr, ndel, 1) - u2dr)
    call calsum2d(mv2d2, mv2d3, mv2d4, mv2d5, mv2d6, mv2d7, mv2d8, mv2d9, mv2d10)

    ! u2dz  - transverse only
    u2dr(1:nx:2,:) =  real( u2dz(1:lx,:) )
    u2dr(2:nx:2,:) = aimag( u2dz(1:lx,:) )

    du2dr = abs(cshift(u2dr, ndel, 2) - u2dr)
    call calsum2d(mu2dz2, mu2dz3, mu2dz4, mu2dz5, mu2dz6, mu2dz7, mu2dz8, mu2dz9, mu2dz10)
    du2dr = abs(cshift(u2dr, ndel, 1) - u2dr)
    call calsum2d(mu2dz2, mu2dz3, mu2dz4, mu2dz5, mu2dz6, mu2dz7, mu2dz8, mu2dz9, mu2dz10)

    ! =================== 3D component : parallel ==============

    ! uz - longitudinal 
    ur(1:nx:2,:,:) =  real( uz(1:lx,:,:) )
    ur(2:nx:2,:,:) = aimag( uz(1:lx,:,:) )

    dur = abs(cshift(ur, ndel, 3) - ur)
    call calsum(mdudz2, mdudz3, mdudz4, mdudz5, mdudz6, mdudz7, mdudz8, mdudz9, mdudz10)

    ! ux - transverse 
    ur(1:nx:2,:,:) =  real( ux(1:lx,:,:) )
    ur(2:nx:2,:,:) = aimag( ux(1:lx,:,:) )

    dur = abs(cshift(ur, ndel, 3) - ur)
    call calsum(mdvdz2, mdvdz3, mdvdz4, mdvdz5, mdvdz6, mdvdz7, mdvdz8, mdvdz9, mdvdz10)

    ! uy - transverse 
    ur(1:nx:2,:,:) =  real( uy(1:lx,:,:) )
    ur(2:nx:2,:,:) = aimag( uy(1:lx,:,:) )

    dur = abs(cshift(ur, ndel, 3) - ur)
    call calsum(mdvdz2, mdvdz3, mdvdz4, mdvdz5, mdvdz6, mdvdz7, mdvdz8, mdvdz9, mdvdz10)


    nfile = nfile + 1

  end do
  close(30)

  ignore_me = 1. / (nx*ny*nz) / nfile 

  mu2 = mu2 * ignore_me / 2
  mu3 = mu3 * ignore_me / 2
  mu4 = mu4 * ignore_me / 2
  mu5 = mu5 * ignore_me / 2
  mu6 = mu6 * ignore_me / 2
  mu7 = mu7 * ignore_me / 2
  mu8 = mu8 * ignore_me / 2
  mu9 = mu9 * ignore_me / 2
  mu10 = mu10 * ignore_me / 2


  mv2 = mv2 * ignore_me / 2
  mv3 = mv3 * ignore_me / 2
  mv4 = mv4 * ignore_me / 2
  mv5 = mv5 * ignore_me / 2
  mv6 = mv6 * ignore_me / 2
  mv7 = mv7 * ignore_me / 2
  mv8 = mv8 * ignore_me / 2
  mv9 = mv9 * ignore_me / 2
  mv10 = mv10 * ignore_me / 2

  mduzdt2  = mduzdt2  * ignore_me / 2
  mduzdt3  = mduzdt3  * ignore_me / 2
  mduzdt4  = mduzdt4  * ignore_me / 2
  mduzdt5  = mduzdt5  * ignore_me / 2
  mduzdt6  = mduzdt6  * ignore_me / 2
  mduzdt7  = mduzdt7  * ignore_me / 2
  mduzdt8  = mduzdt8  * ignore_me / 2
  mduzdt9  = mduzdt9  * ignore_me / 2
  mduzdt10 = mduzdt10 * ignore_me / 2

  mdudz2  = mdudz2  * ignore_me 
  mdudz3  = mdudz3  * ignore_me 
  mdudz4  = mdudz4  * ignore_me 
  mdudz5  = mdudz5  * ignore_me 
  mdudz6  = mdudz6  * ignore_me 
  mdudz7  = mdudz7  * ignore_me 
  mdudz8  = mdudz8  * ignore_me 
  mdudz9  = mdudz9  * ignore_me 
  mdudz10 = mdudz10 * ignore_me 

  mdvdz2  = mdvdz2  * ignore_me  / 2
  mdvdz3  = mdvdz3  * ignore_me  / 2
  mdvdz4  = mdvdz4  * ignore_me  / 2
  mdvdz5  = mdvdz5  * ignore_me  / 2
  mdvdz6  = mdvdz6  * ignore_me  / 2
  mdvdz7  = mdvdz7  * ignore_me  / 2
  mdvdz8  = mdvdz8  * ignore_me  / 2
  mdvdz9  = mdvdz9  * ignore_me  / 2
  mdvdz10 = mdvdz10 * ignore_me  / 2

  mu2d2  = mu2d2  * ignore_me  / 2
  mu2d3  = mu2d3  * ignore_me  / 2
  mu2d4  = mu2d4  * ignore_me  / 2
  mu2d5  = mu2d5  * ignore_me  / 2
  mu2d6  = mu2d6  * ignore_me  / 2
  mu2d7  = mu2d7  * ignore_me  / 2
  mu2d8  = mu2d8  * ignore_me  / 2
  mu2d9  = mu2d9  * ignore_me  / 2
  mu2d10 = mu2d10 * ignore_me  / 2

  mv2d2  = mv2d2  * ignore_me  / 2
  mv2d3  = mv2d3  * ignore_me  / 2
  mv2d4  = mv2d4  * ignore_me  / 2
  mv2d5  = mv2d5  * ignore_me  / 2
  mv2d6  = mv2d6  * ignore_me  / 2
  mv2d7  = mv2d7  * ignore_me  / 2
  mv2d8  = mv2d8  * ignore_me  / 2
  mv2d9  = mv2d9  * ignore_me  / 2
  mv2d10 = mv2d10 * ignore_me  / 2

  mu2dz2  = mu2dz2  * ignore_me  / 2
  mu2dz3  = mu2dz3  * ignore_me  / 2
  mu2dz4  = mu2dz4  * ignore_me  / 2
  mu2dz5  = mu2dz5  * ignore_me  / 2
  mu2dz6  = mu2dz6  * ignore_me  / 2
  mu2dz7  = mu2dz7  * ignore_me  / 2
  mu2dz8  = mu2dz8  * ignore_me  / 2
  mu2dz9  = mu2dz9  * ignore_me  / 2
  mu2dz10 = mu2dz10 * ignore_me  / 2

  path='longvelincm-decmp3d-perp-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)
    write(15,'(I6, 20E12.3)') ndel, mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10
  close(15)
  
  path='transvelincm-decmp3d-perp-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)
    write(15,'(I6, 20E12.3)') ndel, mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, mv10
  close(15)
  
  path='transvelincm-decmp3d-uz-perp-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)
    write(15,'(I6, 20E12.3)') ndel, mduzdt2, mduzdt3, mduzdt4, mduzdt5, mduzdt6, &
                                    mduzdt7, mduzdt8, mduzdt9, mduzdt10
  close(15)
  
  path='longvelincm-decmp3d-para-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)
    write(15,'(I6, 20E12.3)') ndel, mdudz2, mdudz3, mdudz4, mdudz5, mdudz6, mdudz7, mdudz8, mdudz9, mdudz10
  close(15)
  
  path='transvelincm-decmp3d-para-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)
    write(15,'(I6, 20E12.3)') ndel, mdvdz2, mdvdz3, mdvdz4, mdvdz5, mdvdz6, mdvdz7, mdvdz8, mdvdz9, mdvdz10
  close(15)
  
  path='longvelincm-decmp2d-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)
    write(15,'(I6, 20E12.3)') ndel, mu2d2, mu2d3, mu2d4, mu2d5, mu2d6, mu2d7, mu2d8, mu2d9, mu2d10
  close(15)
  
  path='transvelincm-decmp2d-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)
    write(15,'(I6, 20E12.3)') ndel, mv2d2, mv2d3, mv2d4, mv2d5, mv2d6, mv2d7, mv2d8, mv2d9, mv2d10
  close(15)
  
  path='transvelincm-decmp2d-uz-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)
    write(15,'(I6, 20E12.3)') ndel, mu2dz2, mu2dz3, mu2dz4, mu2dz5, mu2dz6, &
                                    mu2dz7, mu2dz8, mu2dz9, mu2dz10
  close(15)
  
  
  deallocate(ux, uy, uz, ur, dur, tmp)
  deallocate(u2dx, u2dy, u2dz, u2dr, du2dr, tmp2d)
  
  call destroyplan3d
  call destroyplan2d

  write(*,*) 'velincm-decmp3d2d.x done.'

contains

  subroutine calsum(m2, m3, m4, m5, m6, m7, m8, m9, m10)
    implicit none

    real(dp) :: m2, m3, m4, m5, m6, m7, m8, m9, m10

    tmp = dur * dur
    m2 = m2 + sum(tmp)
    tmp = tmp * dur
    m3 = m3 + sum(tmp)
    tmp = tmp * dur
    m4 = m4 + sum(tmp)
    tmp = tmp * dur
    m5 = m5 + sum(tmp)
    tmp = tmp * dur
    m6 = m6 + sum(tmp)
    tmp = tmp * dur
    m7 = m7 + sum(tmp)
    tmp = tmp * dur
    m8 = m8 + sum(tmp)
    tmp = tmp * dur
    m9 = m9 + sum(tmp)
    tmp = tmp * dur
    m10 = m10 + sum(tmp)

  end subroutine calsum

  subroutine calsum2d(m2, m3, m4, m5, m6, m7, m8, m9, m10)
    implicit none

    real(dp) :: m2, m3, m4, m5, m6, m7, m8, m9, m10

    tmp2d = du2dr * du2dr
    m2 = m2 + sum(tmp2d)
    tmp2d = tmp2d * du2dr
    m3 = m3 + sum(tmp2d)
    tmp2d = tmp2d * du2dr
    m4 = m4 + sum(tmp2d)
    tmp2d = tmp2d * du2dr
    m5 = m5 + sum(tmp2d)
    tmp2d = tmp2d * du2dr
    m6 = m6 + sum(tmp2d)
    tmp2d = tmp2d * du2dr
    m7 = m7 + sum(tmp2d)
    tmp2d = tmp2d * du2dr
    m8 = m8 + sum(tmp2d)
    tmp2d = tmp2d * du2dr
    m9 = m9 + sum(tmp2d)
    tmp2d = tmp2d * du2dr
    m10 = m10 + sum(tmp2d)

  end subroutine calsum2d

end program velincmdecmp
