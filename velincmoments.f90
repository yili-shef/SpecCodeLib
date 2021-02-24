program velincmoments
  use mconstant
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,lx1,lx,ly,lz,nfile,ndel
  real(sp) :: ignore_me

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz
  real(sp),    allocatable, dimension(:,:,:) :: ur, dur, tmp
  character(80) :: str, flnm, path
  real(dp) :: mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10
  real(dp) :: mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, mv10
  
  nx=iargc()
  if (nx .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./velincmoments.x nx filelist ndel'
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

  allocate( ux(lx1, ly, lz), uy(lx1, ly, lz), uz(lx1, ly, lz) )
  allocate( ur(nx, ny, nz), dur(nx, ny, nz), tmp(nx, ny, nz) )

  write(*,*) 'fftwplan'
  call fftwplan3de(nx, ny, nz)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  mu2  = 0._dp; mv2  = 0._dp
  mu3  = 0._dp; mv3  = 0._dp
  mu4  = 0._dp; mv4  = 0._dp
  mu5  = 0._dp; mv5  = 0._dp
  mu6  = 0._dp; mv6  = 0._dp
  mu7  = 0._dp; mv7  = 0._dp
  mu8  = 0._dp; mv8  = 0._dp
  mu9  = 0._dp; mv9  = 0._dp
  mu10 = 0._dp; mv10 = 0._dp

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
    
    call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)
    
    ! ux
    ur(1:nx:2,:,:) =  real( ux(1:lx,:,:) )
    ur(2:nx:2,:,:) = aimag( ux(1:lx,:,:) )

    ! longitudinal
    dur = abs(cshift(ur, ndel, 1) - ur)
    call calsum(mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10)

    ! transverse
    dur = abs(cshift(ur, ndel, 2) - ur)
    call calsum(mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, mv10)
    dur = abs(cshift(ur, ndel, 3) - ur)
    call calsum(mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, mv10)

    ! uy
    ur(1:nx:2,:,:) =  real( uy(1:lx,:,:) )
    ur(2:nx:2,:,:) = aimag( uy(1:lx,:,:) )

    ! longitudinal
    dur = abs(cshift(ur, ndel, 2) - ur)
    call calsum(mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10)

    ! Transverse
    dur = abs(cshift(ur, ndel, 1) - ur)
    call calsum(mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, mv10)
    dur = abs(cshift(ur, ndel, 3) - ur)
    call calsum(mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, mv10)

    ! uz
    ur(1:nx:2,:,:) =  real( uz(1:lx,:,:) )
    ur(2:nx:2,:,:) = aimag( uz(1:lx,:,:) )

    ! longitudinal
    dur = abs(cshift(ur, ndel, 3) - ur)
    call calsum(mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10)

    ! Transverse
    dur = abs(cshift(ur, ndel, 1) - ur)
    call calsum(mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, mv10)
    dur = abs(cshift(ur, ndel, 2) - ur)
    call calsum(mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, mv10)

    nfile = nfile + 1

  end do
  close(30)

  ignore_me = 1._sp / (nx*ny*nz) / nfile 

  mu2 = mu2 * ignore_me / 3
  mu3 = mu3 * ignore_me / 3
  mu4 = mu4 * ignore_me / 3
  mu5 = mu5 * ignore_me / 3
  mu6 = mu6 * ignore_me / 3
  mu7 = mu7 * ignore_me / 3
  mu8 = mu8 * ignore_me / 3
  mu9 = mu9 * ignore_me / 3
  mu10 = mu10 * ignore_me / 3


  mv2 = mv2 * ignore_me / 6
  mv3 = mv3 * ignore_me / 6
  mv4 = mv4 * ignore_me / 6
  mv5 = mv5 * ignore_me / 6
  mv6 = mv6 * ignore_me / 6
  mv7 = mv7 * ignore_me / 6
  mv8 = mv8 * ignore_me / 6
  mv9 = mv9 * ignore_me / 6
  mv10 = mv10 * ignore_me / 6

  path='longvelincmoments-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)
    write(15,'(I6, 20E12.3)') ndel, mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10
  close(15)
  
  path='transvelincmoments-'//str(1:len_trim(str))//'dx-'//flnm(1:len_trim(flnm))//'.dat'//char(0)
  
  open(15,file=path)
    write(15,'(I6, 20E12.3)') ndel, mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, mv10
  close(15)
  
  
  deallocate(ux, uy, uz, ur, dur, tmp)
  
  call destroyplan3d

  write(*,*) 'velincmoments.x done.'

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

end program velincmoments

    !tmp = cshift(ur, ndel, 1)
    !dur = abs(tmp - ur)
    !tmp = dur * dur
    !mu2 = mu2 + sum(tmp)
    !tmp = tmp * dur
    !mu3 = mu3 + sum(tmp)
    !tmp = tmp * dur
    !mu4 = mu4 + sum(tmp)
    !tmp = tmp * dur
    !mu5 = mu5 + sum(tmp)
    !tmp = tmp * dur
    !mu6 = mu6 + sum(tmp)
    !tmp = tmp * dur
    !mu7 = mu7 + sum(tmp)
    !tmp = tmp * dur
    !mu8 = mu8 + sum(tmp)
    !tmp = tmp * dur
    !mu9 = mu9 + sum(tmp)
    !tmp = tmp * dur
    !mu10 = mu10 + sum(tmp)


    !tmp = dur * dur
    !mv2 = mv2 + sum(tmp)
    !tmp = tmp * dur
    !mv3 = mv3 + sum(tmp)
    !tmp = tmp * dur
    !mv4 = mv4 + sum(tmp)
    !tmp = tmp * dur
    !mv5 = mv5 + sum(tmp)
    !tmp = tmp * dur
    !mv6 = mv6 + sum(tmp)
    !tmp = tmp * dur
    !mv7 = mv7 + sum(tmp)
    !tmp = tmp * dur
    !mv8 = mv8 + sum(tmp)
    !tmp = tmp * dur
    !mv9 = mv9 + sum(tmp)
    !tmp = tmp * dur
    !mv10 = mv10 + sum(tmp)
