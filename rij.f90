! Calculate the velocity correlation tensor rij=<ui(x)uj(x+r)>
program rij
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: lx1,lx,ly,lz,nx,ny,nz,nn,mm,ii,jj,kk
  complex(sp), allocatable, dimension(:,:,:) :: ux,uy,uz
  real(sp), allocatable, dimension(:) :: rxxx,ryyx,rzzx,rxxy,ryyy,rzzy,rxxz,ryyz,rzzz
  real(sp), allocatable, dimension(:,:,:) :: uxr,uyr,uzr
  real(sp) :: const, dx, ignore_me
  character(80) :: str

  write(*,*) 
  write(*,*) '>>> Calculate the correlation functions and integral lengths <<<'
  write(*,*) '>>>         This will take a while...<<<'
  write(*,*) 

  nx=iargc()
  if (nx .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rij.x nx nfile'
          write(*,*) '          nx: resolution of data'
          write(*,*) '          nfile: data file number'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx
  ! file number string
  call getarg(2,str)
  str = adjustl(str)


  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  dx=2*pi/nx

  const=1._sp/(nx*ny*nz)

  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(uxr(nx,ny,nz),uyr(nx,ny,nz),uzr(nx,ny,nz))
  allocate(rxxx(lx),ryyx(lx),rzzx(lx))
  allocate(rxxy(lx),ryyy(lx),rzzy(lx))
  allocate(rxxz(lx),ryyz(lx),rzzz(lx))


  open(15,file='./out/ux'//str(1:len_trim(str))//'.dat',form='unformatted')
  read(15) ux
  close(15)
  open(15,file='./out/uy'//str(1:len_trim(str))//'.dat',form='unformatted')
  read(15) uy
  close(15)
  open(15,file='./out/uz'//str(1:len_trim(str))//'.dat',form='unformatted')
  read(15) uz
  close(15)

  write(*,*) 'fftwplan'
  call fftwplan3de(nx,ny,nz)


  call rfftwnd_f77_one_complex_to_real(c2r3d,ux,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,uy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,uz,ignore_me)

  uxr(1:nx:2,:,:) =  real(ux(1:lx,:,:))
  uxr(2:nx:2,:,:) = aimag(ux(1:lx,:,:))
  uyr(1:nx:2,:,:) =  real(uy(1:lx,:,:))
  uyr(2:nx:2,:,:) = aimag(uy(1:lx,:,:))
  uzr(1:nx:2,:,:) =  real(uz(1:lx,:,:))
  uzr(2:nx:2,:,:) = aimag(uz(1:lx,:,:))


  rxxz = 0.; ryyz = 0.; rzzz = 0.
  rxxy = 0.; ryyy = 0.; rzzy = 0.
  rxxx = 0.; ryyx = 0.; rzzx = 0.
  do nn = 1, lx ! Correlation defined only for half of the domain

    do ii = 1, nx

      mm=modulo(ii+nn-1-1,nx)+1
     
      rxxx(nn) = rxxx(nn) + sum(uxr(ii,:,:)*uxr(mm,:,:))
      ryyx(nn) = ryyx(nn) + sum(uyr(ii,:,:)*uyr(mm,:,:))
      rzzx(nn) = rzzx(nn) + sum(uzr(ii,:,:)*uzr(mm,:,:))
    end do

    do jj = 1, ny
     
      mm=modulo(jj+nn-1-1,ny)+1

      rxxy(nn) = rxxy(nn) + sum(uxr(:,jj,:)*uxr(:,mm,:))
      ryyy(nn) = ryyy(nn) + sum(uyr(:,jj,:)*uyr(:,mm,:))
      rzzy(nn) = rzzy(nn) + sum(uzr(:,jj,:)*uzr(:,mm,:))
    end do

    do kk = 1, nz

      mm=modulo(kk+nn-1-1,nz)+1

      rxxz(nn) = rxxz(nn) + sum(uxr(:,:,kk)*uxr(:,:,mm)) 
      ryyz(nn) = ryyz(nn) + sum(uyr(:,:,kk)*uyr(:,:,mm)) 
      rzzz(nn) = rzzz(nn) + sum(uzr(:,:,kk)*uzr(:,:,mm)) 
     
    end do

    rxxx(nn)=rxxx(nn)*const
    ryyx(nn)=ryyx(nn)*const
    rzzx(nn)=rzzx(nn)*const

    rxxy(nn)=rxxy(nn)*const
    ryyy(nn)=ryyy(nn)*const
    rzzy(nn)=rzzy(nn)*const

    rxxz(nn)=rxxz(nn)*const
    ryyz(nn)=ryyz(nn)*const
    rzzz(nn)=rzzz(nn)*const
  end do

  open(15,file='rij'//str(1:len_trim(str))//'.dat')
    write(15, *) '# rxxx(1) ryyx(1) rzzx(1)'
    write(15, *) '#', rxxx(1), ryyx(1), rzzx(1)
    write(15, *) '# rxxy(1) ryyy(1) rzzy(1)'
    write(15, *) '#', rxxy(1), ryyy(1), rzzy(1)
    write(15, *) '# rxxz(1) ryyz(1) rzzz(1)'
    write(15, *) '#', rxxz(1), ryyz(1), rzzz(1)
    write(15, *) '# Lxxx    Lyyx     Lzzx'
    write(15, *) '#', dx*sum(rxxx), dx*sum(ryyx), dx*sum(rzzx)
    write(15, *) '# Lxxy    Lyyy     Lzzy'
    write(15, *) '#', dx*sum(rxxy), dx*sum(ryyy), dx*sum(rzzy)
    write(15, *) '# Lxxz    Lyyz     Lzzz'
    write(15, *) '#', dx*sum(rxxz), dx*sum(ryyz), dx*sum(rzzz)
    do ii=1,lx
      write(15,'(20f9.4)') (ii-1)*dx,rxxx(ii)/rxxx(1),ryyx(ii)/ryyx(1), &
                   rzzx(ii)/rzzx(1),rxxy(ii)/rxxy(1),ryyy(ii)/ryyy(1), &
                   rzzy(ii)/rzzy(1),rxxz(ii)/rxxz(1),ryyz(ii)/ryyz(1), &
                   rzzz(ii)/rzzz(1)
    end do
  close(15)

  deallocate(rxxx,ryyx,rzzx,rxxy,ryyy,rzzy,rxxz,ryyz,rzzz)
  deallocate(uxr,uyr,uzr, ux, uy, uz)

  call destroyplan3d

  write(*,*) 'done'
end program rij
