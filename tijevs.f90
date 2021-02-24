program tijevs
  use mconstant
  use mfftwplan3d
  use mwavenumber
  implicit none

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz
  complex(sp), allocatable, dimension(:,:,:) :: t11,t12,t13,t22,t23,t33
  real(sp),    allocatable, dimension(:,:,:) :: g
  real(sp),    allocatable, dimension(:)     ::  kx, ky, kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evectors
  integer :: ierr

  real(sp), allocatable, dimension(:,:,:) :: alpha, beta, gmma

  character(80) :: fnm,str,fpath,str1
  integer :: ndel,ii,jj,kk,ll
  integer :: nx,ny,nz,lx,lx1,ly,lz
  real(sp) :: delta_c

  write(*,*) 
  write(*,'(''>>>>>> Eigenvalues and eigenvectors of tij sgs stresses <<<<<<'')')
  write(*,*)

  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./tijevs-s(d)p.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist.list: the list of data file'
          write(*,*) '        ndel: filter scale delta=ndel*dx'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filter parameter
  call getarg(3,str)
  read(str, '(I20)') ndel

  ! file number string
  call getarg(2,fnm)
  fnm = adjustl(fnm)


  ny = nx; nz=nx
  lx = nx/2; lx1=lx+1
  ly = nx; lz=nx

  delta_c=ndel*2*pi/nx

  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(t11(lx1,ly,lz),t12(lx1,ly,lz),t13(lx1,ly,lz))
  allocate(t22(lx1,ly,lz),t23(lx1,ly,lz),t33(lx1,ly,lz))
  allocate(ux(lx1,ly,lz),uy(lx1,ly,lz),uz(lx1,ly,lz))
  allocate(alpha(nx,ny,nz),beta(nx,ny,nz),gmma(nx,ny,nz))
  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20,file=fnm(1:len_trim(fnm))//'.list')

    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

      fpath='./out/ux'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)ux
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)uy
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)uz
      close(10)
      write(*,*) 'after reading data files'
  
      call rsgstauij(ux,ux,t11,g,nx,ny,nz) 
      call rsgstauij(ux,uy,t12,g,nx,ny,nz)
      call rsgstauij(ux,uz,t13,g,nx,ny,nz)
      call rsgstauij(uy,uy,t22,g,nx,ny,nz)
      call rsgstauij(uy,uz,t23,g,nx,ny,nz)
      call rsgstauij(uz,uz,t33,g,nx,ny,nz)
  
      ux = - (t11 + t22 + t33)/3.
  
      t11 = t11 + ux
      t22 = t22 + ux
      t33 = t33 + ux
  
      do kk=1,lz
      do jj=1,ly
      do ii=1,lx
  
        cc(1,1) = -real(t11(ii,jj,kk),sp)
        cc(1,2) = -real(t12(ii,jj,kk),sp)
        cc(1,3) = -real(t13(ii,jj,kk),sp)
        cc(2,1) = cc(1,2)
        cc(2,2) = -real(t22(ii,jj,kk),sp)
        cc(2,3) = -real(t23(ii,jj,kk),sp)
        cc(3,1) = cc(1,3)            
        cc(3,2) = cc(2,3)
        cc(3,3) = -real(t33(ii,jj,kk),sp)
  
        call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
        alpha(2*ii-1,jj,kk)=real(evalues(3),sp)
        beta( 2*ii-1,jj,kk)=real(evalues(2),sp)
        gmma( 2*ii-1,jj,kk)=real(evalues(1),sp)
  
        cc(1,1) = -aimag(t11(ii,jj,kk))
        cc(1,2) = -aimag(t12(ii,jj,kk))
        cc(1,3) = -aimag(t13(ii,jj,kk))
        cc(2,1) = cc(1,2)
        cc(2,2) = -aimag(t22(ii,jj,kk))
        cc(2,3) = -aimag(t23(ii,jj,kk))
        cc(3,1) = cc(1,3)
        cc(3,2) = cc(2,3)
        cc(3,3) = -aimag(t33(ii,jj,kk))
  
        call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
        alpha(2*ii,jj,kk)=real(evalues(3),sp)
        beta( 2*ii,jj,kk)=real(evalues(2),sp)
        gmma( 2*ii,jj,kk)=real(evalues(1),sp)
  
      end do
      end do
      end do

      str1 ='-'//str(1:len_trim(str))//'dx-'//str1(1:len_trim(str1))
   
      open(15, file='./out/tijevs'//str1(1:len_trim(str1)),form='unformatted')
        write(15) alpha
        write(15) beta
        write(15) gmma
      close(15)

    end do
  close(20)

  deallocate(t11,t12,t13,t22,t23,t33,g)
  deallocate(alpha,beta,gmma)
  deallocate(kx,ky,kz,ux,uy,uz)

  call destroyplan3d

  write(*,*) 'finished'

end program tijevs
