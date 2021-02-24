program rijevs
  use mconstant
  use mfftwplan3d
  implicit none


  complex(sp), allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  real(sp), allocatable, dimension(:,:,:) :: g
  real(sp), allocatable, dimension(:) :: kx,ky,kz

  integer, parameter :: matz=5 
  integer :: ierr

  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evectors

  real(sp), allocatable, dimension(:,:,:) :: alpha, beta, gmma

  character(80) :: fnm,str,fpath, str1
  integer :: ndel,ii,jj,kk,ll
  integer :: nx,ny,nz,lx,lx1,ly,lz
  real(sp) :: delta_c, ignore_me

  write(*,*) 
  write(*,'(''>>>>>> Eigenvalues and eigenvectors of rij of vorticity <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rijevs-s(d)p.x nx filelist ndel'
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


  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  delta_c=ndel*2*pi/nx

  call fftwplan3d(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(r11(lx1,ly,lz),r12(lx1,ly,lz),r13(lx1,ly,lz))
  allocate(r22(lx1,ly,lz),r23(lx1,ly,lz),r33(lx1,ly,lz))
  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(alpha(nx,ny,nz),beta(nx,ny,nz),gmma(nx,ny,nz))
  allocate(g(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'
   
  ! Gaussian filter
  g = exp(-g*delta_c**2/24.)

  open(20,file=fnm(1:len_trim(fnm))//'.list')

    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))

      fpath='./out/ux'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)r12
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)r13
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)r23
      close(10)
      write(*,*) 'after reading data files'
   
      r12=r12*g
      r13=r13*g
      r23=r23*g
   
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
        r11(ii,jj,kk)=eye*(ky(jj)*r23(ii,jj,kk)-kz(kk)*r13(ii,jj,kk))
        r22(ii,jj,kk)=eye*(kz(kk)*r12(ii,jj,kk)-kx(ii)*r23(ii,jj,kk))
        r33(ii,jj,kk)=eye*(kx(ii)*r13(ii,jj,kk)-ky(jj)*r12(ii,jj,kk))
    
        r12(ii,jj,kk)=.5*eye*(kx(ii)*r22(ii,jj,kk)+ky(jj)*r11(ii,jj,kk))
        r13(ii,jj,kk)=.5*eye*(kx(ii)*r33(ii,jj,kk)+kz(kk)*r11(ii,jj,kk))
        r23(ii,jj,kk)=.5*eye*(ky(jj)*r33(ii,jj,kk)+kz(kk)*r22(ii,jj,kk))
        r11(ii,jj,kk)=eye*kx(ii)*r11(ii,jj,kk)
        r22(ii,jj,kk)=eye*ky(jj)*r22(ii,jj,kk)
        r33(ii,jj,kk)=eye*kz(kk)*r33(ii,jj,kk)
      end do
      end do
      end do
   
   
      call rfftwnd_f77_one_complex_to_real(c2r3d,r11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,r33,ignore_me)
   
      do kk=1,lz
      do jj=1,ly
      do ii=1,lx
   
        cc(1,1)=real(r11(ii,jj,kk),sp)
        cc(1,2)=real(r12(ii,jj,kk),sp)
        cc(1,3)=real(r13(ii,jj,kk),sp)
        cc(2,1)=cc(1,2)
        cc(2,2)=real(r22(ii,jj,kk),sp)
        cc(2,3)=real(r23(ii,jj,kk),sp)
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
        cc(3,3)=real(r33(ii,jj,kk),sp)
   
        call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
        alpha(2*ii-1,jj,kk)=real(evalues(3),sp)
        beta( 2*ii-1,jj,kk)=real(evalues(2),sp)
        gmma( 2*ii-1,jj,kk)=real(evalues(1),sp)
   
        cc(1,1)=aimag(r11(ii,jj,kk))
        cc(1,2)=aimag(r12(ii,jj,kk))
        cc(1,3)=aimag(r13(ii,jj,kk))
        cc(2,1)=cc(1,2)               
        cc(2,2)=aimag(r22(ii,jj,kk))
        cc(2,3)=aimag(r23(ii,jj,kk))
        cc(3,1)=cc(1,3)               
        cc(3,2)=cc(2,3)               
        cc(3,3)=aimag(r33(ii,jj,kk))
   
        call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
        alpha(2*ii,jj,kk)=real(evalues(3),sp)
        beta( 2*ii,jj,kk)=real(evalues(2),sp)
        gmma( 2*ii,jj,kk)=real(evalues(1),sp)
   
      end do
      end do
      end do
   
      str1 ='-'//str(1:len_trim(str))//'dx-'//str1(1:len_trim(str1))
   
      open(15, file='./out/rijevs'//str1(1:len_trim(str1)),form='unformatted')
        write(15) alpha
        write(15) beta
        write(15) gmma
      close(15)

    end do

  close(20)

  deallocate(r11,r12,r13,r22,r23,r33,kx,ky,kz,g)
  deallocate(alpha,beta,gmma)

  call destroyplan3d

  write(*,*) 
  write(*,*) 'done!'

end program rijevs
