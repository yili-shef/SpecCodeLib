program sijevs
  use mconstant
  use mfftwplan3d
  use mwavenumber
  implicit none

  complex(sp), allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33
  real(sp), allocatable, dimension(:,:,:) :: g
  real(sp), allocatable, dimension(:) :: kx,ky,kz

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: evalues, fv1,fv2
  real(dp), dimension(3,3) :: cc,evectors
  integer :: ierr

  real(sp), allocatable, dimension(:,:,:) :: alpha, beta, gmma

  character(80) :: fnm,str,fpath, str1
  integer :: ndel,ii,jj,kk,ll
  integer :: nx,ny,nz,lx,lx1,ly,lz
  real(sp) :: delta_c, ignore_me

  write(*,*) 
  write(*,'(''>>>>>> Eigenvalues of strain rate <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./sijevs.x nx filelist ndel'
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

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(s11(lx1,ly,lz),s12(lx1,ly,lz),s13(lx1,ly,lz))
  allocate(s22(lx1,ly,lz),s23(lx1,ly,lz),s33(lx1,ly,lz))
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
        read(10)s11
      close(10)
      fpath='./out/uy'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)s22
      close(10)
      fpath='./out/uz'//str1(1:len_trim(str1))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)s33
      close(10)
      write(*,*) 'after reading data files'
  
      ! Gaussian filter
      g=exp(-g*delta_c**2/24.)
      s11=s11*g
      s22=s22*g 
      s33=s33*g
  
      do kk = 1, lz
      do jj = 1, ly
      do ii = 1, lx1
      s12(ii,jj,kk)=.5*eye*(kx(ii)*s22(ii,jj,kk)+ky(jj)*s11(ii,jj,kk))
      s13(ii,jj,kk)=.5*eye*(kx(ii)*s33(ii,jj,kk)+kz(kk)*s11(ii,jj,kk))
      s23(ii,jj,kk)=.5*eye*(ky(jj)*s33(ii,jj,kk)+kz(kk)*s22(ii,jj,kk))
      s11(ii,jj,kk)=eye*kx(ii)*s11(ii,jj,kk)
      s22(ii,jj,kk)=eye*ky(jj)*s22(ii,jj,kk)
      s33(ii,jj,kk)=eye*kz(kk)*s33(ii,jj,kk)
      end do
      end do
      end do
  
  
      call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
  
      do kk=1,lz
      do jj=1,ly
      do ii=1,lx
  
        cc(1,1)=real(s11(ii,jj,kk))
        cc(1,2)=real(s12(ii,jj,kk))
        cc(1,3)=real(s13(ii,jj,kk))
        cc(2,1)=cc(1,2)
        cc(2,2)=real(s22(ii,jj,kk))
        cc(2,3)=real(s23(ii,jj,kk))
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
        cc(3,3)=real(s33(ii,jj,kk))
  
        call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
        alpha(2*ii-1,jj,kk)=real(evalues(3),sp)
        beta( 2*ii-1,jj,kk)=real(evalues(2),sp)
        gmma( 2*ii-1,jj,kk)=real(evalues(1),sp)
  
        cc(1,1)=aimag(s11(ii,jj,kk))
        cc(1,2)=aimag(s12(ii,jj,kk))
        cc(1,3)=aimag(s13(ii,jj,kk))
        cc(2,1)=cc(1,2)
        cc(2,2)=aimag(s22(ii,jj,kk))
        cc(2,3)=aimag(s23(ii,jj,kk))
        cc(3,1)=cc(1,3)
        cc(3,2)=cc(2,3)
        cc(3,3)=aimag(s33(ii,jj,kk))
  
        call rs(3,3,cc,evalues,matz,evectors,fv1,fv2,ierr)
  
        alpha(2*ii,jj,kk)=real(evalues(3),sp)
        beta( 2*ii,jj,kk)=real(evalues(2),sp)
        gmma( 2*ii,jj,kk)=real(evalues(1),sp)
  
      end do
      end do
      end do
  
      str1 ='-'//str(1:len_trim(str))//'dx-'//str1(1:len_trim(str1))
  
      open(15, file='./out/sijevs'//str1(1:len_trim(str1)),form='unformatted')
        write(15) alpha
        write(15) beta
        write(15) gmma
      close(15)
    end do

  close(20)

  deallocate(s11,s22,s33,s12,s13,s23,g,kx,ky,kz,alpha,beta,gmma)

  call destroyplan3d

  write(*,*) 'done'

end program sijevs
