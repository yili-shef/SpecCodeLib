program cndpdforo
  use mconstant
  use mfftwplan3d
  implicit none

  integer :: nx,ny,nz,lx,lx1,ly,lz
  integer :: ndel,ii,jj,kk,ll,mm,nn, nfile

  complex, allocatable, dimension(:,:,:) :: wx,wy,wz
  complex, allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33
  real, allocatable, dimension(:,:,:) :: g, tmp1, tmp2
  real,    allocatable, dimension(:)     :: kx,ky,kz

  integer, parameter :: npnt=160
  real(sp),    parameter :: bound = 10., binw=2.*bound/npnt
  real(dp), dimension(npnt) :: pdfp, pdfm

  real(sp), dimension(3,3) :: rij
  real(dp) :: mean, rms, skew
  real :: delta_c, ignore_me, tmp, const, mean0, rms0, i3r
  character(80) :: fnm,str,fpath
  integer(8) :: nump, numm

  write(*,*) 
  write(*,'(''>>>>>> Conditional PDFs of the Rij*oi*oj (in helical turb) <<<<<<'')')
  write(*,*)
  ll=iargc()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdforo-cndi3r.x nx filelist ndel'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: dns data file list'
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
  ! file list string
  call getarg(2,fnm)
  fnm = adjustl(fnm)


  ny=nx; nz=nx
  lx=nx/2; lx1=lx+1
  ly=nx; lz=nx

  const = 1./(nx*ny*nz)
  delta_c=ndel*2*pi/nx

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan3d'

  allocate(wx(lx1,ly,lz),wy(lx1,ly,lz),wz(lx1,ly,lz))
  allocate(kx(lx1),ky(ly),kz(lz))
  allocate(r11(lx1,ly,lz),r12(lx1,ly,lz),r13(lx1,ly,lz))
  allocate(r22(lx1,ly,lz),r23(lx1,ly,lz),r33(lx1,ly,lz))
  allocate(g(lx1,ly,lz),tmp1(lx1,ly,lz),tmp2(lx1,ly,lz))
  write(*,*) 'arrays allocated'

  call wavenumber(kx,ky,kz,g,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  ! Gaussian filter
  g=exp(-g*delta_c**2/24.)

  open(20, file = fnm(1 : len_trim(fnm))//'.list')

  pdfp = 0._dp
  pdfm = 0._dp
  mean = 0._dp
  rms = 0._dp
  skew=0._dp
  nfile = 0
  nump = 0
  numm = 0
  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) fpath(1 : len_trim(fpath))

    open(10,file='./out/ux'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)r11
    close(10)
    open(10,file='./out/uy'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)r22
    close(10)
    open(10,file='./out/uz'//fpath(1 : len_trim(fpath)),form='unformatted')
      read(10)r33
    close(10)
    write(*,*) 'after reading data files'
 
    r11=r11*g
    r22=r22*g 
    r33=r33*g
 
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      wx(ii,jj,kk)=eye*(ky(jj)*r33(ii,jj,kk)-kz(kk)*r22(ii,jj,kk))
      wy(ii,jj,kk)=eye*(kz(kk)*r11(ii,jj,kk)-kx(ii)*r33(ii,jj,kk))
      wz(ii,jj,kk)=eye*(kx(ii)*r22(ii,jj,kk)-ky(jj)*r11(ii,jj,kk))
      r11(ii,jj,kk)=eye*kx(ii)*wx(ii,jj,kk)
      r22(ii,jj,kk)=eye*ky(jj)*wy(ii,jj,kk)
      r33(ii,jj,kk)=eye*kz(kk)*wz(ii,jj,kk)
      r12(ii,jj,kk)=.5*eye*(kx(ii)*wy(ii,jj,kk)+ky(jj)*wx(ii,jj,kk))
      r13(ii,jj,kk)=.5*eye*(kx(ii)*wz(ii,jj,kk)+kz(kk)*wx(ii,jj,kk))
      r23(ii,jj,kk)=.5*eye*(ky(jj)*wz(ii,jj,kk)+kz(kk)*wy(ii,jj,kk))
    end do
    end do
    end do

 
    call rfftwnd_f77_one_complex_to_real(c2r3d,wx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,wz,ignore_me)
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,r11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,r33,ignore_me)


    tmp1 = real(wx) * real(r11) * real(wx) + real(wy) * real(r22) * real(wy) + &
          real(wz) * real(r33) * real(wz) + 2. * ( real(wx) * real(r12) * real(wy) &
          + real(wx) * real(r13) * real(wz) + real(wy) * real(r23)* real(wz) ) 

    tmp2 = aimag(wx) * aimag(r11) * aimag(wx) + aimag(wy) * aimag(r22) * aimag(wy) + &
          aimag(wz) * aimag(r33) * aimag(wz) + 2. * ( aimag(wx) * aimag(r12) * aimag(wy) &
          + aimag(wx) * aimag(r13) * aimag(wz) + aimag(wy) * aimag(r23)* aimag(wz) ) 

    wx = cmplx( real(wx) * real(wx) + real(wy) * real(wy) + real(wz) * real(wz), &
                aimag(wx) * aimag(wx) + aimag(wy) * aimag(wy) + aimag(wz) * aimag(wz) )
    
    tmp1 = tmp1 / real(wx)
    tmp2 = tmp2 / aimag(wx) 


    if ( nfile .eq. 0 ) then

      ignore_me = sum( tmp1(1:lx,:,:) ) + sum( tmp2(1:lx,:,:) )
      delta_c = sum( tmp1(1:lx,:,:)**2 )+ sum( tmp2(1:lx,:,:)**2 )

      mean0 = ignore_me * const
      delta_c = delta_c * const
      rms0 = sqrt( delta_c - mean0 * mean0)

      write(*,*) 'Estimated mean: ', mean0
      write(*,*) 'Estimated rms: ', rms0
    end if
 
    do kk=1,nz
    do jj=1,ny
    do ii=1,nx

      if ( mod(ii,2) .eq. 1) then
        ll = (ii + 1) / 2

        rij(1,1) = real(r11(ll,jj,kk))
        rij(1,2) = real(r12(ll,jj,kk))
        rij(1,3) = real(r13(ll,jj,kk))
        rij(2,2) = real(r22(ll,jj,kk))
        rij(2,3) = real(r23(ll,jj,kk))
        rij(3,3) = real(r33(ll,jj,kk))

        ignore_me = tmp1(ll,jj,kk) 
      else
        ll = ii / 2

        rij(1,1) = aimag(r11(ll,jj,kk))
        rij(1,2) = aimag(r12(ll,jj,kk))
        rij(1,3) = aimag(r13(ll,jj,kk))
        rij(2,2) = aimag(r22(ll,jj,kk))
        rij(2,3) = aimag(r23(ll,jj,kk))
        rij(3,3) = aimag(r33(ll,jj,kk))

        ignore_me = tmp2(ll,jj,kk) 
      end if
      rij(2,1) = rij(1,2)
      rij(3,1) = rij(1,3)
      rij(3,2) = rij(2,3)

      i3r = - sum( matmul(rij, rij) * rij ) / 3.

      mean = mean + ignore_me 
      rms = rms + ignore_me*ignore_me 
      skew = skew + ignore_me**3

      ignore_me = ignore_me / rms0 + bound
      ll = floor( ignore_me / binw ) + 1

      
      if ( ll .le. npnt .and. ll .ge. 1 ) then 
        if (i3r .ge. 0 ) then 
          pdfp(ll) = pdfp(ll) + 1
          nump = nump + 1
        else
          pdfm(ll) = pdfm(ll) + 1
          numm = numm + 1
        end if
      end if
 
    end do
    end do
    end do

    nfile = nfile + 1
  end do

  close(20)

  tmp = 1._dp/(nx*ny*nz)/nfile

  mean = mean * tmp
  rms = rms * tmp
  skew = skew * tmp

  skew = skew - 3 * rms * mean + 2 * mean * mean * mean

  rms = sqrt( rms - mean*mean) 
  skew = skew / rms**3
  pdfp = pdfp / nump
  pdfm = pdfm / numm

  write(*,*) 'mean is:', mean
  write(*,*) 'rms is:', rms
  write(*,*) 'skew is:', skew
  write(*,*) 'nump = ', nump
  write(*,*) 'numm = ', numm
  write(*,*) 'check pdfp:   ', sum(pdfp)
  write(*,*) 'check pdfm:   ', sum(pdfm)

  pdfp = pdfp / binw
  pdfm = pdfm / binw


  fnm = '-'//str(1:len_trim(str))//'dx-'//fnm(1:len_trim(fnm))//'.dat'

  open(15, file='pdforo-cndi3rp'//fnm(1:len_trim(fnm)))
    do ii=1,npnt
      write(15,'(15E15.6)') (-bound+(ii-.5)*binw)*rms0/rms, pdfp(ii)*rms/rms0
    end do
  close(15)
  open(15, file='pdforo-cndi3rm'//fnm(1:len_trim(fnm)))
    do ii=1,npnt
      write(15,'(15E15.6)') (-bound+(ii-.5)*binw)*rms0/rms, pdfm(ii)*rms/rms0
    end do
  close(15)

  deallocate(wx,wy,wz,g,kx,ky,kz,tmp1,tmp2)
  deallocate(r11,r12,r13,r22,r23,r33)

  call destroyplan3d

  write(*,*) 'done'

end program cndpdforo
