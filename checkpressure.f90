program checkpressure
  use mconstant
  use mfftwplan3d

  integer, parameter :: nx=128, ny=128, nz=128
  integer, parameter :: lx=nx/2, ly=ny, lz=nz, lx1=lx+1

  complex(sp), dimension(lx1,ly,lz) :: p
  real(sp) :: ignore_me
  integer :: ii, jj, kk, ll

  integer, parameter :: npnt = 100
  real(sp), dimension(npnt) :: pdfp = 0.
  real(sp), parameter :: rmsp=1., lowerbnd=-2, upperbnd=2
  real(sp), parameter :: binw=(upperbnd-lowerbnd)/npnt
  real(sp) :: tmp

  call fftwplan3d(nx,ny,nz)

  open(16, file='./out/p0001.dat',form='unformatted')
    read(16) p
  close(16)

  call rfftwnd_f77_one_complex_to_real(C2R3D,p,ignore_me)

  
  open(16, file='./post/p2d.dat')
  write(16,'(''variables = "x", "y", "p"'')')
  write(16, '(''zone T="p", I='', I6, '', J='', I6)') nx, ny
  do jj=1,ly
    do ii=1,lx
      write(16,*) jj, 2*ii-1, real(p(ii,jj, 20))
      write(16,*) jj, 2*ii, aimag(p(ii,jj,20))
    end do
  end do    
  close(16)


  do ii=1,lx
    do jj=1,ly
      do kk=1,lz
        tmp=real(p(ii,jj,kk))/rmsp
        if( tmp .gt. lowerbnd .and. tmp .lt. upperbnd) then
          ll=FLOOR((tmp-lowerbnd)/binw)+1
          pdfp(ll)=pdfp(ll)+1
        end if
        tmp=aimag(p(ii,jj,kk))/rmsp
        if( tmp .gt. lowerbnd .and. tmp .lt. upperbnd) then
          ll=FLOOR((tmp-lowerbnd)/binw)+1
          pdfp(ll)=pdfp(ll)+1
        end if
      end do
    end do 
  end do
  pdfp=pdfp/real(nx*ny*nz)
  write(*,*) 'check pdfp:', sum(pdfp)
  pdfp=pdfp/binw

  open(16, file='./post/pdfp.dat')
    write(16, '(''variables="p", "pdf of p"'')')
    write(16, '(''zone T="pdf of p", I='', I6)') npnt
    do ii=1,npnt
      write(16, *) lowerbnd-binw/2+ii*binw, pdfp(ii)
    end do
  close(16)
  
  call destroyplan3d

end program checkpressure      
