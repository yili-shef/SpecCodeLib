program sijincvar
  use, intrinsic :: iso_c_binding
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,npnt,naverage
  real(sp) :: const, ignore_me
  integer, parameter :: lgap = 8

  complex(sp), allocatable, dimension(:,:,:) :: s11, s12, s13, s22, s23, s33
  real(sp),    pointer,     dimension(:,:,:) :: r11, r12 ,r13, r22, r23, r33
  real(sp), allocatable, dimension(:,:,:) :: kx, ky, kz
  real(dp), allocatable, dimension(:) :: dslll, dsltt, dsllt, dsltn
  type(C_ptr) :: cptr

  character(80) :: str, flnm, fpath
  
  if (iargc() .ne. 2) then
    write(*,*)
    write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
    write(*,*) 
    write(*,*) ' Usage: ./sij-inc-variance.x nx filelist'
    write(*,*) '                     nx: resolution of data'
    write(*,*) '                     filelist: data file list: *.list'
    write(*,*)
    write(*,*) ' Stopped'
    stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! file list 
  call getarg(2,flnm)
  flnm = adjustl(flnm)

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz

  ! The max displacement for correlation is half of the size of the cube
  npnt = nx/2
  allocate( dslll(npnt), dsltt(npnt), dsllt(npnt), dsltn(npnt) )
  allocate( s11(lx1,ly,lz), s12(lx1,ly,lz), s13(lx1,ly,lz) )
  allocate( s22(lx1,ly,lz), s23(lx1,ly,lz) )
  allocate( s33(lx1,ly,lz) )
  allocate( kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz))

  call fftwplan3de(nx,ny,nz)
  call wavenumber(kx,ky,kz,lx1,ly,lz)

  open(20, file = flnm(1 : len_trim(flnm))//'.list')

  naverage = 0
  dslll = 0._dp; dsltt = 0._dp; dsllt = 0._dp; dsltn = 0._dp
  do while ( .not. eof(20) )

    read(20,*) fpath
    write(*,*) "sij-inc-variance.x: ", fpath(1 : len_trim(fpath))

    open(15,file='./out/ux'//trim(fpath),form='unformatted')
      read(15) s11
    close(15)
    open(15,file='./out/uy'//trim(fpath),form='unformatted')
      read(15) s22
    close(15)
    open(15,file='./out/uz'//trim(fpath),form='unformatted')
      read(15) s33
    close(15)
  
    s12=.5_sp*eye*(kx*s22+ky*s11)
    s13=.5_sp*eye*(kx*s33+kz*s11)
    s23=.5_sp*eye*(ky*s33+kz*s22)
    s11=eye*kx*s11
    s22=eye*ky*s22
    s33=eye*kz*s33
 
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)

    cptr = c_loc(s11)
    call c_f_pointer(cptr, r11, [2*lx1,ly,lz])
    cptr = c_loc(s22)
    call c_f_pointer(cptr, r22, [2*lx1,ly,lz])
    cptr = c_loc(s33)
    call c_f_pointer(cptr, r33, [2*lx1,ly,lz])
    cptr = c_loc(s12)
    call c_f_pointer(cptr, r12, [2*lx1,ly,lz])
    cptr = c_loc(s13)
    call c_f_pointer(cptr, r13, [2*lx1,ly,lz])
    cptr = c_loc(s23)
    call c_f_pointer(cptr, r23, [2*lx1,ly,lz])

    do jj = 1, ny, lgap
    do ii = 1, nx, lgap
      do kk = 1, npnt
        dslll(kk) = dslll(kk) + ( r33(ii,jj,kk) - r33(ii,jj,1) )**2
        dsltt(kk) = dsltt(kk) + ( r22(ii,jj,kk) - r22(ii,jj,1) )**2
        dsltt(kk) = dsltt(kk) + ( r11(ii,jj,kk) - r11(ii,jj,1) )**2
        dsllt(kk) = dsllt(kk) + ( r23(ii,jj,kk) - r23(ii,jj,1) )**2
        dsllt(kk) = dsllt(kk) + ( r13(ii,jj,kk) - r13(ii,jj,1) )**2
        dsltn(kk) = dsltn(kk) + ( r12(ii,jj,kk) - r12(ii,jj,1) )**2
      end do
      naverage = naverage + 1
    end do
    end do

    do kk = 1, nz, lgap
    do ii = 1, nx, lgap
      do jj = 1, npnt
        dslll(jj) = dslll(jj) + ( r22(ii,jj,kk) - r22(ii,1,kk) )**2
        dsltt(jj) = dsltt(jj) + ( r33(ii,jj,kk) - r33(ii,1,kk) )**2
        dsltt(jj) = dsltt(jj) + ( r11(ii,jj,kk) - r11(ii,1,kk) )**2
        dsllt(jj) = dsllt(jj) + ( r23(ii,jj,kk) - r23(ii,1,kk) )**2
        dsllt(jj) = dsllt(jj) + ( r12(ii,jj,kk) - r12(ii,1,kk) )**2
        dsltn(jj) = dsltn(jj) + ( r13(ii,jj,kk) - r13(ii,1,kk) )**2
      end do
      naverage = naverage + 1
    end do
    end do


    do kk = 1, nz, lgap
    do jj = 1, ny, lgap
      do ii = 1, npnt
        dslll(ii) = dslll(ii) + ( r11(ii,jj,kk) - r11(1,jj,kk) )**2
        dsltt(ii) = dsltt(ii) + ( r22(ii,jj,kk) - r22(1,jj,kk) )**2
        dsltt(ii) = dsltt(ii) + ( r33(ii,jj,kk) - r33(1,jj,kk) )**2
        dsllt(ii) = dsllt(ii) + ( r12(ii,jj,kk) - r12(1,jj,kk) )**2
        dsllt(ii) = dsllt(ii) + ( r13(ii,jj,kk) - r13(1,jj,kk) )**2
        dsltn(ii) = dsltn(ii) + ( r23(ii,jj,kk) - r23(1,jj,kk) )**2
      end do
      naverage = naverage + 1
    end do
    end do


  end do
  close(20)

  const = 1._dp/naverage
  dslll = dslll * const
  dsltt = dsltt * const / 2 
  dsllt = dsllt * const / 2
  dsltn = dsltn * const

  open(26, file = 'sij-inc-variance-'//trim(flnm)//'.dat')

    write(26,*) 'dx,dslll,dsltt,dsllt,dsltn'
    do ll = 1, npnt
      write(26, '(15E15.3)') real(ll), dslll(ll), dsltt(ll), dsllt(ll), dsltn(ll)
    end do

  close(26)

  call destroyplan3d

  deallocate(s11,s12,s13,s22,s23,s33,kx,ky,kz)
  deallocate(dslll,dsltt,dsllt,dsltn)
  
  write(*,*) 'sij-inc-variance.x done.'

end program sijincvar
