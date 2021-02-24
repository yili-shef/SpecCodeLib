program defmcijcorbw
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer, parameter :: norder = 6
  real(sp) :: y(3), dxyz(3), bg(norder,3)
  integer  :: ix(norder), iy(norder), iz(norder), lhnode(3)
  integer  :: ix1, iy1, iz1
 
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,pp,nfile
  real(sp) :: const, ignore_me

  integer, parameter :: ngroups = 512
  integer :: npntgrp
  integer, allocatable, dimension(:,:) :: xx, yy, zz

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  ! ----------------------------------------------------------------------------------

  real(dp), dimension(evsize, evsize)  :: cc, sij, evtrcc, evtrcc0
  real(dp), dimension(evsize) :: evcc, evcc0, omp, ominc, rotinc, rotinc0, ominc0

  real(dp) :: mcntr1, mcntr2, mcntr3, mcntr4, align, cntr1, cntr2, cntr3, cntr4
  real(dp) :: mcntr1a, mcntr2a, mcntr3a, mcntr4a, cntr1a, cntr2a, cntr3a, cntr4a

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33,omx,omy,omz
  real(sp),    allocatable, dimension(:,:,:) :: r11,r12,r13,r22,r23,r33,rox,roy,roz
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz
  real(sp),    allocatable, dimension(:)     :: xp, yp, zp

  character(80) :: str, flnm, prf, str1, xyzlst, mark
  integer :: nlst, forward
  
  if (iargc() .ne. 6) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-cor-rot-bw.x nx filelist prefix inicoors npntgrp fw'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     prefix: prefix for datafiles for bij'
          write(*,*) '                     inicoors: data file for initial coordinates'
          write(*,*) '                     npntgrp: # of points in each group'
          write(*,*) '                     fw: 1 for forward, 0 for backward'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! file list 
  call getarg(2,flnm)

  ! prefix
  call getarg(3,prf)

  call getarg(5, str)
  read(str, '(I20)') npntgrp

  call getarg(6, str)
  read(str, '(I20)') forward

  ! initial data file
  call getarg(4, str)


  ny=nx; nz=nx; lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz

  !--- generate the prefix for the xyz coordinate data files --------
  xyzlst = flnm( len_trim(flnm) - 3 : len_trim(flnm) ) ! extract the starting file number
  read(xyzlst, "(I20)") nlst
  if (forward .eq. 0) then 
      nlst = nlst - 150 
      mark = 'bw'
  else
      nlst = nlst + 150
      mark = 'fw'
  end if
  write(str1, "(I20)") nlst ! get the last file number
  str1 = adjustl(str1)
  xyzlst = trim(mark)//trim(xyzlst)//"to"//trim(str1)
  !--------------------------------------------------------------


  dxyz(1) = 2 * pi / nx; dxyz(2) = 2 * pi / ny; dxyz(3) = 2 * pi / nz
  const = 1._sp / ((npntgrp-1) * ngroups) ! excluding the center of each group

  allocate( xx(npntgrp,ngroups), yy(npntgrp,ngroups), zz(npntgrp,ngroups) )

  allocate( b11(lx1,ly,lz), b12(lx1,ly,lz), b13(lx1,ly,lz) )
  allocate( b21(lx1,ly,lz), b22(lx1,ly,lz), b23(lx1,ly,lz) )
  allocate( b31(lx1,ly,lz), b32(lx1,ly,lz), b33(lx1,ly,lz) )
  
  allocate( s11(lx1,ly,lz), s12(lx1,ly,lz), s13(lx1,ly,lz) )
  allocate( omx(lx1,ly,lz), s22(lx1,ly,lz), s23(lx1,ly,lz) )
  allocate( omy(lx1,ly,lz), omz(lx1,ly,lz), s33(lx1,ly,lz) )

  allocate( r11(nx,ny,nz), r12(nx,ny,nz), r13(nx,ny,nz) )
  allocate( rox(nx,ny,nz), r22(nx,ny,nz), r23(nx,ny,nz) )
  allocate( roy(nx,ny,nz), roz(nx,ny,nz), r33(nx,ny,nz) )

  allocate( kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz) )
  allocate( xp(nx*ny*nz), yp(nx*ny*nz), zp(nx*ny*nz) )

  call fftwplan3de(nx,ny,nz)
  !write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  !write(*,*) 'after wavenumber'

  open(30, file = './out/'//trim(str)//'.dat', form = 'unformatted')
    read(30) xx, yy, zz
  close(30)

  open(27, file = 'defm-rotcorel-'//trim(str)//'-'//trim(flnm)//'.dat')
  open(30, file = trim(flnm)//'.list')

  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    !write(*,*) 'defm-cor-rot-bw.x: data file ', str1(1:len_trim(str1))

    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b11
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b22
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b33
    close(15)

    s11 = eye * kx * b11
    s22 = eye * ky * b22
    s33 = eye * kz * b33
    s12 = eye * (kx * b22 + ky * b11)*0.5_sp
    s13 = eye * (kx * b33 + kz * b11)*0.5_sp
    s23 = eye * (ky * b33 + kz * b22)*0.5_sp
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)

    omx = eye * (ky * b33 - kz * b22)
    omy = eye * (kz * b11 - kx * b33)
    omz = eye * (kx * b22 - ky * b11)
    call rfftwnd_f77_one_complex_to_real(c2r3d,omx,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,omy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,omz,ignore_me)

    r11(1:nx:2,:,:)=real(s11(1:lx,:,:),sp); r11(2:nx:2,:,:)=aimag(s11(1:lx,:,:))
    r12(1:nx:2,:,:)=real(s12(1:lx,:,:),sp); r12(2:nx:2,:,:)=aimag(s12(1:lx,:,:))
    r13(1:nx:2,:,:)=real(s13(1:lx,:,:),sp); r13(2:nx:2,:,:)=aimag(s13(1:lx,:,:))
                                             
    r22(1:nx:2,:,:)=real(s22(1:lx,:,:),sp); r22(2:nx:2,:,:)=aimag(s22(1:lx,:,:))
    r23(1:nx:2,:,:)=real(s23(1:lx,:,:),sp); r23(2:nx:2,:,:)=aimag(s23(1:lx,:,:))
    r33(1:nx:2,:,:)=real(s33(1:lx,:,:),sp); r33(2:nx:2,:,:)=aimag(s33(1:lx,:,:))

    rox(1:nx:2,:,:)=real(omx(1:lx,:,:),sp); rox(2:nx:2,:,:)=aimag(omx(1:lx,:,:))
    roy(1:nx:2,:,:)=real(omy(1:lx,:,:),sp); roy(2:nx:2,:,:)=aimag(omy(1:lx,:,:))
    roz(1:nx:2,:,:)=real(omz(1:lx,:,:),sp); roz(2:nx:2,:,:)=aimag(omz(1:lx,:,:))


    open(15,file='./out/b11'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b11
    close(15)
    open(15,file='./out/b12'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b12
    close(15)
    open(15,file='./out/b13'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b13
    close(15)
    open(15,file='./out/b21'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b21
    close(15)
    open(15,file='./out/b22'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b22
    close(15)
    open(15,file='./out/b23'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b23
    close(15)
    open(15,file='./out/b31'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b31
    close(15)
    open(15,file='./out/b32'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b32
    close(15)
    open(15,file='./out/b33'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
      read(15) b33
    close(15)

    open(25, file = './out/xyz'//trim(xyzlst)//'-'//trim(str1), form = 'unformatted')
      read(25) xp, yp, zp
    close(25)

    mcntr1 = 0._dp; mcntr2 = 0._dp; mcntr3 = 0._dp; mcntr4 = 0._dp
    mcntr1a = 0._dp; mcntr2a = 0._dp; mcntr3a = 0._dp; mcntr4a = 0._dp
    do ll = 1, ngroups
    do mm = 1, npntgrp

      ii = xx(mm, ll); jj = yy(mm,ll); kk = zz(mm,ll)

      if ( mod(ii, 2) .eq. 0 ) then
        cc(1,1) = aimag( b11(ii/2,jj,kk) )
        cc(1,2) = aimag( b12(ii/2,jj,kk) )
        cc(1,3) = aimag( b13(ii/2,jj,kk) )
        cc(2,1) = aimag( b21(ii/2,jj,kk) )
        cc(2,2) = aimag( b22(ii/2,jj,kk) )
        cc(2,3) = aimag( b23(ii/2,jj,kk) )
        cc(3,1) = aimag( b31(ii/2,jj,kk) )
        cc(3,2) = aimag( b32(ii/2,jj,kk) )
        cc(3,3) = aimag( b33(ii/2,jj,kk) )
      else
        cc(1,1) = real( b11(ii/2+1,jj,kk) )
        cc(1,2) = real( b12(ii/2+1,jj,kk) )
        cc(1,3) = real( b13(ii/2+1,jj,kk) )
        cc(2,1) = real( b21(ii/2+1,jj,kk) )
        cc(2,2) = real( b22(ii/2+1,jj,kk) )
        cc(2,3) = real( b23(ii/2+1,jj,kk) )
        cc(3,1) = real( b31(ii/2+1,jj,kk) )
        cc(3,2) = real( b32(ii/2+1,jj,kk) )
        cc(3,3) = real( b33(ii/2+1,jj,kk) )
      end if
      cc = matmul(cc, transpose(cc))


      pp = (kk-1) * nx * ny + (jj-1) * nx + ii
      y(1) = xp(pp); y(2) = yp(pp); y(3) = zp(pp)
      ! NOTE: ii, jj, kk will be destroyed below. 
 
      call pre_interp(y,dxyz,bg,lhnode)
      call ixiyiz(lhnode,nx,ny,nz,ix,iy,iz)

      sij = 0._dp; omp = 0._dp
      do  kk=1,norder
        iz1 = iz(kk)
        do  jj=1,norder
          iy1 = iy(jj)
          do  ii=1,norder
            ix1 = ix(ii)
      
            ignore_me = bg(ii,1)*bg(jj,2)*bg(kk,3)
      
            sij(1,1)=sij(1,1)+r11(ix1,iy1,iz1)*ignore_me
            sij(1,2)=sij(1,2)+r12(ix1,iy1,iz1)*ignore_me
            sij(1,3)=sij(1,3)+r13(ix1,iy1,iz1)*ignore_me
      
            sij(2,2)=sij(2,2)+r22(ix1,iy1,iz1)*ignore_me
            sij(2,3)=sij(2,3)+r23(ix1,iy1,iz1)*ignore_me
            sij(3,3)=sij(3,3)+r33(ix1,iy1,iz1)*ignore_me
      
            omp(1)=omp(1)+rox(ix1,iy1,iz1)*ignore_me
            omp(2)=omp(2)+roy(ix1,iy1,iz1)*ignore_me
            omp(3)=omp(3)+roz(ix1,iy1,iz1)*ignore_me
          end do
        end do
      end do
      sij(2,1) = sij(1,2); sij(3,1) = sij(1,3); sij(3,2) = sij(2,3)

      ! The eigenvalues are in ascending order: smallest first.
      call dsyevr("V", "A", "U", evsize, cc, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                  nfound, evcc, evtrcc, evsize, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info

      ! define a right-handed frame: gamma = alfa X beta
      evtrcc(1,1) = evtrcc(2,3) - evtrcc(3,2)
      evtrcc(2,1) = evtrcc(3,3) - evtrcc(1,2)
      evtrcc(3,1) = evtrcc(1,3) - evtrcc(2,2)

      ! The eigenvalues are in ascending order: smallest first.

      ominc(1) = dot_product( omp, evtrcc(:,3) ) * .5_dp
      ominc(2) = dot_product( omp, evtrcc(:,2) ) * .5_dp
      ominc(3) = dot_product( omp, evtrcc(:,1) ) * .5_dp

      rotinc(1) = ominc(1) + dot_product(evtrcc(:,2), matmul(sij, evtrcc(:,1))) &
                        * ( evcc(2) + evcc(1) ) / ( evcc(2) - evcc(1) )
      rotinc(2) = ominc(2) + dot_product(evtrcc(:,1), matmul(sij, evtrcc(:,3))) &
                        * ( evcc(3) + evcc(1) ) / ( evcc(1) - evcc(3) )
      rotinc(3) = ominc(3) + dot_product(evtrcc(:,3), matmul(sij, evtrcc(:,2))) &
                        * ( evcc(3) + evcc(2) ) / ( evcc(3) - evcc(2) ) 

      if (mm .eq. 1) then

          rotinc0 = rotinc
          ominc0  = ominc
          evtrcc0 = evtrcc

      else

          align = dot_product( evtrcc0(:,3), evtrcc(:,3) ) 
 
          cntr1 = - rotinc0(2) * dot_product( evtrcc0(:,1), evtrcc(:,3) ) * align / abs(align)
          cntr2 =   rotinc0(3) * dot_product( evtrcc0(:,2), evtrcc(:,3) ) * align / abs(align)
          cntr3 = - rotinc(2)  * dot_product( evtrcc0(:,3), evtrcc(:,1) ) * align / abs(align)
          cntr4 =   rotinc(3)  * dot_product( evtrcc0(:,3), evtrcc(:,2) ) * align / abs(align)
 
          cntr1a = - rotinc0(2) * dot_product( evtrcc0(:,1), evtrcc(:,3) ) * align
          cntr2a =   rotinc0(3) * dot_product( evtrcc0(:,2), evtrcc(:,3) ) * align
          cntr3a = - rotinc(2)  * dot_product( evtrcc0(:,3), evtrcc(:,1) ) * align
          cntr4a =   rotinc(3)  * dot_product( evtrcc0(:,3), evtrcc(:,2) ) * align
 
          mcntr1 = mcntr1 + cntr1
          mcntr2 = mcntr2 + cntr2
          mcntr3 = mcntr3 + cntr3
          mcntr4 = mcntr4 + cntr4
 
          mcntr1a = mcntr1a + cntr1a
          mcntr2a = mcntr2a + cntr2a
          mcntr3a = mcntr3a + cntr3a
          mcntr4a = mcntr4a + cntr4a

      end if

    end do
    end do

    mcntr1 = mcntr1 * const
    mcntr2 = mcntr2 * const
    mcntr3 = mcntr3 * const
    mcntr4 = mcntr4 * const

    mcntr1a = mcntr1a * const
    mcntr2a = mcntr2a * const
    mcntr3a = mcntr3a * const
    mcntr4a = mcntr4a * const

    write(27, '(I6, 15E15.3)') nfile, mcntr1a, mcntr2a, mcntr3a, mcntr4a, &
                                      mcntr1,  mcntr2,  mcntr3,  mcntr4

    nfile = nfile + 1
  end do
  close(30)
  close(27)

  deallocate(b11,b12,b13,b21,b22,b23,b31,b32,b33)
  deallocate(s11,s12,s13,s22,s23,s33,omx,omy,omz)
  deallocate(r11,r12,r13,r22,r23,r33,rox,roy,roz)
  deallocate(kx,ky,kz,xp,yp,zp,xx,yy,zz)
  
  call destroyplan3d

  write(*,*) 'defm-cor-rot-bw.x done.'

end program defmcijcorbw
