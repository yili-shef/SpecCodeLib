program defmcijcorbw
  use, intrinsic :: iso_c_binding
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

  ! Total memory of the machine would be about 80G
  integer, parameter :: batchsize = 20
  real(dp), dimension(batchsize) :: mcntr1, mcntr2

  integer, parameter :: npdf = 40 
  real(sp), parameter :: binw = 1._sp/npdf
  real(dp), dimension(npdf,batchsize) :: pcosth, rotcnd1, rotcnd2

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  ! ----------------------------------------------------------------------------------

  real(dp), dimension(evsize, evsize)  :: cc, sij, evtrcc, evtrcc0
  real(dp), dimension(evsize) :: evcc, evcc0, omp, ominc, rotinc, rotinc0, ominc0, ecross
  real(dp) :: e1, e2, e3

  real(dp) :: mcntr1, mcntr2, mcntr3, mcntr4, align, cntr1, cntr2, cntr3, cntr4
  real(dp) :: mcntr1a, mcntr2a, mcntr3a, mcntr4a, cntr1a, cntr2a, cntr3a, cntr4a

  complex(sp), allocatable, dimension(:,:,:,:) :: ux,uy,uz
  complex(sp), allocatable, dimension(:,:,:,:) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33,omx,omy,omz
  real(sp),    pointer,     dimension(:,:,:) :: r11,r12,r13,r22,r23,r33,rox,roy,roz
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz
  real(sp),    allocatable, dimension(:,:)     :: xp, yp, zp
  type(C_ptr) :: cptr

  character(80) :: str, flnm, prf, str1, xyzlst, mark
  integer :: nlst, forward
  
  if (iargc() .ne. 6) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-rotcndalign.x nx filelist prefix inicoors npntgrp fw'
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

  allocate( b11(lx1,ly,lz,batchsize), b12(lx1,ly,lz,batchsize), b13(lx1,ly,lz,batchsize) )
  allocate( b21(lx1,ly,lz,batchsize), b22(lx1,ly,lz,batchsize), b23(lx1,ly,lz,batchsize) )
  allocate( b31(lx1,ly,lz,batchsize), b32(lx1,ly,lz,batchsize), b33(lx1,ly,lz,batchsize) )
  allocate( ux(lx1,ly,lz,batchsize),  uy(lx1,ly,lz,batchsize),  uz(lx1,ly,lz,batchsize)  )
  
  allocate( s11(lx1,ly,lz), s12(lx1,ly,lz), s13(lx1,ly,lz) )
  allocate( omx(lx1,ly,lz), s22(lx1,ly,lz), s23(lx1,ly,lz) )
  allocate( omy(lx1,ly,lz), omz(lx1,ly,lz), s33(lx1,ly,lz) )

  !allocate( r11(nx,ny,nz), r12(nx,ny,nz), r13(nx,ny,nz) )
  !allocate( rox(nx,ny,nz), r22(nx,ny,nz), r23(nx,ny,nz) )
  !allocate( roy(nx,ny,nz), roz(nx,ny,nz), r33(nx,ny,nz) )

  allocate( kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz) )
  allocate( xp(nx*ny*nz,batchsize), yp(nx*ny*nz,batchsize), zp(nx*ny*nz,batchsize) )

  call fftwplan3de(nx,ny,nz)
  call wavenumber(kx,ky,kz,lx1,ly,lz)

  open(30, file = './out/'//trim(str)//'.dat', form = 'unformatted')
    read(30) xx, yy, zz
  close(30)

  open(27, file = 'defm-rotcorel-mean-'//trim(str)//'-'//trim(flnm)//'.dat')
  open(28, file = 'defm-rotcndalign-'//trim(str)//'-'//trim(flnm)//'.dat')

  open(30, file = trim(flnm)//'.list')

  nfile = 0
  do while ( .not. eof(30) )

    ifile = 1
    do while ( ifile .le. batchsize .and. .not. eof(30) )

      read(30,*) str1
 
      open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
        read(15) ux(:,:,:,ifile)
      close(15)
      open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
        read(15) uy(:,:,:,ifile)
      close(15)
      open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
        read(15) uz(:,:,:,ifile)
      close(15)
 
      open(15,file='./out/b11'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
        read(15) b11(:,:,:,ifile)
      close(15)
      open(15,file='./out/b12'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
        read(15) b12(:,:,:,ifile)
      close(15)
      open(15,file='./out/b13'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
        read(15) b13(:,:,:,ifile)
      close(15)
      open(15,file='./out/b21'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
        read(15) b21(:,:,:,ifile)
      close(15)
      open(15,file='./out/b22'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
        read(15) b22(:,:,:,ifile)
      close(15)
      open(15,file='./out/b23'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
        read(15) b23(:,:,:,ifile)
      close(15)
      open(15,file='./out/b31'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
        read(15) b31(:,:,:,ifile)
      close(15)
      open(15,file='./out/b32'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
        read(15) b32(:,:,:,ifile)
      close(15)
      open(15,file='./out/b33'//prf(1:len_trim(prf))//str1(1:len_trim(str1)),form='unformatted')
        read(15) b33(:,:,:,ifile)
      close(15)
 
      open(25, file = './out/xyz'//trim(xyzlst)//'-'//trim(str1), form = 'unformatted')
        read(25) xp(:,ifile), yp(:,ifile), zp(:,ifile)
      close(25)

      ifile = ifile + 1

    end do
    filecnt = ifile - 1
    write(*,*) 'defm-rotcndalign.x: data file ', str1(1:len_trim(str1))

    mcntr1 = 0._dp; mcntr2 = 0._dp
    pcosth = 0._dp; rotcnd1 = 0._dp; rotcnd2 = 0._dp
    do ifile = 1, filecnt

      s11 = eye * kx * ux(:,:,:,ifile)
      s22 = eye * ky * uy(:,:,:,ifile)
      s33 = eye * kz * uz(:,:,:,ifile)
      s12 = eye * (kx * uy(:,:,:,ifile) + ky * ux(:,:,:,ifile))*0.5_sp
      s13 = eye * (kx * uz(:,:,:,ifile) + kz * ux(:,:,:,ifile))*0.5_sp
      s23 = eye * (ky * uz(:,:,:,ifile) + kz * uy(:,:,:,ifile))*0.5_sp

      call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
 
      omx = eye * (ky * uz(:,:,:,ifile) - kz * uy(:,:,:,ifile))
      omy = eye * (kz * ux(:,:,:,ifile) - kx * uz(:,:,:,ifile))
      omz = eye * (kx * uy(:,:,:,ifile) - ky * ux(:,:,:,ifile))
      call rfftwnd_f77_one_complex_to_real(c2r3d,omx,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,omy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,omz,ignore_me)
 
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

      cptr = c_loc(omx)
      call c_f_pointer(cptr, rox, [2*lx1,ly,lz])
      cptr = c_loc(omy)
      call c_f_pointer(cptr, roy, [2*lx1,ly,lz])
      cptr = c_loc(omz)
      call c_f_pointer(cptr, roz, [2*lx1,ly,lz])

!      r11(1:nx:2,:,:)=real(ux(1:lx,:,:,ifile),sp); r11(2:nx:2,:,:)=aimag(ux(1:lx,:,:,ifile))
!      r12(1:nx:2,:,:)=real(s12(1:lx,:,:),sp); r12(2:nx:2,:,:)=aimag(s12(1:lx,:,:))
!      r13(1:nx:2,:,:)=real(s13(1:lx,:,:),sp); r13(2:nx:2,:,:)=aimag(s13(1:lx,:,:))
                                               
!      r22(1:nx:2,:,:)=real(uy(1:lx,:,:,ifile),sp); r22(2:nx:2,:,:)=aimag(uy(1:lx,:,:,ifile))
!      r23(1:nx:2,:,:)=real(s23(1:lx,:,:),sp); r23(2:nx:2,:,:)=aimag(s23(1:lx,:,:))
!      r33(1:nx:2,:,:)=real(uz(1:lx,:,:,ifile),sp); r33(2:nx:2,:,:)=aimag(uz(1:lx,:,:,ifile))
 
!      rox(1:nx:2,:,:)=real(omx(1:lx,:,:),sp); rox(2:nx:2,:,:)=aimag(omx(1:lx,:,:))
!      roy(1:nx:2,:,:)=real(omy(1:lx,:,:),sp); roy(2:nx:2,:,:)=aimag(omy(1:lx,:,:))
!      roz(1:nx:2,:,:)=real(omz(1:lx,:,:),sp); roz(2:nx:2,:,:)=aimag(omz(1:lx,:,:))

      do ll = 1, ngroups
      do mm = 1, npntgrp
     
        ii = xx(mm, ll); jj = yy(mm,ll); kk = zz(mm,ll)
     
        if ( mod(ii, 2) .eq. 0 ) then
            cc(1,1) = aimag( b11(ii/2,jj,kk,ifile) )
            cc(1,2) = aimag( b12(ii/2,jj,kk,ifile) )
            cc(1,3) = aimag( b13(ii/2,jj,kk,ifile) )
            cc(2,1) = aimag( b21(ii/2,jj,kk,ifile) )
            cc(2,2) = aimag( b22(ii/2,jj,kk,ifile) )
            cc(2,3) = aimag( b23(ii/2,jj,kk,ifile) )
            cc(3,1) = aimag( b31(ii/2,jj,kk,ifile) )
            cc(3,2) = aimag( b32(ii/2,jj,kk,ifile) )
            cc(3,3) = aimag( b33(ii/2,jj,kk,ifile) )
        else                             
            cc(1,1) = real( b11(ii/2+1,jj,kk,ifile) )
            cc(1,2) = real( b12(ii/2+1,jj,kk,ifile) )
            cc(1,3) = real( b13(ii/2+1,jj,kk,ifile) )
            cc(2,1) = real( b21(ii/2+1,jj,kk,ifile) )
            cc(2,2) = real( b22(ii/2+1,jj,kk,ifile) )
            cc(2,3) = real( b23(ii/2+1,jj,kk,ifile) )
            cc(3,1) = real( b31(ii/2+1,jj,kk,ifile) )
            cc(3,2) = real( b32(ii/2+1,jj,kk,ifile) )
            cc(3,3) = real( b33(ii/2+1,jj,kk,ifile) )
        end if
        cc = matmul(cc, transpose(cc))
     
     
        pp = (kk-1) * nx * ny + (jj-1) * nx + ii
        y(1) = xp(pp,ifile); y(2) = yp(pp,ifile); y(3) = zp(pp,ifile)
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
     
        ! Reverse the order of the eigen values and vectors:
        evcc = evcc(evsize:1:-1)
        evtrcc = evtrcc(:,evsize:1:-1)
   
        ! define a right-handed frame: gamma = alfa X beta
        evtrcc(1,3) = evtrcc(2,1) * evtrcc(3,2) - evtrcc(3,1) * evtrcc(2,2)
        evtrcc(2,3) = evtrcc(3,1) * evtrcc(1,2) - evtrcc(1,1) * evtrcc(3,2)
        evtrcc(3,3) = evtrcc(1,1) * evtrcc(2,2) - evtrcc(2,1) * evtrcc(1,2)
       
        ! The eigenvalues are in ascending order: smallest first.
     
        ominc(1) = dot_product( omp, evtrcc(:,1) ) * .5_dp
        ominc(2) = dot_product( omp, evtrcc(:,2) ) * .5_dp
        ominc(3) = dot_product( omp, evtrcc(:,3) ) * .5_dp
     
        rotinc(1) = ominc(1) + dot_product(evtrcc(:,2), matmul(sij, evtrcc(:,3))) &
                          * ( evcc(2) + evcc(1) ) / ( evcc(2) - evcc(3) )
        rotinc(2) = ominc(2) + dot_product(evtrcc(:,3), matmul(sij, evtrcc(:,1))) &
                          * ( evcc(3) + evcc(1) ) / ( evcc(3) - evcc(1) )
        rotinc(3) = ominc(3) + dot_product(evtrcc(:,1), matmul(sij, evtrcc(:,2))) &
                          * ( evcc(1) + evcc(2) ) / ( evcc(1) - evcc(2) ) 
     
        if (mm .eq. 1) then
     
            rotinc0 = rotinc
            ominc0  = ominc
            evtrcc0 = evtrcc
     
        else
     
            align = dot_product( evtrcc0(:,1), evtrcc(:,1) ) 
            costh = abs( align ) 
            align = align / costh
     
            ecross(1) = evtrcc0(2,1) * evtrcc(3,1) - evtrcc0(3,1) * evtrcc(2,1)
            ecross(2) = evtrcc0(3,1) * evtrcc(1,1) - evtrcc0(1,1) * evtrcc(3,1)
            ecross(3) = evtrcc0(1,1) * evtrcc(2,1) - evtrcc0(2,1) * evtrcc(1,1)
     
            e1 = dot_product( ecross, evtrcc0(:,1) )
            e2 = dot_product( ecross, evtrcc0(:,2) )
            e3 = dot_product( ecross, evtrcc0(:,3) )
            cntr1 = rotinc0(1) * e1 + rotinc0(2) * e2 + rotinc0(3) * e3
            cntr1 = cntr1 * align
     
            e1 = -dot_product( ecross, evtrcc(:,1) )
            e2 = -dot_product( ecross, evtrcc(:,2) )
            e3 = -dot_product( ecross, evtrcc(:,3) )
            cntr2 = rotinc(1) * e1 + rotinc(2) * e2 + rotinc(3) * e3
            cntr2 = cntr2 * align
     
            mcntr1(ifile) = mcntr1(ifile) + cntr1 * align/abs(align)
            mcntr2(ifile) = mcntr2(ifile) + cntr2 * align/abs(align)
     
            ii = floor( costh/binw ) + 1
            pcosth (ii,ifile) = pcosth (ii,ifile) + 1
            rotcnd1(ii,ifile) = rotcnd1(ii,ifile) + cntr1
            rotcnd2(ii,ifile) = rotcnd2(ii,ifile) + cntr2
        end if
     
      end do
      end do

    end do

    rotcnd1 = rotcnd1 / (pcosth + mytiny)
    rotcnd2 = rotcnd2 / (pcosth + mytiny)
    pcosth = pcosth * const / binw

    mcntr1 = mcntr1 * const
    mcntr2 = mcntr2 * const
    do ifile = 1, filecnt
      write(27, '(I6, 15E15.3)') nfile, mcntr1(ifile),  mcntr2(ifile)

      write(28, '("# index = ", I6)') nfile
      do jj = 1, npdf
        write(28, '(15E15.3)')  binw * (jj-.5_sp), pcosth(jj,ifile), rotcnd1(jj,ifile), rotcnd2(jj,ifile)
      end do
      write(28,*)
      write(28,*)

      nfile = nfile + 1
    end do

  end do
  close(30)
  close(27)
  close(28)

  deallocate(b11,b12,b13,b21,b22,b23,b31,b32,b33)
  deallocate(s11,s12,s13,s22,s23,s33,omx,omy,omz)
  deallocate(kx,ky,kz,xp,yp,zp,xx,yy,zz)
  deallocate(ux,uy,uz)
  
  call destroyplan3d

  write(*,*) 'defm-rotcndalign.x done.'

end program defmcijcorbw
