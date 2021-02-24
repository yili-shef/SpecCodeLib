program defmcijcorbw
  use mconstant
  use mwavenumber
  implicit none
  
  !------================= WARNING ========================--------
  !----------------------------------------------------------------
  ! This code is deprecated. Use defm-aligncorel-pdf.f90 instead, !
  ! which will handle both forward and backward simulations.      !
  !----------------------------------------------------------------


  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nfile
  real(sp) :: const

  integer, parameter :: ngroups = 512, npntgrp = 81
  integer, dimension(npntgrp, ngroups) :: xx, yy, zz

  integer, parameter :: npdf = 50 
  real(sp), parameter :: binw = 1._sp/npdf
  real(dp), dimension(npdf) :: pcosth

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  ! ----------------------------------------------------------------------------------

  real(dp), dimension(evsize,evsize)  :: cc
  real(dp), dimension(evsize, npntgrp) :: evcc
  real(dp), dimension(evsize, evsize, npntgrp) :: evtrcc

  real(dp) :: costh, mlogal, vlogal, logal, logal1, malcosth, corelal, coeff

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13,b21,b22,b23,b31,b32,b33

  character(80) :: str, flnm, prf, str1
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-aligncorel-pdf-bw.x nx filelist prefix'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     prefix: prefix for datafiles for bij'
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

  ! prefix
  call getarg(3,prf)
  prf = adjustl(prf)

  ny=nx; nz=nx; lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz

  const = 1._sp / (npntgrp * ngroups)

  allocate( b11(lx1,ly,lz), b12(lx1,ly,lz), b13(lx1,ly,lz) )
  allocate( b21(lx1,ly,lz), b22(lx1,ly,lz), b23(lx1,ly,lz) )
  allocate( b31(lx1,ly,lz), b32(lx1,ly,lz), b33(lx1,ly,lz) )

  open(30, file = './out/xyz512grps81pnts.dat', form = 'unformatted')
    read(30) xx, yy, zz
  close(30)

  open(27, file = 'defm-aligncorpdf-'//trim(flnm)//'.dat')

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

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

    pcosth = 0._dp
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

      ! The eigenvalues are in ascending order: smallest first.
      call dsyevr("V", "A", "U", evsize, cc, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                  nfound, evcc(:,mm), evtrcc(:,:,mm), evsize, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info

      costh = abs( dot_product( evtrcc(:,3,1), evtrcc(:,3,mm) ) )

      ii = floor( costh/binw ) + 1
      pcosth(ii) = pcosth(ii) + 1

    end do
    end do
    pcosth = pcosth / (npntgrp * ngroups) / binw


    write(27, '("# index = ", I6)') nfile
    do jj = 1, npdf
      write(27, '(15E15.3)')  binw * (jj-.5_sp), pcosth(jj)
    end do
    write(27,*)
    write(27,*)

    nfile = nfile + 1
  end do
  close(30)
  close(27)

  deallocate(b11,b12,b13,b21,b22,b23,b31,b32,b33)
  
  write(*,*) 'defm-aligncorel-pdf-bw.x done.'

end program defmcijcorbw
