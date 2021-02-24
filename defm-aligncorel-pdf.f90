program defmcijcorbw
  use mconstant
  use mwavenumber
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nfile,ifile,filecnt
  real(sp) :: const

  integer, parameter :: ngroups = 512
  integer :: npntgrp

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  ! ----------------------------------------------------------------------------------

  real(dp), dimension(evsize, evsize)  :: cc
  real(dp), dimension(evsize) :: evcc, evcc0
  real(dp), dimension(evsize, evsize) :: evtrcc, evtrcc0

  ! Total memory of the machine would be about 80G
  integer, parameter :: batchsize = 50
  integer, parameter :: npdf = 50 
  real(sp), parameter :: binw = 1._sp/npdf
  real(dp), dimension(npdf,batchsize) :: palcos, pbecos, pgacos

  complex(sp), allocatable, dimension(:,:,:,:) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
  integer, allocatable, dimension(:, :) :: xx, yy, zz

  integer :: forward
  real(dp) :: costh
  character(80) :: str, flnm, prf, str1
  
  if (iargc() .ne. 6) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-aligncorel-pdf.x nx filelist prefix inifile npntgrp fw'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     prefix: prefix for datafiles for bij'
          write(*,*) '                     initfile: initial data file for particles'
          write(*,*) '                     npntgrp: number of points in each group'
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
  flnm = adjustl(flnm)

  ! prefix
  call getarg(3,prf)
  prf = adjustl(prf)

  call getarg(5, str)
  read(str, '(I20)') npntgrp

  call getarg(6, str)
  read(str, '(I20)') forward

  ! initial data file
  call getarg(4, str)

  ny=nx; nz=nx; lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz

  const = 1._sp / ((npntgrp-1) * ngroups)

  allocate( xx(npntgrp, ngroups), yy(npntgrp, ngroups), zz(npntgrp, ngroups) )
  allocate( b11(lx1,ly,lz,batchsize), b12(lx1,ly,lz,batchsize), b13(lx1,ly,lz,batchsize) )
  allocate( b21(lx1,ly,lz,batchsize), b22(lx1,ly,lz,batchsize), b23(lx1,ly,lz,batchsize) )
  allocate( b31(lx1,ly,lz,batchsize), b32(lx1,ly,lz,batchsize), b33(lx1,ly,lz,batchsize) )

  open(30, file = './out/'//trim(str)//'.dat', form = 'unformatted')
    read(30) xx, yy, zz
  close(30)

  open(27, file = 'defm-aligncorpdf-'//trim(str)//'-'//trim(flnm)//'.dat')
  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  do while ( .not. eof(30) )

    ifile = 1
    do while ( ifile .le. batchsize .and. .not. eof(30) )

      read(30,*) str1
 
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

      ifile = ifile + 1
    end do
    filecnt = ifile - 1
    write(*,*) 'defm-aligncorel-pdf.x: data file ', str1(1:len_trim(str1))

    palcos = 0._dp; pbecos = 0._dp; pgacos = 0._dp
    do ifile = 1, filecnt

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
     
        ! The eigenvalues are in ascending order: smallest first.
        call dsyevr("V", "A", "U", evsize, cc, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                    nfound, evcc, evtrcc, evsize, isuppz, work, lwork, iwork, liwork, info)
        if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info
     

        if ( mm .eq. 1) then
            evtrcc0 = evtrcc
            evcc0 = evcc
        else
          costh = abs( dot_product( evtrcc0(:,3), evtrcc(:,3) ) )
          ii = floor( costh/binw ) + 1
          palcos(ii,ifile) = palcos(ii,ifile) + 1

          costh = abs( dot_product( evtrcc0(:,2), evtrcc(:,2) ) )
          ii = floor( costh/binw ) + 1
          pbecos(ii,ifile) = pbecos(ii,ifile) + 1

          costh = abs( dot_product( evtrcc0(:,1), evtrcc(:,1) ) )
          ii = floor( costh/binw ) + 1
          pgacos(ii,ifile) = pgacos(ii,ifile) + 1
        end if
     
      end do
      end do
      palcos(:,ifile) = palcos(:,ifile) * const / binw
      pbecos(:,ifile) = pbecos(:,ifile) * const / binw
      pgacos(:,ifile) = pgacos(:,ifile) * const / binw
    end do


    do ifile = 1, filecnt

      write(27, '("# index = ", I6)') nfile
      do jj = 1, npdf
        write(27, '(15E15.3)')  binw * (jj-.5_sp), palcos(jj,ifile), pbecos(jj,ifile), pgacos(jj,ifile)
      end do
      write(27,*)
      write(27,*)
 
      nfile = nfile + 1
    end do

  end do
  close(30)
  close(27)

  deallocate(b11,b12,b13,b21,b22,b23,b31,b32,b33)
  deallocate(xx,yy,zz)
  
  write(*,*) 'defm-aligncorel-pdf.x done.'

end program defmcijcorbw
