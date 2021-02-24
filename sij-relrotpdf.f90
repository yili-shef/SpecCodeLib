program defmsijincij
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nfile
  integer :: ifile, filecnt, pp, npnt

  real(sp), parameter :: rnu = 0.006_sp
  integer,  parameter :: npdf = 80
  real(sp), parameter :: bnd = 8._sp, binw = 2*bnd/npdf

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  ! ----------------------------------------------------------------------------------

  real(dp), dimension(evsize,evsize) :: sij, evtrsij, pij, vsij, PijInSij, VijInSij
  real(dp), dimension(evsize,evsize) :: AxesSij
  real(dp), dimension(evsize) :: evsij, vort, PijContSij
  real(dp), dimension(evsize) :: ominsij, VortContSij, VortContSijLab, PijContSijLab
  real(dp), dimension(evsize) :: VijContSij, VijContSijLab, ecross
  real(dp), dimension(evsize) :: AllRot

  ! Total memory of the machine would be about 80G
  integer, parameter :: batchsize = 20

  integer, parameter :: naverage = 500
  real(sp), dimension(naverage) :: xx, yy, zz
  real(dp), dimension(evsize,naverage) :: AllRot0
  real(dp), dimension(evsize,evsize,naverage) :: AxesSij0

  real(dp), allocatable, dimension(:) :: meanrelrot, varrelrot, skewrelrot, allpnt, pntcnt
  real(dp), allocatable, dimension(:,:) :: prelrot

  complex(sp), allocatable, dimension(:,:,:,:) :: s11,s22,s33
  complex(sp), allocatable, dimension(:,:,:,:) :: p33
  complex(sp), allocatable, dimension(:,:,:) :: s12,s13,s23, ox, oy, oz
  complex(sp), allocatable, dimension(:,:,:) :: p11,p22,p12,p13,p23
  complex(sp), allocatable, dimension(:,:,:) :: vs11,vs22,vs33,vs12,vs13,vs23
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz, k2

  character(80) :: str, flnm, str1
  real(sp) :: ignore_me, const, xpp, ypp, zpp
  real(sp) :: dx, dy, dz, ran1
  real(dp) :: relrot, rr, var, skew

  
  if (iargc() .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./sij-relrotpdf.x nx filelist'
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

  ny=nx; nz=nx; lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  const = 1._sp/(nx*ny*nz)
  dx = 2*pi/nx; dy = dx; dz = dx
  npnt = nx/2

  allocate( meanrelrot(npnt), varrelrot(npnt), skewrelrot(npnt) )
  allocate( pntcnt(npnt), allpnt(npnt), prelrot(npdf,npnt) )

  allocate( s11(lx1,ly,lz,batchsize), s22(lx1,ly,lz,batchsize), s33(lx1,ly,lz,batchsize) )
  allocate( s12(lx1,ly,lz), s13(lx1,ly,lz), s23(lx1,ly,lz) )

  allocate( ox(lx1,ly,lz), oy(lx1,ly,lz), oz(lx1,ly,lz) )

  allocate( p11(lx1,ly,lz), p22(lx1,ly,lz), p33(lx1,ly,lz,batchsize) )
  allocate( p12(lx1,ly,lz), p13(lx1,ly,lz), p23(lx1,ly,lz) )

  allocate( vs11(lx1,ly,lz), vs22(lx1,ly,lz), vs33(lx1,ly,lz) )
  allocate( vs12(lx1,ly,lz), vs13(lx1,ly,lz), vs23(lx1,ly,lz) )

  allocate( kx(lx1,ly,lz), ky(lx1,ly,lz), kz(lx1,ly,lz), k2(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  prelrot = 0._dp; meanrelrot = 0._dp; varrelrot = 0._dp; skewrelrot = 0._dp
  pntcnt = 0._dp; allpnt = 0._dp
  do while ( .not. eof(30) )

    ifile = 1
    do while ( ifile .le. batchsize .and. .not. eof(30) )

      read(30,*) str1
 
      open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
        read(15) s11(:,:,:,ifile)
      close(15)
      open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
        read(15) s22(:,:,:,ifile)
      close(15)
      open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
        read(15) s33(:,:,:,ifile)
      close(15)
      open(15,file='./out/p'//str1(1:len_trim(str1)),form='unformatted')
        read(15) p33(:,:,:,ifile)
      close(15)

      ifile = ifile + 1

    end do
    filecnt = ifile - 1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    ! Loop over the batch
    do ifile = 1, filecnt

      ! Vorticity
      ox = eye * (ky * s33(:,:,:,ifile) - kz * s22(:,:,:,ifile))
      oy = eye * (kz * s11(:,:,:,ifile) - kx * s33(:,:,:,ifile))
      oz = eye * (kx * s22(:,:,:,ifile) - ky * s11(:,:,:,ifile))
      call rfftwnd_f77_one_complex_to_real(c2r3d,ox,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,oy,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,oz,ignore_me)
 

      s12 = eye * (kx * s22(:,:,:,ifile) + ky * s11(:,:,:,ifile))*0.5_sp
      s13 = eye * (kx * s33(:,:,:,ifile) + kz * s11(:,:,:,ifile))*0.5_sp
      s23 = eye * (ky * s33(:,:,:,ifile) + kz * s22(:,:,:,ifile))*0.5_sp
      s11(:,:,:,ifile) = eye * kx * s11(:,:,:,ifile)
      s22(:,:,:,ifile) = eye * ky * s22(:,:,:,ifile)
      s33(:,:,:,ifile) = eye * kz * s33(:,:,:,ifile)

      vs11 = - rnu * k2 * s11(:,:,:,ifile)
      vs22 = - rnu * k2 * s22(:,:,:,ifile)
      vs33 = - rnu * k2 * s33(:,:,:,ifile)
      vs12 = - rnu * k2 * s12
      vs13 = - rnu * k2 * s13
      vs23 = - rnu * k2 * s23
      call rfftwnd_f77_one_complex_to_real(c2r3d,s11(:,:,:,ifile),ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s22(:,:,:,ifile),ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s33(:,:,:,ifile),ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
 
      call rfftwnd_f77_one_complex_to_real(c2r3d,vs11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,vs22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,vs33,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,vs12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,vs13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,vs23,ignore_me)
 
      ! pressure Hessian
      p11 = -kx * kx * p33(:,:,:,ifile)
      p22 = -ky * ky * p33(:,:,:,ifile)
      p12 = -kx * ky * p33(:,:,:,ifile)
      p13 = -kx * kz * p33(:,:,:,ifile)
      p23 = -ky * kz * p33(:,:,:,ifile)
      p33(:,:,:,ifile) = - kz * kz * p33(:,:,:,ifile)

      call rfftwnd_f77_one_complex_to_real(c2r3d,p11,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,p12,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,p13,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,p22,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,p23,ignore_me)
      call rfftwnd_f77_one_complex_to_real(c2r3d,p33(:,:,:,ifile),ignore_me)


      ! Calculate the centers of the groups
      ll = -100
      do mm = 1, naverage
  
        ii = floor( ran1(ll) * nx ) + 1
        jj = floor( ran1(ll) * ny ) + 1
        kk = floor( ran1(ll) * nz ) + 1
  
        xx(mm) = (ii-1)*dx; yy(mm) = (jj-1)*dy; zz(mm) = (kk-1)*dz
        call CalRotOnPoint(ii,jj,kk,AllRot0(:,mm),AxesSij0(:,:,mm))
  
      end do
      ! END Calculate the centers of the groups


      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx

        xpp = dx * (ii-1); ypp = dy * (jj-1); zpp = dz * (kk-1)
     
        call CalRotOnPoint(ii,jj,kk,AllRot,AxesSij)

        ! Loop over groups
        do ll = 1, naverage

          ecross(1) = AxesSij0(2,1,ll) * AxesSij(3,1) - AxesSij0(3,1,ll) * AxesSij(2,1)
          ecross(2) = AxesSij0(3,1,ll) * AxesSij(1,1) - AxesSij0(1,1,ll) * AxesSij(3,1)
          ecross(3) = AxesSij0(1,1,ll) * AxesSij(2,1) - AxesSij0(2,1,ll) * AxesSij(1,1)
 
          relrot = dot_product( ( AllRot0(:,ll) - AllRot ), ecross )
          relrot = relrot * sign( 1._dp, dot_product( AxesSij0(:,1,ll), AxesSij(:,1) ) )

          mm = floor( sqrt( (xpp-xx(ll))**2 + (ypp-yy(ll))**2 + (zpp-zz(ll))**2 ) / dx + 0.5_sp )  + 1

          if ( mm .ge. 1 .and. mm .le. npnt ) then

            rr = floor( (relrot + bnd) / binw ) + 1
            if ( rr .ge. 1 .and. rr .le. npdf) then
              prelrot(rr,mm) = prelrot(rr,mm) + 1
              meanrelrot(mm) = meanrelrot(mm) + relrot
              varrelrot(mm) = varrelrot(mm) + relrot * relrot
              skewrelrot(mm) = skewrelrot(mm) + relrot ** 3
              pntcnt(mm) = pntcnt(mm) + 1
            end if
            allpnt(mm) = allpnt(mm) + 1

          end if

        end do
        ! END Loop over groups
     
      end do
      end do
      end do

    end do
    ! END Loop over batch

  end do
  close(30)

  meanrelrot = meanrelrot / pntcnt
  varrelrot = varrelrot / pntcnt
  skewrelrot = skewrelrot / pntcnt

  open(26, file = 'sij-relrotmeans-'//trim(flnm)//'.dat')
  open(27, file = 'sij-relrotpdf-'//trim(flnm)//'.dat')

  do ll = 1, npnt

    var = varrelrot(ll) - meanrelrot(ll)**2
    skew = skewrelrot(ll) - 3 * varrelrot(ll) * meanrelrot(ll) + 2 * meanrelrot(ll) ** 3
    skew = skew / ( var + mytiny )**1.5_dp

    write(26, '(I6, 15E15.3)') ll, meanrelrot(ll), varrelrot(ll), skewrelrot(ll), var, skew

    prelrot(:,ll) = prelrot(:,ll) / allpnt(ll)
    write(*,*) 'sij-relrotpdf.x: normalization of prelrot:', sum(prelrot(:,ll))
    prelrot(:,ll) = prelrot(:,ll) / binw

    write(27, '("# index = ", I6)') ll
    do pp = 1, npdf
      write(27, '(15E15.3)') binw * (pp-.5_sp)-bnd, prelrot(pp,ll)
    end do

  end do

  close(27)
  close(26)

  deallocate(kx,ky,kz,s11,s22,s33,s12,s13,s23)
  deallocate(k2,ox,oy,oz,p11,p12,p13,p22,p23,p33,vs11,vs22,vs33,vs12,vs13,vs23)
  
  call destroyplan3d

  write(*,*) 'sij-relrotpdf.x done.'

contains

  subroutine CalRotOnPoint(ipx,ipy,ipz,TotalRotOnPnt,AxesSijOnPnt)

    integer :: ipx, ipy, ipz, locll
    real(sp), dimension(evsize, evsize) :: AxesSijOnPnt
    real(sp), dimension(evsize) :: TotalRotOnPnt

        if ( mod(ipx, 2) .eq. 1 ) then
            locll = ipx/2 + 1

            sij(1,1) = real( s11(locll,ipy,ipz,ifile),sp )
            sij(2,2) = real( s22(locll,ipy,ipz,ifile),sp )
            sij(3,3) = real( s33(locll,ipy,ipz,ifile),sp )
            sij(1,2) = real( s12(locll,ipy,ipz),sp )
            sij(1,3) = real( s13(locll,ipy,ipz),sp )
            sij(2,3) = real( s23(locll,ipy,ipz),sp )
     
            vsij(1,1) = real( vs11(locll,ipy,ipz),sp )
            vsij(2,2) = real( vs22(locll,ipy,ipz),sp )
            vsij(3,3) = real( vs33(locll,ipy,ipz),sp )
            vsij(1,2) = real( vs12(locll,ipy,ipz),sp )
            vsij(1,3) = real( vs13(locll,ipy,ipz),sp )
            vsij(2,3) = real( vs23(locll,ipy,ipz),sp )
     
            ! Pressure Hessian
            pij(1,1) = real( p11(locll,ipy,ipz),sp )
            pij(1,2) = real( p12(locll,ipy,ipz),sp )
            pij(1,3) = real( p13(locll,ipy,ipz),sp )
            pij(2,2) = real( p22(locll,ipy,ipz),sp )
            pij(2,3) = real( p23(locll,ipy,ipz),sp )
            pij(3,3) = real( p33(locll,ipy,ipz, ifile),sp )

            vort(1) = real( ox(locll,ipy,ipz),sp )
            vort(2) = real( oy(locll,ipy,ipz),sp )
            vort(3) = real( oz(locll,ipy,ipz),sp )
        else
            locll = ipx/2

            sij(1,1) = aimag( s11(locll,ipy,ipz,ifile) )
            sij(2,2) = aimag( s22(locll,ipy,ipz,ifile) )
            sij(3,3) = aimag( s33(locll,ipy,ipz,ifile) )
            sij(1,2) = aimag( s12(locll,ipy,ipz) )
            sij(1,3) = aimag( s13(locll,ipy,ipz) )
            sij(2,3) = aimag( s23(locll,ipy,ipz) )
     
            vsij(1,1) = aimag( vs11(locll,ipy,ipz) )
            vsij(2,2) = aimag( vs22(locll,ipy,ipz) )
            vsij(3,3) = aimag( vs33(locll,ipy,ipz) )
            vsij(1,2) = aimag( vs12(locll,ipy,ipz) )
            vsij(1,3) = aimag( vs13(locll,ipy,ipz) )
            vsij(2,3) = aimag( vs23(locll,ipy,ipz) )
     
            ! Pressure Hessian
            pij(1,1) = aimag( p11(locll,ipy,ipz) )
            pij(1,2) = aimag( p12(locll,ipy,ipz) )
            pij(1,3) = aimag( p13(locll,ipy,ipz) )
            pij(2,2) = aimag( p22(locll,ipy,ipz) )
            pij(2,3) = aimag( p23(locll,ipy,ipz) )
            pij(3,3) = aimag( p33(locll,ipy,ipz, ifile) )

            vort(1) = aimag( ox(locll,ipy,ipz) )
            vort(2) = aimag( oy(locll,ipy,ipz) )
            vort(3) = aimag( oz(locll,ipy,ipz) )
        endif
        sij(2,1) = sij(1,2); sij(3,1) = sij(1,3); sij(3,2) = sij(2,3)
        pij(2,1) = pij(1,2); pij(3,1) = pij(1,3); pij(3,2) = pij(2,3)
        vsij(2,1) = vsij(1,2); vsij(3,1) = vsij(1,3); vsij(3,2) = vsij(2,3)


        call dsyevr("V", "A", "U", evsize, sij, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                    nfound, evsij, evtrsij, evsize, isuppz, work, lwork, iwork, liwork, info)
        ! sij will be destroyed on return
        if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev sij. info = ', info
        ! The eigenvalues are in ascending order: smallest first.

        ! Reverse the order of the eigen values and vectors:
        evsij = evsij(evsize:1:-1)
        evtrsij = evtrsij(:,evsize:1:-1)

        ! Define a right-handed frame: gamma = alfa X beta
        evtrsij(1,3) = evtrsij(2,1) * evtrsij(3,2) - evtrsij(3,1) * evtrsij(2,2)
        evtrsij(2,3) = evtrsij(3,1) * evtrsij(1,2) - evtrsij(1,1) * evtrsij(3,2)
        evtrsij(3,3) = evtrsij(1,1) * evtrsij(2,2) - evtrsij(2,1) * evtrsij(1,2)
       
        ominsij = matmul( transpose(evtrsij), vort )

        ! Vorticity induced rotation
        VortContSij(1) = -0.25_dp * ominsij(2)*ominsij(3)/(evsij(2) - evsij(3)) 
        VortContSij(2) = -0.25_dp * ominsij(3)*ominsij(1)/(evsij(3) - evsij(1)) 
        VortContSij(3) = -0.25_dp * ominsij(1)*ominsij(2)/(evsij(1) - evsij(2)) 

        VortContSijLab = matmul(evtrsij, VortContSij)

        ! Pressure Hessian induced rotation
        PijInSij = matmul( transpose(evtrsij), matmul(pij, evtrsij) )

        PijContSij(1) = -PijInSij(2,3)/(evsij(2) - evsij(3)) 
        PijContSij(2) = -PijInSij(3,1)/(evsij(3) - evsij(1)) 
        PijContSij(3) = -PijInSij(1,2)/(evsij(1) - evsij(2)) 
        
        PijContSijLab = matmul(evtrsij, PijContSij)

        ! Viscous diffusion induced rotation
        VijInSij = matmul( transpose(evtrsij), matmul(vsij, evtrsij) )

        VijContSij(1) = VijInSij(2,3)/(evsij(2) - evsij(3)) 
        VijContSij(2) = VijInSij(3,1)/(evsij(3) - evsij(1)) 
        VijContSij(3) = VijInSij(1,2)/(evsij(1) - evsij(2)) 
        
        VijContSijLab = matmul(evtrsij, VijContSij)

        AxesSijOnPnt = evtrsij
        TotalRotOnPnt = VortContSijLab + PijContSijLab + VijContSijLab

  end subroutine CalRotOnPoint

end program defmsijincij
