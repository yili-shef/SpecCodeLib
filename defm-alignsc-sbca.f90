program defmsijincij
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nfile
  integer :: ifile, filecnt, pp
  real(sp) :: ignore_me, const

  real(sp), parameter :: rnu = 0.006_sp

  integer,  parameter :: npnt = 50
  real(sp), parameter :: bnd = 1._sp, binw = bnd/npnt

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  ! ----------------------------------------------------------------------------------

  real(dp), dimension(evsize,evsize)  :: cij, sij, evtrcij, evtrsij, pij, vsij, PijInSij, VijInSij
  real(dp), dimension(evsize) :: evcij, evsij, vort, OmInCijLab, straincij, PijContSij
  real(dp), dimension(evsize) :: SijInCijLab, ominsij, VortContSij, VortContSijLab, PijContSijLab
  real(dp), dimension(evsize) :: VijContSij, VijContSijLab, evcross

  real(dp) :: scalign, signalign, srot1, srot2, srot3, crot1, crot2

  ! Total memory of the machine would be about 80G
  integer, parameter :: batchsize = 20

  real(dp), dimension(batchsize) :: meanscalign, meansrot1, meansrot2, meansrot3
  real(dp), dimension(batchsize) :: meancrot1, meancrot2

  real(dp), dimension(npnt, batchsize) :: pdfalign, cndsrot1, cndsrot2, cndsrot3, cndcrot1, cndcrot2

  complex(sp), allocatable, dimension(:,:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:,:) :: s11,s22,s33
  complex(sp), allocatable, dimension(:,:,:,:) :: p33
  complex(sp), allocatable, dimension(:,:,:) :: s12,s13,s23, ox, oy, oz
  complex(sp), allocatable, dimension(:,:,:) :: p11,p22,p12,p13,p23
  complex(sp), allocatable, dimension(:,:,:) :: vs11,vs22,vs33,vs12,vs13,vs23
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz, k2

  character(80) :: str, flnm, prf, str1

  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-alignsc-sbca.x nx filelist prefix'
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

  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  const = 1._sp/(nx*ny*nz)

  allocate( b11(lx1,ly,lz,batchsize), b12(lx1,ly,lz,batchsize), b13(lx1,ly,lz,batchsize) )
  allocate( b21(lx1,ly,lz,batchsize), b22(lx1,ly,lz,batchsize), b23(lx1,ly,lz,batchsize) )
  allocate( b31(lx1,ly,lz,batchsize), b32(lx1,ly,lz,batchsize), b33(lx1,ly,lz,batchsize) )

  allocate( s11(lx1,ly,lz,batchsize), s22(lx1,ly,lz,batchsize), s33(lx1,ly,lz,batchsize) )
  allocate( s12(lx1,ly,lz), s13(lx1,ly,lz), s23(lx1,ly,lz) )

  allocate( ox(lx1,ly,lz), oy(lx1,ly,lz), oz(lx1,ly,lz) )

  allocate( p11(lx1,ly,lz), p22(lx1,ly,lz), p33(lx1,ly,lz,batchsize) )
  allocate( p12(lx1,ly,lz), p13(lx1,ly,lz), p23(lx1,ly,lz) )

  allocate( vs11(lx1,ly,lz), vs22(lx1,ly,lz), vs33(lx1,ly,lz) )
  allocate( vs12(lx1,ly,lz), vs13(lx1,ly,lz), vs23(lx1,ly,lz) )

  allocate(  kx(lx1,ly,lz),  ky(lx1,ly,lz),  kz(lx1,ly,lz), k2(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)

  open(27, file = 'defm-alignsc-sbca-pdf-'//trim(flnm)//'.dat')
  open(26, file = 'defm-alignsc-sbca-mean-'//trim(flnm)//'.dat')

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
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
    write(*,*) 'data file ', str1(1:len_trim(str1))

    meanscalign = 0._dp; meansrot1 = 0._dp; meansrot2 = 0._dp; meansrot3 = 0._dp
    meancrot1 = 0._dp; meancrot2 = 0._dp
    pdfalign = 0._dp; cndsrot1 = 0._dp; cndsrot2 = 0._dp; cndsrot3 = 0._dp
    cndcrot1 = 0._dp; cndcrot2 = 0._dp
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

      do kk = 1, nz
      do jj = 1, ny
      do ii = 1, nx
     
        if ( mod(ii, 2) .eq. 1 ) then
            ll = ii/2 + 1

            sij(1,1) = real( s11(ll,jj,kk,ifile),sp )
            sij(2,2) = real( s22(ll,jj,kk,ifile),sp )
            sij(3,3) = real( s33(ll,jj,kk,ifile),sp )
            sij(1,2) = real( s12(ll,jj,kk),sp )
            sij(1,3) = real( s13(ll,jj,kk),sp )
            sij(2,3) = real( s23(ll,jj,kk),sp )
     
            vsij(1,1) = real( vs11(ll,jj,kk),sp )
            vsij(2,2) = real( vs22(ll,jj,kk),sp )
            vsij(3,3) = real( vs33(ll,jj,kk),sp )
            vsij(1,2) = real( vs12(ll,jj,kk),sp )
            vsij(1,3) = real( vs13(ll,jj,kk),sp )
            vsij(2,3) = real( vs23(ll,jj,kk),sp )
     
            ! Pressure Hessian
            pij(1,1) = real( p11(ll,jj,kk),sp )
            pij(1,2) = real( p12(ll,jj,kk),sp )
            pij(1,3) = real( p13(ll,jj,kk),sp )
            pij(2,2) = real( p22(ll,jj,kk),sp )
            pij(2,3) = real( p23(ll,jj,kk),sp )
            pij(3,3) = real( p33(ll,jj,kk, ifile),sp )

            vort(1) = real( ox(ll,jj,kk),sp )
            vort(2) = real( oy(ll,jj,kk),sp )
            vort(3) = real( oz(ll,jj,kk),sp )

            cij(1,1) = real( b11(ll,jj,kk,ifile),sp )
            cij(1,2) = real( b12(ll,jj,kk,ifile),sp )
            cij(1,3) = real( b13(ll,jj,kk,ifile),sp )
            cij(2,1) = real( b21(ll,jj,kk,ifile),sp )
            cij(2,2) = real( b22(ll,jj,kk,ifile),sp )
            cij(2,3) = real( b23(ll,jj,kk,ifile),sp )
            cij(3,1) = real( b31(ll,jj,kk,ifile),sp )
            cij(3,2) = real( b32(ll,jj,kk,ifile),sp )
            cij(3,3) = real( b33(ll,jj,kk,ifile),sp )
        else
            ll = ii/2

            sij(1,1) = aimag( s11(ll,jj,kk,ifile) )
            sij(2,2) = aimag( s22(ll,jj,kk,ifile) )
            sij(3,3) = aimag( s33(ll,jj,kk,ifile) )
            sij(1,2) = aimag( s12(ll,jj,kk) )
            sij(1,3) = aimag( s13(ll,jj,kk) )
            sij(2,3) = aimag( s23(ll,jj,kk) )
     
            vsij(1,1) = aimag( vs11(ll,jj,kk) )
            vsij(2,2) = aimag( vs22(ll,jj,kk) )
            vsij(3,3) = aimag( vs33(ll,jj,kk) )
            vsij(1,2) = aimag( vs12(ll,jj,kk) )
            vsij(1,3) = aimag( vs13(ll,jj,kk) )
            vsij(2,3) = aimag( vs23(ll,jj,kk) )
     
            ! Pressure Hessian
            pij(1,1) = aimag( p11(ll,jj,kk) )
            pij(1,2) = aimag( p12(ll,jj,kk) )
            pij(1,3) = aimag( p13(ll,jj,kk) )
            pij(2,2) = aimag( p22(ll,jj,kk) )
            pij(2,3) = aimag( p23(ll,jj,kk) )
            pij(3,3) = aimag( p33(ll,jj,kk, ifile) )

            vort(1) = aimag( ox(ll,jj,kk) )
            vort(2) = aimag( oy(ll,jj,kk) )
            vort(3) = aimag( oz(ll,jj,kk) )

            cij(1,1) = aimag( b11(ll,jj,kk,ifile) )
            cij(1,2) = aimag( b12(ll,jj,kk,ifile) )
            cij(1,3) = aimag( b13(ll,jj,kk,ifile) )
            cij(2,1) = aimag( b21(ll,jj,kk,ifile) )
            cij(2,2) = aimag( b22(ll,jj,kk,ifile) )
            cij(2,3) = aimag( b23(ll,jj,kk,ifile) )
            cij(3,1) = aimag( b31(ll,jj,kk,ifile) )
            cij(3,2) = aimag( b32(ll,jj,kk,ifile) )
            cij(3,3) = aimag( b33(ll,jj,kk,ifile) )
        endif
        sij(2,1) = sij(1,2); sij(3,1) = sij(1,3); sij(3,2) = sij(2,3)
        pij(2,1) = pij(1,2); pij(3,1) = pij(1,3); pij(3,2) = pij(2,3)
        vsij(2,1) = vsij(1,2); vsij(3,1) = vsij(1,3); vsij(3,2) = vsij(2,3)
        cij = matmul(cij, transpose(cij))

       
        call dsyevr("V", "A", "U", evsize, cij, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                    nfound, evcij, evtrcij, evsize, isuppz, work, lwork, iwork, liwork, info)
                ! cij will be destroyed on return
        if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cij. info = ', info
     
        ! Reverse the order of the eigen values and vectors:
        evcij = evcij(evsize:1:-1)
        evtrcij = evtrcij(:,evsize:1:-1)

        ! define a right-handed frame: gamma = alfa X beta
        evtrcij(1,3) = evtrcij(2,1) * evtrcij(3,2) - evtrcij(3,1) * evtrcij(2,2)
        evtrcij(2,3) = evtrcij(3,1) * evtrcij(1,2) - evtrcij(1,1) * evtrcij(3,2)
        evtrcij(3,3) = evtrcij(1,1) * evtrcij(2,2) - evtrcij(2,1) * evtrcij(1,2)
       
        ! vorticity contribution in Lab frame
        OmInCijLab = vort * .5_sp
 
        ! strain contribution in C-frame
        straincij(1) = dot_product(evtrcij(:,2), matmul(sij, evtrcij(:,3))) &
                      * ( evcij(2) + evcij(3) ) / ( evcij(2) - evcij(3) )
        straincij(2) = dot_product(evtrcij(:,1), matmul(sij, evtrcij(:,3))) &
                      * ( evcij(3) + evcij(1) ) / ( evcij(3) - evcij(1) )
        straincij(3) = dot_product(evtrcij(:,1), matmul(sij, evtrcij(:,2))) &
                      * ( evcij(1) + evcij(2) ) / ( evcij(1) - evcij(2) ) 

        ! strain contribution in Lab-frame
        SijInCijLab = matmul(evtrcij, straincij)


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

        ! evtrsij(:,2) cross evtrcij(:,1)
        evcross(1) = evtrsij(2,2) * evtrcij(3,1) - evtrsij(3,2) * evtrcij(2,1)
        evcross(2) = evtrsij(3,2) * evtrcij(1,1) - evtrsij(1,2) * evtrcij(3,1)
        evcross(3) = evtrsij(1,2) * evtrcij(2,1) - evtrsij(2,2) * evtrcij(1,1)

        scalign = dot_product(evtrsij(:,2), evtrcij(:,1))
        signalign = scalign / abs(scalign)
        scalign = abs(scalign)

        srot1 = dot_product(VortContSijLab, evcross) * signalign
        srot2 = dot_product(PijContSijLab,  evcross) * signalign
        srot3 = dot_product(VijContSijLab,  evcross) * signalign

        crot1 = - dot_product(OmInCijLab,  evcross) * signalign
        crot2 = - dot_product(SijInCijLab, evcross) * signalign

        meanscalign(ifile) = meanscalign(ifile) + scalign
        meansrot1(ifile) = meansrot1(ifile) + srot1
        meansrot2(ifile) = meansrot2(ifile) + srot2
        meansrot3(ifile) = meansrot3(ifile) + srot3
        meancrot1(ifile) = meancrot1(ifile) + crot1
        meancrot2(ifile) = meancrot2(ifile) + crot2

        pp = floor( scalign / binw ) + 1
        pdfalign(pp,ifile) = pdfalign(pp,ifile) + 1
        cndsrot1(pp,ifile) = cndsrot1(pp,ifile) + srot1
        cndsrot2(pp,ifile) = cndsrot2(pp,ifile) + srot2
        cndsrot3(pp,ifile) = cndsrot3(pp,ifile) + srot3
        cndcrot1(pp,ifile) = cndcrot1(pp,ifile) + crot1
        cndcrot2(pp,ifile) = cndcrot2(pp,ifile) + crot2
      end do
      end do
      end do

      meanscalign(ifile) = meanscalign(ifile) * const
      meansrot1(ifile)   = meansrot1(ifile)   * const
      meansrot2(ifile)   = meansrot2(ifile)   * const
      meansrot3(ifile)   = meansrot3(ifile)   * const
      meancrot1(ifile)   = meancrot1(ifile)   * const
      meancrot2(ifile)   = meancrot2(ifile)   * const

      cndsrot1(:,ifile) = cndsrot1(:,ifile) / (pdfalign(:,ifile) + mytiny)
      cndsrot2(:,ifile) = cndsrot2(:,ifile) / (pdfalign(:,ifile) + mytiny)
      cndsrot3(:,ifile) = cndsrot3(:,ifile) / (pdfalign(:,ifile) + mytiny)
      cndcrot1(:,ifile) = cndcrot1(:,ifile) / (pdfalign(:,ifile) + mytiny)
      cndcrot2(:,ifile) = cndcrot2(:,ifile) / (pdfalign(:,ifile) + mytiny)

      pdfalign(:,ifile) = pdfalign(:,ifile) * const / binw

    end do

    do ifile = 1, filecnt

      write(26, '(I6, 15E15.3)') nfile, meanscalign(ifile), meansrot1(ifile), meansrot2(ifile), &
                               meansrot3(ifile), meancrot1(ifile), meancrot2(ifile)

      write(27, '("# index = ", I6)') nfile
      do pp = 1, npnt
        write(27, '(15E15.3)') binw * (pp-.5_sp), pdfalign(pp,ifile), cndsrot1(pp,ifile), &
            cndsrot2(pp,ifile), cndsrot3(pp,ifile), cndcrot1(pp,ifile), cndcrot2(pp,ifile)
      end do
      write(27,*)
      write(27,*)

      nfile = nfile + 1
    end do

  end do
  close(30)
  close(27)
  close(26)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,s11,s22,s33,s12,s13,s23)
  deallocate(k2,ox,oy,oz,p11,p12,p13,p22,p23,p33,vs11,vs22,vs33,vs12,vs13,vs23)
  
  call destroyplan3d

  write(*,*) 'defm-alignsc.x done.'

end program defmsijincij
