module mfields
  use mconstant
  include 'fftw3.f'

  integer :: nx,ny,nz,lx1,lx,ly,lz

  complex(sp), allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33,omx,omy,omz,visx,visy,visz
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz, k2

  integer(8) :: IFTs11, IFTs12, IFTs13, IFTs22, IFTs23, IFTs33, IFTomx, IFTomy, IFTomz
  integer(8) :: IFTvisx, IFTvisy, IFTvisz

contains

  subroutine initfields

    allocate(s11(lx1,ly,lz), s12(lx1,ly,lz), s13(lx1,ly,lz), s22(lx1,ly,lz), s23(lx1,ly,lz), s33(lx1,ly,lz))
    allocate(omx(lx1,ly,lz), omy(lx1,ly,lz), omz(lx1,ly,lz))
    allocate(visx(lx1,ly,lz), visy(lx1,ly,lz), visz(lx1,ly,lz))
    allocate(  kx(lx1,ly,lz),  ky(lx1,ly,lz),  kz(lx1,ly,lz), k2(lx1,ly,lz) )


    call dfftw_plan_dft_c2r_3d(IFTs11,nx,ny,nz,s11,s11,FFTW_MEASURE+FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_3d(IFTs12,nx,ny,nz,s12,s12,FFTW_MEASURE+FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_3d(IFTs13,nx,ny,nz,s13,s13,FFTW_MEASURE+FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_3d(IFTs22,nx,ny,nz,s22,s22,FFTW_MEASURE+FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_3d(IFTs23,nx,ny,nz,s23,s23,FFTW_MEASURE+FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_3d(IFTs33,nx,ny,nz,s33,s33,FFTW_MEASURE+FFTW_UNALIGNED)

    call dfftw_plan_dft_c2r_3d(IFTomx,nx,ny,nz,omx,omx,FFTW_MEASURE+FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_3d(IFTomy,nx,ny,nz,omy,omy,FFTW_MEASURE+FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_3d(IFTomz,nx,ny,nz,omz,omz,FFTW_MEASURE+FFTW_UNALIGNED)

    call dfftw_plan_dft_c2r_3d(IFTvisx,nx,ny,nz,visx,visx,FFTW_MEASURE+FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_3d(IFTvisy,nx,ny,nz,visy,visy,FFTW_MEASURE+FFTW_UNALIGNED)
    call dfftw_plan_dft_c2r_3d(IFTvisz,nx,ny,nz,visz,visz,FFTW_MEASURE+FFTW_UNALIGNED)

  end subroutine initfields

  subroutine fftwsij
    call dfftw_execute_dft_c2r(IFTs11,s11,s11)
    call dfftw_execute_dft_c2r(IFTs12,s12,s12)
    call dfftw_execute_dft_c2r(IFTs13,s13,s13)
    call dfftw_execute_dft_c2r(IFTs22,s22,s22)
    call dfftw_execute_dft_c2r(IFTs23,s23,s23)
    call dfftw_execute_dft_c2r(IFTs33,s33,s33)
  end subroutine fftwsij

  subroutine fftwom
    call dfftw_execute_dft_c2r(IFTomx,omx,omx)
    call dfftw_execute_dft_c2r(IFTomy,omy,omy)
    call dfftw_execute_dft_c2r(IFTomz,omz,omz)
  end subroutine fftwom

  subroutine fftwvis
    call dfftw_execute_dft_c2r(IFTvisx,visx,visx)
    call dfftw_execute_dft_c2r(IFTvisy,visy,visy)
    call dfftw_execute_dft_c2r(IFTvisz,visz,visz)
  end subroutine fftwvis

  subroutine cleanupfield
    call dfftw_destroy_plan(IFTs11)
    call dfftw_destroy_plan(IFTs12)
    call dfftw_destroy_plan(IFTs13)
    call dfftw_destroy_plan(IFTs22)
    call dfftw_destroy_plan(IFTs23)
    call dfftw_destroy_plan(IFTs33)

    call dfftw_destroy_plan(IFTomx)
    call dfftw_destroy_plan(IFTomy)
    call dfftw_destroy_plan(IFTomz)

    call dfftw_destroy_plan(IFTvisx)
    call dfftw_destroy_plan(IFTvisy)
    call dfftw_destroy_plan(IFTvisz)

    deallocate(s11,s12,s13,s22,s23,s33,omx,omy,omz)
    deallocate(kx,ky,kz,k2,visx,visy,visz)
  end subroutine cleanupfield

end module mfields

program vortalignpdfeulerian
  use, intrinsic :: iso_c_binding
  use mconstant
  use mfields
  use mwavenumber
  implicit none
  
  integer :: ii,jj,kk,ll,mm,nn,pp,npnt,rr
  real(sp) :: const 

  real(sp), parameter :: rnu = 0.006

  integer, parameter :: naverage = 500
  integer, dimension(naverage) :: xx, yy, zz
  real(sp) :: xxll, yyll, zzll, xpp, ypp, zpp, xxrr, yyrr, zzrr, dx, dy, dz
  integer  :: iiloc, jjloc, kkloc

  ! For finding the PDF of relative rotation rate
  integer,  parameter :: npdf = 80
  real(sp), parameter :: bnd =8._sp, binw = 2*bnd/npdf

  real(sp) :: ran1
  external ran1

  integer, parameter :: evsize = 3
  real(dp), dimension(evsize) :: rotvis, rotsst, rottot, omp, ecross
  real(dp), dimension(evsize, naverage) :: rotvis0, rotsst0, rottot0, omp0

  real(sp),    pointer,     dimension(:,:,:) :: r11,r12,r13,r22,r23,r33,rox,roy,roz,rvx,rvy,rvz
  type(C_ptr) :: cptr

  integer,  allocatable, dimension(:) :: pntcnt, allpnt
  real(dp), allocatable, dimension(:) :: mean_relrottot, vari_relrottot, skew_relrottot
  real(dp), allocatable, dimension(:) :: mean_relrotvis, vari_relrotvis, skew_relrotvis
  real(dp), allocatable, dimension(:) :: mean_relrotsst, vari_relrotsst, skew_relrotsst
  real(dp), allocatable, dimension(:,:) :: prelrottot, prelrotvis, prelrotsst
  real(dp) :: relrottot, relrotvis, relrotsst 
  real(dp) :: vari_tot, skew_tot, vari_vis, skew_vis, vari_sst, skew_sst

  character(80) :: str, flnm, strll
  
  if (iargc() .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./vort-relrotpdf-eulerian.x nx filelist'
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
  dx = 2 * pi / nx; dy = 2 * pi / ny; dz = 2 * pi / nz
  const = 1._sp/(nx*ny*nz)

  ! The max displacement for correlation is half of the size of the cube
  npnt = nx/2

  allocate( pntcnt(npnt), allpnt(npnt), prelrottot(npdf,npnt), prelrotvis(npdf,npnt), prelrotsst(npdf,npnt) )
  allocate( mean_relrottot(npnt), vari_relrottot(npnt), skew_relrottot(npnt) )
  allocate( mean_relrotvis(npnt), vari_relrotvis(npnt), skew_relrotvis(npnt) )
  allocate( mean_relrotsst(npnt), vari_relrotsst(npnt), skew_relrotsst(npnt) )

  call initfields
  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)

  ! centers of the groups
  ! will calculate the displacement from these centers
  ll = -100
  do mm = 1, naverage
 
    ii = floor( ran1(ll) * nx ) + 1
    jj = floor( ran1(ll) * ny ) + 1
    kk = floor( ran1(ll) * nz ) + 1
    xx(mm) = ii; yy(mm) = jj; zz(mm) = kk

  end do

  prelrottot = 0._dp; prelrotvis = 0._dp; prelrotsst = 0._dp
  mean_relrottot = 0._dp; mean_relrotvis = 0._dp; mean_relrotsst = 0._dp
  vari_relrottot = 0._dp; vari_relrotvis = 0._dp; vari_relrotsst = 0._dp
  skew_relrottot = 0._dp; skew_relrotvis = 0._dp; skew_relrotsst = 0._dp
  pntcnt = 0
  allpnt = 0
  open(30, file = trim(flnm)//'.list')
  do while ( .not. eof(30) )

    read(30,*) strll
    strll = adjustl(strll)
    write(*,*) 'reading file ', strll

    open(15,file='./out/ux'//trim(strll),form='unformatted')
      read(15) s11
    close(15)
    open(15,file='./out/uy'//trim(strll),form='unformatted')
      read(15) s22
    close(15)
    open(15,file='./out/uz'//trim(strll),form='unformatted')
      read(15) s33
    close(15)
 
    ! vorticity and strain rate tensor
    omx = eye * (ky * s33 - kz * s22)
    omy = eye * (kz * s11 - kx * s33)
    omz = eye * (kx * s22 - ky * s11)

    visx = - k2 * rnu * omx; visy = - k2 * rnu * omy; visz = -k2 * rnu * omz
 
    call fftwom
    call fftwvis

    ! magnitude of vorticity
    s12 = cmplx( sqrt(real(omx,sp) * real(omx,sp) + real(omy,sp) * real(omy,sp) + real(omz,sp) * real(omz,sp)), &
                 sqrt(aimag(omx) * aimag(omx) + aimag(omy) * aimag(omy) + aimag(omz) * aimag(omz)) )

    ! normalize the viscous diffusion term
    visx = cmplx( real(visx,sp) / real(s12,sp), aimag(visx) / aimag(s12) )
    visy = cmplx( real(visy,sp) / real(s12,sp), aimag(visy) / aimag(s12) )
    visz = cmplx( real(visz,sp) / real(s12,sp), aimag(visz) / aimag(s12) )

    ! normalize the vorticity: omx is now the x-component of the direction
    ! vector
    omx  = cmplx( real(omx,sp) / real(s12,sp), aimag(omx) / aimag(s12) )
    omy  = cmplx( real(omy,sp) / real(s12,sp), aimag(omy) / aimag(s12) )
    omz  = cmplx( real(omz,sp) / real(s12,sp), aimag(omz) / aimag(s12) )
  
    ! strain rate tensor in Fourier space
    s12 = eye * (kx * s22 + ky * s11)*0.5_sp
    s13 = eye * (kx * s33 + kz * s11)*0.5_sp
    s23 = eye * (ky * s33 + kz * s22)*0.5_sp
    s11 = eye * kx * s11
    s22 = eye * ky * s22
    s33 = eye * kz * s33
 
    call fftwsij

    ! strain rate tensor sij in real space
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
 
    cptr = c_loc(visx)
    call c_f_pointer(cptr, rvx, [2*lx1,ly,lz])
    cptr = c_loc(visy)
    call c_f_pointer(cptr, rvy, [2*lx1,ly,lz])
    cptr = c_loc(visz)
    call c_f_pointer(cptr, rvz, [2*lx1,ly,lz])
 
    ! RHS of the equation of the direction vector for the vorticity

    ! vortex-stretching term 
    r11 =  r11 * rox + r12 * roy + r13 * roz  
    r22 =  r12 * rox + r22 * roy + r23 * roz  
    r33 =  r13 * rox + r23 * roy + r33 * roz  

    ! projection onto a plane perpendicular to vorticity: viscous diffusion term 
    r12 = rvx * rox + rvy * roy + rvz * roz

    rvx = rvx - r12 * rox
    rvy = rvy - r12 * roy
    rvz = rvz - r12 * roz
    ! (rvx, rvy, rvz) is the viscous diffusion contribution for the rate of change of (rox, roy, roz)

    ! projection: vortex stretching term
    r12 = r11 * rox + r22 * roy + r33 * roz

    r11 = r11 - r12 * rox
    r22 = r22 - r12 * roy
    r33 = r33 - r12 * roz
    ! (r11, r22, r33) is the vortex-stretching contribution for the rate of change of (rox, roy, roz)

 
    do mm = 1, naverage
      call findRotOnPnt(xx(mm),yy(mm),zz(mm),rotvis0(:,mm),rotsst0(:,mm),rottot0(:,mm),omp0(:,mm))
    end do
 
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx
 
      xpp = dx * (ii-1)
      ypp = dy * (jj-1)
      zpp = dz * (kk-1)
 
      call findRotOnPnt(ii,jj,kk,rotvis,rotsst,rottot,omp)
 
      ! loop over different groups
      do ll = 1, naverage
 
        ecross(1) = omp0(2,ll) * omp(3) - omp0(3,ll) * omp(2)
        ecross(2) = omp0(3,ll) * omp(1) - omp0(1,ll) * omp(3)
        ecross(3) = omp0(1,ll) * omp(2) - omp0(2,ll) * omp(1)
       
        relrotvis = dot_product( ( rotvis0(:,ll) - rotvis ), ecross )
        relrotvis = relrotvis * sign( 1._dp, dot_product( omp0(:,ll), omp ) )
    
        relrotsst = dot_product( ( rotsst0(:,ll) - rotsst ), ecross )
        relrotsst = relrotsst * sign( 1._dp, dot_product( omp0(:,ll), omp ) )
    
        relrottot = dot_product( ( rottot0(:,ll) - rottot ), ecross )
        relrottot = relrottot * sign( 1._dp, dot_product( omp0(:,ll), omp ) )
    
        
        xxrr = (xx(ll) - 1) * dx; yyrr = (yy(ll) - 1) * dx; zzrr = (zz(ll) - 1) * dx
 
        xxrr = modulo(xxrr, 2*pi)
        yyrr = modulo(yyrr, 2*pi)
        zzrr = modulo(zzrr, 2*pi)
 
        do kkloc = -1, 1 ! checking the images of the flow field
        do jjloc = -1, 1
        do iiloc = -1, 1
 
          xxll = xxrr + iiloc * 2 * pi
          yyll = yyrr + jjloc * 2 * pi
          zzll = zzrr + kkloc * 2 * pi
 
          mm = floor( sqrt( (xpp-xxll)**2 + (ypp-yyll)**2 + (zpp-zzll)**2 ) / dx + 0.5_sp )  + 1
  
          if ( mm .ge. 1 .and. mm .le. npnt) then
 
              rr = floor( (relrottot + bnd) / binw ) + 1
              nn = floor( (relrotvis + bnd) / binw ) + 1
              pp = floor( (relrotsst + bnd) / binw ) + 1
              
              if ( (rr .ge. 1 .and. rr .le. npdf) .and. &
                   (nn .ge. 1 .and. nn .le. npdf) .and. & 
                   (pp .ge. 1 .and. pp .le. npdf) ) then

                prelrottot(rr,mm) = prelrottot(rr,mm) + 1
                prelrotvis(nn,mm) = prelrotvis(nn,mm) + 1
                prelrotsst(pp,mm) = prelrotsst(pp,mm) + 1

                mean_relrottot(mm) = mean_relrottot(mm) + relrottot
                vari_relrottot(mm) = vari_relrottot(mm) + relrottot * relrottot
                skew_relrottot(mm) = skew_relrottot(mm) + relrottot ** 3

                mean_relrotvis(mm) = mean_relrotvis(mm) + relrotvis
                vari_relrotvis(mm) = vari_relrotvis(mm) + relrotvis * relrotvis
                skew_relrotvis(mm) = skew_relrotvis(mm) + relrotvis ** 3

                mean_relrotsst(mm) = mean_relrotsst(mm) + relrotsst
                vari_relrotsst(mm) = vari_relrotsst(mm) + relrotsst * relrotsst
                skew_relrotsst(mm) = skew_relrotsst(mm) + relrotsst ** 3

                pntcnt(mm) = pntcnt(mm) + 1
              end if
              allpnt(mm) = allpnt(mm) + 1
            
          end if
 
        end do
        end do
        end do
  
      end do
 
    end do
    end do
    end do
  end do
  close(30)

  open(25, file = 'vort-relrotmean-eulerian-'//trim(flnm)//'.dat')
  open(26, file = 'vort-relrotpdf-eulerian-'//trim(flnm)//'.dat')

    do mm = 1, npnt

      prelrottot(:,mm) = prelrottot(:,mm) / (allpnt(mm) + mytiny) 
      prelrotvis(:,mm) = prelrotvis(:,mm) / (allpnt(mm) + mytiny) 
      prelrotsst(:,mm) = prelrotsst(:,mm) / (allpnt(mm) + mytiny) 
      write(*,*) 'vort-relrotpdf-eulerian.x: prelrottot normalization: ', sum(prelrottot(:,mm))
      write(*,*) 'vort-relrotpdf-eulerian.x: prelrotvis normalization: ', sum(prelrotvis(:,mm))
      write(*,*) 'vort-relrotpdf-eulerian.x: prelrotsst normalization: ', sum(prelrotsst(:,mm))

      prelrottot(:,mm) = prelrottot(:,mm) / binw
      prelrotvis(:,mm) = prelrotvis(:,mm) / binw
      prelrotsst(:,mm) = prelrotsst(:,mm) / binw

      mean_relrottot(mm) = mean_relrottot(mm) / pntcnt(mm)
      vari_relrottot(mm) = vari_relrottot(mm) / pntcnt(mm)
      skew_relrottot(mm) = skew_relrottot(mm) / pntcnt(mm)

      vari_tot = vari_relrottot(mm) - mean_relrottot(mm)**2
      skew_tot = skew_relrottot(mm) - 3 * vari_relrottot(mm) * mean_relrottot(mm) + &
             2 * mean_relrottot(mm) ** 3
      skew_tot = skew_tot / ( vari_tot + mytiny )**1.5_dp

      mean_relrotvis(mm) = mean_relrotvis(mm) / pntcnt(mm)
      vari_relrotvis(mm) = vari_relrotvis(mm) / pntcnt(mm)
      skew_relrotvis(mm) = skew_relrotvis(mm) / pntcnt(mm)

      vari_vis = vari_relrotvis(mm) - mean_relrotvis(mm)**2
      skew_vis = skew_relrotvis(mm) - 3 * vari_relrotvis(mm) * mean_relrotvis(mm) + &
             2 * mean_relrotvis(mm) ** 3
      skew_vis = skew_vis / ( vari_vis + mytiny )**1.5_dp

      mean_relrotsst(mm) = mean_relrotsst(mm) / pntcnt(mm)
      vari_relrotsst(mm) = vari_relrotsst(mm) / pntcnt(mm)
      skew_relrotsst(mm) = skew_relrotsst(mm) / pntcnt(mm)

      vari_sst = vari_relrotsst(mm) - mean_relrotsst(mm)**2
      skew_sst = skew_relrotsst(mm) - 3 * vari_relrotsst(mm) * mean_relrotsst(mm) + &
             2 * mean_relrotsst(mm) ** 3
      skew_sst = skew_sst / ( vari_sst + mytiny )**1.5_dp

      write(25, '(15E13.3)') real(mm), mean_relrottot(mm), vari_tot, skew_tot, &
                                       mean_relrotvis(mm), vari_vis, skew_vis, &
                                       mean_relrotsst(mm), vari_sst, skew_sst

      write(26, '( "#index = ", I6)') mm
      do rr = 1, npdf
        write(26, '(15E15.3)') -bnd + (rr-0.5_sp)*binw, prelrottot(rr,mm), prelrotvis(rr,mm), &
        prelrotsst(rr,mm)
      end do
      write(26, *)
      write(26, *)

    end do

  close(26)
  close(25)

  deallocate(prelrottot, prelrotvis, prelrotsst, pntcnt, allpnt) 
  deallocate(mean_relrottot, vari_relrottot, skew_relrottot)
  deallocate(mean_relrotvis, vari_relrotvis, skew_relrotvis)
  deallocate(mean_relrotsst, vari_relrotsst, skew_relrotsst)

  call cleanupfield
  
  write(*,*) 'vort-relrotpdf-eulerian.x done.'

contains

  subroutine findRotOnPnt(ipntx, ipnty, ipntz, rotvis, rotsst,rottot,omp)

    integer, intent(in) :: ipntx, ipnty, ipntz
    real(dp), dimension(evsize), intent(out) :: rotvis, rotsst, rottot, omp

    real(dp), dimension(evsize) :: vis, sst

    vis(1) = rvx(ipntx,ipnty,ipntz)
    vis(2) = rvy(ipntx,ipnty,ipntz)
    vis(3) = rvz(ipntx,ipnty,ipntz)

    sst(1) = r11(ipntx,ipnty,ipntz)
    sst(2) = r22(ipntx,ipnty,ipntz)
    sst(3) = r33(ipntx,ipnty,ipntz)
    
    omp(1) = rox(ipntx,ipnty,ipntz)
    omp(2) = roy(ipntx,ipnty,ipntz)
    omp(3) = roz(ipntx,ipnty,ipntz)

    rotvis(1) = omp(2) * vis(3) - omp(3) * vis(2)
    rotvis(2) = omp(3) * vis(1) - omp(1) * vis(3)
    rotvis(3) = omp(1) * vis(2) - omp(2) * vis(1)

    rotsst(1) = omp(2) * sst(3) - omp(3) * sst(2)
    rotsst(2) = omp(3) * sst(1) - omp(1) * sst(3)
    rotsst(3) = omp(1) * sst(2) - omp(2) * sst(1)

    rottot = rotvis + rotsst

  end subroutine findRotOnPnt

end program vortalignpdfeulerian
