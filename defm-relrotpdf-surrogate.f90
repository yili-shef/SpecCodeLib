program defmalignpdfsurrogate
  use, intrinsic :: iso_c_binding
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,npnt,rr
  real(sp) :: const, ignore_me

  real(dp), parameter :: walfa = 1._dp, wbeta = -1._dp, wgmma = 1._dp

  ! For interpolation
  integer, parameter :: norder = 6
  real(sp) :: dxyz(3), bg(norder,3)
  integer  :: ix(norder), iy(norder), iz(norder), lhnode(3)

  integer, parameter :: naverage = 400
  integer, dimension(naverage) :: xx, yy, zz
  real(sp) :: xxll, yyll, zzll, xpp, ypp, zpp, xxrr, yyrr, zzrr, dx
  integer  :: iiloc, jjloc, kkloc

  ! For finding the PDF of relative rotation rate
  integer, parameter :: npdf = 60
  real(sp), parameter :: bnd =15._sp, binw = 2*bnd/npdf

  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  ! ----------------------------------------------------------------------------------
  real(sp) :: ran1
  external ran1

  real(dp), dimension(evsize) :: rotInLab, ecross
  real(dp), dimension(evsize, evsize) :: axescij
  real(dp), dimension(evsize, naverage) :: rotInLab0
  real(dp), dimension(evsize, evsize, naverage) :: axescij0

  integer :: noffset, forward, nlst

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33,omx,omy,omz
  real(sp),    pointer,     dimension(:,:,:) :: r11,r12,r13,r22,r23,r33,rox,roy,roz
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz
  type(C_ptr) :: cptr

  real(sp), allocatable, dimension(:) :: xp, yp, zp
  integer,  allocatable, dimension(:) :: pntcnt, allpnt
  real(dp), allocatable, dimension(:,:) :: prelrot
  real(dp) :: relrot, var, skew
  real(dp), allocatable, dimension(:) :: meanrelrot, varrelrot, skewrelrot

  character(80) :: str, flnm, prf, str1, strll, mark, xyzlst
  
  if (iargc() .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-relrotpdf-surrogate.x nx filelist prefix noffset forward'
          write(*,*) '                     nx: resolution of data'
          write(*,*) '                     filelist: data file list: *.list'
          write(*,*) '                     prefix: prefix for datafiles for bij'
          write(*,*) '                     noffset: offset in file number in the list'
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

  call getarg(4,str)
  read(str, '(I20)') noffset

  call getarg(5,mark)
  read(mark, '(I20)') forward


  ny=nx; nz=nx
  lx=nx/2; lx1=nx/2+1; ly=ny; lz=nz
  const = 1._sp/(nx*ny*nz)

  dx = 2*pi/nx
  dxyz = dx

  ! The max displacement for correlation is half of the size of the cube
  npnt = nx/2

  allocate( pntcnt(npnt), prelrot(npdf,npnt) )
  allocate( allpnt(npnt) )
  allocate( meanrelrot(npnt), varrelrot(npnt), skewrelrot(npnt) )

  allocate( b11(lx1,ly,lz), b12(lx1,ly,lz), b13(lx1,ly,lz) )
  allocate( b21(lx1,ly,lz), b22(lx1,ly,lz), b23(lx1,ly,lz) )
  allocate( b31(lx1,ly,lz), b32(lx1,ly,lz), b33(lx1,ly,lz) )

  allocate( s11(lx1,ly,lz), s12(lx1,ly,lz), s13(lx1,ly,lz) )
  allocate( omx(lx1,ly,lz), s22(lx1,ly,lz), s23(lx1,ly,lz) )
  allocate( omy(lx1,ly,lz), omz(lx1,ly,lz), s33(lx1,ly,lz) )
  allocate(  kx(lx1,ly,lz),  ky(lx1,ly,lz),  kz(lx1,ly,lz) )

  allocate( xp(nx*ny*nz), yp(nx*ny*nz), zp(nx*ny*nz) )

  call fftwplan3de(nx,ny,nz)
  call wavenumber(kx,ky,kz,lx1,ly,lz)

  !-------------------------------------------------------------------------------------
  ! Will calculate the pdf for the "noffset" data file in each list only.
  ! For this we need to generate the file names instead of reading them from the lists. 
  ! Will use the smartrunner(.py) to run the code, so the interface is kept the same.

  xyzlst = flnm( len_trim(flnm) - 3 : len_trim(flnm) ) ! extract the starting file number
  read(xyzlst, "(I20)") nlst

  if (forward .eq. 0) then 
      ll = nlst - noffset
      nlst = nlst - 150 
      mark = 'bw'
  else
      ll = nlst + noffset
      nlst = nlst + 150
      mark = 'fw'
  end if

  write(str1, "(I20)") nlst ! get the last file number
  str1 = adjustl(str1)

  write(strll, "(I20)") ll
  strll = adjustl(strll)


  xyzlst = trim(mark)//trim(xyzlst)//"to"//trim(str1)//'-'//trim(strll)//'.dat'
  flnm = trim(flnm)//'-'//trim(strll)//'.dat'
  !--------------------------------------------------------------

  open(15,file='./out/ux'//trim(strll)//'.dat',form='unformatted')
    read(15) s11
  close(15)
  open(15,file='./out/uy'//trim(strll)//'.dat',form='unformatted')
    read(15) s22
  close(15)
  open(15,file='./out/uz'//trim(strll)//'.dat',form='unformatted')
    read(15) s33
  close(15)

  open(15,file='./out/b11'//trim(flnm),form='unformatted')
    read(15) b11
  close(15)
  open(15,file='./out/b12'//trim(flnm),form='unformatted')
    read(15) b12
  close(15)
  open(15,file='./out/b13'//trim(flnm),form='unformatted')
    read(15) b13
  close(15)
  open(15,file='./out/b21'//trim(flnm),form='unformatted')
    read(15) b21
  close(15)
  open(15,file='./out/b22'//trim(flnm),form='unformatted')
    read(15) b22
  close(15)
  open(15,file='./out/b23'//trim(flnm),form='unformatted')
    read(15) b23
  close(15)
  open(15,file='./out/b31'//trim(flnm),form='unformatted')
    read(15) b31
  close(15)
  open(15,file='./out/b32'//trim(flnm),form='unformatted')
    read(15) b32
  close(15)
  open(15,file='./out/b33'//trim(flnm),form='unformatted')
    read(15) b33
  close(15)
  open(15,file='./out/xyz'//trim(xyzlst),form='unformatted')
    read(15) xp, yp, zp
  close(15)

  omx = eye * (ky * s33 - kz * s22)
  omy = eye * (kz * s11 - kx * s33)
  omz = eye * (kx * s22 - ky * s11)
  call rfftwnd_f77_one_complex_to_real(c2r3d,omx,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,omy,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,omz,ignore_me)
 
  s12 = eye * (kx * s22 + ky * s11)*0.5_sp
  s13 = eye * (kx * s33 + kz * s11)*0.5_sp
  s23 = eye * (ky * s33 + kz * s22)*0.5_sp
  s11 = eye * kx * s11
  s22 = eye * ky * s22
  s33 = eye * kz * s33

  call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
  call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
 
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


  ll = -100
  do mm = 1, naverage

    ii = floor( ran1(ll) * nx ) + 1
    jj = floor( ran1(ll) * ny ) + 1
    kk = floor( ran1(ll) * nz ) + 1

    xx(mm) = ii; yy(mm) = jj; zz(mm) = kk

    call findRotOnPnt(ii,jj,kk,rotInLab0(:,mm), axescij0(:,:,mm))

  end do


  prelrot = 0._dp
  meanrelrot = 0._dp
  varrelrot = 0._dp
  skewrelrot = 0._dp
  pntcnt = 0
  allpnt = 0
  do kk = 1, nz
  do jj = 1, ny
  do ii = 1, nx

    mm = (kk-1) * nx * ny + (jj-1) * nx + ii
    xpp = xp(mm); ypp = yp(mm); zpp = zp(mm)
    xpp = modulo(xpp, 2*pi)
    ypp = modulo(ypp, 2*pi)
    zpp = modulo(zpp, 2*pi)

    call findRotOnPnt(ii,jj,kk,rotInLab,axescij)

    ! loop over different groups
    do ll = 1, naverage

      ecross(1) = axescij0(2,1,ll) * axescij(3,1) - axescij0(3,1,ll) * axescij(2,1)
      ecross(2) = axescij0(3,1,ll) * axescij(1,1) - axescij0(1,1,ll) * axescij(3,1)
      ecross(3) = axescij0(1,1,ll) * axescij(2,1) - axescij0(2,1,ll) * axescij(1,1)
     
      relrot = dot_product( ( rotInLab0(:,ll) - rotInLab ), ecross )
      relrot = relrot * sign( 1._dp, dot_product( axescij0(:,1,ll), axescij(:,1) ) )
  
      rr = ( zz(ll) - 1 ) * nx * ny + ( yy(ll) - 1 ) * nx + xx(ll)
      xxrr = xp(rr); yyrr = yp(rr); zzrr = zp(rr)

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
      end do
      end do
 
    end do

  end do
  end do
  end do

  open(25, file = 'defm-relrotmean-surrogate-'//trim(flnm))
  open(26, file = 'defm-relrotpdf-surrogate-'//trim(flnm))

    do mm = 1, npnt

      prelrot(:,mm) = prelrot(:,mm) / (allpnt(mm) + mytiny) 
      write(*,*) 'defm-relrotpdf-surrogate.x: prelrot normalization: ', sum(prelrot(:,mm))

      prelrot(:,mm) = prelrot(:,mm) / binw
      meanrelrot(mm) = meanrelrot(mm) / pntcnt(mm)
      varrelrot(mm)  = varrelrot(mm)  / pntcnt(mm)
      skewrelrot(mm) = skewrelrot(mm) / pntcnt(mm)

      var = varrelrot(mm) - meanrelrot(mm)**2
      skew = skewrelrot(mm) - 3 * varrelrot(mm) * meanrelrot(mm) + 2 * meanrelrot(mm) ** 3
      skew = skew / ( var + mytiny )**1.5_dp

      write(25, '(15E13.3)') real(mm), real(allpnt(mm)), meanrelrot(mm), &
        varrelrot(mm), skewrelrot(mm), var, skew

      write(26, '( "#index = ", I6)') mm
      do rr = 1, npdf
        write(26, '(15E15.3)') -bnd + (rr-0.5_sp)*binw, prelrot(rr,mm)
      end do
      write(26, *)
      write(26, *)

    end do

  close(26)
  close(25)

  deallocate(b11,b12,b13,b21,b22,b23,b31,b32,b33)
  deallocate(s11,s12,s13,s22,s23,s33,omx,omy,omz)
  deallocate(prelrot, pntcnt, allpnt, meanrelrot, varrelrot, skewrelrot)
  deallocate(xp, yp, zp, kx, ky, kz)

  call destroyplan3d
  
  write(*,*) 'defm-relrotpdf-surrogate.x done.'

contains

  subroutine findRotOnPnt(ipntx, ipnty, ipntz, rotInLabLocal, evtrcc)

    integer, intent(in) :: ipntx, ipnty, ipntz
    real(dp), dimension(evsize), intent(out) :: rotInLabLocal
    real(dp), dimension(evsize, evsize), intent(out) :: evtrcc

    real(dp), dimension(evsize) :: ominc, rotinc, evcc, y, omp
    real(dp), dimension(evsize, evsize) :: cc, sij
    integer :: pp, ix1, iy1, iz1, iiloc, jjloc, kkloc

    if ( mod(ipntx, 2) .eq. 0 ) then
        cc(1,1) = aimag( b11(ipntx/2,ipnty,ipntz) )
        cc(1,2) = aimag( b12(ipntx/2,ipnty,ipntz) )
        cc(1,3) = aimag( b13(ipntx/2,ipnty,ipntz) )
        cc(2,1) = aimag( b21(ipntx/2,ipnty,ipntz) )
        cc(2,2) = aimag( b22(ipntx/2,ipnty,ipntz) )
        cc(2,3) = aimag( b23(ipntx/2,ipnty,ipntz) )
        cc(3,1) = aimag( b31(ipntx/2,ipnty,ipntz) )
        cc(3,2) = aimag( b32(ipntx/2,ipnty,ipntz) )
        cc(3,3) = aimag( b33(ipntx/2,ipnty,ipntz) )
    else                             
        cc(1,1) = real( b11(ipntx/2+1,ipnty,ipntz) )
        cc(1,2) = real( b12(ipntx/2+1,ipnty,ipntz) )
        cc(1,3) = real( b13(ipntx/2+1,ipnty,ipntz) )
        cc(2,1) = real( b21(ipntx/2+1,ipnty,ipntz) )
        cc(2,2) = real( b22(ipntx/2+1,ipnty,ipntz) )
        cc(2,3) = real( b23(ipntx/2+1,ipnty,ipntz) )
        cc(3,1) = real( b31(ipntx/2+1,ipnty,ipntz) )
        cc(3,2) = real( b32(ipntx/2+1,ipnty,ipntz) )
        cc(3,3) = real( b33(ipntx/2+1,ipnty,ipntz) )
    end if
    cc = matmul(cc, transpose(cc))
 
    pp = (ipntz-1) * nx * ny + (ipnty-1) * nx + ipntx
    y(1) = xp(pp); y(2) = yp(pp); y(3) = zp(pp)
 
    call pre_interp(y,dxyz,bg,lhnode)
    call ixiyiz(lhnode,nx,ny,nz,ix,iy,iz)
 
    sij = 0._dp; omp = 0._dp
    do  kkloc = 1, norder
      iz1 = iz(kkloc)
      do  jjloc = 1, norder
        iy1 = iy(jjloc)
        do  iiloc = 1, norder
          ix1 = ix(iiloc)
    
          ignore_me = bg(iiloc,1)*bg(jjloc,2)*bg(kkloc,3)
    
          sij(1,1) = sij(1,1) + r11(ix1,iy1,iz1) * ignore_me
          sij(1,2) = sij(1,2) + r12(ix1,iy1,iz1) * ignore_me
          sij(1,3) = sij(1,3) + r13(ix1,iy1,iz1) * ignore_me
    
          sij(2,2) = sij(2,2) + r22(ix1,iy1,iz1) * ignore_me
          sij(2,3) = sij(2,3) + r23(ix1,iy1,iz1) * ignore_me
          sij(3,3) = sij(3,3) + r33(ix1,iy1,iz1) * ignore_me
    
          omp(1) = omp(1) + rox(ix1,iy1,iz1) * ignore_me
          omp(2) = omp(2) + roy(ix1,iy1,iz1) * ignore_me
          omp(3) = omp(3) + roz(ix1,iy1,iz1) * ignore_me
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
   
    ominc(1) = dot_product( omp, evtrcc(:,1) ) * .5_dp
    ominc(2) = dot_product( omp, evtrcc(:,2) ) * .5_dp
    ominc(3) = dot_product( omp, evtrcc(:,3) ) * .5_dp
 
    rotinc(1) = ominc(1) + dot_product(evtrcc(:,2), matmul(sij, evtrcc(:,3))) * walfa
    rotinc(2) = ominc(2) + dot_product(evtrcc(:,3), matmul(sij, evtrcc(:,1))) * wbeta
    rotinc(3) = ominc(3) + dot_product(evtrcc(:,1), matmul(sij, evtrcc(:,2))) * wgmma
     
    rotInLabLocal = matmul(evtrcc, rotinc)

  end subroutine findRotOnPnt

end program defmalignpdfsurrogate
