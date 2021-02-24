program defmpdfome
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt = 20
  real(sp), parameter :: bnd = 1._sp, binw = 2._sp* bnd/npnt

  integer,  parameter :: npome = 150
  real(sp), parameter :: bndome = 30._sp, bwome = bndome/npome


  real(dp), dimension(npnt) :: pdfalign, pdfalignomemom
  real(dp), dimension(npome) :: pdfome, pdfovec
  real(dp), dimension(3,3)  :: cc, aa, bb
  real(dp), dimension(3)    :: iangmom, ovec, ome
  real(dp) :: detbij, alignomemom
  real(dp) :: ome2norm, ovec2norm, align
  integer :: ierr, pointcnt

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: a11,a12,a13
  complex(sp), allocatable, dimension(:,:,:) :: a21,a22,a23
  complex(sp), allocatable, dimension(:,:,:) :: a31,a32,a33
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1

  external :: dpotrf, dpotrs ! subroutines from lapack in MKL
  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-pdf-ome.x nx filelist prefix'
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

  allocate( b11(lx1,ly,lz), b12(lx1,ly,lz), b13(lx1,ly,lz) )
  allocate( b21(lx1,ly,lz), b22(lx1,ly,lz), b23(lx1,ly,lz) )
  allocate( b31(lx1,ly,lz), b32(lx1,ly,lz), b33(lx1,ly,lz) )
  allocate(  kx(lx1,ly,lz),  ky(lx1,ly,lz),  kz(lx1,ly,lz) )
  allocate( a11(lx1,ly,lz), a12(lx1,ly,lz), a13(lx1,ly,lz) )
  allocate( a21(lx1,ly,lz), a22(lx1,ly,lz), a23(lx1,ly,lz) )
  allocate( a31(lx1,ly,lz), a32(lx1,ly,lz), a33(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(27, file = 'pdf-ome-align-'//trim(flnm)//'.dat')
  write(27, '("# cosine  pdfalignomeomega  pdfalignomemom")') 
  open(28, file = 'pdf-ome-size-'//trim(flnm)//'.dat')
  write(28, '("# ome/ovec  pdfome  pdfovec")') 

  open(30, file = flnm(1:len_trim(flnm))//'.list')

  nfile = 0
  do while ( .not. eof(30) )
    read(30,*) str1
    write(*,*) 'data file ', str1(1:len_trim(str1))

    open(15,file='./out/ux'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b11
    close(15)
    open(15,file='./out/uy'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b22
    close(15)
    open(15,file='./out/uz'//str1(1:len_trim(str1)),form='unformatted')
      read(15) b33
    close(15)
    write(*,*) 'finishing reading ux uy uz'

    ! velocity gradients
    a11 = eye * kx * b11; a12 = eye * ky * b11; a13 = eye * kz * b11
    a21 = eye * kx * b22; a22 = eye * ky * b22; a23 = eye * kz * b22
    a31 = eye * kx * b33; a32 = eye * ky * b33; a33 = eye * kz * b33 

    call rfftwnd_f77_one_complex_to_real(c2r3d,a11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a21,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a31,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a32,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,a33,ignore_me)

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
    write(*,*) 'finishing reading bij'

    pdfome = 0._dp
    pdfovec = 0._dp
    pdfalign = 0._dp
    pdfalignomemom = 0._dp
    pointcnt = 0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1

          aa(1,1) = real( a11(ll,jj,kk),sp )
          aa(1,2) = real( a12(ll,jj,kk),sp )
          aa(1,3) = real( a13(ll,jj,kk),sp )
          aa(2,1) = real( a21(ll,jj,kk),sp )
          aa(2,2) = real( a22(ll,jj,kk),sp )
          aa(2,3) = real( a23(ll,jj,kk),sp )
          aa(3,1) = real( a31(ll,jj,kk),sp )
          aa(3,2) = real( a32(ll,jj,kk),sp )
          aa(3,3) = real( a33(ll,jj,kk),sp )

          cc(1,1) = real( b11(ll,jj,kk),sp )
          cc(1,2) = real( b12(ll,jj,kk),sp )
          cc(1,3) = real( b13(ll,jj,kk),sp )
          cc(2,1) = real( b21(ll,jj,kk),sp )
          cc(2,2) = real( b22(ll,jj,kk),sp )
          cc(2,3) = real( b23(ll,jj,kk),sp )
          cc(3,1) = real( b31(ll,jj,kk),sp )
          cc(3,2) = real( b32(ll,jj,kk),sp )
          cc(3,3) = real( b33(ll,jj,kk),sp )
      else
          ll = ii/2

          aa(1,1) = aimag( a11(ll,jj,kk) )
          aa(1,2) = aimag( a12(ll,jj,kk) )
          aa(1,3) = aimag( a13(ll,jj,kk) )
          aa(2,1) = aimag( a21(ll,jj,kk) )
          aa(2,2) = aimag( a22(ll,jj,kk) )
          aa(2,3) = aimag( a23(ll,jj,kk) )
          aa(3,1) = aimag( a31(ll,jj,kk) )
          aa(3,2) = aimag( a32(ll,jj,kk) )
          aa(3,3) = aimag( a33(ll,jj,kk) )

          cc(1,1) = aimag( b11(ll,jj,kk) )
          cc(1,2) = aimag( b12(ll,jj,kk) )
          cc(1,3) = aimag( b13(ll,jj,kk) )
          cc(2,1) = aimag( b21(ll,jj,kk) )
          cc(2,2) = aimag( b22(ll,jj,kk) )
          cc(2,3) = aimag( b23(ll,jj,kk) )
          cc(3,1) = aimag( b31(ll,jj,kk) )
          cc(3,2) = aimag( b32(ll,jj,kk) )
          cc(3,3) = aimag( b33(ll,jj,kk) )
      endif

      detbij =  cc(1,1) * cc(2,2) * cc(3,3) + cc(1,2) * cc(2,3) * cc(3,1) &
              + cc(1,3) * cc(2,1) * cc(3,2) - cc(1,2) * cc(2,1) * cc(3,3) &
              - cc(1,1) * cc(2,3) * cc(3,2) - cc(1,3) * cc(2,2) * cc(3,1)
      if ( abs(detbij - 1._sp) .le. 0.01) then 
          pointcnt = pointcnt + 1
      else
          cycle
      end if

      ! Vorticity
      ovec(1) = aa(3,2) - aa(2,3)
      ovec(2) = aa(1,3) - aa(3,1)
      ovec(3) = aa(2,1) - aa(1,2)

      ! Lagrangian velocity gradient
      bb = matmul(aa, cc)

      ! Angular momentum
      bb = matmul(cc, transpose(bb))
      iangmom(1) = bb(2,3) - bb(3,2)
      iangmom(2) = bb(3,1) - bb(1,3)
      iangmom(3) = bb(1,2) - bb(2,1)

      ! CG tensor
      cc = matmul(cc, transpose(cc))

      ! normalization
      ome2norm = 1._dp / ( cc(1,1) + cc(2,2) + cc(3,3) ) ! temporary use of ome2norm
      cc = cc * ome2norm
      ome = iangmom * ome2norm

      ! ----- Here we need to solve equation for omega_e ---
      ! omega_e is the equivalent rotation angular velocity
      ! But this is not true probably. This code is in dormant
      call dpotrf("U", 3, cc, 3, ierr) ! cc will be overwritten
      if ( ierr .eq. 0 ) then
          call dpotrs("U", 3, 1, cc, 3, ome, 3, ierr)
      else 
          write(*,*) 'Solution for ome fails. ierr = ', ierr
      end if

      ! Now iangmom = omega_e
      ome2norm = dot_product(ome, ome)
      ovec2norm = dot_product(ovec, ovec)

      align = dot_product(ome, ovec) / sqrt(ome2norm * ovec2norm)
      alignomemom = dot_product(ome, iangmom) / sqrt(ome2norm * dot_product(iangmom, iangmom))

      ll = floor( ome2norm / bwome ) + 1
      if (ll .ge. 1 .and. ll .le. npome) pdfome(ll) = pdfome(ll) + 1

      ll = floor( ovec2norm / bwome ) + 1
      if (ll .ge. 1 .and. ll .le. npome) pdfovec(ll) = pdfovec(ll) + 1

      ll = floor( (align + bnd) / binw ) + 1
      if (ll .ge. 1 .and. ll .le. npnt) pdfalign(ll) = pdfalign(ll) + 1

      ll = floor( (alignomemom + bnd) / binw ) + 1
      if (ll .ge. 1 .and. ll .le. npnt) pdfalignomemom(ll) = pdfalignomemom(ll) + 1

    end do
    end do
    end do
    pdfome = pdfome * const
    pdfovec = pdfovec * const
    
    write(*,*) 'normalization check pdfovec: ', sum(pdfovec)
    write(*,*) 'normalization check pdfome: ', sum(pdfome)

    write(*,*) 'fraction of points: ', pointcnt * const

    pdfome = pdfome / bwome
    pdfovec = pdfovec / bwome

    pdfalign = pdfalign * const / binw
    pdfalignomemom = pdfalignomemom * const / binw

    write(27, '( "#index = ", I6)') nfile
    do ll = 1, npnt
      write(27, '(15E15.3)') -bnd + (ll - 0.5_sp) * binw, pdfalign(ll), pdfalignomemom(ll)
    end do
    write(27, *)
    write(27, *)

    write(28, '( "#index = ", I6)') nfile
    do ll = 1, npome
      write(28, '(15E15.3)') (ll - 0.5_sp) * bwome, pdfome(ll), pdfovec(ll)
    end do
    write(28, *)
    write(28, *)
    nfile = nfile + 1
  end do
  close(30)
  close(27)
  close(28)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33)
  deallocate(a11,a12,a13,a21,a22,a23,a31,a32,a33)
  
  call destroyplan3d

  write(*,*) 'defm-pdf-ome.x done.'

end program defmpdfome
