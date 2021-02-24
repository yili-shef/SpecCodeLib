program defmrotincij
  use mconstant
  use mwavenumber
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,ii,jj,kk,lx1,lx,ly,lz,ll,mm,nfile
  real(sp) :: ignore_me, const

  integer,  parameter :: npnt = 80
  real(sp), parameter :: bnd = 8._sp, binw = bnd/npnt
  real(dp), dimension(npnt) :: pdf_omalfa, pdf_ombeta, pdf_omgmma
  real(dp), dimension(npnt) :: pdf_oma1, pdf_omb2, pdf_omg3
  real(dp), dimension(npnt) :: pdf_oms1, pdf_oms2, pdf_oms3


  ! ----------------- For eigenvalue subroutines from MKL ----------------------------
  external :: dsyevr, dlamch
  integer, parameter :: evsize = 3, lwork = 26 * evsize, liwork = 10 * evsize
  real(dp) :: work(lwork), dlamch
  integer  :: iwork(liwork), isuppz(2*evsize), iignore, nfound, info
  real(dp), dimension(evsize,evsize)  :: cc, sij, evtrcc
  real(dp), dimension(evsize) :: evcc, omega
  real(dp) :: meanalfa, meanbeta, meangmma, meanom1, meanom2, meanom3
  real(dp) :: means1, means2, means3, omalfa, ombeta, omgmma, oma1, omb2, omg3
  ! ----------------------------------------------------------------------------------
  real(dp), parameter :: ugly = 16._dp

  integer :: pointcnt

  complex(sp), allocatable, dimension(:,:,:) :: b11,b12,b13
  complex(sp), allocatable, dimension(:,:,:) :: b21,b22,b23
  complex(sp), allocatable, dimension(:,:,:) :: b31,b32,b33
  complex(sp), allocatable, dimension(:,:,:) :: s11,s22,s33,s12,s13,s23,ox,oy,oz
  real(sp),    allocatable, dimension(:,:,:) :: kx, ky, kz

  character(80) :: str, flnm, prf, str1

  
  if (iargc() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./defm-rot-in-cij.x nx filelist prefix'
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
  allocate( s11(lx1,ly,lz), s22(lx1,ly,lz), s33(lx1,ly,lz) )
  allocate( s12(lx1,ly,lz), s13(lx1,ly,lz), s23(lx1,ly,lz) )
  allocate(  ox(lx1,ly,lz),  oy(lx1,ly,lz),  oz(lx1,ly,lz) )

  call fftwplan3de(nx,ny,nz)
  write(*,*) 'after fftwplan'

  call wavenumber(kx,ky,kz,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(27, file = 'pdf-rot-in-cij-'//trim(flnm)//'.dat')
  open(28, file = 'mean-rot-in-cij-'//trim(flnm)//'.dat')

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

    s11 = eye * kx * b11
    s22 = eye * ky * b22
    s33 = eye * kz * b33
    s12 = eye * (kx * b22 + ky * b11)*0.5_sp
    s13 = eye * (kx * b33 + kz * b11)*0.5_sp
    s23 = eye * (ky * b33 + kz * b22)*0.5_sp
    ox  = eye * (ky * b33 - kz * b22)
    oy  = eye * (kz * b11 - kx * b33)
    oz  = eye * (kx * b22 - ky * b11)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,ox,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,oy,ignore_me)
    call rfftwnd_f77_one_complex_to_real(c2r3d,oz,ignore_me)


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


    pdf_omalfa = 0._dp; pdf_ombeta = 0._dp; pdf_omgmma = 0._dp
    pdf_oma1 = 0._dp; pdf_omb2 = 0._dp; pdf_omg3 = 0._dp
    pdf_oms1 = 0._dp; pdf_oms2 = 0._dp; pdf_oms3 = 0._dp
    meanalfa = 0._dp; meanbeta = 0._dp; meangmma = 0._dp
    meanom1 = 0._dp; meanom2 = 0._dp; meanom3 = 0._dp
    means1 = 0._dp; means2 = 0._dp; means3 = 0._dp
    pointcnt = 0
    do kk = 1, nz
    do jj = 1, ny
    do ii = 1, nx

      if ( mod(ii, 2) .eq. 1 ) then
          ll = ii/2 + 1
          sij(1,1) = real( s11(ll,jj,kk),sp )
          sij(1,2) = real( s12(ll,jj,kk),sp )
          sij(1,3) = real( s13(ll,jj,kk),sp )
          sij(2,2) = real( s22(ll,jj,kk),sp )
          sij(2,3) = real( s23(ll,jj,kk),sp )
          sij(3,3) = real( s33(ll,jj,kk),sp )

          cc(1,1) = real( b11(ll,jj,kk),sp )
          cc(1,2) = real( b12(ll,jj,kk),sp )
          cc(1,3) = real( b13(ll,jj,kk),sp )
          cc(2,1) = real( b21(ll,jj,kk),sp )
          cc(2,2) = real( b22(ll,jj,kk),sp )
          cc(2,3) = real( b23(ll,jj,kk),sp )
          cc(3,1) = real( b31(ll,jj,kk),sp )
          cc(3,2) = real( b32(ll,jj,kk),sp )
          cc(3,3) = real( b33(ll,jj,kk),sp )
          omega(1) = real( ox(ll,jj,kk),sp )
          omega(2) = real( oy(ll,jj,kk),sp )
          omega(3) = real( oz(ll,jj,kk),sp )
      else
          ll = ii/2
          sij(1,1) = aimag( s11(ll,jj,kk) )
          sij(1,2) = aimag( s12(ll,jj,kk) )
          sij(1,3) = aimag( s13(ll,jj,kk) )
          sij(2,2) = aimag( s22(ll,jj,kk) )
          sij(2,3) = aimag( s23(ll,jj,kk) )
          sij(3,3) = aimag( s33(ll,jj,kk) )

          cc(1,1) = aimag( b11(ll,jj,kk) )
          cc(1,2) = aimag( b12(ll,jj,kk) )
          cc(1,3) = aimag( b13(ll,jj,kk) )
          cc(2,1) = aimag( b21(ll,jj,kk) )
          cc(2,2) = aimag( b22(ll,jj,kk) )
          cc(2,3) = aimag( b23(ll,jj,kk) )
          cc(3,1) = aimag( b31(ll,jj,kk) )
          cc(3,2) = aimag( b32(ll,jj,kk) )
          cc(3,3) = aimag( b33(ll,jj,kk) )
          omega(1) = aimag( ox(ll,jj,kk) )
          omega(2) = aimag( oy(ll,jj,kk) )
          omega(3) = aimag( oz(ll,jj,kk) )
      endif
      cc = matmul( cc, transpose(cc) )
      sij(2,1) = sij(1,2); sij(3,1) = sij(1,3); sij(3,2) = sij(2,3)

      ! The interface of subroutine dsyevr
      ! call dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work,
      ! lwork, iwork, liwork, info)

      call dsyevr("V", "A", "U", evsize, cc, evsize, iignore, iignore, iignore, iignore, dlamch('S'), &
                  nfound, evcc, evtrcc, evsize, isuppz, work, lwork, iwork, liwork, info)
      if ( .not. (info .eq. 0) ) write(*,*) 'Something wrong with ev cc. info = ', info

      if (      abs(evcc(1) - evcc(2)) .le. ugly * dlamch('e') * abs(evcc(2)) &
           .or. abs(evcc(2) - evcc(3)) .le. ugly * dlamch('e') * abs(evcc(3)) &
           .or. abs(evcc(3) - evcc(1)) .le. ugly * dlamch('e') * abs(evcc(3)) ) then 
          cycle
      else
          pointcnt = pointcnt + 1
      end if


      ! define a right-handed frame: gamma = alfa X beta
      evtrcc(1,1) = evtrcc(2,3) - evtrcc(3,2)
      evtrcc(2,1) = evtrcc(3,3) - evtrcc(1,2)
      evtrcc(3,1) = evtrcc(1,3) - evtrcc(2,2)
     
      ! The eigenvalues are in ascending order: smallest first.

      oma1 = dot_product( omega, evtrcc(:,3) ) * .5_sp
      omb2 = dot_product( omega, evtrcc(:,2) ) * .5_sp
      omg3 = dot_product( omega, evtrcc(:,1) ) * .5_sp

      omalfa = oma1 + dot_product(evtrcc(:,2), matmul(sij, evtrcc(:,1))) &
                    * ( evcc(2) + evcc(1) ) / ( evcc(2) - evcc(1) )
      ombeta = omb2 + dot_product(evtrcc(:,1), matmul(sij, evtrcc(:,3))) &
                    * ( evcc(3) + evcc(1) ) / ( evcc(1) - evcc(3) )
      omgmma = omg3 + dot_product(evtrcc(:,3), matmul(sij, evtrcc(:,2))) &
                    * ( evcc(3) + evcc(2) ) / ( evcc(3) - evcc(2) ) 

      meanalfa = meanalfa + abs(omalfa)**2
      meanbeta = meanbeta + abs(ombeta)**2
      meangmma = meangmma + abs(omgmma)**2
      meanom1 = meanom1 + abs(oma1)**2
      meanom2 = meanom2 + abs(omb2)**2
      meanom3 = meanom3 + abs(omg3)**2
      means1 = means1 + abs(omalfa - oma1)**2
      means2 = means2 + abs(ombeta - omb2)**2
      means3 = means3 + abs(omgmma - omg3)**2

      ll = floor(abs(omalfa)/binw) + 1
      if (ll .ge. 1 .and. ll .le. npnt) pdf_omalfa(ll)=pdf_omalfa(ll)+1
      ll = floor(abs(ombeta)/binw) + 1
      if (ll .ge. 1 .and. ll .le. npnt) pdf_ombeta(ll)=pdf_ombeta(ll)+1
      ll = floor(abs(omgmma)/binw) + 1
      if (ll .ge. 1 .and. ll .le. npnt) pdf_omgmma(ll)=pdf_omgmma(ll)+1

      ll = floor(abs(oma1)/binw) + 1
      if (ll .ge. 1 .and. ll .le. npnt) pdf_oma1(ll)=pdf_oma1(ll)+1
      ll = floor(abs(omb2)/binw) + 1
      if (ll .ge. 1 .and. ll .le. npnt) pdf_omb2(ll)=pdf_omb2(ll)+1
      ll = floor(abs(omg3)/binw) + 1
      if (ll .ge. 1 .and. ll .le. npnt) pdf_omg3(ll)=pdf_omg3(ll)+1

      ll = floor(abs(omalfa-oma1)/binw) + 1
      if (ll .ge. 1 .and. ll .le. npnt) pdf_oms1(ll)=pdf_oms1(ll)+1
      ll = floor(abs(ombeta-omb2)/binw) + 1
      if (ll .ge. 1 .and. ll .le. npnt) pdf_oms2(ll)=pdf_oms2(ll)+1
      ll = floor(abs(omgmma-omg3)/binw) + 1
      if (ll .ge. 1 .and. ll .le. npnt) pdf_oms3(ll)=pdf_oms3(ll)+1

    end do
    end do
    end do
    pdf_omalfa = pdf_omalfa * const / binw
    pdf_ombeta = pdf_ombeta * const / binw
    pdf_omgmma = pdf_omgmma * const / binw

    pdf_oma1 = pdf_oma1 * const / binw
    pdf_omb2 = pdf_omb2 * const / binw
    pdf_omg3 = pdf_omg3 * const / binw

    pdf_oms1 = pdf_oms1 * const / binw
    pdf_oms2 = pdf_oms2 * const / binw
    pdf_oms3 = pdf_oms3 * const / binw

    write(*,*) 'fraction of points: ', pointcnt * const

    write(*,*) 'normalization pdf_omalfa: ', sum(pdf_omalfa)*binw
    write(*,*) 'normalization pdf_ombeta: ', sum(pdf_ombeta)*binw
    write(*,*) 'normalization pdf_omgmma: ', sum(pdf_omgmma)*binw

    meanalfa = meanalfa * const
    meanbeta = meanbeta * const
    meangmma = meangmma * const
    meanom1  = meanom1  * const
    meanom2  = meanom2  * const
    meanom3  = meanom3  * const
    means1   = means1   * const
    means2   = means2   * const
    means3   = means3   * const

    write(27, '( "#index = ", I6)') nfile
    do ll = 1, npnt
      write(27, '(15E15.3)') (ll - 0.5_sp) * binw, pdf_omalfa(ll), pdf_ombeta(ll), pdf_omgmma(ll), &
      pdf_oma1(ll), pdf_omb2(ll), pdf_omg3(ll), pdf_oms1(ll), pdf_oms2(ll), pdf_oms3(ll)
    end do
    write(27, *)
    write(27, *)

    write(28, '(I4, 15E15.3)') nfile, meanalfa, meanbeta, meangmma, meanom1, meanom2, meanom3, &
                                      means1, means2, means3

    nfile = nfile + 1
  end do
  close(30)
  close(28)
  close(27)

  deallocate(kx,ky,kz,b11,b12,b13,b21,b22,b23,b31,b32,b33,s11,s22,s33,s12,s13,s23,ox,oy,oz)
  
  call destroyplan3d

  write(*,*) 'defm-rot-in-cij.x done.'

end program defmrotincij
