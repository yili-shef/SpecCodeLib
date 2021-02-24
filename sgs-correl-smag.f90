program correlation
  use mconstant
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz, vx, vy, vz
  complex(sp), allocatable, dimension(:,:,:) :: tij
  complex(sp), allocatable, dimension(:,:,:) :: mt11, mt12, mt13, mt22, mt23, mt33
  real(sp), allocatable, dimension(:,:,:) :: k2, g
  real(sp), allocatable, dimension(:) :: kx, ky, kz

  real(sp) :: delta_c, const, ignore_me

  real(dp) :: mdt11, mdt12, mdt13, mdt22, mdt23, mdt33
  real(dp) :: vdt11, vdt12, vdt13, vdt22, vdt23, vdt33
  real(dp) :: mmt11, mmt12, mmt13, mmt22, mmt23, mmt33
  real(dp) :: vmt11, vmt12, vmt13, vmt22, vmt23, vmt33
  real(dp) :: mdtii, mdtij, mmtii, mmtij, vdtii, vdtij, vmtii, vmtij
  real(dp) :: covt11, covt12, covt13, covt22, covt23, covt33, covtii, covtij
  real(dp) :: rhot11, rhot12, rhot13, rhot22, rhot23, rhot33, rhotii, rhotij

  character(80) :: str, fnm, ndellist, str1, fpath

  write(*,*)
  write(*,'('' >>> correlation coefficients of the sgs stress componenets for smag model <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./sgs-correl-smag.x nx filelist ndellist filtertype '
          write(*,*) '        nx: resolution'
          write(*,*) '        filelist: the list of data files, *.list'
          write(*,*) '        ndellist: list of filter scales ndel: Delta = ndel * dx'
          write(*,*) '        filtertype: 0 for Gaussian 1 for cutoff'
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! filelist string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ! filter scale
  call getarg(3,ndellist)
  ndellist = adjustl(ndellist)

  ! filter type
  call getarg(4,str)
  read(str, '(I20)') mm

  str='-smag-'//fnm(1:len_trim(fnm))//'.dat'
  fnm = fnm(1:len_trim(fnm))//'.list'

  ny = nx; nz = nx
  lx = nx / 2; ly = ny; lz = nz; lx1 = lx + 1
  const = 1./ (nx * ny * nz)


  call fftwplan3de(nx,ny,nz) ! fftwplan3de: estimated plan
  write(*,*) 'after fftwplan3d'


  allocate( ux(lx1,ly,lz),   uy(lx1,ly,lz),   uz(lx1,ly,lz)  )
  allocate( vx(lx1,ly,lz),   vy(lx1,ly,lz),   vz(lx1,ly,lz)  )
  allocate( mt11(lx1,ly,lz), mt12(lx1,ly,lz), mt13(lx1,ly,lz))
  allocate( mt22(lx1,ly,lz), mt23(lx1,ly,lz), mt33(lx1,ly,lz))
  allocate( tij(lx1,ly,lz),  k2(lx1,ly,lz),   g(lx1,ly,lz)   )
  allocate( kx(lx1), ky(ly), kz(lz) )

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(28, file = 'mean-modtij'//str( 1:len_trim(str) ) )
  open(29, file = 'var-modtij'//str( 1:len_trim(str) ) )
  open(30, file = 'mean-dnstij'//str( 1:len_trim(str) ) )
  open(31, file = 'var-dnstij'//str( 1:len_trim(str) ) )
  open(32, file = 'correltij'//str( 1:len_trim(str) ) )
  open(33, file = 'rhotij'//str( 1:len_trim(str) ) )

  open(25, file = ndellist(1:len_trim(ndellist))//'.list')

    do while ( .not. eof(25) )

      read(25,*) str1
      read(str1, '(I20)') ndel
      write(*,*) 'ndel = ', str1(1:len_trim(str1))
  
      delta_c=ndel*2*pi/nx
      write(*,*) 'delta_c = ', delta_c
  
      if ( mm .eq. 1 ) then
        where (k2 .ge. (pi/delta_c) * (pi/delta_c))
                g = 0.
        elsewhere
                g = 1.
        endwhere
      else if ( mm .eq. 0 ) then 
        g=exp(-k2*delta_c**2/24.)
      else
        write(*,*) 'Unknown filter type! Stopping'
        stop
      end if
  
      mdt11  = 0.d0; mdt12  = 0.d0; mdt13  = 0.d0; mdt22  = 0.d0; mdt23  = 0.d0; mdt33  = 0.d0
      vdt11  = 0.d0; vdt12  = 0.d0; vdt13  = 0.d0; vdt22  = 0.d0; vdt23  = 0.d0; vdt33  = 0.d0
      mmt11  = 0.d0; mmt12  = 0.d0; mmt13  = 0.d0; mmt22  = 0.d0; mmt23  = 0.d0; mmt33  = 0.d0
      vmt11  = 0.d0; vmt12  = 0.d0; vmt13  = 0.d0; vmt22  = 0.d0; vmt23  = 0.d0; vmt33  = 0.d0
      covt11 = 0.d0; covt12 = 0.d0; covt13 = 0.d0; covt22 = 0.d0; covt23 = 0.d0; covt33 = 0.d0
  
      open(20,file=fnm(1:len_trim(fnm)))
  
        nfile = 0
        do while ( .not. eof(20)) 
          read(20,*) str1
          write(*,*) str1(1:len_trim(str1))
  
          fpath='./out/ux'//str1(1:len_trim(str1))
          open(10,file=fpath,status='unknown',form='unformatted')
            read(10)ux
          close(10)
          fpath='./out/uy'//str1(1:len_trim(str1))
          open(10,file=fpath,status='unknown',form='unformatted')
            read(10)uy
          close(10)
          fpath='./out/uz'//str1(1:len_trim(str1))
          open(10,file=fpath,status='unknown',form='unformatted')
            read(10)uz
          close(10)
          write(*,*) 'after reading data files'
  
          vx = ux * g
          vy = uy * g
          vz = uz * g
  
          call rsmag(vx,vy,vz,mt11,mt12,mt13,mt22,mt23,mt33,kx,ky,kz,lx1,ly,lz,nx,ny,nz,delta_c)
          call rfftwnd_f77_one_real_to_complex(r2c3d,mt11,ignore_me)
          call rfftwnd_f77_one_real_to_complex(r2c3d,mt12,ignore_me)
          call rfftwnd_f77_one_real_to_complex(r2c3d,mt13,ignore_me)
          call rfftwnd_f77_one_real_to_complex(r2c3d,mt22,ignore_me)
          call rfftwnd_f77_one_real_to_complex(r2c3d,mt23,ignore_me)
          call rfftwnd_f77_one_real_to_complex(r2c3d,mt33,ignore_me)
          mt11 = mt11 * g * const
          mt12 = mt12 * g * const
          mt13 = mt13 * g * const
          mt22 = mt22 * g * const
          mt23 = mt23 * g * const
          mt33 = mt33 * g * const
          call rfftwnd_f77_one_complex_to_real(c2r3d,mt11,ignore_me)
          call rfftwnd_f77_one_complex_to_real(c2r3d,mt12,ignore_me)
          call rfftwnd_f77_one_complex_to_real(c2r3d,mt13,ignore_me)
          call rfftwnd_f77_one_complex_to_real(c2r3d,mt22,ignore_me)
          call rfftwnd_f77_one_complex_to_real(c2r3d,mt23,ignore_me)
          call rfftwnd_f77_one_complex_to_real(c2r3d,mt33,ignore_me)

          mmt11 = mmt11 + sum( real( mt11(1:lx,ly,lz) ) ) + sum( aimag( mt11(1:lx,ly,lz) ) )
          mmt12 = mmt12 + sum( real( mt12(1:lx,ly,lz) ) ) + sum( aimag( mt12(1:lx,ly,lz) ) )
          mmt13 = mmt13 + sum( real( mt13(1:lx,ly,lz) ) ) + sum( aimag( mt13(1:lx,ly,lz) ) )
          mmt22 = mmt22 + sum( real( mt22(1:lx,ly,lz) ) ) + sum( aimag( mt22(1:lx,ly,lz) ) )
          mmt23 = mmt23 + sum( real( mt23(1:lx,ly,lz) ) ) + sum( aimag( mt23(1:lx,ly,lz) ) )
          mmt33 = mmt33 + sum( real( mt33(1:lx,ly,lz) ) ) + sum( aimag( mt33(1:lx,ly,lz) ) )
  
          vmt11 = vmt11 + sum( real( mt11(1:lx,ly,lz) * conjg( mt11(1:lx,ly,lz) ) ) )
          vmt12 = vmt12 + sum( real( mt12(1:lx,ly,lz) * conjg( mt12(1:lx,ly,lz) ) ) )
          vmt13 = vmt13 + sum( real( mt13(1:lx,ly,lz) * conjg( mt13(1:lx,ly,lz) ) ) )
          vmt22 = vmt22 + sum( real( mt22(1:lx,ly,lz) * conjg( mt22(1:lx,ly,lz) ) ) )
          vmt23 = vmt23 + sum( real( mt23(1:lx,ly,lz) * conjg( mt23(1:lx,ly,lz) ) ) )
          vmt33 = vmt33 + sum( real( mt33(1:lx,ly,lz) * conjg( mt33(1:lx,ly,lz) ) ) )
  

          call rsgstauij(ux,ux,vx,g,nx,ny,nz) 
          call rsgstauij(uy,uy,vy,g,nx,ny,nz) 
          call rsgstauij(uz,uz,vz,g,nx,ny,nz) 
          tij = - (vx + vy + vz) /3.
          vx = vx + tij; vy = vy + tij; vz = vz + tij

          mdt11 = mdt11 + sum( real( vx(1:lx,ly,lz) ) ) + sum( aimag( vx(1:lx,ly,lz) ) )
          vdt11 = vdt11 + real( sum( vx(1:lx,ly,lz) * conjg( vx(1:lx,ly,lz) ) ) )
  
          covt11 = covt11 + sum( real( mt11(1:lx,ly,lz) ) * real( vx(1:lx,ly,lz) ) ) + &
                         sum( aimag( mt11(1:lx,ly,lz) ) * aimag( vx(1:lx,ly,lz) ) )
  
          mdt22 = mdt22 + sum( real( vy(1:lx,ly,lz) ) ) + sum( aimag( vy(1:lx,ly,lz) ) )
          vdt22 = vdt22 + real( sum( vy(1:lx,ly,lz) * conjg( vy(1:lx,ly,lz) ) ) )
  
          covt22 = covt22 + sum( real( mt22(1:lx,ly,lz) ) * real( vy(1:lx,ly,lz) ) ) + &
                         sum( aimag( mt22(1:lx,ly,lz) ) * aimag( vy(1:lx,ly,lz) ) )
  
          mdt33 = mdt33 + sum( real( vz(1:lx,ly,lz) ) ) + sum( aimag( vz(1:lx,ly,lz) ) )
          vdt33 = vdt33 + real( sum( vz(1:lx,ly,lz) * conjg( vz(1:lx,ly,lz) ) ) )
  
          covt33 = covt33 + sum( real( mt33(1:lx,ly,lz) ) * real( vz(1:lx,ly,lz) ) ) + &
                         sum( aimag( mt33(1:lx,ly,lz) ) * aimag( vz(1:lx,ly,lz) ) )
  
          call rsgstauij(ux,uy,tij,g,nx,ny,nz) 
          mdt12 = mdt12 + sum( real( tij(1:lx,ly,lz) ) ) + sum( aimag( tij(1:lx,ly,lz) ) )
          vdt12 = vdt12 + real( sum( tij(1:lx,ly,lz) * conjg( tij(1:lx,ly,lz) ) ) )
  
          covt12 = covt12 + sum( real( mt12(1:lx,ly,lz) ) * real( tij(1:lx,ly,lz) ) ) + &
                         sum( aimag( mt12(1:lx,ly,lz) ) * aimag( tij(1:lx,ly,lz) ) )
  
          call rsgstauij(ux,uz,tij,g,nx,ny,nz) 
          mdt13 = mdt13 + sum( real( tij(1:lx,ly,lz) ) ) + sum( aimag( tij(1:lx,ly,lz) ) )
          vdt13 = vdt13 + real( sum( tij(1:lx,ly,lz) * conjg( tij(1:lx,ly,lz) ) ) )
  
          covt13 = covt13 + sum( real( mt13(1:lx,ly,lz) ) * real( tij(1:lx,ly,lz) ) ) + &
                         sum( aimag( mt13(1:lx,ly,lz) ) * aimag( tij(1:lx,ly,lz) ) )
  
          call rsgstauij(uy,uz,tij,g,nx,ny,nz) 
          mdt23 = mdt23 + sum( real( tij(1:lx,ly,lz) ) ) + sum( aimag( tij(1:lx,ly,lz) ) )
          vdt23 = vdt23 + real( sum( tij(1:lx,ly,lz) * conjg( tij(1:lx,ly,lz) ) ) )
  
          covt23 = covt23 + sum( real( mt23(1:lx,ly,lz) ) * real( tij(1:lx,ly,lz) ) ) + &
                         sum( aimag( mt23(1:lx,ly,lz) ) * aimag( tij(1:lx,ly,lz) ) )
  
          nfile = nfile + 1
        end do
  
      close(20)

      ignore_me = const / nfile
  
      mdt11 = mdt11 * ignore_me
      mdt12 = mdt12 * ignore_me
      mdt13 = mdt13 * ignore_me
      mdt22 = mdt22 * ignore_me
      mdt23 = mdt23 * ignore_me
      mdt33 = mdt33 * ignore_me
  
      vdt11 = vdt11 * ignore_me
      vdt12 = vdt12 * ignore_me
      vdt13 = vdt13 * ignore_me
      vdt22 = vdt22 * ignore_me
      vdt23 = vdt23 * ignore_me
      vdt33 = vdt33 * ignore_me
  
      vdt11 = vdt11 - mdt11 * mdt11
      vdt12 = vdt12 - mdt12 * mdt12
      vdt13 = vdt13 - mdt13 * mdt13
      vdt22 = vdt22 - mdt22 * mdt22
      vdt23 = vdt23 - mdt23 * mdt23
      vdt33 = vdt33 - mdt33 * mdt33
  
      mmt11 = mmt11 * ignore_me
      mmt12 = mmt12 * ignore_me
      mmt13 = mmt13 * ignore_me
      mmt22 = mmt22 * ignore_me
      mmt23 = mmt23 * ignore_me
      mmt33 = mmt33 * ignore_me
  
      vmt11 = vmt11 * ignore_me
      vmt12 = vmt12 * ignore_me
      vmt13 = vmt13 * ignore_me
      vmt22 = vmt22 * ignore_me
      vmt23 = vmt23 * ignore_me
      vmt33 = vmt33 * ignore_me
  
      vmt11 = vmt11 - mmt11 * mmt11
      vmt12 = vmt12 - mmt12 * mmt12
      vmt13 = vmt13 - mmt13 * mmt13
      vmt22 = vmt22 - mmt22 * mmt22
      vmt23 = vmt23 - mmt23 * mmt23
      vmt33 = vmt33 - mmt33 * mmt33
  
      covt11 = covt11 * ignore_me
      covt12 = covt12 * ignore_me
      covt13 = covt13 * ignore_me
      covt22 = covt22 * ignore_me
      covt23 = covt23 * ignore_me
      covt33 = covt33 * ignore_me
  
      covt11 = covt11 - mmt11 * mdt11
      covt12 = covt12 - mmt12 * mdt12
      covt13 = covt13 - mmt13 * mdt13
      covt22 = covt22 - mmt22 * mdt22
      covt23 = covt23 - mmt23 * mdt23
      covt33 = covt33 - mmt33 * mdt33
  

      rhot11 = covt11 / sqrt( vdt11 ) / sqrt( vmt11 )
      rhot22 = covt22 / sqrt( vdt22 ) / sqrt( vmt22 )
      rhot33 = covt33 / sqrt( vdt33 ) / sqrt( vmt33 )
      rhot12 = covt12 / sqrt( vdt12 ) / sqrt( vmt12 )
      rhot13 = covt13 / sqrt( vdt13 ) / sqrt( vmt13 )
      rhot23 = covt23 / sqrt( vdt23 ) / sqrt( vmt23 )

      covtij = covt11 + covt22 + covt22 + 2. * (covt12 + covt13 + covt23)
      vmtij = vmt11 + vmt22 + vmt33 + 2. * (vmt12 + vmt13 + vmt23)
      vdtij = vdt11 + vdt22 + vdt33 + 2. * (vdt12 + vdt13 + vdt23)
  
      rhotij = covtij / sqrt( vdtij ) / sqrt( vmtij )
!      covtij = covtij - mdt11 * mmt11 - mdt22 * mmt22 - mdt33 * mmt33 - 2. * &
!               (mdt12 * mmt12 + mdt13 * mmt13 + mdt23 * mmt23)
!      vmtij = vmtij - mmt11 * mmt11 - mmt22 * mmt22 - mmt33 * mmt33 - 2. *  &
!              (mmt12 * mmt12 + mmt13 * mmt13 + mmt23 * mmt23)
!      vdtij = vdtij - mdt11 * mdt11 - mdt22 * mdt22 - mdt33 * mdt33 - 2. * &
!              (mdt12 * mdt12 + mdt13 * mdt13 + mdt23 * mdt23)

!      rhotii = (covtii - mdtii * mmtii) / sqrt( vdtii - mdtii*mdtii ) / sqrt( vmtii - mmtii*mmtii )
!      rhotij = (covtij - mdtij * mmtij) / sqrt( vdtij - mdtij*mdtij ) / sqrt( vmtij - mmtij*mmtij )


      write(28, '(I4, 30E15.4)') ndel, mmt11, mmt12, mmt13, mmt22, mmt23, mmt33, mmtij
      write(29, '(I4, 30E15.4)') ndel, vmt11, vmt12, vmt13, vmt22, vmt23, vmt33, vmtij
      write(30, '(I4, 30E15.4)') ndel, mdt11, mdt12, mdt13, mdt22, mdt23, mdt33, mdtij
      write(31, '(I4, 30E15.4)') ndel, vdt11, vdt12, vdt13, vdt22, vdt23, vdt33, vdtij
      write(32, '(I4, 30E15.4)') ndel, covt11, covt12, covt13, covt22, covt23, &
        covt33, covtij
      write(33, '(I4, 30E15.4)') ndel, rhot11, rhot12, rhot13, rhot22, rhot23, &
        rhot33, rhotij

  
    end do

  close(25)

  close(28)
  close(29)
  close(30)
  close(31)
  close(32)
  close(33)

  call destroyplan3d

  deallocate( ux, uy, uz, vx, vy, vz )
  deallocate( mt11, mt12, mt13, mt22, mt23, mt33 )
  deallocate( tij, k2, g, kx, ky, kz )

  write(*,*) 'finished'

end program correlation 
