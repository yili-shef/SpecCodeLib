program tijcospec
  use mconstant
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz, vx, vy, vz
  complex(sp), allocatable, dimension(:,:,:) :: tij
  complex(sp), allocatable, dimension(:,:,:) :: mt11, mt12, mt13, mt22, mt23, mt33
  real(sp), allocatable, dimension(:,:,:) :: k2, g
  real(sp), allocatable, dimension(:) :: kx, ky, kz 
  real(dp), allocatable, dimension(:) :: cospec, dnstijspec, modtijspec

  real(sp) :: delta_c, const, ignore_me, kxii, kyjj, kzkk
  real(dp) :: rmsdnstij, rmsmodtij

  character(80) :: str, fnm, ndelstr, fpath, filtstr

  write(*,*)
  write(*,'('' >>> co-spectrum of the sgs stress componenets between dns and several models <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./tij-cospec.x nx filelist ndel filtertype modeltype'
          write(*,*) '        nx: resolution'
          write(*,*) '        filelist: the list of data files, *.list'
          write(*,*) '        ndel: filter scales ndel: Delta = ndel * dx'
          write(*,*) '        filtertype: 0 for Gaussian 1 for cutoff'
          write(*,*) '        modeltype: 1 Smag; 2 similarity; 3 nonlinear; 4 mat exp'
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  ! resolution
  call getarg(1,str)
  read(str, '(I10)') nx

  ! filelist string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ! filter scale
  call getarg(3,ndelstr)
  read(ndelstr, '(I10)') ndel
  ndelstr=adjustl(ndelstr)

  ! filter type
  call getarg(4,str)
  read(str, '(I10)') mm

  call getarg(5,str)
  read(str, '(I10)') nn

  ny = nx; nz = nx
  lx = nx / 2; ly = ny; lz = nz; lx1 = lx + 1
  const = 1./ (nx * ny * nz)

  call fftwplan3d(nx,ny,nz) ! fftwplan3de: estimated plan
  write(*,*) 'after fftwplan3d'

  allocate( ux(lx1,ly,lz),   uy(lx1,ly,lz),   uz(lx1,ly,lz)  )
  allocate( vx(lx1,ly,lz),   vy(lx1,ly,lz),   vz(lx1,ly,lz)  )
  allocate( mt11(lx1,ly,lz), mt12(lx1,ly,lz), mt13(lx1,ly,lz))
  allocate( mt22(lx1,ly,lz), mt23(lx1,ly,lz), mt33(lx1,ly,lz))
  allocate( tij(lx1,ly,lz),  k2(lx1,ly,lz),   g(lx1,ly,lz)   )
  allocate( kx(lx1), ky(ly), kz(lz) )
  allocate( cospec(lx), dnstijspec(lx), modtijspec(lx) )

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  delta_c=ndel*2*pi/nx
  write(*,*) 'delta_c = ', delta_c

  if ( mm .eq. 1 ) then
    where (k2 .ge. (pi/delta_c) * (pi/delta_c))
      g = 0.
    elsewhere
      g = 1.
    endwhere
    filtstr = 'cutoff'
  else if ( mm .eq. 0 ) then 
    g=exp(-k2*delta_c**2/24.)
    filtstr = 'gaussian'
  else
    write(*,*) 'Unknown filter type! Stopping'
    stop
  end if
  
  open(25, file = fnm(1:len_trim(fnm))//'.list')

    cospec = 0.d0
    dnstijspec = 0.d0
    modtijspec = 0.d0
    nfile = 0
    do while ( .not. eof(25) )

      read(25,*) str
      write(*,*) str(1:len_trim(str))
  
      fpath='./out/ux'//str(1:len_trim(str))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)ux
      close(10)
      fpath='./out/uy'//str(1:len_trim(str))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)uy
      close(10)
      fpath='./out/uz'//str(1:len_trim(str))
      open(10,file=fpath,status='unknown',form='unformatted')
        read(10)uz
      close(10)
      write(*,*) 'after reading data files'
  
      vx = ux * g
      vy = uy * g
      vz = uz * g
      select case (nn)
      case (1)
        call rsmag(vx,vy,vz,mt11,mt12,mt13,mt22,mt23,mt33,kx,ky,kz,lx1,ly,lz,nx,ny,nz,delta_c)
      case (2)
        call rsgstauij(vx,vx,mt11,g,nx,ny,nz)
        call rsgstauij(vx,vy,mt12,g,nx,ny,nz)
        call rsgstauij(vx,vz,mt13,g,nx,ny,nz)
        call rsgstauij(vy,vy,mt22,g,nx,ny,nz)
        call rsgstauij(vy,vz,mt23,g,nx,ny,nz)
        call rsgstauij(vz,vz,mt33,g,nx,ny,nz)
        vx = - (mt11 + mt22 + mt33)/3.
        mt11 = mt11 + vx; mt22 = mt22 + vx; mt33 = mt33 + vx
      case (3)
        write(*,*) 'not implemented yet'
        stop
      case (4)
        call rmatexp(vx,vy,vz,mt11,mt12,mt13,mt22,mt23,mt33,kx,ky,kz,lx1,ly,lz,nx,ny,nz,delta_c)
      end select
      call rfftwnd_f77_one_real_to_complex(r2c3d,mt11,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,mt12,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,mt13,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,mt22,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,mt23,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,mt33,ignore_me)
      mt11 = mt11 * const
      mt12 = mt12 * const
      mt13 = mt13 * const
      mt22 = mt22 * const
      mt23 = mt23 * const
      mt33 = mt33 * const

      call rsgstauij(ux,ux,vx,g,nx,ny,nz) 
      call rsgstauij(uy,uy,vy,g,nx,ny,nz) 
      call rsgstauij(uz,uz,vz,g,nx,ny,nz) 
      tij = - (vx + vy + vz) /3.
      vx = vx + tij; vy = vy + tij; vz = vz + tij
      call rfftwnd_f77_one_real_to_complex(r2c3d,vx,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,vy,ignore_me)
      call rfftwnd_f77_one_real_to_complex(r2c3d,vz,ignore_me)
      vx = vx * const
      vy = vy * const
      vz = vz * const

      mt11 = mt11 * conjg(mt11) + mt22 * conjg(mt22) + mt33 * conjg(mt33)
      mt11 = mt11 + 2. * ( mt12 * conjg(mt12) + mt13 * conjg(mt13) + mt23 * conjg(mt23) )      
      mt11 = 2. * mt11
      mt11(1,:,:) = .5 * mt11(1,:,:)

      mt22 =  vx * conjg(vx) + vy * conjg(vy) + vz * conjg(vz)
      tij = mt11 * conjg(vx) + mt22 * conjg(vy) + mt33 * conjg(vz)

      call rsgstauij(ux,uy,vx,g,nx,ny,nz) 
      call rsgstauij(ux,uz,vy,g,nx,ny,nz) 
      call rsgstauij(uy,uz,vz,g,nx,ny,nz) 
      vx = vx * const
      vy = vy * const
      vz = vz * const

      tij = tij + 2. * ( mt12 * conjg(vx) + mt13 * conjg(vy) + mt23 * conjg(vz) ) 
      tij = 2. * tij
      tij(1,:,:) = .5 * tij(1,:,:)
  
      mt22 = mt22 + 2. * ( vx * conjg(vx) + vy * conjg(vy) + vz * conjg(vz) )
      mt22 = 2. * mt22
      mt22(1,:,:) = .5 * mt22(1,:,:)

      do kk=1,lz
        kzkk=kz(kk)
        do jj=1,ly
          kyjj=ky(jj)
          do ii=1,lx1
            kxii=kx(ii)
      
            ll=floor(sqrt(kxii*kxii+kyjj*kyjj+kzkk*kzkk)+.5)
            if (ll .ge. 1 .and. ll .le. lx) then
              cospec(ll) = cospec(ll) + real( tij(ii,jj,kk) )
              modtijspec(ll) = modtijspec(ll) + real( mt11(ii,jj,kk) )
              dnstijspec(ll) = dnstijspec(ll) + real( mt22(ii,jj,kk) )
            end if
      
          end do
        end do
      end do
      nfile = nfile + 1

    end do

    cospec = cospec / nfile
    modtijspec = modtijspec / nfile
    dnstijspec = dnstijspec / nfile
    rmsmodtij = sqrt( sum( modtijspec ) )
    rmsdnstij = sqrt( sum( dnstijspec ) )

  close(25)

  select case(nn)
  case (1)
    str = 'smag-ndel'//ndelstr(1:len_trim(ndelstr))//'-'//fnm(1:len_trim(fnm))//'-'//filtstr(1:len_trim(filtstr))//'.dat'
  case (2)
    str = 'sim-ndel'//ndelstr(1:len_trim(ndelstr))//'-'//fnm(1:len_trim(fnm))//'-'//filtstr(1:len_trim(filtstr))//'.dat'
  case (3)
    str = 'nl-ndel'//ndelstr(1:len_trim(ndelstr))//'-'//fnm(1:len_trim(fnm))//'-'//filtstr(1:len_trim(filtstr))//'.dat'
  case (4)
    str = 'matexp-ndel'//ndelstr(1:len_trim(ndelstr))//'-'//fnm(1:len_trim(fnm))//'-'//filtstr(1:len_trim(filtstr))//'.dat'
  case default
    write(*,*) 'Wrong model type. Stopping'
    stop
  end select

  open(33, file = 'tij-cospec-'//str( 1:len_trim(str) ) )
    do ii = 1, lx
      write(33, '(I5, 30E13.4)') ii, cospec(ii) / rmsmodtij / rmsdnstij, &
           dnstijspec(ii) / rmsdnstij**2,  modtijspec(ii) / rmsmodtij**2
    end do
  close(33)

  deallocate(ux, uy, uz, vx, vy, vz)
  deallocate(mt11, mt12, mt13, mt22, mt23, mt33)
  deallocate(tij, kx, ky, kz, k2, cospec, dnstijspec, modtijspec, g)

  call destroyplan3d

  write(*,*) 'done'

end program tijcospec
