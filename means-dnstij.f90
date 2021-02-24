program meansdnstij
  use mconstant
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz, t11, t22, t33
  complex(sp), allocatable, dimension(:,:,:) :: tij
  real(sp), allocatable, dimension(:,:,:) :: k2, g
  real(sp), allocatable, dimension(:) :: kx, ky, kz

  real(sp) :: delta_c, const, ignore_me

  real(dp) :: mdt11, mdt12, mdt13, mdt22, mdt23, mdt33
  real(dp) :: vdt11, vdt12, vdt13, vdt22, vdt23, vdt33
  real(dp) :: mdtii, mdtij, vdtii, vdtij

  character(80) :: str, fnm, ndellist, str1, fpath

  write(*,*)
  write(*,'('' >>> mean statistics of the sgs stress componenets from dns <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./means-dnstij.x nx filelist ndellist filtertype '
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

  str='-'//fnm(1:len_trim(fnm))//'.dat'
  fnm = fnm(1:len_trim(fnm))//'.list'

  ny = nx; nz = nx
  lx = nx / 2; ly = ny; lz = nz; lx1 = lx + 1
  const = 1./ (nx * ny * nz)


  call fftwplan3de(nx,ny,nz) ! fftwplan3de: estimated plan
  write(*,*) 'after fftwplan3d'


  allocate( ux(lx1,ly,lz),   uy(lx1,ly,lz),   uz(lx1,ly,lz)  )
  allocate( t11(lx1,ly,lz),  t22(lx1,ly,lz),  t33(lx1,ly,lz)  )
  allocate( tij(lx1,ly,lz),  k2(lx1,ly,lz),   g(lx1,ly,lz)   )
  allocate( kx(lx1), ky(ly), kz(lz) )

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(30, file = 'mean-dnstij'//str( 1:len_trim(str) ) )
  open(31, file = 'var-dnstij'//str( 1:len_trim(str) ) )

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
  

          call rsgstauij(ux,ux,t11,g,nx,ny,nz) 
          call rsgstauij(uy,uy,t22,g,nx,ny,nz) 
          call rsgstauij(uz,uz,t33,g,nx,ny,nz) 
          tij = - (t11 + t22 + t33) /3.
          t11 = t11 + tij; t22 = t22 + tij; t33 = t33 + tij

          mdt11 = mdt11 + sum( real( t11(1:lx,ly,lz) ) ) + sum( aimag( t11(1:lx,ly,lz) ) )
          vdt11 = vdt11 + real( sum( t11(1:lx,ly,lz) * conjg( t11(1:lx,ly,lz) ) ) )
  
          mdt22 = mdt22 + sum( real( t22(1:lx,ly,lz) ) ) + sum( aimag( t22(1:lx,ly,lz) ) )
          vdt22 = vdt22 + real( sum( t22(1:lx,ly,lz) * conjg( t22(1:lx,ly,lz) ) ) )
  
          mdt33 = mdt33 + sum( real( t33(1:lx,ly,lz) ) ) + sum( aimag( t33(1:lx,ly,lz) ) )
          vdt33 = vdt33 + real( sum( t33(1:lx,ly,lz) * conjg( t33(1:lx,ly,lz) ) ) )
  
          call rsgstauij(ux,uy,tij,g,nx,ny,nz) 
          mdt12 = mdt12 + sum( real( tij(1:lx,ly,lz) ) ) + sum( aimag( tij(1:lx,ly,lz) ) )
          vdt12 = vdt12 + real( sum( tij(1:lx,ly,lz) * conjg( tij(1:lx,ly,lz) ) ) )
  
          call rsgstauij(ux,uz,tij,g,nx,ny,nz) 
          mdt13 = mdt13 + sum( real( tij(1:lx,ly,lz) ) ) + sum( aimag( tij(1:lx,ly,lz) ) )
          vdt13 = vdt13 + real( sum( tij(1:lx,ly,lz) * conjg( tij(1:lx,ly,lz) ) ) )
  
          call rsgstauij(uy,uz,tij,g,nx,ny,nz) 
          mdt23 = mdt23 + sum( real( tij(1:lx,ly,lz) ) ) + sum( aimag( tij(1:lx,ly,lz) ) )
          vdt23 = vdt23 + real( sum( tij(1:lx,ly,lz) * conjg( tij(1:lx,ly,lz) ) ) )
  
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

      write(30, '(I4, 30E15.4)') ndel, mdt11, mdt12, mdt13, mdt22, mdt23, mdt33, mdtij
      write(31, '(I4, 30E15.4)') ndel, vdt11, vdt12, vdt13, vdt22, vdt23, vdt33, vdtij

  
    end do

  close(25)

  close(30)
  close(31)

  call destroyplan3d

  deallocate( ux, uy, uz, t11, t22, t33 )
  deallocate( tij, k2, g, kx, ky, kz )

  write(*,*) 'finished'

end program meansdnstij 
