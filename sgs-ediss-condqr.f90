program condenerdiss
  use mconstant
  use mfftwplan3d
  implicit none
  
  integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel

  complex(sp), allocatable, dimension(:,:,:) :: ux, uy, uz, vx, vy, vz
  complex(sp), allocatable, dimension(:,:,:) :: tij
  complex(sp), allocatable, dimension(:,:,:) :: t11, t12, t13, t22, t23, t33
  real(sp), allocatable, dimension(:,:,:) :: k2, g
  real(sp), allocatable, dimension(:) :: kx, ky, kz

  real(sp) :: delta_c, const, ignore_me

  real(dp) :: mdt11, mdt12, mdt13, mdt22, mdt23, mdt33
  real(dp) :: vdt11, vdt12, vdt13, vdt22, vdt23, vdt33
  real(dp) :: t11, t12, t13, t22, t23, t33
  real(dp) :: vt11, vt12, vt13, vt22, vt23, vt33
  real(dp) :: mdtii, mdtij, mmtii, mmtij, vdtii, vdtij, vmtii, vmtij
  real(dp) :: covt11, covt12, covt13, covt22, covt23, covt33, covtii, covtij
  real(dp) :: rhot11, rhot12, rhot13, rhot22, rhot23, rhot33, rhotii, rhotij

  character(80) :: str, fnm, ndellist, str1, fpath

  write(*,*)
  write(*,'('' >>> sgs ener dissipation conditioned on q r plane for several model <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 5) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./sgs-ediss-condqr.x nx filelist ndellist filtertype modeltype'
          write(*,*) '        nx: resolution'
          write(*,*) '        filelist: the list of data files, *.list'
          write(*,*) '        ndellist: list of filter scales ndel: Delta = ndel * dx'
          write(*,*) '        filtertype: 0 for Gaussian 1 for cutoff'
          write(*,*) '        modeltype: 0 DNS; 1 Smag; 2 similarity; 3 nonlinear; 4 mat exp'
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  ! ============= Getting command line arguments =================

  ! resolution
  call getarg(1,str)
  read(str, '(I10)') nx

  ! filelist string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ! filter scale
  call getarg(3,ndellist)
  ndellist = adjustl(ndellist)

  ! filter type
  call getarg(4,str)
  read(str, '(I10)') mm

  ! model type
  call getarg(5,str)
  read(str, '(I10)') nn

  ! ===============================================================


  ! ================= Strings for output =======================

  ! Define strings for output for chosen model.
  select case(nn)
  case (1)
    str = '-smag-'//fnm(1:len_trim(fnm))//'.dat'
  case (2)
    str = '-sim-'//fnm(1:len_trim(fnm))//'.dat'
  case (3)
    str = '-nl-'//fnm(1:len_trim(fnm))//'.dat'
  case (4)
    str = '-matexp-'//fnm(1:len_trim(fnm))//'.dat'
  case default
    write(*,*) 'Wrong model type. Stopping'
    stop
  end select

  ! file name for file list
  fnm = fnm(1:len_trim(fnm))//'.list'

  ! ============================================================

  ny = nx; nz = nx
  lx = nx / 2; ly = ny; lz = nz; lx1 = lx + 1
  const = 1./ (nx * ny * nz)

  call fftwplan3d(nx,ny,nz) ! fftwplan3de: estimated plan
  write(*,*) 'after fftwplan3d'


  allocate( ux(lx1,ly,lz),   uy(lx1,ly,lz),   uz(lx1,ly,lz)  )
  allocate( vx(lx1,ly,lz),   vy(lx1,ly,lz),   vz(lx1,ly,lz)  )
  allocate( t11(lx1,ly,lz), t12(lx1,ly,lz), t13(lx1,ly,lz))
  allocate( t22(lx1,ly,lz), t23(lx1,ly,lz), t33(lx1,ly,lz))
  allocate( tij(lx1,ly,lz),  k2(lx1,ly,lz),   g(lx1,ly,lz)   )
  allocate( kx(lx1), ky(ly), kz(lz) )

  call wavenumber(kx,ky,kz,k2,lx1,ly,lz)
  write(*,*) 'after wavenumber'

  open(32, file = 'sgs-ediss-condqr-'//str( 1:len_trim(str) ) )

  open(25, file = ndellist(1:len_trim(ndellist))//'.list')

    do while ( .not. eof(25) )

      read(25,*) str1
      read(str1, '(I20)') ndel
      write(*,*) 'ndel = ', str1(1:len_trim(str1))
  
      delta_c=ndel*2*pi/nx
      write(*,*) 'delta_c = ', delta_c
  
      ! =============== Defining filters ===================
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
      ! ====================================================
  
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

          ! ==================== Calculate subgrid-scale stress ===========================
          if ( nn .eq. 0) then
            ! DNS sgs stress
            call rsgstauij(ux,ux,t11,g,nx,ny,nz) 
            call rsgstauij(uy,uy,t22,g,nx,ny,nz) 
            call rsgstauij(uz,uz,t33,g,nx,ny,nz) 
            call rsgstauij(ux,uy,t12,g,nx,ny,nz) 
            call rsgstauij(ux,uz,t13,g,nx,ny,nz) 
            call rsgstauij(uy,uz,t23,g,nx,ny,nz) 
          end if
          ux = ux * g
          uy = uy * g
          uz = uz * g

          select case (nn)
          case (1)
            ! Smagorinsky model
            call rsmag(ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,lx1,ly,lz,nx,ny,nz,delta_c)
          case (2)
            ! Similarity model
            call rsgstauij(ux,ux,t11,g,nx,ny,nz)
            call rsgstauij(ux,uy,t12,g,nx,ny,nz)
            call rsgstauij(ux,uz,t13,g,nx,ny,nz)
            call rsgstauij(uy,uy,t22,g,nx,ny,nz)
            call rsgstauij(uy,uz,t23,g,nx,ny,nz)
            call rsgstauij(uz,uz,t33,g,nx,ny,nz)
            s12 = - (t11 + t22 + t33)/3.
            t11 = t11 + s12; t22 = t22 + s12; t33 = t33 + s12
          case (3)
            ! Nonlinear model
            write(*,*) 'not implemented yet'
            stop
          case (4)
            ! Matrix exponential model
            call rmatexp(ux,uy,uz,t11,t12,t13,t22,t23,t33,kx,ky,kz,lx1,ly,lz,nx,ny,nz,delta_c)
          end select

          s12 = eye * .5 * ( kx * uy + ky * ux )
          call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
          diss = - s12

  
          nfile = nfile + 1
        end do
  
      close(20)

      ignore_me = const / nfile
  
      write(28, '(I4, 30E15.4)') ndel, t11, t12, t13, t22, t23, t33, mmtij
      write(29, '(I4, 30E15.4)') ndel, vt11, vt12, vt13, vt22, vt23, vt33, vmtij
      write(32, '(I4, 30E15.4)') ndel, covt11, covt12, covt13, covt22, covt23, &
        covt33, covtij
      write(33, '(I4, 30E15.4)') ndel, rhot11, rhot12, rhot13, rhot22, rhot23, &
        rhot33, rhotij

  
    end do

  close(25)

  close(28)
  close(29)
  close(32)
  close(33)

  call destroyplan3d

  deallocate( ux, uy, uz, vx, vy, vz )
  deallocate( t11, t12, t13, t22, t23, t33 )
  deallocate( tij, k2, g, kx, ky, kz )

  write(*,*) 'finished'

end program condenerdiss 
