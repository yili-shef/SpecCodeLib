program horiuti
    use mconstant
    use mfftwplan3d
    use mwavenumber

    integer :: nx,ny,nz,lx,lx1,ly,lz,ii,jj,kk,ll,mm,nn,nfile,ndel

    complex(sp), allocatable, dimension(:,:,:) :: s11,s12,s13,s22,s23,s33, wx, wy, wz
    real(sp),    allocatable, dimension(:,:,:) :: g
    real(sp),    allocatable, dimension(:)     :: kx,ky,kz

    integer, parameter :: matz=5 
    real(dp), dimension(3)   :: evalues, fv1, fv2, omega, odotev
    real(dp), dimension(3,3) :: cc, evcij
    integer :: ierr, mloc
    integer(8) :: ncurvesheet, ntube, nflatsheet
    character(80) :: str1, str, fnm, fpath
    real(sp) :: ignore_me, delta_c, const

    write(*,*) 
    write(*,'(''>>> Horiuti classification <<<'')')
    write(*,*)
    ll=iargc()
    if (ll .ne. 3) then
        write(*,*)
        write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
        write(*,*) 
        write(*,*) ' Usage: ./horiuticlassification.x nx filelist ndel'
        write(*,*) '        nx: resolution of data'
        write(*,*) '        filelist: the list of the data files, *.list'
        write(*,*) '        ndel: filter scale Delta=ndel*dx'
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
    call getarg(3,str)
    read(str,'(I20)') ndel
    str=adjustl(str)

    ny=nx; nz=nx
    lx=nx/2; ly=ny;lz=nz;lx1=lx+1
    const=1./(nx*ny*nz)
 
    delta_c=ndel*2*pi/nx
 
    call fftwplan3de(nx,ny,nz)
    write(*,*) 'after fftwplan3d'
 
    allocate( s11(lx1,ly,lz), s12(lx1,ly,lz), s13(lx1,ly,lz) )
    allocate( s22(lx1,ly,lz), s23(lx1,ly,lz), s33(lx1,ly,lz) )
    allocate(  wx(lx1,ly,lz),  wy(lx1,ly,lz),  wz(lx1,ly,lz) )
    allocate( kx(lx1), ky(ly), kz(lz) )
    allocate( g(lx1,ly,lz) )
    write(*,*) 'arrays allocated'

    call wavenumber(kx,ky,kz,g,lx1,ly,lz)
    write(*,*) 'after wavenumber'

    ! Gaussian filter
    g=exp(-g*delta_c**2/24.)

    nfile = 0
    ncurvesheet = 0
    nflatsheet  = 0
    ntube       = 0
    open(20, file = fnm(1:len_trim(fnm))//'.list')

        do while ( .not. eof(20)) 

            read(20,*) str1
            write(*,*) str1(1:len_trim(str1))
           
            ! s11 s22 s33 are velocities
            fpath='./out/ux'//str1(1:len_trim(str1))
            open(10,file=fpath,status='unknown',form='unformatted')
              read(10) s11
            close(10)
            fpath='./out/uy'//str1(1:len_trim(str1))
            open(10,file=fpath,status='unknown',form='unformatted')
              read(10) s22
            close(10)
            fpath='./out/uz'//str1(1:len_trim(str1))
            open(10,file=fpath,status='unknown',form='unformatted')
              read(10) s33
            close(10)
            write(*,*) 'after reading data files'

            s11 = s11 * g
            s22 = s22 * g
            s33 = s33 * g

            do kk = 1, lz
            do jj = 1, ly
            do ii = 1, lx1
   
              s12(ii,jj,kk) = .5*eye*(kx(ii)*s22(ii,jj,kk)+ky(jj)*s11(ii,jj,kk))
              s13(ii,jj,kk) = .5*eye*(kx(ii)*s33(ii,jj,kk)+kz(kk)*s11(ii,jj,kk))
              s23(ii,jj,kk) = .5*eye*(ky(jj)*s33(ii,jj,kk)+kz(kk)*s22(ii,jj,kk))

              wx(ii,jj,kk) = eye * ( ky(jj) * s33(ii,jj,kk) - kz(kk) * s22(ii,jj,kk) ) 
              wy(ii,jj,kk) = eye * ( kz(kk) * s11(ii,jj,kk) - kx(ii) * s33(ii,jj,kk) )
              wz(ii,jj,kk) = eye * ( kx(ii) * s22(ii,jj,kk) - ky(jj) * s11(ii,jj,kk) )

              s11(ii,jj,kk) = eye*kx(ii)*s11(ii,jj,kk)
              s22(ii,jj,kk) = eye*ky(jj)*s22(ii,jj,kk)
              s33(ii,jj,kk) = eye*kz(kk)*s33(ii,jj,kk)

            end do
            end do
            end do
  
            call rfftwnd_f77_one_complex_to_real(c2r3d,s11,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d,s12,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d,s13,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d,s22,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d,s23,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d,s33,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d, wx,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d, wy,ignore_me)
            call rfftwnd_f77_one_complex_to_real(c2r3d, wz,ignore_me)


            do kk=1,nz
            do jj=1,ny
            do ii=1,nx
 
              if ( mod(ii,2) .eq. 1) then
                  ll = (ii + 1)/2
         
                  cc(1,1) = real(s11(ll,jj,kk))
                  cc(1,2) = real(s12(ll,jj,kk))
                  cc(1,3) = real(s13(ll,jj,kk))
                  cc(2,2) = real(s22(ll,jj,kk))
                  cc(2,3) = real(s23(ll,jj,kk))
                  cc(3,3) = real(s33(ll,jj,kk))

                  omega(1) = real( wx(ll,jj,kk) )
                  omega(2) = real( wy(ll,jj,kk) )
                  omega(3) = real( wz(ll,jj,kk) )
         
              else
                  ll = ii/2

                  cc(1,1) = aimag(s11(ll,jj,kk))
                  cc(1,2) = aimag(s12(ll,jj,kk))
                  cc(1,3) = aimag(s13(ll,jj,kk))
                  cc(2,2) = aimag(s22(ll,jj,kk))
                  cc(2,3) = aimag(s23(ll,jj,kk))
                  cc(3,3) = aimag(s33(ll,jj,kk))
          
                  omega(1) = aimag( wx(ll,jj,kk) )
                  omega(2) = aimag( wy(ll,jj,kk) )
                  omega(3) = aimag( wz(ll,jj,kk) )
         
              end if
              cc(2,1) = cc(1,2)
              cc(3,1) = cc(1,3)
              cc(3,2) = cc(2,3)

              cc = matmul(cc, cc) 
              do mm = 1, 3
              do nn = 1, 3
                cc(mm,nn) = cc(mm,nn) + .25 * omega(mm) * omega(nn)
              end do
              end do

              cc(1,1) = cc(1,1) - sum(omega*omega) / 4
              cc(2,2) = cc(2,2) - sum(omega*omega) / 4
              cc(3,3) = cc(3,3) - sum(omega*omega) / 4

              call rs(3,3,cc,evalues,matz,evcij,fv1,fv2,ierr)
              do ll = 1, 3
                evcij(:,ll) = evcij(:,ll) / sqrt( sum( evcij(:,ll)**2 ) )
                odotev(ll) = abs( dot_product( evcij(:,ll), omega ) )
              end do
              mloc = maxloc( odotev, 1 )

              select case (mloc)
              case (1)
                  if ( evalues(2) .gt. 0 .and. evalues(3) .gt. 0 ) then
                      ncurvesheet = ncurvesheet + 1
                  else if ( evalues(2) .lt. 0 .and. evalues(3) .lt. 0 ) then 
                      ntube = ntube + 1
                  else
                      nflatsheet = nflatsheet + 1
                  end if
              case (2)
                  if ( evalues(1) .gt. 0 .and. evalues(3) .gt. 0 ) then
                      ncurvesheet = ncurvesheet + 1
                  else if ( evalues(1) .lt. 0 .and. evalues(3) .lt. 0 ) then 
                      ntube = ntube + 1
                  else
                      nflatsheet = nflatsheet + 1
                  end if
              case (3)
                  if ( evalues(1) .gt. 0 .and. evalues(2) .gt. 0 ) then
                      ncurvesheet = ncurvesheet + 1
                  else if ( evalues(1) .lt. 0 .and. evalues(2) .lt. 0 ) then 
                      ntube = ntube + 1
                  else
                      nflatsheet = nflatsheet + 1
                  end if
              case default
                  stop 'Something wrong'
              end select
        
            end do
            end do
            end do

            nfile = nfile + 1

        end do

    close(20)

    write(*,*) '% flat sheet:  ', real(nflatsheet,  dp) / (nfile * nx * ny * nz)
    write(*,*) '% curve sheet: ', real(ncurvesheet, dp) / (nfile * nx * ny * nz)
    write(*,*) '% tube:        ', real(ntube,       dp) / (nfile * nx * ny * nz)


    call destroyplan3d

    deallocate(s11, s12, s13, s22, s23, s33, wx, wy, wz, g, kx, ky, kz)

    write(*,*) 'Finished'

end program horiuti      
