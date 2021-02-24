program cndsgshd
  use mconstant
  implicit none

  integer :: nx,ny,nz,ii,jj,kk,ll,mm,nn,lx,lx1,ly,lz,nfile
  real, allocatable, dimension(:,:,:) :: alpha, beta, gmma
  complex, allocatable, dimension(:,:,:) :: dissrate

  integer, parameter :: npnt=100
  real,    parameter :: bnd=8., binw1=2./npnt, binw2=bnd/npnt
  real(dp), dimension(npnt) :: pbe, psat, cndhd, cndhdsat

  real(dp) :: rmsr, const, mdissrate
  real :: rmsr0, mdissrate0, atmp, btmp, gtmp, ss, disstmp
  character(80) :: fnm,str,str1,disslist

  write(*,*) 
  write(*,'(''>>>>>> mean SGS helicity dissipation conditioned on eigenvalues of rij <<<<<<'')')
  write(*,*)

  ll=iargc()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./cndsgshd.x nx evlist disslist normdissrate'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        evlist: the list of the eigenvalue data files'
          write(*,*) '        disslist: dissipation data file list'
          write(*,*) '        normdissrate: normalization factor for dissipation, normally the mean'
          write(*,*)
          write(*,*) ' Stopped'
          stop
  end if
  ! resolution
  call getarg(1,str)
  read(str, '(I20)') nx

  ! evlist string
  call getarg(2,fnm)
  fnm = adjustl(fnm)

  ! disslist 
  call getarg(3,disslist)
  disslist = adjustl(disslist)

  call getarg(4,str)
  read(str, '(F10.6)') mdissrate0

  str='-'//fnm(1:len_trim(fnm))//'.dat'
  disslist = disslist( 1 : len_trim( disslist ) )//'.list'

  lx = nx / 2; lx1 = lx + 1; ly = nx; lz = nx
  ny = nx; nz = nx
  const = 1.d0 / (nx * ny * nz)

  allocate(alpha(nx,ny,nz),beta(nx,ny,nz),gmma(nx,ny,nz))
  allocate(dissrate(lx1,ly,lz))

  ! Start calculation
  open(20, file = fnm( 1 : len_trim(fnm))//'.list')
  open(21, file = disslist( 1 : len_trim( disslist ) ) )

    mdissrate = 0.d0
    rmsr=0.d0
    psat=0.d0
    pbe=0.d0 
    cndhd = 0.d0
    cndhdsat = 0.d0

    nfile = 0
    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))
 
      open(15,file='./out/'//str1(1:len_trim(str1)),form='unformatted')
        read(15) alpha
        read(15) beta
        read(15) gmma
      close(15)

      read(21,*) str1
      write(*,*) str1(1:len_trim(str1))
 
      open(15,file='./out/'//str1(1:len_trim(str1)),form='unformatted')
        read(15) dissrate
      close(15)
  
      write(*,*) 'after reading data files'

      if (nfile .eq. 0) then
        ! Estimate RijRij
        rmsr0 = sum( alpha*alpha + beta*beta + gmma*gmma ) * const
        rmsr0 = sqrt( rmsr0 )
        write(*,*) 'Estimate of rmsr = ', rmsr0
      end if

      do kk=1,nz
      do jj=1,ny
      do ii=1,nx
        if ( mod(ii, 2) .eq. 1) then
              disstmp = real( dissrate( (ii+1)/2, jj, kk ) )
        else 
              disstmp = aimag( dissrate( ii/2, jj,kk ) )
        end if

        mdissrate = mdissrate + disstmp

        atmp = alpha(ii,jj,kk)
        btmp = beta(ii,jj,kk)
        gtmp = gmma(ii,jj,kk)
        
        ss =  atmp * atmp + btmp * btmp + gtmp * gtmp 
        rmsr = rmsr + ss 
         
        ss = -3 * sqrt(6.) * atmp * btmp * gtmp / ss**1.5  ! s^* defined by Lund and Rogers
       
        ll = floor( (ss + 1) / binw1 ) + 1
        if (ll .ge. 1 .and. ll .le. npnt) then 
          psat(ll) = psat(ll) + 1
          cndhdsat(ll) = cndhdsat(ll) + disstmp
        end if
       
       
        btmp = btmp / rmsr0
        mm = floor( (btmp + bnd/2) / binw2 ) + 1                       ! beta in [-bnd/2,bnd/2]
        if (mm .ge. 1 .and. mm .le. npnt) then
          pbe(mm) = pbe(mm) + 1
          cndhd(mm) = cndhd(mm) + disstmp
        end if

      end do
      end do
      end do

      nfile = nfile + 1
 
    end do

  close(20)
  close(21)

  write(*,*) 'finished loops'

  const = const / nfile

  cndhd = cndhd / (pbe + tiny)
  cndhdsat = cndhdsat / (psat + tiny)

  cndhd = cndhd / mdissrate0
  cndhdsat = cndhdsat / mdissrate0
  mdissrate = mdissrate * const
  write(*,*) 'input normalization dissrate: ', mdissrate0
  write(*,*) 'check mead dissrate: ', mdissrate
  
  psat = psat * const 
  pbe = pbe * const

  write(*,*) 'check psat: ', sum(psat)
  write(*,*) 'check pbe: ', sum(pbe)

  psat=psat/binw1
  pbe=pbe/binw2

  rmsr = rmsr * const
  rmsr = sqrt(rmsr)

  write(*,*) 'rmsr: ', rmsr

  open(15,file='cndhd'//str(1:len_trim(str)))
    do ii=1,npnt
      write(15,'(15E15.5)') (-bnd*.5+(ii-.5)*binw2)*rmsr0/rmsr, pbe(ii)*rmsr/rmsr0, cndhd(ii), &
        -1 + (ii-.5) * binw1, psat(ii), cndhdsat(ii)
    end do
  close(15)

  deallocate(alpha,beta,gmma)
  deallocate(dissrate)

  write(*,*) 'finished'

end program cndsgshd
