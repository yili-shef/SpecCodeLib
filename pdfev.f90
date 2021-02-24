program pdfev
  use mconstant
  implicit none

  integer :: nx,ny,nz,ii,jj,kk,ll,mm,nn
  real(sp), allocatable, dimension(:,:,:) :: alpha, beta, gmma

  integer,  parameter :: npnt=100
  real(sp), parameter :: bnd=10., binw1=2./npnt, binw2=bnd/npnt
  real(dp), dimension(npnt) :: pal, pbe, pga, psat
  real(dp), dimension(npnt,npnt) :: pab, pag, pbg

  integer  :: nfile
  real(dp) :: rmsr, malpha, mbeta, mgamma, ralpha, rbeta, rgamma
  real(dp) :: sbeta
  real(sp) :: rmsr0, atmp, btmp, gtmp, ss
  real(dp) :: const
  character(80) :: fnm,str,str1

  write(*,*) 
  write(*,'(''>>>>>> Joint PDFs of the eigenvalues of given data file<<<<<<'')')
  write(*,*)

  ll=iargc()
  if (ll .ne. 2) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./pdfev.x nx filelist'
          write(*,*) '        nx: resolution of data'
          write(*,*) '        filelist: the list of the data files'
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


  str='-'//fnm(1:len_trim(fnm))//'.dat'

  ny=nx; nz=nx
  const=1.0_dp/(nx*ny*nz)
  allocate(alpha(nx,ny,nz),beta(nx,ny,nz),gmma(nx,ny,nz))

  
  ! Start calculation
  open(20,file=fnm(1:len_trim(fnm))//'.list')

    malpha = 0.0_dp; mbeta = 0.0_dp; mgamma = 0.0_dp
    ralpha = 0.0_dp; rbeta = 0.0_dp; rgamma = 0.0_dp
    rmsr = 0.0_dp
    sbeta = 0.0_dp; 
    psat = 0.0_dp
    pal = 0.0_dp; pbe=0.0_dp; pga=0.0_dp
    pab=0.d0; pag=0.d0; pbg=0.d0

    nfile=0
    do while ( .not. eof(20)) 
      read(20,*) str1
      write(*,*) str1(1:len_trim(str1))
 
      open(15,file='./out/'//str1(1:len_trim(str1)),form='unformatted')
        read(15) alpha
        read(15) beta
        read(15) gmma
      close(15)

      if (nfile .eq. 0) then
        ! Estimate RijRij
        rmsr0 = sum( alpha*alpha + beta*beta + gmma*gmma ) * const
        rmsr0 = sqrt( rmsr0 )
        write(*,*) 'Estimate of rmsr =', rmsr0
      end if

  
      do kk=1,nz
      do jj=1,ny
      do ii=1,nx
        atmp = alpha(ii,jj,kk)
        btmp = beta(ii,jj,kk)
        gtmp = gmma(ii,jj,kk)

        ss =  atmp * atmp + btmp * btmp + gtmp * gtmp 

        !!!! a rms of the eigenvalues
        rmsr = rmsr + ss 

        ss = -3 * sqrt(6.) * atmp * btmp * gtmp / ss**1.5  ! s^* defined by Lund and Rogers
  
        ll = floor( (ss + 1) / binw1 ) + 1
        if (ll .ge. 1 .and. ll .le. npnt) psat(ll) = psat(ll) + 1


        malpha = malpha + atmp
        mbeta = mbeta + btmp
        mgamma = mgamma + gtmp
       
        ralpha = ralpha + atmp * atmp
        rbeta = rbeta + btmp * btmp
        rgamma = rgamma + gtmp * gtmp

        sbeta = sbeta + btmp * btmp * btmp
        
        atmp = atmp / rmsr0
        btmp = btmp / rmsr0
        gtmp = gtmp / rmsr0
  
        ll = floor( atmp / binw2 ) + 1                           ! alpha in [0,bnd]
        mm = floor( (btmp + bnd/2) / binw2 ) + 1                       ! beta in [-bnd/2,bnd/2]
        nn = floor( (gtmp + bnd) / binw2 ) + 1                        ! gmma in [-bnd,0]


        if (ll .ge. 1 .and. ll .le. npnt .and. & 
            mm .ge. 1 .and. mm .le. npnt .and. &
            nn .ge. 1 .and. nn .le. npnt) then
          pal(ll) = pal(ll) + 1
          pbe(mm) = pbe(mm) + 1
          pga(nn) = pga(nn) + 1
          pab(ll,mm) = pab(ll,mm) + 1
          pag(ll,nn) = pag(ll,nn) + 1
          pbg(mm,nn) = pbg(mm,nn) + 1
        end if

      end do
      end do
      end do
 
      nfile = nfile + 1

    end do

  close(20)

  write(*,*) 'finished loops'

  psat = psat * const / nfile
  pal = pal * const / nfile

  pbe = pbe * const / nfile
  pga=pga*const/nfile
  pab=pab*const/nfile
  pag=pag*const/nfile
  pbg=pbg*const/nfile
  write(*,*) 'check psat: ', sum(psat)
  write(*,*) 'check pal: ', sum(pal)
  write(*,*) 'check pbe: ', sum(pbe)
  write(*,*) 'check pga: ', sum(pga)
  write(*,*) 'check pab: ', sum(pab)
  write(*,*) 'check pag: ', sum(pag)
  write(*,*) 'check pbg: ', sum(pbg)
  psat=psat/binw1
  pal=pal/binw2
  pbe=pbe/binw2
  pga=pga/binw2
  pab=pab/binw2/binw2
  pag=pag/binw2/binw2
  pbg=pbg/binw2/binw2

  rmsr=rmsr*const/nfile
  rmsr=sqrt(rmsr)

  write(*,*) 'rmsr: ', rmsr

  malpha = malpha * const / nfile
  mbeta = mbeta * const / nfile
  mgamma = mgamma * const / nfile

  ralpha = ralpha * const / nfile
  rbeta = rbeta * const / nfile
  rgamma = rgamma * const / nfile

  sbeta = sbeta * const / nfile

  sbeta = sbeta - 3 * rbeta * mbeta + 2 * mbeta * mbeta * mbeta

  ralpha = sqrt( ralpha - malpha * malpha )
  rbeta = sqrt( rbeta - mbeta * mbeta )
  rgamma = sqrt( rgamma - mgamma * mgamma )

  sbeta = sbeta / rbeta**3

  open(15, file = 'ms'//str(1:len_trim(str)) )
    write(15,*) 'malpha mbeta mgamma ralpha rbeta rgamma rmsr sbeta'
    write(15,*) malpha, mbeta, mgamma, ralpha, rbeta, rgamma, rmsr, sbeta
  close(15) 
  
  open(15,file='psat'//str(1:len_trim(str)))
    do ii=1,npnt
      write(15,*) -1.+(ii-.5)*binw1, psat(ii)
    end do
  close(15)

  open(15,file='pev'//str(1:len_trim(str)))
    do ii=1,npnt
      write(15,'(15E15.5)') (ii-.5)*binw2*rmsr0/rmsr, pal(ii)*rmsr/rmsr0, &
                  (-bnd*.5+(ii-.5)*binw2)*rmsr0/rmsr, pbe(ii)*rmsr/rmsr0, &
                  (-bnd+(ii-.5)*binw2)*rmsr0/rmsr, pga(ii)*rmsr/rmsr0
    end do
  close(15)

  open(15,file='jpev'//str(1:len_trim(str)))
    write(15,*) 'Zone T= "(alpha,beta)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,*) (ii-.5)*binw2*rmsr0/rmsr, (-bnd*.5+(jj-.5)*binw2)*rmsr0/rmsr, &
                  pab(ii,jj)*rmsr*rmsr/rmsr0/rmsr0
    end do
    end do

    write(15,*) 'Zone T= "(alpha,gmma)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,*) (ii-.5)*binw2*rmsr0/rmsr, (-bnd+(jj-.5)*binw2)*rmsr0/rmsr, &
                  pag(ii,jj)*rmsr*rmsr/rmsr0/rmsr0
    end do
    end do

    write(15,*) 'Zone T= "(beta,gmma)", i=', npnt, ', j=', npnt, ', F=point'
    do jj=1,npnt
    do ii=1,npnt
      write(15,*) (-bnd*.5+(ii-.5)*binw2)*rmsr0/rmsr, (-bnd+(jj-.5)*binw2)*rmsr0/rmsr, &
                  pbg(ii,jj)*rmsr*rmsr/rmsr0/rmsr0
    end do
    end do
  close(15)
 
  deallocate(alpha,beta,gmma)

  write(*,*) 'pdfev.x finished'

end program pdfev
