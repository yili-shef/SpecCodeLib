program rcmalignos
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,three,npnt) :: aijall
  real(dp), dimension(three,three) :: sij, evsij
  real(dp), dimension(three) :: evalues

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: fv1,fv2
  integer :: ierr

  integer, parameter :: npdf = 40 
  real(sp), parameter :: bnd = 1., binw = 2.*bnd/npdf
  real(dp), dimension(npdf) :: palpha, pbeta, pgamma, psstar

  integer :: ii, jj, kk, ll, nreal, ifirst
  real(dp) :: tmp, omx, omy, omz, osalpha, osbeta, osgamma, sstar
  character(80) :: str, strro, strreal

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rcm-alignos.x filelist nreal ifirst'
          write(*,*) '        filelist: filelist.list is the list of data files'
          write(*,*) '        nreal: number of realizations'
          write(*,*) '        ifirst: the index of the first realization'
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  call getarg(1,strro)
  strro = adjustl(strro)
  call getarg(2,strreal)
  read(strreal, '(I20)') nreal
  call getarg(3,str)
  read(str, '(I20)') ifirst

  str = 'alignos-'//strro(1:len_trim(strro))//'.dat'
  open(20, file = str(1:len_trim(str)) )

  open(22, file = strro(1:len_trim(strro))//'.list' )
  jj = 0
  do while ( .not. eof(22))

    read(22,*) str
    write(*, *) str
    str = adjustl(str)
 
    palpha = 0._dp
    pbeta = 0._dp
    pgamma = 0._dp
    psstar = 0._dp

    do kk = ifirst, nreal + ifirst -1
      write(*,*) 'realization ', kk
      write(strreal, '(I20)') kk
      strreal = adjustl(strreal)
      strreal = str(1:len_trim(str))//strreal(1:len_trim(strreal))//'.data'
 
      open(10, file = strreal(1:len_trim(strreal)), form = 'binary')
        read(10) aijall
      close(10)
    
      do ii = 1, npnt

        sij = aijall(:,:,ii)
        omx = sij(3,2) - sij(2,3)
        omy = sij(1,3) - sij(3,1)
        omz = sij(2,1) - sij(1,2)
        tmp = sqrt( omx * omx + omy * omy + omz * omz )
        omx = omx / tmp
        omy = omy / tmp
        omz = omz / tmp

        sij = ( sij + transpose(sij) ) * 0.5
     
        call rs(3,3,sij,evalues,matz,evsij,fv1,fv2,ierr)
        do ll=1,3
          evsij(:,ll)=evsij(:,ll)/sqrt(sum(evsij(:,ll)**2))
        end do

        osalpha = omx * evsij(1,3) + omy * evsij(2,3) + omz * evsij(3,3)
        osbeta  = omx * evsij(1,2) + omy * evsij(2,2) + omz * evsij(3,2)
        osgamma = omx * evsij(1,1) + omy * evsij(2,1) + omz * evsij(3,1)
     
        sstar = -3. *sqrt(6.)*evalues(1)*evalues(2)*evalues(3) &
                /sqrt( ( evalues(1)**2 + evalues(2)**2 + evalues(3)**2 )**3 )

        ll = floor( (osalpha+1) / binw ) + 1
        if ( ll .ge. 1 .and. ll .le. npnt) palpha(ll) = palpha(ll) + 1

        ll = floor( (osbeta + 1) / binw ) + 1
        if (ll .ge. 1 .and. ll .le. npnt) pbeta(ll) = pbeta(ll) + 1

        ll = floor( (osgamma + 1) / binw ) + 1
        if (ll .ge. 1 .and. ll .le. npnt) pgamma(ll) = pgamma(ll) + 1

        ll = floor( ( sstar + 1 ) / binw ) + 1
        if (ll .ge. 1 .and. ll .le. npnt) psstar(ll) = psstar(ll) + 1

      end do

    end do

    tmp = 1. / nreal / npnt

    palpha = palpha * tmp
    pbeta  = pbeta  * tmp
    pgamma = pgamma * tmp
    psstar = psstar * tmp

    write(*,*) 'Check normalization palpha: ', sum(palpha)
    write(*,*) 'Check normalization pbeta:  ', sum(pbeta)
    write(*,*) 'Check normalization pgamma: ', sum(pgamma)
    write(*,*) 'Check normalization psstar: ', sum(psstar)

    palpha = palpha / binw
    pbeta  = pbeta  / binw
    pgamma = pgamma / binw
    psstar = psstar / binw
 
    write(20,'(''# zone T = " step '',I4, ''", I='', I4, '' F = point'')') jj, npdf
    do ii = 1, npdf
      tmp = -bnd + (ii-.5) * binw
      write(20,'(15E15.4)') tmp, palpha(ii), pbeta(ii), pgamma(ii), psstar(ii)
    end do
    write(20, *)
 
    jj = jj + 1

  end do
  close(20)
  close(22)

  write(*,*) 'finished'

end program rcmalignos 
