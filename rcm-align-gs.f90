program rcmaligngs
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,three,npnt) :: aijall
  real(sp), dimension(three,npnt) :: giall
  real(dp), dimension(three,three) :: sij, evsij
  real(dp), dimension(three) :: evalues, gi

  integer, parameter :: matz=5 
  real(dp), dimension(3) :: fv1,fv2
  integer :: ierr

  integer, parameter :: npdf = 40 
  real(sp), parameter :: bnd = 1., binw = 2.*bnd/npdf
  real(dp), dimension(npdf) :: palpha, pbeta, pgamma

  integer :: ii, jj, kk, ll, nreal, ifirst
  real(dp) :: tmp, gsalpha, gsbeta, gsgamma
  character(80) :: str, straij, strgi, strreal

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rcm-aligngs.x aijfilelist gifilelist nreal ifirst'
          write(*,*) '        aijfilelist: aijfilelist.list is the list of aij data files'
          write(*,*) '        gifilelist: gifilelist.list is the list of gi data files'
          write(*,*) '        nreal: number of realizations'
          write(*,*) '        ifirst: the index of the first realization'
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  call getarg(1,straij)
  straij = adjustl(straij)
  call getarg(2,strgi)
  strgi = adjustl(strgi)
  call getarg(3,strreal)
  read(strreal, '(I20)') nreal
  call getarg(4,str)
  read(str, '(I20)') ifirst

  str = 'aligngs-'//strgi(1:len_trim(strgi))//'.dat'
  open(20, file = str(1:len_trim(str)) )

  open(22, file = straij(1:len_trim(straij))//'.list' )
  open(23, file = strgi (1:len_trim(strgi ))//'.list' )

  jj = 0
  do while ( .not. eof(22) )

    read(22,*) straij
    write(*,*) straij
    straij = adjustl(straij)

    read(23,*) strgi
    write(*,*) strgi
    strgi = adjustl(strgi)
 
    palpha = 0._dp
    pbeta = 0._dp
    pgamma = 0._dp
    do kk = ifirst, nreal + ifirst -1

      write(*,*) 'realization ', kk
      write(strreal, '(I20)') kk
      strreal = adjustl(strreal)

      str = straij(1:len_trim(straij))//strreal(1:len_trim(strreal))//'.data'
      open(10, file = str(1:len_trim(str)), form = 'binary')
        read(10) aijall
      close(10)
    
      str = strgi (1:len_trim(strgi ))//strreal(1:len_trim(strreal))//'.data'
      open(10, file = str(1:len_trim(str)), form = 'binary')
        read(10) giall
      close(10)
    
      do ii = 1, npnt

        gi = giall(:,ii)
        gi = gi/ sqrt( sum(gi * gi) )

        sij = aijall(:,:,ii)
        sij = ( sij + transpose(sij) ) * 0.5
     
        call rs(3,3,sij,evalues,matz,evsij,fv1,fv2,ierr)
        do ll=1,3
          evsij(:,ll)=evsij(:,ll)/sqrt(sum(evsij(:,ll)**2))
        end do

        gsalpha = gi(1) * evsij(1,3) + gi(2) * evsij(2,3) + gi(3) * evsij(3,3)
        gsbeta  = gi(1) * evsij(1,2) + gi(2) * evsij(2,2) + gi(3) * evsij(3,2)
        gsgamma = gi(1) * evsij(1,1) + gi(2) * evsij(2,1) + gi(3) * evsij(3,1)
     
        ll = floor( (gsalpha + 1) / binw ) + 1
        if ( ll .ge. 1 .and. ll .le. npnt) palpha(ll) = palpha(ll) + 1

        ll = floor( (gsbeta + 1) / binw ) + 1
        if (ll .ge. 1 .and. ll .le. npnt) pbeta(ll) = pbeta(ll) + 1

        ll = floor( (gsgamma + 1) / binw ) + 1
        if (ll .ge. 1 .and. ll .le. npnt) pgamma(ll) = pgamma(ll) + 1

      end do

    end do

    tmp = 1. / nreal / npnt

    palpha = palpha * tmp
    pbeta  = pbeta  * tmp
    pgamma = pgamma * tmp

    write(*,*) 'Check normalization palpha: ', sum(palpha)
    write(*,*) 'Check normalization pbeta:  ', sum(pbeta)
    write(*,*) 'Check normalization pgamma: ', sum(pgamma)

    palpha = palpha / binw
    pbeta  = pbeta  / binw
    pgamma = pgamma / binw
 
    write(20,'(''# zone T = " step '',I4, ''", I='', I4, '' F = point'')') jj, npdf
    do ii = 1, npdf
      tmp = -bnd + (ii-.5) * binw
      write(20,'(15E15.4)') tmp, palpha(ii), pbeta(ii), pgamma(ii)
    end do
    write(20, *)
 
    jj = jj + 1

  end do
  close(20)
  close(22)
  close(23)

  write(*,*) 'finished'

end program rcmaligngs 
