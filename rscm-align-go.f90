program rscmaligngo
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,three,npnt) :: aijall
  real(sp), dimension(three,npnt) :: giall
  real(sp), dimension(three,three) :: aij
  real(dp), dimension(three) :: wi, gi

  integer, parameter :: npdf = 40 
  real(sp), parameter :: bnd = 1., binw = bnd/npdf
  real(dp), dimension(npdf) :: pgo, pge3, pwe3

  integer :: ii, jj, kk, ll, nreal, ifirst
  real(dp) :: tmp, go 
  character(80) :: str, straij, strgi, strreal, strout

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./rscm-align-go.x aijfilelist gifilelist nreal ifirst'
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

  strout = strgi(1:len_trim(strgi))

  open(22, file = straij(1:len_trim(straij))//'.list' )
  open(23, file = strgi (1:len_trim(strgi ))//'.list' )

  jj = 0
  pgo = 0._dp
  pge3 = 0._dp
  pwe3 = 0._dp
  do while ( (.not. eof(22) ) .and. (.not. eof(23) ) )

    read(22,*) straij
    write(*,*) straij
    straij = adjustl(straij)

    read(23,*) strgi
    write(*,*) strgi
    strgi = adjustl(strgi)
 
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

        aij = aijall(:,:,ii)
        wi(1) = aij(2,3) - aij(3,2)
        wi(2) = aij(3,1) - aij(1,3)
        wi(3) = aij(1,2) - aij(2,1)
        wi = wi / sqrt( sum(wi * wi) )
     

        go = abs(sum(gi * wi))
     
        ll = floor( go / binw ) + 1
        if ( ll .ge. 1 .and. ll .le. npnt) pgo(ll) = pgo(ll) + 1

        ll = floor( abs(gi(3)) / binw ) + 1
        if (ll .ge. 1 .and. ll .le. npnt) pge3(ll) = pge3(ll) + 1

        ll = floor( abs(wi(3)) / binw ) + 1
        if (ll .ge. 1 .and. ll .le. npnt) pwe3(ll) = pwe3(ll) + 1

      end do

    end do

    jj = jj + 1

  end do
  close(22)
  close(23)

  tmp = 1. / nreal / npnt / jj

  pgo = pgo * tmp
  pge3  = pge3  * tmp
  pwe3 = pwe3 * tmp

  write(*,*) 'Check normalization pgo: ', sum(pgo)
  write(*,*) 'Check normalization pge3:  ', sum(pge3)
  write(*,*) 'Check normalization pwe3: ', sum(pwe3)

  pgo = pgo / binw
  pge3  = pge3  / binw
  pwe3 = pwe3 / binw
 
  str = 'align-go-'//strout(1:len_trim(strout))//'.dat'
  open(20, file = str(1:len_trim(str)) )

  write(20,'(''# x pdfgo pdfge3 pdfwe3'')') 
  do ii = 1, npdf
    tmp = (ii-.5) * binw
    write(20,'(15E15.4)') tmp, pgo(ii), pge3(ii), pwe3(ii)
  end do
 
  close(20)

  write(*,*) 'finished'

end program rscmaligngo 
