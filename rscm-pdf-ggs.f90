program rscmpdfggs
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,npnt) :: giall
  real(sp), dimension(three) :: gi
  real(sp), dimension(three,three,npnt) :: aijall
  real(sp), dimension(three,three) :: aij

  integer,  parameter :: npdf = 80
  real(sp), parameter :: bndggs = 15., bwggs = 2.*bndggs/npdf
  real(dp), dimension(npdf) :: pdfggs

  real(dp) :: meanggs, meanggs0, rmsggs, rmsggs0

  integer :: ii, jj, ll, nfile, nreal, ifirst
  real(sp) :: tmp, ggs
  character(80) :: str, strgilist, straijlist, strreal, strgi, straij

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
      write(*,*)
      write(*,*) ' >>>>>> Wrong number of arguments <<<<<< '
      write(*,*) 
      write(*,*) ' Usage: ./rscm-pdf-ggs-s(d)p.x gifilelist aijfilelist nreal ifirst'
      write(*,*) '        gifilelist: gifilelist.list is the list of gi data files'
      write(*,*) '        aijfilelist: aijfilelist.list is the list of aij data files'
      write(*,*) '        nreal: number of realizations'
      write(*,*) '        ifirst: the index of the first realization'
      write(*,*)
      write(*,*) ' Stopped'
      stop 
  end if

  call getarg(1,strgilist)
  strgilist = adjustl(strgilist)

  call getarg(2,straijlist)
  straijlist = adjustl(straijlist)

  call getarg(3,strreal)
  read(strreal, '(I20)') nreal

  call getarg(4,str)
  read(str, '(I20)') ifirst

  open(22, file =  strgilist(1:len_trim( strgilist))//'.list')
  open(30, file = straijlist(1:len_trim(straijlist))//'.list')

  meanggs = 0.0_dp
  rmsggs = 0.0_dp
  nfile = 0
  pdfggs = 0._dp
  do while ( .not. eof(22) )

    read (22,*) strgi
    write( *,*) strgi
    strgi = adjustl(strgi)
 
    read (30,*) straij
    write( *,*) straij
    straij = adjustl(straij)
  
    do ll = ifirst, nreal + ifirst -1
 
 
      write(*,*) 'realization ', ll
      write(strreal, '(I20)') ll
      strreal = adjustl(strreal)
  
      str  = strgi (1:len_trim( strgi))//strreal(1:len_trim(strreal))//'.data'
      open(10, file = str(1:len_trim( str)), form = 'binary')
        read(10) giall
      close(10)

      str = straij(1:len_trim(straij))//strreal(1:len_trim(strreal))//'.data'
      open(10, file = str(1:len_trim(str)), form = 'binary')
        read(10) aijall
      close(10)
 
      if ( nfile .eq. 0 ) then
          meanggs0 = 0.d0
          do ii = 1, npnt
            gi  =  giall(:,ii)
            aij = aijall(:,:,ii)
            aij = .5 * (aij + transpose(aij))
            ggs = - dot_product( gi , matmul( aij, gi ) )

            meanggs0 = meanggs0 + ggs
            rmsggs0 = rmsggs0 + ggs * ggs 
          end do
          meanggs0 = meanggs0 / npnt
          rmsggs0 = rmsggs0 / npnt
          rmsggs0 = sqrt( rmsggs0 - meanggs0 * meanggs0 )

          write(*,*) 'estimate of mean ggs: ', meanggs0
          write(*,*) 'estimate of rms ggs: ', rmsggs0
      end if
    
      do ii = 1, npnt
        gi  =  giall(:,ii)
        aij = aijall(:,:,ii)

        aij = ( aij + transpose(aij) ) * .5 
 
        ggs = - dot_product( gi, matmul( aij, gi ) ) 
 
        meanggs = meanggs + ggs
        rmsggs = rmsggs + ggs * ggs
        ggs = (ggs - meanggs0) / rmsggs0
 
        jj = floor( ( ggs + bndggs) / bwggs ) + 1
 
        if ( jj .ge. 1 .and. jj .le. npdf ) then 
            pdfggs(jj) = pdfggs(jj) + 1
        end if
      end do

      nfile = nfile + 1

    end do

  end do
 
  close(22)
  close(30)

  tmp = 1. / nfile / npnt

  meanggs = meanggs * tmp
  rmsggs = rmsggs * tmp
  rmsggs = sqrt(rmsggs - meanggs * meanggs)

  pdfggs = pdfggs * tmp
 
  write(*,*) 'check the normalization of pdfggs: ', sum(pdfggs)

  write(*,*) 'meanggs: ', meanggs
  write(*,*) 'rmsggs: ', rmsggs

  pdfggs = pdfggs / bwggs

  str = 'rcm-pdf-ggs-'//strgilist(1:len_trim(strgilist))//'.dat'
  open(20, file = str(1:len_trim(str)) )

    write(20,'(''# zone T = " pdf ", I='', I4, '', F = point'')') npdf
    do jj = 1, npdf
        tmp = ( ( -bndggs + (jj-.5) * bwggs ) * rmsggs0 + meanggs0 - meanggs ) / rmsggs
        write(20,'(15E15.4)') tmp, pdfggs(jj) * rmsggs / rmsggs0
    end do

  close(20)

  write(*,*) 'finished'

end program rscmpdfggs
