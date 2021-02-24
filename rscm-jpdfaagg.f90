program rscmjpdf
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,npnt) :: giall
  real(sp), dimension(three) :: gi
  real(sp), dimension(three,three,npnt) :: aijall
  real(sp), dimension(three,three) :: aij

  integer,  parameter :: npdfgg = 60, npdfaa = 40
  real(sp), parameter :: bndgg = 15., bwgg = bndgg/npdfgg, bndaa = 12., bwaa = bndaa/npdfaa
  real(dp), dimension(npdfaa,npdfgg) :: jpdf

  real(dp) :: meanaa, meangg, meanaa0, meangg0

  integer :: ii, jj, kk, ll, nfile, nreal, ifirst
  real(sp) :: tmp, aa, gg
  character(80) :: str, strgilist, straijlist, strreal, strgi, straij

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
      write(*,*)
      write(*,*) ' >>>>>> Wrong number of arguments <<<<<< '
      write(*,*) 
      write(*,*) ' Usage: ./rscm-jpdfaagg-s(d)p.x gifilelist aijfilelist nreal ifirst'
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

  meangg = 0.0_dp
  meanaa = 0.0_dp
  nfile = 0
  jpdf = 0._dp
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
          gg = 0.d0
          aa = 0.d0
          do ii = 1, npnt
            gi  =  giall(:,ii)
            aij = aijall(:,:,ii)
            gg = gg + sum(gi * gi) 
            aa = aa + sum(aij * aij)
          end do
          meangg0 = gg / npnt
          meanaa0 = aa / npnt

          write(*,*) 'mean gg estimate: ', meangg0
          write(*,*) 'mean aa estimate: ', meanaa0
      end if
    
      do ii = 1, npnt
        gi  =  giall(:,ii)
        aij = aijall(:,:,ii)
 
        gg = sum(gi * gi) 
        aa = sum(aij * aij) 
 
        meangg = meangg + gg
        meanaa = meanaa + aa
        gg = gg / meangg0
        aa = aa / meanaa0
 
        jj = floor( aa / bwaa ) + 1
        kk = floor( gg / bwgg ) + 1
 
        if ( jj .ge. 1 .and. jj .le. npdfaa .and. & 
             kk .ge. 1 .and. kk .le. npdfgg ) then
            jpdf(jj,kk) = jpdf(jj,kk) + 1
        end if
      end do

      nfile = nfile + 1

    end do

  end do
 
  close(22)
  close(30)

  tmp = 1. / nfile / npnt

  meangg = meangg * tmp
  meanaa = meanaa * tmp
  jpdf = jpdf * tmp
 
  write(*,*) 'check normalization jpdf: ', sum(jpdf)

  write(*,*) 'meangg: ', meangg
  write(*,*) 'meanaa: ', meanaa

  jpdf = jpdf / bwgg / bwaa

  str = 'jpdfaagg-'//strgilist(1:len_trim(strgilist))//'.dat'
  open(20, file = str(1:len_trim(str)) )

    write(20,'(''# meangg'', E13.4, ''meanaa'', E13.4)') meangg, meanaa
    write(20,'(''# zone T = " jpdf ", I='', I4,'', J='', I4, '', F = point'')') npdfaa, npdfgg
    do jj = 1, npdfaa
 
      tmp =  (jj-.5) * bwaa * meanaa0 / meanaa
      do kk = 1, npdfgg
        write(20,'(15E15.4)') tmp, (kk-.5)*bwgg * meangg0 / meangg, &
                              jpdf(jj,kk) * meanaa/meanaa0 * meangg/meangg0
      end do
      write(20,*)
 
    end do

  close(20)

  write(*,*) 'finished'

end program rscmjpdf
