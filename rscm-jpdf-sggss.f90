program rscmjpdfsggss
  use mconstant
  implicit none

  integer, parameter :: npnt = 100000, three = 3
  real(sp), parameter :: dt = 0.0001

  real(sp), dimension(three,npnt) :: giall
  real(sp), dimension(three) :: gi
  real(sp), dimension(three,three,npnt) :: aijall
  real(sp), dimension(three,three) :: aij

  integer,  parameter :: npdfss = 60, npdfsgg = 40
  real(sp), parameter :: bndss = 10., bwss = bndss/npdfss, bndsgg = 4., bwsgg = 2*bndsgg/npdfsgg
  real(dp), dimension(npdfsgg,npdfss) :: jpdf

  real(dp) :: meansgg, rmssgg, meanss, meanss0, rmsss, corr

  integer :: ii, jj, kk, ll, nfile, nreal, ifirst
  real(sp) :: tmp, sgg, ss, gg
  character(80) :: str, strgilist, straijlist, strreal, strgi, straij

  write(*,*)
  write(*,'('' >>> Postprocessing Chevillard-Meneveau model for velocity gradient <<< '')')
  write(*,*) 
  ll = iarg()
  if (ll .ne. 4) then
      write(*,*)
      write(*,*) ' >>>>>> Wrong number of arguments <<<<<< '
      write(*,*) 
      write(*,*) ' Usage: ./rscm-jpdfsggss-s(d)p.x gifilelist aijfilelist nreal ifirst'
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

  meanss = 0.0_dp
  meansgg = 0.0_dp
  rmssgg = 0.0_dp
  rmsss = 0.0_dp
  corr = 0.0_dp
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
          meanss0 = 0.d0
          do ii = 1, npnt
            gi  =  giall(:,ii)
            aij = aijall(:,:,ii)
            aij = .5 * (aij + transpose(aij))

            meanss0 = meanss0 + sum(aij * aij)
          end do
          meanss0 = meanss0 / npnt

          write(*,*) 'mean ss estimate: ', meanss0
      end if
    
      do ii = 1, npnt
        gi  =  giall(:,ii)
        aij = aijall(:,:,ii)
        aij = .5 * (aij + transpose(aij))
 
        ss = sum(aij * aij) 

        gg = sum(gi * gi)
        sgg = - dot_product( gi, matmul(aij, gi) ) /gg
 
        meanss = meanss + ss
        rmsss = rmsss + ss * ss

        meansgg = meansgg + sgg
        rmssgg = rmssgg + sgg * sgg

        corr = corr + ss * sgg

        ss = ss / meanss0 
        sgg = sgg / sqrt(meanss0)
 
        jj = floor( (sgg + bndsgg) / bwsgg ) + 1
        kk = floor( ss / bwss ) + 1
 
        if ( jj .ge. 1 .and. jj .le. npdfsgg .and. & 
             kk .ge. 1 .and. kk .le. npdfss ) then
            jpdf(jj,kk) = jpdf(jj,kk) + 1
        end if
      end do

      nfile = nfile + 1

    end do

  end do
 
  close(22)
  close(30)

  tmp = 1. / nfile / npnt

  meanss = meanss * tmp
  rmsss = rmsss * tmp
  rmsss = sqrt(rmsss - meanss * meanss)

  meansgg = meansgg * tmp
  rmssgg = rmssgg * tmp
  rmssgg = sqrt(rmssgg - meansgg * meansgg)

  corr = corr * tmp
  corr = (corr - meanss * meansgg) / ( rmsss * rmssgg )

  jpdf = jpdf * tmp
 
  write(*,*) 'check normalization jpdf: ', sum(jpdf)

  write(*,*) 'meanss: ', meanss
  write(*,*) 'rmsss: ', rmsss
  write(*,*) 'meansgg: ', meansgg
  write(*,*) 'rmssgg: ', rmssgg

  jpdf = jpdf / bwss / bwsgg

  str = 'jpdfsggss-'//strgilist(1:len_trim(strgilist))//'.dat'
  open(20, file = str(1:len_trim(str)) )

    write(20,'(''# meanss'', E13.4, ''rmsss'', E13.4)') meanss, rmsss
    write(20,'(''# meansgg'', E13.4, ''rmssgg'', E13.4)') meansgg, rmssgg
    write(20,'(''# correlation'', E13.4)') corr
    write(20,'(''# zone T = " jpdf ", I='', I4,'', J='', I4, '', F = point'')') npdfsgg, npdfss
    do jj = 1, npdfsgg
 
      tmp = ( (-bndsgg + (jj-.5) * bwsgg) * sqrt(meanss0) ) / sqrt(meanss)
      do kk = 1, npdfss
        write(20,'(15E15.4)') tmp, (kk-.5)*bwss * meanss0 / meanss, &
                              jpdf(jj,kk) * sqrt(meanss/meanss0) * meanss/meanss0
      end do
      write(20,*)
 
    end do

  close(20)

  write(*,*) 'finished'

end program rscmjpdfsggss
