program getgroups
    implicit none

    integer, parameter :: nskip = 16, nx = 128, ny = nx, nz = nx
    integer, parameter :: radius1 = 10, radius2 = 11, ngroups = (nx / nskip)**3
    integer, parameter :: npntgrp = (2*radius2+1)**3 - (2*radius1 + 1)**3

    integer :: ii, jj, kk, ll, mm, iiloc, jjloc, kkloc, ipnt
    integer, dimension(npntgrp, ngroups) :: xx, yy, zz
    integer :: xxll, yyll, zzll
    character(80) :: strgrp, strpnt, fdata


    write(*,*) 'Get the centers'
    ll = 0
    do kk = nskip, nz, nskip
    do jj = nskip, ny, nskip
    do ii = nskip, nx, nskip
        ll = ll + 1
        xx(1, ll) = ii; yy(1, ll) = jj; zz(1, ll) = kk ! NOTE: the center of each group
    end do
    end do
    end do

    write(*,*) 'll = ', ll

    write(*,*) 'Get the groups'
    do ll = 1, ngroups

        ipnt = 1
          
        do kkloc = -1, 1 ! checking the images of the flow field
        do jjloc = -1, 1
        do iiloc = -1, 1
        
            xxll = xx(1,ll) + iiloc * nx
            yyll = yy(1,ll) + jjloc * ny
            zzll = zz(1,ll) + kkloc * nz
            
            do kk = 1, nz
            do jj = 1, ny
            do ii = 1, nx
 
                    mm = floor( sqrt( real(ii-xxll)**2 + (jj-yyll)**2 + (kk-zzll)**2 ) + 0.5d0 ) + 1
 
                    if ( mm .gt. radius1 .and. mm .le. radius2 ) then

                        ipnt = ipnt + 1

                        if (ipnt .gt. npntgrp ) then 
                            write(*,*) 'ipnt gt npntgrp', ipnt
                            stop
                        end if

                        xx(ipnt,ll) = ii; yy(ipnt, ll) = jj; zz(ipnt, ll) = kk
                    end if

            end do
            end do
            end do
 
        end do
        end do
        end do
 
    end do

 
    write(*,*) 'ipnt = ', ipnt
    write(*,*) 'ngroups = ', ngroups

    write(strgrp, '(I20)') ngroups
    strgrp = adjustl(strgrp)
    write(strpnt, '(I20)') ipnt
    strpnt = adjustl(strpnt)

    fdata = 'xyz10ring'//strgrp(1:len_trim(strgrp))//'grps'//strpnt(1:len_trim(strpnt))//'pnts.dat'

    open(30, file = fdata(1:len_trim(fdata)), form = 'unformatted')
        write(30) xx(1:ipnt,1:ngroups), yy(1:ipnt,1:ngroups), zz(1:ipnt,1:ngroups)
    close(30)
 
end program getgroups
