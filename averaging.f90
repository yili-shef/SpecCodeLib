program averaging
    implicit none

    integer, parameter :: numfile = 13, numstep = 100

    ! for the output of pt-vort-dp.x
    real(8), dimension(numstep) :: time, ww, meanalpha, meancurv, mmeanalfcurv, meanalfww
    real(8) :: time1, ww1, meanalpha1, meancurv1, mmeanalfcurv1, meanalfww1
    ! ------------------------

    character(80) :: flnm, flsuffix, flprefix

    integer :: ii, jj, nprtcle
    integer(8) :: allprt

    ! for the output of pt-vort-dp.x
    flprefix = "ww-new-"
    flprefix = adjustl(flprefix)
    ! ------------------------

    flsuffix = "highs-fw.dat"
    flsuffix = adjustl(flsuffix)

    time = 0.


    ! for the output of pt-vort-dp.x
    ww = 0.
    meanalpha = 0.
    meancurv = 0.
    mmeanalfcurv = 0.
    meanalfww = 0.
    ! ---------------

    allprt = 0
    open(20, file = 'nprtcle-highs.dat')

    do ii = 1, numfile
        read(20,*) nprtcle
        allprt = allprt + nprtcle

        write(flnm, '(I20)') ii
        flnm = adjustl(flnm)

        flnm = flprefix(1:len_trim(flprefix))//flnm(1:len_trim(flnm)) &
               //flsuffix(1:len_trim(flsuffix))

        open(15, file = flnm)
            do jj = 1, numstep

            ! for the output of pt-vort-dp.x

                read(15,*) time1, ww1, meanalpha1, meancurv1, mmeanalfcurv1, meanalfww1

                time(jj) = time1
                ww(jj) = ww(jj) + ww1 * real(nprtcle, 8)
                meanalpha(jj) = meanalpha(jj) + meanalpha1  * real(nprtcle, 8)
                meancurv(jj) = meancurv(jj) + meancurv1 * real(nprtcle, 8)
                mmeanalfcurv(jj) = mmeanalfcurv(jj) + mmeanalfcurv1 * real(nprtcle, 8)
                meanalfww(jj) = meanalfww(jj) + meanalfww1 * real(nprtcle, 8)
            ! --------------------------------

            end do
        close(15) 

    end do
    close(20)

    ! for the output of pt-vort-dp.x
    ww = ww / real(allprt, 8)
    meanalpha = meanalpha / real(allprt, 8)
    meancurv = meancurv / real(allprt, 8)
    mmeanalfcurv = mmeanalfcurv / real(allprt, 8)
    meanalfww = meanalfww / real(allprt, 8)
    ! ----------------------------------

    flnm = "mean-"//flprefix(1:len_trim(flprefix))//flsuffix(1:len_trim(flsuffix))
    open(15, file = flnm)
        do jj = 1, numstep

        ! for the output of pt-vort-dp.x
            write(15,'(20E12.4)') time(jj), ww(jj), meanalpha(jj), meancurv(jj), mmeanalfcurv(jj), &
                                  meanalfww(jj)
        ! ----------------------------------

        end do
    close(15)
    
end program averaging      
