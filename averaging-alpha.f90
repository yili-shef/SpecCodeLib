program averaging
    implicit none

    integer, parameter :: numfile = 9, numstep = 100

    ! for the output of pt-vort-dp.x
    real(8), dimension(numstep) :: time, meanss, meanvs, meandiffs
    real(8), dimension(numstep) :: meanphloc, meanphnonloc, meanph, meandalphadt
    real :: time1, meanss1, meanvs1, meandiffs1, meanphloc1, meanphnonloc1, meanph1, meandalphadt1
    ! ------------------------

    character(80) :: flnm, flsuffix, flprefix

    integer :: ii, jj, nprtcle
    integer(8) :: allprt

    ! for the output of pt-vort-dp.x
    flprefix = "alpha-"
    flprefix = adjustl(flprefix)
    ! ------------------------

    flsuffix = "higho-bw.dat"
    flsuffix = adjustl(flsuffix)

    time = 0.

    meanss = 0.
    meanvs = 0.
    meandiffs = 0.
    meanphloc = 0.
    meanphnonloc = 0.
    meanph = 0.
    meandalphadt = 0.

    allprt = 0
    open(20, file = 'nprtcle.dat')

    do ii = 1, numfile
        read(20,*) nprtcle
        allprt = allprt + nprtcle

        write(flnm, '(I20)') ii
        flnm = adjustl(flnm)

        flnm = flprefix(1:len_trim(flprefix))//flnm(1:len_trim(flnm)) &
               //flsuffix(1:len_trim(flsuffix))

        open(15, file = flnm)
            do jj = 1, numstep

                read(15,*) time1, meanss1, meanvs1, meandiffs1, meanphloc1, meanphnonloc1, &
                            meanph1, meandalphadt1

                time(jj) = time1
                meanss(jj) = meanss(jj) + meanss1 * real(nprtcle, 8)
                meanvs(jj) = meanvs(jj) + meanvs1  * real(nprtcle, 8)
                meandiffs(jj) = meandiffs(jj) + meandiffs1 * real(nprtcle, 8)
                meanphloc(jj) = meanphloc(jj) + meanphloc1 * real(nprtcle, 8)
                meanphnonloc(jj) = meanphnonloc(jj) + meanphnonloc1 * real(nprtcle, 8)
                meanph(jj) = meanph(jj) + meanph1 * real(nprtcle,8)
                meandalphadt(jj) = meandalphadt(jj) + meandalphadt1 * real(nprtcle,8)

            end do
        close(15) 

    end do
    close(20)

    meanss = meanss / real(allprt, 8)
    meanvs = meanvs / real(allprt, 8)
    meandiffs = meandiffs / real(allprt, 8)
    meanphloc = meanphloc / real(allprt, 8)
    meanphnonloc = meanphnonloc / real(allprt, 8)
    meanph = meanph / real(allprt, 8)
    meandalphadt = meandalphadt / real(allprt, 8)

    flnm = "mean-"//flprefix(1:len_trim(flprefix))//flsuffix(1:len_trim(flsuffix))
    open(15, file = flnm)
        do jj = 1, numstep

            write(15,'(20E12.4)') time(jj), meanss(jj), meanvs(jj), meandiffs(jj), meanphloc(jj), &
                                  meanphnonloc(jj), meanph(jj), meandalphadt(jj)

        end do
    close(15)
    
end program averaging      
