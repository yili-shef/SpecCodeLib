program genrand
    implicit none

    integer, parameter :: nx = 256, np = 64*64*64
    real, parameter :: pi = 3.1415926
    real, dimension(np) :: xp, yp, zp

    real :: ran1

    integer :: i, idum

    idum = -356

    DO i=1,np
        xp(i)=2.*pi*ran1(idum)
        yp(i)=2.*pi*ran1(idum)
        zp(i)=2.*pi*ran1(idum)
    END DO

    open(16, file = 'random13.dat', form = 'unformatted')
        write(16) xp, yp, zp
    close(16)

end program genrand      
