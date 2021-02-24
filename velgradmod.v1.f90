
! The stochastic model for velocity gradient coupled with its rate of change
module mprmtr
  use mconstant
  implicit none

  integer, parameter :: three = 3, npnt = 100

  ! Gaussian random noise forcing term.
  real(sp), parameter :: af = (3 + 15**.5) / (10**.5 + 6**.5) / 3.
  real(sp), parameter :: bf = -(10**.5 + 6**.5)/4.
  real(sp), parameter :: cf = 1. / (10**.5 + 6**.5)

  ! controls
  real(sp), parameter :: dt = 2.e-4, dtwrite = 0.1 
  ! Integral time scale is 1.

  ! parameters for models
  real(sp), parameter :: gmma = 0.1

  integer, dimension(three, three), parameter :: delij = (/(/1,0,0/),(/0,1,0/),(/0,0,1/)/)

end module mprmtr

program velgradmod
  use mconstant
  use mprmtr
  implicit none

  real(sp), dimension(three, three, npnt) :: aijall
  real(sp), dimension(three, three) :: aij, dwij, driftaij, driftaij1, phiaij
  real(sp), dimension(three, three) :: tmpvalmn1, tmpvalmn2, tmpvalmn3, tmpvalmn4
  real(sp), dimension(three, three, three, three) :: bijmn, tmpval
  real(sp), dimension(three, three, three, three) :: v, upaijmn, umaijmn
  real(sp), dimension(three, three, three, three) :: rpaijmn, rmaijmn

  integer :: idum, ii, jj, kk, ll, mm, nn, pp, qq, rr, nini 

  real(sp) :: twrite, tmax, time, gasdev, ran1, tmp1
  character(80) :: str


  write(*,*)
  write(*,'('' >>> Velocity Gradient Model <<< '')')
  write(*,*) 
  ll = iarg()
  if (iarg() .ne. 3) then
          write(*,*)
          write(*,*) ' >>>>>> Wrong number of arguments <<<<<<'
          write(*,*) 
          write(*,*) ' Usage: ./velgradmod.x nini tmax idum'
          write(*,*) '        nini: 0 initialized with Gaussian; otherwise the file number'
          write(*,*) '        tmax: time to run '
          write(*,*) '        idum: seed for random number'
          write(*,*)
          write(*,*) ' Stopped'
          stop 
  end if

  call getarg(1,str)
  read(str, '(I20)') nini

  call getarg(2,str)
  read(str, '(F10.6)') tmax

  call getarg(3, str)
  read(str, '(I20)') idum


  write(*,*) 'tmax = ', tmax
  write(*,*) 'dtwrite = ', dtwrite

  if ( nini .eq. 0 ) then

    idum = - idum
    write(*,*) 'idum = ', idum
    do ll = 1, npnt

      do jj = 1, three
      do ii = 1, three
        aij(ii,jj)=gasdev(idum)
      end do
      end do
      
      aij = af * delij * sum(aij * delij) + bf * aij + cf * transpose(aij)
      aall(:,:,ll) = aij

    end do

  else

    write(str, '(I20)') nini
    str = adjustl(str)
    str = 'velgrad-aij-'//str(1:len_trim(str))//'.data' 
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) aall
    close(10)

  end if

  time = 0.
  twrite = 0.
  ll = nini
  do while ( time .le. tmax + dt)

    if ( abs(time - twrite) .le. .5 * dt ) then

      write(str, '(I20)') ll
      str = adjustl(str)
      str = 'velgrad-aij-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) aall
      close(10)

      ll = ll + 1
      twrite = twrite + dtwrite

      write(*,*) 'time', time, 'aii', sum(aall(1,1,:) + aall(2,2,:)+aall(3,3,:))/npnt

      ! Remove the trace... Be careful about this... it's better to do it after the output
      do rr = 1, npnt
        aij = aall(:,:,rr)
        tmp1 = ( aij(1,1) + aij(2,2) + aij(3,3) )/3
        aall(:,:,rr) = aij - tmp1 * delij
      end do

    end if

    do rr = 1, npnt

      ! White noise force modeling part of the effects of the pressure Hessian
      ! TODO: this enforces different rms for off-diagonal and on-diagonal terms, 
      !       do we need this? given now dwij is a model for pressure Hessian.
      do kk = 1, three
      do jj = 1, three
        aij(jj,kk)=gasdev(idum) * sqrt(2.*dt)
      end do
      end do
      dwij = af * delij * sum(aij * delij) + bf * aij + cf * transpose(aij)

      ! Random variable modelling the integral of white noise
      do jj = 1, three * three
      do ii = 1, three * three

        ! The mapping between 1D and 2D
        mm = ( jj - 1 ) / 3 + 1    ! jj -> (mm,nn)
        nn = jj - 3 * (mm - 1)

        pp = ( ii - 1 ) / 3 + 1    ! ii -> (pp, qq)
        qq = ii - 3 * (pp - 1)

        if ( jj .eq. ii ) then 

            v(pp,qq,mm,nn) = - dt

        else if ( jj .lt. ii ) then

            if ( ran1(idum) .le. 0.5 ) then
                v(pp,qq,mm,nn) = - dt
            else
                v(pp,qq,mm,nn) =   dt
            end if

        else 

            v(pp,qq,mm,nn) = - v(mm,nn,pp,qq) 

        end if

      end do
      end do

      aij = aall(:,:,rr)

      ! order 2 weak convergent predictor-corrector method

      call drift(aij, driftaij)
      call diffusion(aij, bijmn)

      ! aij1 
      aij1 = aij + driftaij * dt

      do nn = 1, three
      do mm = 1, three

        upaijmn(:,:,mm,nn) = aij + bijmn(:,:,mm,nn) * sqrt(dt)
        umaijmn(:,:,mm,nn) = aij - bijmn(:,:,mm,nn) * sqrt(dt)
 
        rpaijmn(:,:,mm,nn) = aij + driftaij * dt + bijmn(:,:,mm,nn) * sqrt(dt)
        rmaijmn(:,:,mm,nn) = aij + driftaij * dt - bijmn(:,:,mm,nn) * sqrt(dt)

        aij1 = aij1 + bijmn(:,:,mm,nn) * dwij(mm,nn) !Diffusion term 

      end do
      end do

      ! phi
      phiaij = 0
      do nn = 1, three
      do mm = 1, three

        tmpvalmn1 = 0
        tmpvalmn2 = 0
        do qq = 1, three
        do pp = 1, three

          if ( (pp .eq. mm) .and. (qq .eq. nn) ) cycle

          call diffusion(upaijmn(:,:,pp,qq), tmpval)
          tmpvalmn3 = tmpval(:,:,mm,nn)

          call diffusion(umaijmn(:,:,pp,qq), tmpval)
          tmpvalmn4 = tmpval(:,:,mm,nn)

          tmpvalmn1 = tmpvalmn1 + tmpvalmn3 + tmpvalmn4 - 2 * bijmn(:,:,mm,nn)

          tmpvalmn3 = tmpvalmn3 - tmpvalmn4
          tmpvalmn3 = tmpvalmn3 * ( dwij(mm,nn) * dwij(pp,qq) - v(pp,qq,mm,nn) )
          tmpvalmn2 = tmpvalmn2 + tmpvalmn3

        end do
        end do

        call diffusion(rpaijmn(:,:,mm,nn), tmpval)
        tmpvalmn3 = tmpval(:,:,mm,nn)

        call diffusion(rmaijmn(:,:,mm,nn), tmpval)
        tmpvalmn4 = tmpval(:,:,mm,nn)

        phiaij = phiaij + ( tmpvalmn3 + tmpvalmn4 + 2 * bijmn(:,:,mm,nn) + tmpvalmn1 / sqrt(dt) ) &
                          * dwij(mm,nn)  &
                        + ( tmpvalmn3 - tmpvalmn4 ) * ( dwij(mm,nn) ** 2 - dt ) / sqrt(dt)        &
                        + tmpvalmn2 / sqrt(dt)
      end do
      end do
      phiaij = phiaij / 4

      call drift(aij1, driftaij1)

      aij1 = aij + (driftaij1 + driftaij) * dt / 2 + phiaij

      call drift(aij1, driftaij1)

      aij = aij + (driftaij1 + driftaij) * dt / 2 + phiaij


      aall(:,:,rr) = aij
 
    end do
 
    time = time + dt
  end do

  write(*,*) 'finished'
  
end program velgradmod      

subroutine drift(aij, driftaij)
  use mconstant
  use mprmtr
  implicit none

  real(sp), dimension(three, three), intent(in)  :: aij
  real(sp), dimension(three, three), intent(out) :: driftaij

  real(sp), dimension(three, three) :: a2ij, cinv, dinv
  real(sp) :: tra2, trcinv

  ! drift term for the aij equation

  ! Matrix exponential 
  dinv = -gmma * aij
  tra2 = sqrt( sum ( dinv * dinv ) )
  call taylorch(dinv, tra2, three)
  cinv = matmul(transpose(dinv), dinv)
  trcinv = cinv(1,1) + cinv(2,2) + cinv(3,3)

  a2ij = matmul(aij, aij)
  tra2 = a2ij(1,1) + a2ij(2,2) + a2ij(3,3)
  driftaij =  -a2ij + tra2 / trcinv * cinv - (trcinv / 3) * aij 


end subroutine drift

subroutine diffusion (aij, diffaij)
  use mconstant
  use mprmtr
  implicit none

  real(sp), dimension(three, three), intent(in)  :: aij
  real(sp), dimension(three, three, three, three), intent(out) :: diffaij
  integer :: ii, jj, pp, qq

  real(sp), dimension(three, three) :: cinv, dinv
  real(sp) :: trcinv

  ! Matrix exponential 
  dinv = - gmma * aij
  trcinv = sqrt( sum ( dinv * dinv ) )
  call taylorch(dinv, trcinv, three)
  cinv = matmul(transpose(dinv), dinv)
  trcinv = cinv(1,1) + cinv(2,2) + cinv(3,3)

  do qq = 1, three
  do pp = 1, three
  do jj = 1, three
  do ii = 1, three
    diffaij(ii,jj,pp,qq) = .5 * dinv(qq,ii) * delij(jj,pp) + .5 * dinv(qq,jj) * delij(ii,pp) &
                          - cinv(ii,jj) * dinv(qq,pp) / trcinv
  end do
  end do
  end do
  end do

end subroutine diffusion
