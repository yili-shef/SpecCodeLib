
! The stochastic model for velocity gradient coupled with its rate of change
module mprmtr
  use mconstant
  implicit none

  integer, parameter :: three = 3, npnt = 100

  ! Gaussian random noise forcing term.
  real(sp), parameter :: af = (3 + 15**.5_sp) / (10**.5_sp + 6**.5_sp) / 3
  real(sp), parameter :: bf = -(10**.5_sp + 6**.5_sp)/4
  real(sp), parameter :: cf = 1 / (10**.5_sp + 6**.5_sp)

  ! controls
  real(sp), parameter :: dt = 1.e-5_sp, dtwrite = 0.01_sp
  ! Integral time scale is 1.

  ! parameters for models
  real(sp), parameter :: gmma = 0.1_sp
  !real(sp), parameter :: gmma = 1.

  integer, dimension(three, three), parameter :: delij = (/(/1,0,0/),(/0,1,0/),(/0,0,1/)/)

end module mprmtr

program velgradmod
  use mconstant
  use mprmtr
  implicit none

  real(sp), dimension(three, three, npnt) :: aall, gall
  real(sp), dimension(three, three) :: aij, gij, aij1, gij1
  real(sp), dimension(three, three) :: dwij, driftgij, driftgij1, phigij
  !real(sp), dimension(three, three) :: tmpvalmn1, tmpvalmn2, tmpvalmn3, tmpvalmn4
  real(sp), dimension(three, three, three, three) :: bijmn, tmpval1, tmpval2
  !real(sp), dimension(three, three, three, three) :: v
  !real(sp), dimension(three, three, three, three) :: rpgijmn, rmgijmn, upgijmn, umgijmn
  real(sp), dimension(three, three) :: uaijmn, raijmn

  integer :: idum, ii, jj, kk, ll, mm, nn, pp, qq, rr, nini 

  real(sp) :: twrite, tmax, time, gasdev, ran1, tmp1, tmp2
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


      gij = - gmma * aij
      tmp2 = sqrt( sum ( gij * gij ) )
      call taylorch(gij, tmp2, three)   ! exp(-gmma * aij)
  
      dwij = matmul(transpose(gij), gij)  ! dwij = C^{-1} = exp(-gamma * aij^T) exp(-gamma * aij)
      tmp2 = dwij(1,1) + dwij(2,2) + dwij(3,3)   ! Tr(C^{-1})

      gij = matmul(aij, aij)
      tmp1 = gij(1,1) + gij(2,2) + gij(3,3)     ! Tr(A^2)

      ! Linearized pressure model
      !gij = - ( gij - tmp1 * delij / 3 ) &
      !      - gmma * tmp1 * ( aij + transpose(aij) ) / 3  &
      !      - tmp2 * aij/3

      ! Full version of pressure
      gij =  - gij + tmp1 / tmp2 * dwij - (tmp2 / 3) * aij 

      gall(:,:,ll) = gij

    end do

  else

    write(str, '(I20)') nini
    str = adjustl(str)
    str = 'velgrad-aij-'//str(1:len_trim(str))//'.data' 
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) aall
    close(10)

    write(str, '(I20)') nini
    str = adjustl(str)
    str = 'velgrad-gij-'//str(1:len_trim(str))//'.data' 
    open(10, file = str(1:len_trim(str)), form = 'binary')
      read(10) gall
    close(10)

  end if

  time = 0.
  twrite = 0.
  ll = nini

  open(15, file = 'aijseries.dat')
  do while ( time .le. tmax + dt)

    if ( abs(time - twrite) .le. .5 * dt ) then

      write(str, '(I20)') ll
      str = adjustl(str)
      str = 'velgrad-aij-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) aall
      close(10)

      write(str, '(I20)') ll
      str = adjustl(str)
      str = 'velgrad-gij-'//str(1:len_trim(str)) 
      str = str(1:len_trim(str))//'.data'
 
      open(10, file = str(1:len_trim(str)), form = 'binary')
        write(10) gall
      close(10)

      ll = ll + 1
      twrite = twrite + dtwrite

      write(*,*) 'time', time, 'aii', sum(aall(1,1,:) + aall(2,2,:)+aall(3,3,:))/npnt
      write(*,*) 'time', time, 'gii', sum(gall(1,1,:) + gall(2,2,:)+gall(3,3,:))/npnt
      write(*,*) 'time', time, 'max aij', maxval(aall(1,1,:) + aall(2,2,:)+aall(3,3,:))
      write(*,*) 'time', time, 'max gij', maxval(gall(1,1,:) + gall(2,2,:)+gall(3,3,:))
      write(*,*)

      write(15,'(10E12.3)') time, aall(1,1,1), aall(2,2,1), aall(3,3,1), aall(1,2,1)

      ! Remove the trace... Be careful about this... it should be done after the output
      !do rr = 1, npnt
      !  aij = aall(:,:,rr)
      !  tmp1 = ( aij(1,1) + aij(2,2) + aij(3,3) )/3
      !  aall(:,:,rr) = aij - tmp1 * delij
      !
      !  gij = gall(:,:,rr)
      !  tmp1 = ( gij(1,1) + gij(2,2) + gij(3,3) )/3
      !  gall(:,:,rr) = gij - tmp1 * delij
      !end do

    end if

    do rr = 1, npnt

      ! White noise
      do kk = 1, three
      do jj = 1, three
        dwij(jj,kk)=gasdev(idum) * sqrt(2._sp*dt)
      end do
      end do
      ! This part is removed. The projection has been absorbed into the definition of bijmn
      ! dwij = af * delij * sum(aij * delij) + bf * aij + cf * transpose(aij)
      ! dwij = ( dwij + transpose(dwij) ) / 2


      ! Random variable modelling the integral of white noise
      ! --- For our model this part does not appear ---
      !do jj = 1, three * three
      !do ii = 1, three * three

      !  ! The mapping between 1D and 2D
      !  mm = ( jj - 1 ) / 3 + 1    ! jj -> (mm,nn)
      !  nn = jj - 3 * (mm - 1)

      !  pp = ( ii - 1 ) / 3 + 1    ! ii -> (pp, qq)
      !  qq = ii - 3 * (pp - 1)
      !
      !  if ( jj .eq. ii ) then 
      !
      !      v(pp,qq,mm,nn) = - dt

      !  else if ( jj .lt. ii ) then
      !
      !      if ( ran1(idum) .le. 0.5 ) then
      !          v(pp,qq,mm,nn) = - dt
      !      else
      !          v(pp,qq,mm,nn) =   dt
      !      end if
      !
      !  else 

      !      v(pp,qq,mm,nn) = - v(mm,nn,pp,qq) 
      !
      !  end if

      !end do
      !end do
      !------------------------------------------------------------

      aij = aall(:,:,rr)
      gij = gall(:,:,rr)

      ! order 2 weak convergent predictor-corrector method

      call drift(aij, gij, driftgij)
      call diffusion(aij, bijmn)

      ! aij1 and gij1
      aij1 = aij + gij * dt        ! Drift term for aij. The diffusion term for aij is zero
      gij1 = gij + driftgij * dt

      do nn = 1, three
      do mm = 1, three

        ! U and R for aij. Do not depend on mm and nn. Thus taken out of the loop.

        !upaijmn(:,:,mm,nn) = aij
        !umaijmn(:,:,mm,nn) = aij
        !rpaijmn(:,:,mm,nn) = aij + gij * dt 
        !rmaijmn(:,:,mm,nn) = aij + gij * dt 


        ! U and R are needed to calculate bijmn.
        ! But bijmn does not depend on gij. Therefore upgijmn and umgijmn 
        ! are not relevant, can be omitted.

        !upgijmn(:,:,mm,nn) = gij + bijmn(:,:,mm,nn) * sqrt(dt)
        !umgijmn(:,:,mm,nn) = gij - bijmn(:,:,mm,nn) * sqrt(dt)
        !rpgijmn(:,:,mm,nn) = gij + driftgij * dt + bijmn(:,:,mm,nn) * sqrt(dt)
        !rmgijmn(:,:,mm,nn) = gij + driftgij * dt - bijmn(:,:,mm,nn) * sqrt(dt)

        gij1 = gij1 + bijmn(:,:,mm,nn) * dwij(mm,nn) !Diffusion term for gij

      end do
      end do

      ! Calculate phi
      ! ---- For our problem the general formulas can be simplied to the following ----

      uaijmn = aij
      raijmn = aij + gij * dt

      call diffusion(uaijmn, tmpval1)
      tmpval1 = (three * three - 1) * 2 * ( tmpval1 - bijmn ) / sqrt(dt)

      call diffusion(raijmn, tmpval2)

      phigij = 0
      do nn = 1, three
      do mm = 1, three
        phigij = phigij + ( 2 * tmpval2(:,:,mm,nn) + 2 * bijmn(:,:,mm,nn) + tmpval1(:,:,mm,nn) ) &
                          * dwij(mm,nn) 
      end do
      end do
      phigij = phigij / 4
      ! ------------------ For the full formula see end of the program --------------------
      ! -----------------------------------------------------------------------------------


      call drift(aij1, gij1, driftgij1)

      aij1 = aij + (gij1 + gij ) * dt / 2
      gij1 = gij + (driftgij1 + driftgij) * dt / 2 + phigij

      call drift(aij1, gij1, driftgij1)

      aij = aij + (gij1 + gij ) * dt / 2
      gij = gij + (driftgij1 + driftgij) * dt / 2 + phigij


      aall(:,:,rr) = aij
      gall(:,:,rr) = gij
 
    end do
 
    time = time + dt
  end do
  close(15)

  write(*,*) 'finished'
  
end program velgradmod      

      ! The following general formula can be simplified in current problem

      ! phi
      !phigij = 0
      !do nn = 1, three
      !do mm = 1, three


        !tmpvalmn1 = 0
        !tmpvalmn2 = 0
        !do qq = 1, three
        !do pp = 1, three
        !
        !  if ( (pp .eq. mm) .and. (qq .eq. nn) ) cycle
        !
        !  !call diffusion(upaijmn(:,:,pp,qq), tmpval)
        !  !tmpvalmn3 = tmpval(:,:,mm,nn)
        !
        !  !call diffusion(umaijmn(:,:,pp,qq), tmpval)
        !  !tmpvalmn4 = tmpval(:,:,mm,nn)

        !  tmpvalmn1 = tmpvalmn1 + tmpvalmn3 + tmpvalmn4 - 2 * bijmn(:,:,mm,nn)
        !
        !  tmpvalmn3 = tmpvalmn3 - tmpvalmn4
        !  tmpvalmn3 = tmpvalmn3 * ( dwij(mm,nn) * dwij(pp,qq) - v(pp,qq,mm,nn) )
        !  tmpvalmn2 = tmpvalmn2 + tmpvalmn3
        !
        !end do
        !end do

        !call diffusion(rpaijmn, tmpval)
        !tmpvalmn3 = tmpval(:,:,mm,nn)

        !call diffusion(rmaijmn, tmpval)
        !tmpvalmn4 = tmpval(:,:,mm,nn)

        !phigij = phigij + ( tmpvalmn3 + tmpvalmn4 + 2 * bijmn(:,:,mm,nn) + tmpvalmn1 / sqrt(dt) ) &
        !                  * dwij(mm,nn)  &
        !                + ( tmpvalmn3 - tmpvalmn4 ) * ( dwij(mm,nn) ** 2 - dt ) / sqrt(dt)        &
        !                + tmpvalmn2 / sqrt(dt)
      !end do
      !end do
      !phigij = phigij / 4


subroutine drift(aij, gij, driftgij)
  use mconstant
  use mprmtr
  implicit none

  real(sp), dimension(three, three), intent(in)  :: aij, gij
  real(sp), dimension(three, three), intent(out) :: driftgij

  real(sp), dimension(three, three) :: a2ij, sij, dinv, cinv
  real(sp) :: tra2, trcinv, trcinvg, trcinva, trag, aii, gii

  ! drift term for the gij equation

  aii = aij(1,1) + aij(2,2) + aij(3,3)
  gii = gij(1,1) + gij(2,2) + gij(3,3)

  dinv = - gmma * aij
  tra2 = max( sqrt( sum ( dinv * dinv ) ) , 1.)
  call taylorch(dinv, tra2, three)

  cinv = matmul(transpose(dinv), dinv)
  trcinv = cinv(1,1) + cinv(2,2) + cinv(3,3)

  a2ij = matmul(aij, aij)
  tra2 = a2ij(1,1) + a2ij(2,2) + a2ij(3,3)

  dinv = matmul(cinv, gij)
  trcinvg = dinv(1,1) + dinv(2,2) + dinv(3,3)
  !trcinvg = sum(transpose(cinv) * gij)

  dinv = matmul(cinv, aij)
  trcinva =  dinv(1,1) + dinv(2,2) + dinv(3,3)
  !trcinva = sum(transpose(cinv) * aij)
  
  trag = sum( gij * transpose(aij) ) 

  ! sij = .5 * ( aij + transpose(aij) )

  ! What is this?
  ! In this version pressure is totally white noise. Viscous diffusion is modelled with CM model
  ! This one is not working
  !driftgij = tra2 * aij / 3 &
             !- gmma * tra2 * (gij + transpose(gij)) / 3            & 
             !- 4 * gmma * sij * sum( gij * transpose(aij) ) / 3                   &
  !           - trcinv * gij / 3 + 2 * gmma * trcinvg * aij / 3  &  ! From viscous diffusion term
             !+ 2 * gmma * tra2 / 3 * ( matmul(aij, sij) + matmul(sij, aij)   &
             !- 2 * delij * sum(sij * transpose(aij)) / 3 )                        &
  !           + 2 * trcinv / 3 * ( a2ij - tra2 * delij / 3 ) 


  ! linearized pressure model + full viscous diffusion 
  !driftgij = tra2 * aij / 3 &
  !           - gmma * tra2 * (gij + transpose(gij)) / 3            & 
  !           - 4 * gmma * sij * sum( gij * transpose(aij) ) / 3                   &
  !           - trcinv * gij / 3 + 2 * gmma * trcinvg * aij / 3  &  ! From viscous diffusion term
  !           + 2 * gmma * tra2 / 3 * ( matmul(aij, sij) + matmul(sij, aij)   &
  !           - 2 * delij * sum(sij * transpose(aij)) / 3 )                        &
  !           + 2 * trcinv / 3 * ( a2ij - tra2 * delij / 3 ) 


  ! full pressure model + full viscous diffusion 
  driftgij = tra2 * (aij - aii * delij/3) / 3 &
             - tra2 / trcinv * ( matmul(aij,cinv) + matmul(cinv,aij) - 2 * trcinva / 3 * delij ) &
             + (aij - aii * delij / 3 ) * 2 * tra2 / 3 &
             + 2 * trag * ( cinv / trcinv -  delij / 3 )                        &
             + gmma * tra2 / trcinv *  2 * trcinvg * cinv / trcinv   &
             - gmma * tra2 / trcinv * ( matmul(transpose(gij),cinv) + matmul(cinv,gij) ) &
             - trcinv * (gij - gii * delij/3) / 3  & ! From viscous term 
             + 2 * gmma * trcinvg * (aij - aii * delij/3) / 3  &  ! From viscous diffusion term
             + 2 * trcinv / 3 * ( a2ij - tra2 * delij / 3 )  ! viscous damping term

end subroutine drift

subroutine diffusion(aij, bijmn)
  use mconstant
  use mprmtr
  implicit none

  real(sp), dimension(three,three,three,three), intent(out) :: bijmn
  real(sp), dimension(three,three), intent(in) :: aij

  integer :: ii, jj, mm, nn
  real(sp) :: aii

  !aii = aij(1,1) + aij(2,2) + aij(3,3)
  !
  !do nn = 1, three
  !do mm = 1, three
  !do jj = 1, three
  !do ii = 1, three

  !  This version erranously took the forcing term in CM model as the pure white noise.

  !  bijmn(ii,jj,mm,nn) = - delij(jj,nn) * aij(ii,mm) - delij(ii,mm) * aij(nn,jj) &
  !                       + 2 * delij(ii,jj) * aij(nn,mm) / 3


    ! The following should be correct now.
  !  bijmn(ii,jj,mm,nn) = sqrt(6.) * ( - aij(ii,jj) * delij(mm,nn) / 3 &
  !                                    + delij(ii,jj) * delij(mm,nn) * aii / 9  &
  !                      + ( delij(jj,nn) * aij(ii,mm) + delij(jj,mm) * aij(ii,nn) ) / 4 &
  !                      + ( delij(ii,mm) * aij(nn,jj) + delij(ii,nn) * aij(mm,jj) ) / 4 &
  !                      - ( aij(nn,mm) + aij(mm,nn) ) * delij(ii,jj) / 6 )
  !end do
  !end do
  !end do
  !end do
  !
  ! TESTING TODO
   bijmn = 0._sp
  
  
end subroutine diffusion 
