! updated on 02/23/03---Single precision!
!!!hange forcing scheme-- fix enstrophy input rate eta! 
!... add force into equation.
! Add hypoviscosity everywhere. a = rnu_i * k ^(-16), ki=2
! Add hyperviscosity everywhere. a = rnu_u * k ^(16), ku=740
! sigma = rnu_u ^ 1/16, nu = run_i ^ 1/16
! Add forcing at 3 shells :  k1=4, k2=5, k3=6
! energy,dissipation,flux and enstophy in k-space are divides by 4 !!
! to keep consistent with x-space values.
! reduce large-damping order to k^(-4).
!*************************************
      program spectral_mpi
      include 'include.h'

      complex,allocatable,dimension(:, :)::ux,uy,vx,vy,wz,  &
	                                   owz
      complex,allocatable,dimension(:, :)::wt,uxt,vxt
      real,allocatable,dimension(:, :)::kx,ky,k2,tmp,tmp1,ftmp,ftmp1
      real, allocatable,dimension(:)::wx,wy 
      integer,allocatable,dimension(:)::ipx,ipy         
      integer,allocatable,dimension(:)::iseed

      integer new,id, numb 
      integer nstep,itout, ieout,iseed_m0
      integer jj,jjj, i, j, k, irr,istep,idp,idp1,iii
      real dt, u0, pi, dt_h, time, ak0
      real etak, rnu_i,rnu_u,kf,sigma,alpha,nu,eta,scale
      real ff, Re_lmda, lmda,turb_u, beta
      real timer1,timer2,timer3,timer4,timer5,timer6
      real timer7,timer8,timer9,timer10,timer11,timer12
      real ek, e_t, tt, ttt, ent, ek2, ek4 ,ek1
      real umax,max1,max2,umax_t,xmax,ymax,dt_max 
      real flatness,vor(2),vort(2),sum1
      real (8) num 
      character(46) fnm2
      character(5) ncase
      common/casenum/ncase
!  setup MPP environment

      ncase='case1'
      call mpi_init(ierror)
      call mpi_comm_size(mpi_comm_world,nproc,ierror)
      call mpi_comm_rank(mpi_comm_world,id,ierror)
      call mpi_barrier(mpi_comm_world,ierror)
      nallgrp = mpi_comm_world

!   read parameters from node id = 0

      if (id.eq.0)  then              
        open(1,file='parameter.d2',status='old') 
           read(1,*) my
           write(*,*) '  my=', my
           read(1,*) nproc
           write(*,*)'  nproc=', nproc
           read(1,*) iseed_m0
           write(*,*) '   seed=', iseed_m0
           read(1,*) ak0
           write(*,*) '   ak0=', ak0
           read(1,*) nstep
           write(*,*) '   nstep=', nstep
           read(1,*) itout
           write(*,*) '   itout=', itout
           read(1,*) ieout
           write(*,*) '   ieout=', ieout
           read(1,*) dt
           write(*,*) '   dt=', dt
           read(1,*)  nu
           write(*,*) '   nu=', nu
           read(1,*)  sigma
           write(*,*) '   sigma=', sigma
           read(1,*) u0
           write(*,*) '   u0=', u0
           read(1,*) new
           write(*,*) '   new=', new
           read(1,*) time
           write(*,*) '   beginning time', time
           read(1,*) idp
           write(*,*) '   idp(for rerun)=', idp
           read(1,*) eta
           write(*,*) '   force=', eta 
           read(1,*) kf 
           write(*,*) '   kf=', kf
         
         close(1)
      endif

! Broadcast inputs across processors, since just read into id=0.
      call mpi_bcast(my,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(nproc,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(iseed_m0,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(ak0,1,MPI_REAL,0,nallgrp,ierror)
      call mpi_bcast(nstep,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(itout,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(ieout,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(dt,1,MPI_REAL,0,nallgrp,ierror)
      call mpi_bcast(nu,1,MPI_REAL,0,nallgrp,ierror)
      call mpi_bcast(sigma,1,MPI_REAL,0,nallgrp,ierror)
      call mpi_bcast(u0,1,MPI_REAL,0,nallgrp,ierror)
      call mpi_bcast(new,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(time,1,MPI_REAL,0,nallgrp,ierror)
      call mpi_bcast(idp,1,MPI_INTEGER,0,nallgrp,ierror)
      call mpi_bcast(eta,1,MPI_REAL,0,nallgrp,ierror)
      call mpi_bcast(kf,1,MPI_REAL,0,nallgrp,ierror) 
!-------------------------------
      mx2=my/2
      mx=my
      mmx2=mx2/nproc
      mmy=my/nproc
      pi = 4.0 * atan(1.0)
      nek = int(sqrt(2.0)*my/3.0)  
      etak = eta/104.            
      if (id.eq.0) print*, 'etak=',etak
      scale = 1.0/my/mx        
!---------------------------------
      jj = 1
      jjj = 821
      iii = 821
      dt_h = 0.50 * dt        
      idp1 = idp        
!----------------------------------
!... allocate memory...................................................

        allocate (vx(my,mmx2) )
        allocate (vy(my,mmx2) )
        allocate (ux(my,mmx2) )
        allocate (uy(my,mmx2) )
        allocate (wz(my,mmx2) )
        allocate (owz(my,mmx2) )
        allocate (kx(my,mmx2) )
        allocate (ky(my,mmx2) )
        allocate (k2(my,mmx2) )
        allocate (tmp(my,mmx2) )
        allocate (tmp1(my,mmx2) )
        allocate (ftmp(my,mmx2) )
	allocate (ftmp1(my,mmx2) )
        allocate (wt(mx2,mmy) )
        allocate (uxt(mx2,mmy) ) 
	allocate (vxt(mx2,mmy) )    
        allocate (wx(0:(mx2-1)) )
        allocate (wy(0:(my-1)) )
        allocate (ipx(0:mx2) )  
        allocate (ipy(0:my) )
        allocate (iseed(nproc) )
        if (id.eq.0) print *, '  memory allocated '

!***********************************
        if (id.eq.0) then
         fnm2='data/2048/hypervis/'//ncase//'/spectrum.dat'
         open(20,file=fnm2,status='unknown')
         fnm2='data/2048/hypervis/'//ncase//'/dissipat.dat'
         open(2,file=fnm2,status='unknown')
         fnm2='data/2048/hypervis/'//ncase//'/inispect.dat'
         open(70,file=fnm2,status='unknown')
         fnm2='data/2048/hypervis/'//ncase//'/enstflux.dat'
         open(82,file=fnm2,status='unknown')
        endif

!...INITIAL CONDITIONS
      if (new.eq.1) then
        if (id.eq.0) then
          numb = irand (iseed_m0)
          do i = 1,nproc
            iseed(i) = irand(0)
          enddo
        endif
        call mpi_bcast(iseed,nproc,MPI_INTEGER,0,nallgrp,ierror)
        irr = iseed(id+1)
        num =  irand(irr)

         wz = (0.0, 0.0)
         call gaussian(vx,pi)
         call gaussian(vy,pi)
         call wavenumber(kx,ky,k2,id)

!...following expression because tmp, kx, ky,k2 are real arrays
         tmp = (kx*real(vx) + ky*real(vy))/k2
         vx = cmplx(real(vx) - kx*tmp, aimag(vx))
         vy = cmplx(real(vy) - ky*tmp, aimag(vy))
         tmp = (kx*aimag(vx) + ky*aimag(vy))/k2
         vx = cmplx(real(vx), aimag(vx) - kx*tmp)
         vy = cmplx(real(vy), aimag(vy) - ky*tmp)

         tmp1 = sqrt(k2)+0.50
         tmp = u0 / (tmp1**2) * exp (-0.5*tmp1/ak0)
         vx = vx * tmp
         vy = vy * tmp
         call symmetrize(vx,id)
         call symmetrize(vy,id)
         
        wz = (0.,1.) * (kx*vy - ky*vx)
        call symmetrize(wz,id)

!--calculating initial spectrum------
        ! initial enstrophy spectrum
        tmp= 0.25*wz*conjg(wz) / k2
        if (id.eq.0)  tmp(:, 1)=0.50*tmp(:, 1)  
        do i=1,nek
            ek=0.0
            e_t=0.0
            ek=sum(tmp,mask=((sqrt(k2)-i).lt.0.5.and.  &
                  (sqrt(k2)-i).ge.-0.5))
            call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM, &
                             0,mpi_comm_world,ierror)
           if (id.eq.0) then
              write(70,1002) real(i),  e_t    
           endif
        enddo
        if (id.eq.0) close(70)

      else
!        call input1 (wz,idp,id,nallgrp)    
         call input (wz,idp,id,nallgrp)
         if (id.eq.0) then
           wz(1,1)=(0.0, 0.0)
         end if     
      endif

      call symmetrize (wz, id)
      owz = wz

!---Prepare for call mpifft----------

         ipx(0) = 0  
         ipy(0) = 0

!...*********** MAIN LOOP *************

      do istep = 0, nstep
         timer1 = mpi_wtime( )
         call wavenumber(kx,ky,k2,id)

!...WRITE ENERGY SPECTRUM
         if (mod(istep,ieout).eq.0) then
            tmp= 0.25*wz*conjg(wz) / k2
            if (id.eq.0) then
              tmp(:, 1)=0.50*tmp(:, 1)  
            end if
            write(20,*) jjj

            beta = 0.0  
            do i=1,nek           
                 ek=0.0
                 e_t=0.0
                 ek=sum(tmp,mask=((sqrt(k2)-i).lt.0.5.and.  &
                        (sqrt(k2)-i).ge.-0.5))
                 call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0, &
                                            nallgrp,ierror)
                 if (id.eq.0) then
                       write(20, 1006) real(i), e_t
                 endif    
             enddo
             jjj = jjj + 1
         endif

! !!!WARNING!!! need times 2 !!!
!...Calculate hyper-dissipation (more precisely)..........................
        if (mod(istep,ieout).eq.0) then
           tmp = 0.25*2*(sigma*sqrt(k2))**16*wz*conjg(wz) 
           if (id.eq.0)  tmp(:, 1) = .50 * tmp(:, 1)  
           ek = sum(tmp)        
           call mpi_reduce(ek,beta,1,MPI_REAL,MPI_SUM,0, &
                                            nallgrp,ierror)
            if (id.eq.0) then
              write(2, 1000)  time,beta
            endif
       endif
   
!...STORE VORTICITY TEMPORARILY
         tmp = real(wz)
         tmp1 = aimag(wz)

!...output K-space vorticity if it's time
         if ( (mod(istep,itout).eq.0).and.(istep.ne.0) ) then
             idp1=idp1+1
             call output (wz,idp1,id,nallgrp) 
!  for back-up.
             if (istep.eq.nstep) then
               call output (wz,idp1,id,nallgrp) 
             endif
         endif

!         timer4=mpi_wtime( )
!           if (id.eq.0)  print *, ' timer4: ',timer4-timer1

!-------------NONLINEAR TERM----------

!...TRANSFER VELOCITY AND VORTICITY TO X-SPACE
         vx = (0.0,1.0) * ky * wz / k2
         call symmetrize (vx, id)
         isign= -1      
         call newfft (vx,uxt,isign,ipx,ipy,wx,wy,id,nallgrp) 
         call newfft (wz,wt,isign,ipx,ipy,wx,wy,id,nallgrp)
         
         vxt=uxt ! vxt now is the real space vx

!...FORM THE PRODUCT V*W IN X-SPACE
         uxt=cmplx( real(uxt)*real(wt), aimag(uxt) )  
         uxt=cmplx( real(uxt), aimag(uxt)*aimag(wt) )   
         ! uxt now is the real space product vx wz

!...TRANSFORM BACK TO K-SPACE
         isign = 1       
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp)   
         ! Now vx is the k-space vx*wz
         call symmetrize ( vx, id)
         call dealiasing (vx, k2) 

         vy = vx * (0.0,-0.50)*kx 
         ! minus one half of the derivative of vx*wz?

!____________________________________________________________
!------------------------------------------------------------

!===============================================
!     now do on y-component of velocity
!     vort. in x-space is already saved in wt
!     wz also changed 

         vx=-(0.0,1.0)* kx*cmplx(tmp,tmp1)/k2   
         ! cmplx(tmp, tmp1) = wz, vorticity
         ! This vx is actually k-space vy

         isign = -1
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp)

	 vxt=uxt ! now vxt is the realspace vy.

!-------------------------------------------------------
         uxt=cmplx( real(uxt)*real(wt), aimag(uxt) )
         uxt=cmplx( real(uxt), aimag(uxt)*aimag(wt) )
         ! uxt is now real-space vy*wz
         
         isign = 1
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp) 
         call symmetrize (vx, id)
         call dealiasing (vx, k2)
         ! vx is now k-space vy*wz
   
         vy = vy + (0.0,-0.50) * ky * vx
         ! the second term is the one half of the derivative of vy*wz
         ! The sum is minus one half of the convection term
         ! the other half seems to come from phase shift

!*******************************************************         

!-Recover vort.,do phase shift dealiasing on x-comp of velo.
         wz = cmplx (tmp, tmp1)
         vx = (0.0,1.0) * ky * wz / k2

!...PHASE SHIFT
      vx=vx*cmplx(cos(pi/my*(kx+ky)),sin(pi/my*(kx+ky)))
      wz=wz*cmplx(cos(pi/my*(kx+ky)),sin(pi/my*(kx+ky)))

         call symmetrize(vx,id)
         call symmetrize(wz,id)


!...TRANSFER VELOCITY AND VORTICITY TO X-SPACE
         isign = -1
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp) 
         call newfft (wz,wt, isign,ipx,ipy,wx,wy,id,nallgrp) 
	 vxt=uxt 

!...FORM THE PRODUCT V*W IN X-SPACE
         uxt=cmplx( real(uxt)*real(wt), aimag(uxt) )
         uxt=cmplx( real(uxt), aimag(uxt)*aimag(wt) )
        
!...TRANSFORM TO K-SPACE
         isign = 1
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp) 
!...PHASE SHIFT
         vx=vx*cmplx(cos(pi/my*(kx+ky)),-sin(pi/my*(kx+ky)) )
         call symmetrize(vx, id)

         vy = vy + vx*(0.0,-0.50)*kx

!###################################################

!--- do phase shift dealiasing for y-comp
!--- phase-shifted wz in real space already saved in wt

         vx=-(0.0,1.0)* kx*cmplx(tmp, tmp1)/k2

         vx=vx*cmplx(cos(pi/my*(kx+ky)),sin(pi/my*(kx+ky)))

         call symmetrize(vx,id)

         isign = -1
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp) 

         vxt=uxt

         uxt=cmplx( real(uxt)*real(wt), aimag(uxt) )
         uxt=cmplx( real(uxt), aimag(uxt)*aimag(wt) )

         isign = 1
         call newfft (vx,uxt, isign,ipx,ipy,wx,wy,id,nallgrp) 
         vx=vx*cmplx(cos(pi/my*(kx+ky)),-sin(pi/my*(kx+ky)))

         call symmetrize(vx, id)

         vy = vy + vx*(0.0,-0.50)*ky

!*********************************************************

!-----------END CONVOLUTION PLUS DEALISING-----------
!-----now, vy = -i k.fft(Vw)
!...RECOVER VORTICITY
           wz = cmplx( tmp,tmp1)
!...CALCULATE ENSTROPHY FLUX (it must be done HERE)
!...enstrophy flux is
!   2 Imag [ w^*(k) k.F(uw) ]  (F(.) means fourier transform)
         if (mod(istep,ieout).eq.0.and.istep.ne.0) then
            vx = 0.25* 2.* (0.0, 1.0)*vy
            kx = real(wz)*aimag(vx) - aimag(wz)*real(vx)
            if (id.eq.0)  kx (:,1)=0.50*kx (:,1)  
            if (id.eq.0)  write(82, *) iii
         
            sum1=0.0 
            do i = 1,nek          
               a = i - 1 + 0.5
               tt =sum (kx,mask=(sqrt(k2).ge.a.and.sqrt(k2).lt.a+1) )
               ttt = 0.0
               call mpi_reduce(tt, ttt, 1, MPI_REAL,MPI_SUM,0, &
                nallgrp,ierror)
               if (id.eq.0) then
                   sum1=sum1+ttt
                   write(82, 1002) real(i), -sum1/beta 
               endif
            end do
            iii = iii + 1

         endif

!-------------END NONLINEAR TERM----------
         call dealiasing (vy, k2)
         call symmetrize (vy, id)
         call dealiasing (uy, k2)
         call symmetrize (uy, id)
!---If first-step then use modified Euler method
         if (istep.eq.0) then
            where ( sqrt(k2).ge.kf.and.sqrt(k2).le.(kf+3.) )   
              wz =wz * (1.0-((sqrt(k2)*sigma)**16 +  &
                             (sqrt(k2)/nu)**(-4))*dt_h ) +  &
                        dt_h*( vy + etak/conjg(wz) )
            elsewhere 
              wz = wz * (1.0-((sqrt(k2)*sigma)**16 +  &
                         (sqrt(k2)/nu)**(-4))*dt_h ) + dt_h*vy
            endwhere
            call output(vy,0,id,nallgrp)
            time=time+dt_h
         
         else if (istep.eq.1) then
             where ( sqrt(k2).ge.kf.and.sqrt(k2).le.(kf+3.) )   
               wz = owz + dt *( vy + etak/conjg(wz)) -    &
                          ((sqrt(k2)*sigma)**16 +  &
                          (sqrt(k2)/nu)**(-4)) * dt * wz 
             elsewhere                                         
               wz = owz + dt *vy-((sqrt(k2)*sigma)**16 +  &
                            (sqrt(k2)/nu)**(-4)) * dt * wz
             endwhere
             call input(owz,0,id,nallgrp)
             time=time+dt_h

!---Adams-banshford--------
          else
            where (k2.gt.0.9) 
              kx=exp (-((sqrt(k2)*sigma)**16 +  &
                     (sqrt(k2)/nu)**(-4)) *dt)		       
            endwhere  
            if (id.eq.0) kx(1,1)=1.0
			
            where ( sqrt(k2).ge.kf.and.sqrt(k2).le.(kf+3.))   
              wz = wz + dt_h * ( 3.0* (vy + etak/conjg(wz))  &
                      - kx * (owz + etak/conjg(owz)))
             elsewhere
               wz = wz + dt_h * (3.0* vy-kx * owz)
             endwhere
             wz = wz * kx
             owz=vy
             time=time+dt
          end if
       
        if (mod(istep,ieout).eq.0) then
           timer2=mpi_wtime( )
           if (id.eq.0)  print *, ' time for one step: ',timer2-timer1
        endif

        enddo                     

!===================================================
!...Output k-space vorticity for next run        
!       call output (wz,idp1+1,id,nallgrp)         

      if (id.eq.0) then
         close(20)
         close(2)
         close(70)
         close(82)
      endif

1004    format( i8,' ',f14.7,' ',f14.8)
1002    format( f10.1,' ',e14.8)
1000    format( f11.6,'  ', e14.6)
1006    format( f12.6,'  ', 2e16.8)

      call MPI_FINALIZE(ierror)
      stop
      end

!--------------------------------------------

      subroutine dealiasing(vx,k2)
!...8/9 rule for dealiasing convolution in FFT

      include 'include.h'
      complex, dimension(my, mmx2)::vx
      real,dimension (my, mmx2)::k2
      real ass

      ass = 2.0/9.0*(my*my)
      where(k2.ge.ass)
        vx =  (0.0, 0.0)
      endwhere
      return
      end
!----------------------------------------------
      subroutine wavenumber(kx,ky,k2,id)

      include 'include.h'
      real,dimension (my, mmx2)::kx,ky,k2
      integer id
      integer i, j, j1

      do j=1,mmx2
         j1 = id * mmx2 + j
         kx (:, j) = float( j1-1 )
      enddo

      do i=1,my
         ky (i, :) = float(mod( i-1+my/2, my)-my/2)
      enddo

      k2 = kx*kx + ky*ky
      if (id.eq.0) then
         k2(1,1) = 0.5         
      endif

      return
      end

!------------------------------------------------
      subroutine gaussian (u,pi)
      include 'include.h'
      complex, dimension(my, mmx2) :: u
      real pi,t1,t2

      u  = (0.0, 0.0)
      do i = 1,my
         do j = 1,mmx2
           call RANDOM_NUMBER(t1)
           call RANDOM_NUMBER(t2)
!           t1 = drand(0)
!           t2 = drand(0)
            if (t1.le.1.e-10) t1 = 1.e-10
            if (t2.le.1.e-10) t2 = 1.e-10
            t2 = 2.0*pi*t2
            u (i, j) = sqrt(-2.0*log(t1))*cmplx(cos(t2),sin(t2)) 
         enddo
      enddo

      return
      end

!-----------------------------------------------------

      subroutine symmetrize(c,id)
      include 'include.h'
      complex,dimension(my, mmx2)::c

        c(my/2+1,:) = (0.0, 0.0)

        if (id.eq.0) then
          c(1,1) = (0.0, 0.0)
          do iy = 2,my/2-1
            c(iy,1) = .50*(c(iy,1)+ conjg(c(mod(my+1-iy,my)+1,1)))  
          enddo

           do iy = 2, my/2-1
             iy2 = mod(my+1-iy,my) + 1
             c(iy2, 1) = conjg (c(iy, 1) )
            enddo
           
         endif

      return
      end
! ---------------------------------------------------
      subroutine output (ux,idp,id,nallgrp)

      include 'include.h'

      complex,dimension(my, mmx2)::ux
      integer i1d,i2d,i3d,i4d,i5d,i6d,idp,id
      character(70) fnm1
      character(5) ncase
      common/casenum/ncase
      i1d=int(id/100)
      i2d=int((id-i1d*100)/10)
      i3d=id-i1d*100-i2d*10
  
      i4d=int(idp/100)
      i5d=int((idp-i4d*100)/10)
      i6d=idp-i4d*100-i5d*10
  
      fnm1='data/2048/hypervis/'//ncase//'/vort'//   &          
!     fnm1='vort'//  &
     
         char(i4d+48)//char(i5d+48)//char(i6d+48)//'.'// &
         char(i1d+48)//char(i2d+48)//char(i3d+48)

      open(10,file=fnm1,status='unknown', form='unformatted')
         write(10) ux 
      close(10)

      return
      end
! ---------------------------------------------------
      subroutine input (ux,idp,id,nallgrp)

      include 'include.h'

      complex,dimension(my, mmx2)::ux
      integer i1d,i2d,i3d,i4d,i5d,i6d,idp,id
      character(70) fnm1
      character(5) ncase
      common/casenum/ncase
      i1d=int(id/100)
      i2d=int((id-i1d*100)/10)
      i3d=id-i1d*100-i2d*10

      i4d=int(idp/100)
      i5d=int((idp-i4d*100)/10)
      i6d=idp-i4d*100-i5d*10       

      fnm1='data/2048/hypervis/'//ncase//'/vort'//  & 
!     fnm1='vort'//  &
         char(i4d+48)//char(i5d+48)//char(i6d+48)//'.'//  &
         char(i1d+48)//char(i2d+48)//char(i3d+48)

      open(10,file=fnm1,status='unknown',form='unformatted')
         read (10) ux
      close(10)

      return
      end

!******************************************************
! --------------------------------------------------
        subroutine newfft (ux,uxt, isign,ipx,ipy,wx,wy,id,nallgrp)
               
        include 'include.h'
        complex, dimension (mx2, mmy):: uxt
        complex, dimension (my, mmx2):: ux
        real  ux1(0:mx-1)
        complex  uy1(0:my-1)
       
        real  wx (0:(mx2-1) ) 
        real  wy (0:(my-1) )
        integer ipx (0:mx2)         
        integer ipy (0:my)
        integer isign,i, j
        real time1,time2

       if (isign.eq.1) then     
!-- rc fft in x-dir-------------
        do j=1, mmy
           do i=1, mx2
              ux1(2*i-2) = real(uxt (i, j) )
              ux1(2*i-1) = aimag(uxt (i, j) )
           enddo

           call rdft (mx, isign, ux1, ipx, wx)
           ux1=ux1/float(mx2)
           ux1(1)=0.0

           do i=1, mx2
              uxt(i, j)=cmplx( ux1(2*i-2), ux1(2*i-1) )
           enddo
        enddo
        
        call transpose_xtok (uxt, ux, id, nallgrp)

!-- cc fft in y-dir-------------       
        do i=1, mmx2
            uy1=ux(:, i)
            call cdft (my*2, isign, uy1, ipy, wy)
            ux(:, i)= uy1/float(my)
         enddo    
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
        else if (isign.eq. -1)  then  
         do i=1, mmx2
            uy1=ux(:, i)
            call cdft (my*2, isign, uy1, ipy, wy)
            ux(:, i)= uy1
         enddo

        call transpose_ktox (ux, uxt, id, nallgrp)
        do j=1, mmy
           do i=1, mx2
              ux1(2*i-2) = real(uxt (i, j) )
              ux1(2*i-1) = aimag(uxt (i, j) )
            enddo
            ux1(1)=0.0
            call rdft (mx, isign, ux1, ipx, wx)
            do i=1, mx2
               uxt(i, j)=cmplx( ux1(2*i-2), ux1(2*i-1))
            enddo
        enddo     
       endif

       return
       end
! -------------------------------------------------
      subroutine transpose_ktox (ux,uxt,id,nallgrp)
!     transpose ux to uxt so can do x-dir fft on uxt

      include 'include.h'
      complex,dimension(my,mmx2)::ux
      complex,dimension(mx2,mmy)::uxt
      complex,dimension(mmy,mmx2)::tmp1,tmp2  
      complex,dimension(mmy,mx2)::tmp
      integer isize,nzm,nzp,status(MPI_STATUS_SIZE,2),req(2)
      integer i,j,k,js,j1,ks,k1, l
 
      if (nproc.gt.1) then
        isize = mmy*mmx2
        do i = 1,nproc-1
         nzp=mod(id+i,nproc)
         nzm=mod(id-i+nproc,nproc)
         js = nzp*mmy

         do j = 1,mmy
              j1=js+j
              tmp1( j, : ) = ux ( j1, :)
         enddo
         call mpi_isend(tmp1, isize, MPI_COMPLEX, nzp, i, &
             nallgrp,req(1),ierror)
         call mpi_irecv(tmp2, isize, MPI_COMPLEX, nzm, i, &
             nallgrp,req(2),ierror)
         call mpi_waitall (2,req,status,ierror)
       
         ks = nzm*mmx2
         do k = 1, mmx2
              k1 = ks+k
              tmp( :, k1)= tmp2 ( :, k)
         enddo
       enddo

!     does (id,id) spot from ux to tmp so can transpose to uxt

       ks = id*mmx2
       js = id*mmy
       do k = 1,mmx2
         k1 = ks + k
         do j = 1,mmy
            j1 = js + j
            tmp( j, k1 ) = ux ( j1, k )
         enddo
       enddo

         do k=1, mmy
            do j=1, mx2
                uxt ( j, k ) = tmp(k, j )
            enddo
         enddo
       
       else if (nproc.eq.1) then
         do j=1, mx2
           do k=1, my
                uxt ( j, k)= ux (k, j ) 
           enddo
        enddo
       end if
       
      return
      end

!------------------------------------------------
      subroutine transpose_xtok (uxt,ux,id,nallgrp)
!     transpose uxt to ux so can do y-dir ifft on ux

      include 'include.h'
      complex,dimension(my,mmx2)::ux
      complex,dimension(mx2,mmy)::uxt
      complex,dimension(mmx2, mmy)::tmp1,tmp2
      complex,dimension(mmx2,my)::tmp
      integer isize,nzm,nzp,status(MPI_STATUS_SIZE,2),req(2)
      integer i,j,k,js,j1,ks,k1, m

      if (nproc.gt.1) then
       isize = mmy*mmx2
       do i = 1,nproc-1
         nzp=mod(id+i,nproc)
         nzm=mod(id-i+nproc,nproc)
         js = nzp*mmx2
         do j = 1,mmx2
              j1 = js+j
              do m=1, mmy
                    tmp1(j, m) = uxt( j1, m) 
              enddo
         enddo
          call mpi_isend(tmp1, isize, MPI_COMPLEX, nzp, i, &
             nallgrp,req(1),ierror)
         call mpi_irecv(tmp2, isize, MPI_COMPLEX, nzm, i, &
             nallgrp,req(2),ierror)
         call mpi_waitall(2,req,status,ierror)

         ks = nzm*mmy
         do k = 1,mmy
              k1 = ks+k
              tmp( :, k1) = tmp2( :, k)
         enddo
       enddo

!     does the (id,id) spot from uxt to tmp so can transpose to ux
        ks = id*mmy
        js = id*mmx2
        do k = 1,mmy
            k1 = ks+k
           do j = 1,mmx2
              j1 = js+j
              tmp( j, k1) = uxt ( j1, k)
           enddo
       enddo
!---important  Transpose here!!
       do k = 1,mmx2
         do j = 1,my
            ux( j, k) = tmp( k, j)
         enddo
       enddo

      else if (nproc.eq.1) then
        do j=1, my  
          do k=1, mmx2  
             ux (j, k)= uxt (k, j) 
          enddo
        enddo
      end if

      return
      end
!----------------------------------
!      include 'fft1d.f'

