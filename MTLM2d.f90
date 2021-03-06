MODULE global
    real(8), parameter :: pi = 3.141592653589793238462643d0, two_Pi = 2.d0 * pi
 
    !========  parameter for velocity field   =========
    ! For 3D turbulence

    !real(8), parameter :: urms = 1.000d00
    !real(8), parameter :: Ldom = 6.283185308d0  ! = 2 pi
    !real(8), parameter :: phi = pi
    !real(8), parameter :: cl = 9.811747074127197d-01   
    !real(8), parameter :: Ck = 1.50d+00
    !real(8), parameter :: l = 2.073500000000000d+00
    !real(8), parameter :: p0 = 4.d0

    ! For 2D inverse cascade
    real(8), parameter :: diss = 0.1d0
    real(8), parameter :: cl = 1.d0  !! Use as an approximate to the 3d value. 
    real(8), parameter :: p0 = 2.d0  !! Kind of random guess. 
    real(8), parameter :: ck = 6.d0  !! Coefficient taken from Boffettaetal00 
    real(8), parameter :: l = 0.3667d0   !! The Ekman frictional scale calculated.
      ! Note that l is also the integral length scale, used to estimate epsilon. 
      ! In 2D it is to be interpreted as the Ekman friction length scale. 
    
    !========  parameter for scalar field   =========
    real(8), parameter :: theta_rms = 1.000d00
    real(8), parameter :: qdiss = 2.893658066071859d-01
    real(8), parameter :: gamma_diff = 1.955936002032939d-03
    real(8), parameter :: Cq =  6.800000d-01
    real(8), parameter :: clq = 2.562338113784790d-01
    real(8), parameter :: etaq = 1.116071428717157d-02
    real(8), parameter :: lq = 1.430000000000000d+00
    
    !----------------------------------------------
    integer, parameter :: nt = 256
    integer, parameter :: nx = nt, ny = nt
    integer, parameter :: ld = nx + 2
    integer, parameter :: lh = ld/2
    
    
    complex(8), dimension(nt/2 + 1,nt) :: vx,vy,theta
 
 
    real(8) E(lh-2),Et(lh-2),eta,alpha
    real(8) hx,hy,CFL,max_k2,resc_factor,epsilon,epsl,tiny_real
    real(8) urms0,thetarms0,Int0kmax,Int0kmaxscalar
    integer nxh,nyh,nyhp,lh2,lhm,nth,nbytes
    ! logical output

END MODULE
  
  


PROGRAM MTLM_uT
    use global
    implicit none
    
    !....... definitions in fftw3.f ......................................
    integer,parameter::FFTW_R2HC = 0   , FFTW_HC2R = 1    , FFTW_DHT = 2
    integer,parameter::FFTW_REDFT00 = 3, FFTW_REDFT01 = 4 , FFTW_REDFT10 = 5 
    integer,parameter::FFTW_REDFT11 = 6, FFTW_RODFT00 = 7 , FFTW_RODFT01 = 8 
    integer,parameter::FFTW_RODFT10 = 9, FFTW_RODFT11 = 10, FFTW_FORWARD = -1
    
    integer,parameter::FFTW_BACKWARD = +1          , FFTW_MEASURE = 0  
    integer,parameter::FFTW_DESTROY_INPUT = 1      , FFTW_UNALIGNED = 2
    integer,parameter::FFTW_CONSERVE_MEMORY = 4    , FFTW_EXHAUSTIVE = 8
    integer,parameter::FFTW_PRESERVE_INPUT = 16    , FFTW_PATIENT = 32
    integer,parameter::FFTW_ESTIMATE = 64          , FFTW_ESTIMATE_PATIENT = 128
    integer,parameter::FFTW_BELIEVE_PCOST = 256    , FFTW_DFT_R2HC_ICKY = 512
    integer,parameter::FFTW_NONTHREADED_ICKY = 1024, FFTW_NO_BUFFERING = 2048
    integer,parameter::FFTW_NO_INDIRECT_OP = 4096  , FFTW_ALLOW_LARGE_GENERIC = 8192
    integer,parameter::FFTW_NO_RANK_SPLITS = 16384 , FFTW_NO_VRANK_SPLITS = 32768
    integer,parameter::FFTW_NO_VRECURSE = 65536    , FFTW_NO_SIMD = 131072
    !......................................................................
    
    !...for reduced velocity field
    integer, parameter :: nbnd = 8

    complex(8), dimension(nbnd/2*16 + 1, nbnd*16 ) :: vx128,vy128
    complex(8), dimension(nbnd/2*8  + 1, nbnd*8  ) :: vx64 ,vy64
    complex(8), dimension(nbnd/2*4  + 1, nbnd*4  ) :: vx32 ,vy32
    complex(8), dimension(nbnd/2*2  + 1, nbnd*2  ) :: vx16 ,vy16
    complex(8), dimension(nbnd/2    + 1, nbnd    ) :: vx08 ,vy08
    
    !...for reduced scalar field
    complex(8), dimension(nbnd/2*16 + 1, nbnd*16 ) :: theta128
    complex(8), dimension(nbnd/2*8  + 1, nbnd*8  ) :: theta64
    complex(8), dimension(nbnd/2*4  + 1, nbnd*4  ) :: theta32
    complex(8), dimension(nbnd/2*2  + 1, nbnd*2  ) :: theta16
    complex(8), dimension(nbnd/2    + 1, nbnd    ) :: theta08
    
    real(8), dimension(lh-2) :: Eloc
    
    integer(8) FTvx,FTvy,IFTvx,IFTvy
    integer(8) FTvx128,FTvy128,IFTvx128,IFTvy128
    integer(8) FTvx64,FTvy64,IFTvx64,IFTvy64
    integer(8) FTvx32,FTvy32,IFTvx32,IFTvy32
    integer(8) FTvx16,FTvy16,IFTvx16,IFTvy16
    integer(8) FTvx08,FTvy08,IFTvx08,IFTvy08
    
    integer(8) FTth,IFTth
    integer(8) FTth128,IFTth128
    integer(8) FTth64,IFTth64
    integer(8) FTth32,IFTth32
    integer(8) FTth16,IFTth16
    integer(8) FTth08,IFTth08
    
    complex(8) cn1,cn2,cn3,cp1,cp2,cp3
    complex(8) a,b,c,t,vxx,vyy
    
    real(8) d,mag,ratio,wksp,rn,wn,t2,fact,const,cns,ur,u1,u2,v1,v2,tadv
    real(8) hb,spec2d,k2r,kx,k22,norm_factor,s2,upl,scalar_spec2d,theta_rmsl
    real(8) Rex,Imx,Rey,Imy,Rez,Imz,ss,tn,time_scale,total_adv,next_total
    real(8) thetapl
    
    
    real ran2
    
    integer num(nt+1),ind,jx,jy,idum,ans,ky,r,rs,nadv,ad
    integer kyky,bl,nnr
    integer jxr,jyr,jx2,jy2,i,j,k,ii
    integer n,nit,kxi
    !logical last_adv
    character(24) cmd
    character(4) gen
    
    
    ind = iarg()
    if ( ind .ne. 2 ) stop "need two arguments: ./MTLM2d.x -idum nfile"
    
    call getarg(1,cmd)
    read(cmd, '(I6)') idum
    idum = -idum

    call getarg(2,cmd) ! cmd is the output file number.
    cmd = adjustl(cmd)

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (mod(ld,2) /= 0) stop ' ld must be even'

    ! all in global
    nxh = nx/2
    nyh = ny/2
    nyhp = nyh + 1
    lh2 = lh*lh
    lhm = lh - 1
    nth = nt/2
    
    max_k2 = real((nt/2 + 1)**2,kind=8)
    
    hx = two_Pi/nx
    hy = two_Pi/ny
    
    tiny_real =tiny(d)
    t2 = dble(ran2(idum))    !...initialize random generator
    

    ! Calculate the total KTE and urms from the specified analytic form of the spectrum
    s2 = 0.5d0 * spec2d(1.0d0)
    do jx = 2, lh-2
      s2 = s2 + 0.5d0*(spec2d(dble(jx-1)) + spec2d(dble(jx)))
    end do
    Int0kmax = s2
    urms0 = sqrt(Int0kmax)
    ! The spectrum and urms0 are to be maintained through the mapping. 

    ! The Ekman frictional coefficient.
    alpha = diss / Int0kmax / 2.
    
    write(*,*) 'Parameters: diss, l, urms0, alpha, nx, ny'
    write(*,*) diss, l, urms0, alpha, nx, ny

    
    write(*,"(/,' initialising FFT...')")
    nnr = nt
    
    print *,nnr
    call dfftw_plan_dft_c2r_2d(IFTvx,nnr,nnr,vx,vx,FFTW_MEASURE)
    call dfftw_plan_dft_c2r_2d(IFTvy,nnr,nnr,vy,vy,FFTW_MEASURE)
    
    call dfftw_plan_dft_r2c_2d(FTvx,nnr,nnr,vx,vx,FFTW_MEASURE)
    call dfftw_plan_dft_r2c_2d(FTvy,nnr,nnr,vy,vy,FFTW_MEASURE)

    nnr = nbnd*16
    print *,nnr
    call dfftw_plan_dft_c2r_2d(IFTvx128,nnr,nnr,vx128,vx128,FFTW_MEASURE)
    call dfftw_plan_dft_c2r_2d(IFTvy128,nnr,nnr,vy128,vy128,FFTW_MEASURE)
    
    call dfftw_plan_dft_r2c_2d(FTvx128,nnr,nnr,vx128,vx128,FFTW_MEASURE)
    call dfftw_plan_dft_r2c_2d(FTvy128,nnr,nnr,vy128,vy128,FFTW_MEASURE)
    
    nnr = nbnd*8
    print *,nnr
    call dfftw_plan_dft_c2r_2d(IFTvx64,nnr,nnr,vx64,vx64,FFTW_MEASURE)
    call dfftw_plan_dft_c2r_2d(IFTvy64,nnr,nnr,vy64,vy64,FFTW_MEASURE)
    
    call dfftw_plan_dft_r2c_2d(FTvx64,nnr,nnr,vx64,vx64,FFTW_MEASURE)
    call dfftw_plan_dft_r2c_2d(FTvy64,nnr,nnr,vy64,vy64,FFTW_MEASURE)
    
    nnr = nbnd*4
    print *,nnr
    call dfftw_plan_dft_c2r_2d(IFTvx32,nnr,nnr,vx32,vx32,FFTW_MEASURE)
    call dfftw_plan_dft_c2r_2d(IFTvy32,nnr,nnr,vy32,vy32,FFTW_MEASURE)
    
    call dfftw_plan_dft_r2c_2d(FTvx32,nnr,nnr,vx32,vx32,FFTW_MEASURE)
    call dfftw_plan_dft_r2c_2d(FTvy32,nnr,nnr,vy32,vy32,FFTW_MEASURE)
    
    nnr = nbnd*2
    print *,nnr
    call dfftw_plan_dft_c2r_2d(IFTvx16,nnr,nnr,vx16,vx16,FFTW_MEASURE)
    call dfftw_plan_dft_c2r_2d(IFTvy16,nnr,nnr,vy16,vy16,FFTW_MEASURE)
    
    call dfftw_plan_dft_r2c_2d(FTvx16,nnr,nnr,vx16,vx16,FFTW_MEASURE)
    call dfftw_plan_dft_r2c_2d(FTvy16,nnr,nnr,vy16,vy16,FFTW_MEASURE)
    
    nnr = nbnd
    print *,nnr
    call dfftw_plan_dft_c2r_2d(IFTvx08,nnr,nnr,vx08,vx08,FFTW_MEASURE)
    call dfftw_plan_dft_c2r_2d(IFTvy08,nnr,nnr,vy08,vy08,FFTW_MEASURE)
    
    call dfftw_plan_dft_r2c_2d(FTvx08,nnr,nnr,vx08,vx08,FFTW_MEASURE)
    call dfftw_plan_dft_r2c_2d(FTvy08,nnr,nnr,vy08,vy08,FFTW_MEASURE)
    
    
    !................... end of initialisation .....................
    
    
    ! calculate num(k)
    ! this is the number of wavenumbers in shell k
    num = 0
    do jy = 1,ny
      if (jy == nyhp) cycle
      ky = jy - 1
      if (jy > nyh) ky = ky - ny
      kyky = ky * ky   

      do jx = 1,lhm
        kx = real(jx - 1)
        k22 = kx*kx + kyky 
        if ((k22 < 1.d-20) .or. (k22 >= max_k2)) cycle
        ind = int(sqrt(k22) + 0.5d0) 
        if (kx < 1.d-20) then
          num(ind) = num(ind) + 1
        else
          num(ind) = num(ind) + 2
        end if
      end do

    end do
    
    !...Calculate random velocity field.  abs(v(k))=sqrt(2E(k)/N(k))
    do jy = 1,ny
      if (jy == nyhp) cycle
      ky = jy - 1
      if (jy > nyh) ky = ky - ny
      kyky = ky * ky   
 
      do jx = 1,nth
        kx = dble(jx - 1)
        k22 = kx*kx + kyky 
        if (k22 < 1.d-20) cycle
    
        if(k22 < max_k2) then
          wn = sqrt(k22)
          ind = int(wn + 0.5d0)
          if (ind < 1) ind = 1
    
          t2 = dble(ran2(idum))
          if (t2 <= 1.d-10) t2 = 1.d-10
          t2 = two_Pi*t2
          a = cmplx(dcos(t2),dsin(t2),kind=8)
    
          t2 = dble(ran2(idum))
          if (t2 <= 1.d-10) t2 = 1.d-10
          t2 = two_Pi*t2
          b = cmplx(dcos(t2),dsin(t2),kind=8)
 
          d = sqrt(2.d0*spec2d(wn)/dble(num(ind)))
          t = (a*kx + b*ky )/k22
          vxx = a - t*kx
          vyy = b - t*ky
          mag = sqrt( vxx*conjg(vxx) + vyy*conjg(vyy) )
    
          vx(jx,jy) = vxx/mag*d     !u_i k_i=0
          vy(jx,jy) = vyy/mag*d
        end if

      end do
    end do
    
    
    !...first spectrum, before re-scaling
    call spectrum(8,'spec001.dat', .false.)
    ! subroutine spectrum postprocesses and saves initial data.
    ! E(k) is calculated in subroutine spectrum
    
    Eloc = 0.d0
    !...re-scale spectrum to improve fitting at low wave numbers
    do jy = 1,ny

      if (jy == nyhp) cycle  ! nyhp = ny halved plus 1
      ky = jy - 1
      if (jy > nyh) ky = ky - ny
      kyky = ky * ky
    
      do jx = 1, nth

        kx = dble(jx - 1)
        k22 = kx*kx + kyky

        if (k22 < 1.d-20) cycle
        if (k22 < max_k2) then
          wn = sqrt(k22)
          ind = int(wn + 0.5d0)
    
          if (ind <= 1) then
            ind = 1
            wn = 1.0d0
          end if  
    
          if (ind < nth) then ! nth = half of nt
            ! E(k) is calculated in subroutine spectrum
            if (E(ind) < 2.d-38) then
               ratio = 0.0d0
            else  
               ratio = sqrt(spec2d(wn)/E(ind))
            end if    
            vx(jx,jy) = vx(jx,jy)*ratio
            vy(jx,jy) = vy(jx,jy)*ratio
      
            !...compute contribution to the new energy spectrum
            Rex =  real(vx(jx,jy))
            Imx = aimag(vx(jx,jy))
            Rey =  real(vy(jx,jy))
            Imy = aimag(vy(jx,jy))
            ss = Rex*Rex + Imx*Imx +  Rey*Rey + Imy*Imy

            if (kx < 1.d-20) ss = ss * 0.5d0
            if (ss < tiny_real) ss = 0.d0
            Eloc(ind) = Eloc(ind) +  ss
          else
            vx(jx,jy) = (0.d0,0.d0)
            vy(jx,jy) = (0.d0,0.d0)
          end if

        else
          vx(jx,jy) = (0.d0,0.d0)
          vy(jx,jy) = (0.d0,0.d0)
        end if
    
      end do
    end do
    
    !...second rescale to match the kinetic energy
    s2 = 0.5d0 * Eloc(1)
    do jx = 2,lh-2
      s2 = s2 + 0.5d0*(Eloc(jx-1) + Eloc(jx))
    end do
    upl = sqrt( s2  )
    resc_factor = urms0/upl  !!! TODO: Why is this urms instead of urms0?
    
    vx = vx * resc_factor
    vy = vy * resc_factor
    
    call spectrum(1,'spec002.dat', .true.)
    epsilon = epsl 
    !TODO: Check this. epsl is calculated in spectrum, the dissipation in the initial field

    ! Initial velocity field 
    open(unit=21, file='./out/ux0.dat', form='unformatted')
      write(21) vx
    close(unit=21)
    open(unit=21, file='./out/uy0.dat', form='unformatted')
      write(21) vy
    close(unit=21)
    
    !============================================================================
    
    !open(unit=35, file='control.out')
    
    !...apply the MTLM procedure
    call apply_MTLM(vx08,vy08,theta08,nbnd + 2,                         &
                    IFTvx08,IFTvy08,IFTth08,FTvx08,FTvy08,FTth08)
    
    call apply_MTLM(vx16,vy16,theta16,nbnd*2+ 2,                         &
                    IFTvx16,IFTvy16,IFTth16,FTvx16,FTvy16,FTth16)
    
    call apply_MTLM(vx32,vy32,theta32, nbnd*4+ 2,                          &
                    IFTvx32,IFTvy32,IFTth32,FTvx32,FTvy32,FTth32)
    
    call apply_MTLM(vx64,vy64,theta64, nbnd*8 + 2,                          &
                    IFTvx64,IFTvy64,IFTth64,FTvx64,FTvy64,FTth64)
    
    call apply_MTLM(vx128,vy128,theta128, nbnd*16 + 2,IFTvx128,IFTvy128,    &
                    IFTth128,FTvx128,FTvy128,FTth128)
    
    !...the last scale is direct (no filtered fields):
    
    write(*,"('using epsilon value =',ES25.15)") epsilon
    

    write(*,"('52/',$)") 
    call dfftw_execute(IFTvx)
    call dfftw_execute(IFTvy)
    !!!!! call dfftw_execute(IFTth)
    write(*,"('53/',$)") 
    
    
    ! compute urms  
    ur = 0.0d0
    do jy = 1,ny
    do jx = 1,nxh
       u1 =  abs(real(vx(jx,jy)))
       u2 = abs(aimag(vx(jx,jy)))
    
       v1 =  abs(real(vy(jx,jy)))
       v2 = abs(aimag(vy(jx,jy)))
    
       ur = ur + u1*u1 + v1*v1 + u2*u2 + v2*v2 
    end do
    end do
    ur = sqrt(ur/(2.0d0*nx*ny))
    tn = hx/ur
    time_scale = (hx*hx/epsilon)**(1.d0/3.d0)
    
    write(*,"(/,'    urms,n = ',ES12.5)") ur
    write(*,"('        tn = ',ES12.5)") tn
    write(*,"('time_scale = ',ES12.5)") time_scale
    write(911,"(/,'    urms,n = ',ES24.14)") ur
    write(911,"('        tn = ',ES24.14)") tn
    write(911,"('time_scale = ',ES24.14)") time_scale
    
    
    nadv = nint(time_scale/tn)
    tadv = time_scale/nadv
    CFL = tadv/tn
    if (CFL > 1.05) then
      nadv = nadv + 1
      tadv = time_scale/nadv
      CFL = tadv/tn
    end if  
    write(*,"('       CFL = ',ES12.5)") CFL
    write(911,"('       CFL = ',ES12.5)") CFL  
    

    total_adv = 0.d0
    ! last_adv = .false.
    for_time_scale: DO ad = 1,nadv
       next_total = total_adv + tadv
    
       !... Apply MTLM on disk
       call MTLM(vx,vy,theta,ld,tadv,urms0,thetarms0)
    
       !... Back to Fourier space
       write(*,"('54/',$)") 
       call dfftw_execute(FTvx)
       call dfftw_execute(FTvy)
       write(*,"('55/',$)") 
    
       norm_factor = 1.d0/dble(nx*ny)
       vx = vx*norm_factor
       vy = vy*norm_factor
    
       vx(lh,1:ny) = (0.d0,0.d0)      !Oddballs
       vy(lh,1:ny) = (0.d0,0.d0)
       vx(1:lh,nyhp) = (0.d0,0.d0)
       vy(1:lh,nyhp) = (0.d0,0.d0)
       vx(1,1) = (0.d0,0.d0)   !<u_i>=0
       vy(1,1) = (0.d0,0.d0)
    
    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !.............. project on to divergence-free part ......................
    
       write(*,"('Main projection')")
    
       do jy = 1,ny

         if (jy == nyhp) cycle
         ky = jy - 1
         if (jy > nyh) ky = ky - ny
         kyky = ky * ky

         do jx = 1,lhm
           kx = dble(jx - 1)
           k22 = kx*kx + kyky 
           if (k22 < 1.d-20) cycle
     
           if (k22 < max_k2) then
             cn1 = vx(jx,jy)
             cn2 = vy(jx,jy)
             t = (cn1 * kx + cn2 * ky )/k22
             vxx = cn1 - t*kx
             vyy = cn2 - t*ky
             mag = sqrt(vxx*conjg(vxx) + vyy*conjg(vyy) )
             d = sqrt(cn1 * conjg(cn1) + cn2 * conjg(cn2) )

             vx(jx,jy) = vxx/mag*d    
             vy(jx,jy) = vyy/mag*d
           else
             vx(jx,jy) = (0.d0,0.d0)
             vy(jx,jy) = (0.d0,0.d0)
           end if
      
         end do
       end do
    
       total_adv = total_adv + tadv
       write(*,"(/'advance = ',ES12.5,'    acumulated = ',ES12.5)") tadv, total_adv   
       write(911,"(/'advance = ',ES12.5,'    acumulated = ',ES12.5)") tadv,total_adv 
    
       if (ad /= nadv) then    
    !.........................................................................
    !................     IFT for the next advancement .......................
          write(*,"('sIFT/',$)") 
          call dfftw_execute(IFTvx)
          call dfftw_execute(IFTvy)
          write(*,"('eIFT',$)") 
       end if
    END DO for_time_scale
    
    
    !...calculate new spectrum for rescaling 
    write(gen,"(I4.4)") nx
    call spectrum(19,'afterspec'//gen//'.dat', .true.)
    ! E(k) is calculated in Subroutine spectrum
    
    write(*,"('Final Rescaling (1/2)')")
    Eloc = 0.d0

    do jy = 1,ny

      if (jy == nyhp) cycle
      ky = jy - 1
      if (jy > nyh) ky = ky - ny
      kyky = ky * ky

      do jx = 1,lhm

        kx = dble(jx - 1)
        k22 = kx*kx + kyky 

        if (k22 < 1.d-20) cycle
        if (k22 < max_k2) then
          wn = sqrt(k22)
          ind = int(wn + 0.5d0)
          if (ind == 1) wn = 1.0d0
          if (ind < lh - 1) then

           ! E(k) is calculated in Subroutine spectrum
           if (E(ind) < 2.d-38) then
              ratio = 0.0d0
           else  
              ratio = sqrt(spec2d(wn)/E(ind))
           end if    
           vx(jx,jy) = vx(jx,jy)*ratio
           vy(jx,jy) = vy(jx,jy)*ratio
      
            !....compute contribution to the new energy spectrum
            Rex =  real(vx(jx,jy))
            Imx = aimag(vx(jx,jy))
            Rey =  real(vy(jx,jy))
            Imy = aimag(vy(jx,jy))
            ss = Rex*Rex + Imx*Imx +  Rey*Rey + Imy*Imy 

            if (kx < 1.d-20) ss = ss * 0.5d0
            if (ss < tiny_real) ss = 0.d0
            Eloc(ind) = Eloc(ind) +  ss
      
          else
            vx(jx,jy) = (0.d0,0.d0)
            vy(jx,jy) = (0.d0,0.d0)
          end if
        else
            vx(jx,jy) = (0.d0,0.d0)
            vy(jx,jy) = (0.d0,0.d0)
        end if
      end do
    end do
    
    
    !...second rescale to match the kinetic energy
    write(*,"('Final energy spectrum rescaling (2/2)')")

    s2 = 0.5d0 * Eloc(1)
    do jx = 2,lh-2
      s2 = s2 + 0.5d0*(Eloc(jx-1) + Eloc(jx))
    end do
    upl = sqrt( s2 )
    resc_factor = urms0/upl

    print *,'resc= ',resc_factor
    
    vx = vx * resc_factor
    vy = vy * resc_factor
    
    !................  final spectrum  ....................
    
    !close(35)
    
    
    call spectrum(7,'spec_final.dat', .true.)
    
    
    open(unit=21, file='./out/ux'//cmd(1:len_trim(cmd))//'.dat', form='unformatted')
      write(21) vx
    close(unit=21)
    open(unit=21, file='./out/uy'//cmd(1:len_trim(cmd))//'.dat', form='unformatted')
      write(21) vy
    close(unit=21)
    
    
    WRITE(*,"(//,' MTLM0256 HAS FINISHED NORMALLY')")


END PROGRAM


!======================================================================
subroutine spectrum(iu,name, output)
    ! calculates the energy spectrum and related Fourier-space quantities
    !======================================================================
    use global
    implicit none
    integer, intent(in) :: iu
    logical, intent(in) :: output
    character(*), intent(in) :: name

    integer jx,jy,ind,ky,kyky,n,jp,jxs,jxe,bl
    real(8) dm,s2,s3,spec2d,kx,k22,upl
    real(8) Rex,Imx,Rey,Imy,ss

    
    
    write(*,"(/,'calculating energy spectrum...',/)")

    E = 0.0d0
    do jy = 1,ny
      if (jy == nyhp) cycle
      ky = jy - 1
      if (jy > nyh) ky = ky - ny
      kyky = ky * ky
      do jx = 1,nth ! nth = nt halved
        kx = dble(jx - 1)
        k22 = kx*kx + kyky 
        if (k22 < 1.d-20) cycle
        ind = int(sqrt(k22) + 0.5d0)
        if (ind < 1) ind = 1
        if (ind < nth) then
          Rex =  real(vx(jx,jy))
          Imx = aimag(vx(jx,jy))
          Rey =  real(vy(jx,jy))
          Imy = aimag(vy(jx,jy))
          ss = Rex*Rex + Imx*Imx+  Rey*Rey + Imy*Imy 
          if (kx < 1.d-20) ss = ss * 0.5d0
          if (ss < tiny_real) ss = 0.d0
          E(ind) = E(ind) +  ss
        end if   
      end do
    end do
    
    s2 = 0.5d0 * E(1)
    do jx=2,lh-2
      s2 = s2 + 0.5d0*(E(jx-1) + E(jx))
    end do
    upl = sqrt(s2)
    resc_factor = urms0/upl   
    
    if (iu > 0) then

      open(unit = iu, file = name, status = 'unknown')

        write(iu, '('' # Statistics of the generated velocity field'')')
        write(iu,"('#      k      E(k) discrete   E(k) analytical')")

        do jx = 1,nth - 1
          write(iu,"(f11.5,2(es16.6))") real(jx), E(jx), spec2d(dble(jx))    
        end do
      
        epsl = 2.d0*alpha*s2 ! Ekman dissipation
        
        write(iu,"(//,'# (Fourier) dissipation    =',ES13.6)") epsl
        write(iu,"('# (Fourier) kinetic energy =',ES13.6)") s2
        write(iu,"('# (Fourier) rms velocity   =',ES13.6)") upl
      close(iu)

    end if

    if(output) then

      write(*,"(//,'# (Fourier) dissipation    =',ES13.6)") epsl
      write(*,"('# (Fourier) kinetic energy =',ES13.6)") s2
      write(*,"('# (Fourier) rms velocity   =',ES13.6)") upl
      ! Commentednby YL
      ! output = .false. 

    end if

END subroutine spectrum



real(8) function spec2d(rk)
!****** 2-D spectrum function ********
  use global
  implicit none

  real(8), intent(in) :: rk

  real(8) :: kl
  integer :: jx, ii

  kl = rk * l
  spec2d = Ck * diss**(2.d0/3.d0) * rk**(-5.d0/3.d0)       &
        * (kl/sqrt(kl*kl + cl))**(5.d0/3.d0 + p0)          

end function spec2d




Subroutine apply_MTLM(vxred,vyred,thetared,ldr,IFTvxred,IFTvyred,   &
                      IFTthred,FTvxred,FTvyred,FTthred)
    !...Apply the MMLM with:
    !        reduced arrays in RAM memory.
    !        main field is disk
    !        3-D FFT for reduced arrays.        
    
    use global
    implicit none

    integer, intent(in) :: ldr
    complex(8), dimension(ldr/2,ldr-2) :: vxred,vyred,thetared
    real(8), dimension(lh-2) :: Eloc,Etloc
    complex(8) t,vxx,vyy,vxa,vya
    integer(8) IFTvxred,IFTvyred,FTvxred,FTvyred,IFTthred,FTthred
    real(8) hxr,hyr,cns,wn,spec2d,ratio,mag,d,k2r,kx,k22,ur,tn,time_scale
    real(8) u1,u2,v1,v2,total_adv,tadv,next_total,Int0kmaxr,urb,s2,thetarb
    real(8) Rex,Imx,Rey,Imy,ss,dm,upl,thetapl,scalar_spec2d
    integer lhr,nxr,nyr,nxrh,nyrh,dny,jxr,jyr,ind,jx,jy,r
    integer ky,kyky,Sup_shell,lhrm,n,jp,bl
    integer i,j,k,ad,nadv,i1,i2
    !logical last_adv
    character(4) suf
    
    lhr = ldr/2
    nxr = ldr - 2
    nyr = (ldr - 2)
    nxrh = nxr/2
    nyrh = nyr/2
    dny = ny - nyr             ! difference in size in y between original array and reduced one
    Sup_shell = (ldr - 2)/2
    lhrm = lhr-1
    
    vxred = (0.d0,0.d0)
    vyred = (0.d0,0.d0)
    
    !......apply cut-off spectral filter 
    ! What is being done in the following four filtering pieces is to pick up the lowpass filtered trunk
    ! and store it into the small array vxred, vyred, vzred. Because of the arrangement of wavenumber
    ! components, it has to be done in four parts. (similar to padd and unpadd.f90)
    
    lhrm = lhr-1
    
    write(*,"('filter 1')")
    
    loop_jyr: do jyr = 1,nyrh
      ky = jyr - 1
      kyky = ky * ky

      do jxr = 1,lhrm
        kx = dble(jxr - 1)
        k22 = kx*kx + kyky 
        if (k22 < 1.d-20) cycle
        ind = int(sqrt(k22) + 0.5d0)
        if (ind < Sup_shell) then
          vxred(jxr,jyr) = vx(jxr,jyr)
          vyred(jxr,jyr) = vy(jxr,jyr)
        end if
      end do
    end do loop_jyr
    
    
    write(*,"('filter 2')")
    
    loop_jyr2: do jyr = nyrh+1,nyr  
      jy = jyr + dny
      if (jy == nyhp) cycle
      ky = jy - 1
      if (jy > nyh) ky = ky - ny
      kyky = ky * ky

      do jxr = 1,lhrm
        kx = dble(jxr - 1)
        k22 = kx*kx + kyky 
        if (k22 < 1.d-20) cycle
        ind = int(sqrt(k22) + 0.5d0)
        if (ind < Sup_shell) then
          vxred(jxr,jyr) = vx(jxr,jy)
          vyred(jxr,jyr) = vy(jxr,jy)
        end if
      end do
    end do loop_jyr2
    
    
    !...compute rms velocity in the band, for rescaling
    s2 = 0.5d0 * spec2d(1.0d0)   
    ! why is the wavenumber 1.0d0? I see: each band of velocity includes all the wavenumber
    ! back to the lowest one. The band is the lowpass filtered velocity field, not the band-pass
    ! filtered one. 
    do jx=2,lhr-1
      s2 = s2 + 0.5d0*(spec2d(dble(jx-1)) + spec2d(dble(jx)))
    end do
    Int0kmaxr = s2
    urb = sqrt(Int0kmaxr/Int0kmax) * urms0
    ! The urb is taken as a portion of the total rms. The portion is calculated from the energy ratio.
   
    
    !...first IFT
    write(*,"('52/',$)") 
    call dfftw_execute(IFTvxred)
    call dfftw_execute(IFTvyred)
    write(*,"('53/',$)") 
    
    
    hxr = two_Pi/nxr
    
    ur = 0.0d0
    do j = 1,nyr
    do i = 1,nxr/2
       i2 = 2*i
       i1 = i2-1
       u1 =  abs(real(vxred(i,j)))
       u2 = abs(aimag(vxred(i,j)))
       
       v1 =  abs(real(vyred(i,j)))
       v2 = abs(aimag(vyred(i,j)))
    
       ur = ur + u1*u1 + v1*v1 + u2*u2 + v2*v2 
    end do
    end do
    ur = sqrt(ur/(2.0d0*nxr*nyr)) ! ur ~ rms velocity of banded velocity field
    tn = hxr/ur ! hxr = grid size, tn = time scale by CFL number
    time_scale = (hxr*hxr/epsilon)**(1.d0/3.d0) ! smallest time scale defined by gridsize and eps
    
    write(*,"(/,'    urms,n = ',ES12.5)") ur
    write(*,"('        tn = ',ES12.5)") tn
    write(*,"('time_scale = ',ES12.5)") time_scale
    write(911,"(//,60('='),/,'  LEVEL = ',I5)") nxr
    write(911,"(/,'    urms,n = ',ES12.5)") ur
    write(911,"('        tn = ',ES12.5)") tn
    write(911,"('time_scale = ',ES12.5)") time_scale
    
    total_adv = 0.d0
    ! last_adv = .false. 
    
    nadv = nint(time_scale/tn)
    tadv = time_scale/nadv  !!! time-step-size
    CFL = tadv/tn
    if (CFL > 1.05) then !!! limit the CFL number not too large
      nadv = nadv + 1
      tadv = time_scale/nadv
      CFL = tadv/tn
    end if  
    write(*,"('       CFL = ',ES12.5)") CFL
    write(911,"('       CFL = ',ES12.5)") CFL  
    
    DO ad = 1, nadv  !!! loop over number of steps
       next_total = total_adv + tadv  !!! next_total is not used????
    
       ! advance one step 
       call MTLM(vxred,vyred,thetared,ldr,tadv,urb,thetarb)
    
       cns = dble(nxr*nyr)

       vxred(lhr,1:nyr) = (0.d0,0.d0)
       vyred(lhr,1:nyr) = (0.d0,0.d0)
    
       !... back to Fourier space
       write(*,"('54/',$)") 
       call dfftw_execute(FTvxred)
       call dfftw_execute(FTvyred)
       write(*,"('55/',$)") 
    
       vxred = vxred/cns
       vyred = vyred/cns
    
       !....velocity field
       vxred(lhr,1:nyr) = (0.d0,0.d0)      !Oddballs
       vyred(lhr,1:nyr) = (0.d0,0.d0)
       vxred(1:lhr,nyr/2+1) = (0.d0,0.d0)
       vyred(1:lhr,nyr/2+1) = (0.d0,0.d0)
    
       vxred(1,1) = (0.d0,0.d0)   !<u_i>=0
       vyred(1,1) = (0.d0,0.d0)
    
       !....perform projection in Fourier space
       ! YL: it seems not necessary to do it in parts. Can do it over small array and corresponding
       ! wavenumbers, except that Sup_shell need to be changed.
       do jyr = 1,nyrh
         ky = jyr - 1
         kyky = ky * ky
         do jxr = 1,lhr-1
           kx = dble(jxr - 1)
           k22 = kx*kx + kyky 
           if (k22 < 1.d-20) cycle
           ind = int(sqrt(k22) + 0.5d0)
           if (ind < Sup_shell) then
             vxa = vxred(jxr,jyr)
             vya = vyred(jxr,jyr)
             t = ( vxa*kx + vya*ky )/k22
             vxx = vxa - t*kx
             vyy = vya - t*ky
             mag = sqrt(vxx*conjg(vxx) + vyy*conjg(vyy) )
             d = sqrt(vxa * conjg(vxa) + vya * conjg(vya) )/mag
             vxred(jxr,jyr) = vxx*d    
             vyred(jxr,jyr) = vyy*d
           end if
         end do
       end do
    
    
       do jyr = nyrh+1,nyr
         jy = jyr + dny
         if (jy == nyhp) cycle
         ky = jy - 1
         if (jy > nyh) ky = ky - ny
         kyky = ky * ky

         do jxr = 1,lhr-1
           kx = dble(jxr - 1)
           k22 = kx*kx + kyky 
           if (k22 < 1.d-20) cycle
           ind = int(sqrt(k22) + 0.5d0)
           if (ind < Sup_shell) then
              vxa = vxred(jxr,jyr)
              vya = vyred(jxr,jyr)
              t = (vxa*kx + vya*ky )/k22
              vxx = vxa - t*kx
              vyy = vya - t*ky
              mag = sqrt(vxx*conjg(vxx) + vyy*conjg(vyy) )
              d = sqrt(vxa * conjg(vxa) + vya * conjg(vya) )/mag
              vxred(jxr,jyr) = vxx*d    
              vyred(jxr,jyr) = vyy*d
           end if   
         end do
       end do
    
    
       vxred(1,1) = 0.d0   !Mean velocity=0
       vyred(1,1) = 0.d0
       write(*,"('57/',$)") 
 
       total_adv = total_adv + tadv
       write(*,"(/'advance = ',ES12.5,'    acumulated = ',ES12.5)") tadv, total_adv   
       write(911,"(/'advance = ',ES12.5,'    acumulated = ',ES12.5)") tadv,total_adv
       write(911,"(/,30('+-'),/)")
    
       if (ad /= nadv) then  ! YL: skip IFT for the last round of the loop
          !...IFT for the next advancement
          write(*,"('sIFT/',$)") 
          call dfftw_execute(IFTvxred)
          call dfftw_execute(IFTvyred)
          write(*,"('eIFT',$)") 
       end if

    END DO
    
    
    !...put back the modified shell into the whole field
    
    write(*,"('write-back 1')")
    loop_jyr5: do jyr = 1,nyrh
       ky = jyr - 1
       kyky = ky * ky
       do jxr = 1,lhrm
          kx = dble(jxr - 1)
          k22 = kx*kx + kyky 
          if (k22 < 1.d-20) cycle
          ind = int(sqrt(k22) + 0.5d0)
          if (ind < Sup_shell) then                
             !!! This condition here has some impact, which makes
             !!! calculating k22 necessary. Can we find a way to avoid calculating k22?
             vx(jxr,jyr) = vxred(jxr,jyr)
             vy(jxr,jyr) = vyred(jxr,jyr)
          end if
       end do  
    end do loop_jyr5
    
    
    
    write(*,"('write-back 2')")
    loop_jyr6: do jyr = nyrh+1,nyr
       jy = jyr + dny
       if (jy == nyhp) cycle
       ky = jy - 1
       if (jy > nyh) ky = ky - ny
       kyky = ky * ky 
       do jxr = 1,lhrm
          kx = dble(jxr - 1)
          k22 = kx*kx + kyky 
          if (k22 < 1.d-20) cycle
          ind = int(sqrt(k22) + 0.5d0)
          if (ind < Sup_shell) then
            vx(jxr,jy) = vxred(jxr,jyr)
            vy(jxr,jy) = vyred(jxr,jyr)
          end if
       end do
    end do loop_jyr6
    
    
    !...calculate new spectra for rescaling (velocity and scalar fields)
    
    write(suf,"(I4.4)") nxr
    call spectrum(19,'afterspec'//suf//'.dat', .false.) 


    !...re-scale the spectra
    write(*,"('Rescaling the spectra (1/2)')")
    
    Eloc = 0.0d0
    do jy = 1,ny
      if (jy == nyhp) cycle
      ky = jy - 1
      if (jy > nyh) ky = ky - ny
      kyky = ky * ky

      do jx = 1,lhm
        kx = dble(jx - 1)
        k22 = kx*kx + kyky 
        if (k22 < 1.d-20) cycle
 
        if (k22 < max_k2) then
 
          wn = sqrt(k22)
          ind = int(wn + 0.5d0)
          if (ind <= 1) then
            ind = 1
            wn = 1.0d0
          end if  
 
          if (ind < nth) then

            if (E(ind) < 2.d-38) then  !!! E has been updated in subroutine spectrum
               ratio = 0.0d0
            else  
               ratio = sqrt(spec2d(wn)/E(ind))
            end if    
            vx(jx,jy) = vx(jx,jy)*ratio
            vy(jx,jy) = vy(jx,jy)*ratio
    
            !...compute contribution to the new energy spectrum
            Rex =  real(vx(jx,jy))
            Imx = aimag(vx(jx,jy))
            Rey =  real(vy(jx,jy))
            Imy = aimag(vy(jx,jy))
            ss = Rex*Rex + Imx*Imx+  Rey*Rey + Imy*Imy 
            if (kx < 1.d-20) ss = ss * 0.5d0
            if (ss < tiny_real) ss = 0.d0
            Eloc(ind) = Eloc(ind) +  ss
    
          else
            vx(jx,jy) = (0.d0,0.d0)
            vy(jx,jy) = (0.d0,0.d0)
 
          end if
 
        else
            vx(jx,jy) = (0.d0,0.d0)
            vy(jx,jy) = (0.d0,0.d0)
        end if
      end do
    end do
    
    !...second rescale to match the kinetic energy
    write(*,"('Rescaling the energy spectrum (2/2)')")
    
    s2 = 0.5d0 * Eloc(1)
    do jx = 2,lh-2
      s2 = s2 + 0.5d0*(Eloc(jx-1) + Eloc(jx))
    end do
    upl = sqrt(s2*2.d0/2.d0)
    resc_factor = urms0/upl
    print *, 'urms0 = ', urms0
    print *, 'upl = ', upl
    print *,'resc= ',resc_factor
    
    vx = vx * resc_factor
    vy = vy * resc_factor
    
End subroutine apply_MTLM
    
    
    
    
!============================================================================
    
    
SUBROUTINE MTLM(vxred,vyred,thetared,ldr,tadv,urmsb,thetarmsb)
    !...field distortion using 'MTLM' displacement. Velocities are recovered
    !...by gridding on the irregular distorted mesh.
    !...Full RAM mapping.
    use global
    implicit none
    
    integer, intent(in) :: ldr
    complex(8), dimension(ldr/2,ldr-2) :: vxred,vyred,thetared
    real(8), dimension(ldr-2,ldr-2) :: ugrid,vgrid,Swf,thetagrid
    
    integer, dimension(ldr-2) :: is,ie
    integer, dimension(ldr-2) :: js,je
    integer, dimension(ldr-2) :: ks,ke
    
    real(8) dt,maxd,wf,hxr,hyr,hzr,vel1,vel2,f1,f2,f3,R,minwf,summat,tadv,urmsb,ur
    real(8) ua1,va1,wa1,ua2,va2,wa2,ir,jr,kr,dirl2,dirh2,djrl2,djrh2,dkrl2,dkrh2
    integer il,ih,jl,jh,kl,kh,ia,ja,ka,ia1,ia2,i,j,k,nxr,nyr,nzr,i1,i2
    real(8) swfs,Swfu,Swfv,R2,Ra,u1,u2,v1,v2,rfac,thetarmsb
    real(8) thetaa1,thetaa2,rfacth,Swftheta,thetamean,th,thetar
    real(8) dx,dy,dist2,ip,jp,nxr8,nyr8
    integer :: in,jn,isearch,jsearch,ksearch,inh,ini
    integer voids
    
    
    
    write(*,"('90/',$)") 

    nxr = ldr - 2
    nyr = nxr
    hxr = two_Pi/nxr
    hyr = two_Pi/nyr
    
    nxr8 = dble(nxr) 
    nyr8 = dble(nyr) 

    
    !...define sphere for neighbors searching
    R = 1.01d0
    minwf = 1.d0/R

    !...apply 'MMLM' displacements
    dt = tadv/hxr  ! normalized 
    ugrid = 0.d0
    vgrid = 0.d0
    thetagrid = 0.d0
    Swf = 0.d0
    
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    
    do ja = 1,nyr
    do ia = 1,nxr/2
       ia2 = 2*ia
       ia1 = ia2 - 1
    
       ua1 = real(vxred(ia,ja))
       va1 = real(vyred(ia,ja))
    
       ua2 = aimag(vxred(ia,ja))
       va2 = aimag(vyred(ia,ja))
    
    !****(1) 
       ir = dble(ia1) + dt*ua1  !!! dt has been normalized by dx
       jr = dble(ja) + dt*va1
       
       il = floor(ir)
       ih = il + 1
       jl = floor(jr)
       jh = jl + 1
    
       dirl2 = ir - dble(il)
       dirh2 = ir - dble(ih)   
       djrl2 = jr - dble(jl)                   
       djrh2 = jr - dble(jh)
    
       dirl2 = dirl2 * dirl2
       dirh2 = dirh2 * dirh2
       djrl2 = djrl2 * djrl2
       djrh2 = djrh2 * djrh2
    
       if (il > nxr) then
         il = il - nxr
       else if (il < 1) then
         il = il + nxr
       end if  
    
       if (ih > nxr) then
         ih = ih - nxr
       else if (ih < 1) then
         ih = ih + nxr
       end if  
     
       if (jl > nyr) then
         jl = jl - nyr
       else if (jl < 1) then
         jl = jl + nyr
       end if  

       if (jh > nyr) then
         jh = jh - nyr
       else if (jh < 1) then
         jh = jh + nyr
       end if  

   
       !...00
       wf = 1.d0/sqrt(dirl2 + djrl2 )
       if (wf > minwf) then
         ugrid(il,jl) = ugrid(il,jl) + wf*ua1
         vgrid(il,jl) = vgrid(il,jl) + wf*va1
         Swf(il,jl) = Swf(il,jl) + wf
       end if 
       
       !...01
       wf = 1.d0/sqrt(dirl2 + djrh2 )
       if (wf > minwf) then
         ugrid(il,jh) = ugrid(il,jh) + wf*ua1
         vgrid(il,jh) = vgrid(il,jh) + wf*va1
         Swf(il,jh) = Swf(il,jh) + wf
       end if
    
       !...10
       wf = 1.d0/sqrt(dirh2 + djrl2 )
       if (wf > minwf) then
         ugrid(ih,jl) = ugrid(ih,jl) + wf*ua1
         vgrid(ih,jl) = vgrid(ih,jl) + wf*va1
         Swf(ih,jl) = Swf(ih,jl) + wf
       end if
    
        !...11
       wf = 1.d0/sqrt(dirh2 + djrh2 )
       if (wf > minwf) then
         ugrid(ih,jh) = ugrid(ih,jh) + wf*ua1
         vgrid(ih,jh) = vgrid(ih,jh) + wf*va1
         Swf(ih,jh) = Swf(ih,jh) + wf
       end if
    
    !****(2) 
       ir = dble(ia2) + dt*ua2
       jr = dble(ja) + dt*va2
    
       il = floor(ir)
       ih = il + 1
       jl = floor(jr)
       jh = jl + 1
    
       dirl2 = ir - dble(il)
       dirh2 = ir - dble(ih)   
       djrl2 = jr - dble(jl)                   
       djrh2 = jr - dble(jh)
    
       dirl2 = dirl2 * dirl2
       dirh2 = dirh2 * dirh2
       djrl2 = djrl2 * djrl2
       djrh2 = djrh2 * djrh2

    
       if (il > nxr) then
         il = il - nxr
       else if (il < 1) then
         il = il + nxr
       end if  
    
       if (ih > nxr) then
         ih = ih - nxr
       else if (ih < 1) then
         ih = ih + nxr
       end if  
         
       if (jl > nyr) then
         jl = jl - nyr
       else if (jl < 1) then
         jl = jl + nyr
       end if  
    
       if (jh > nyr) then
         jh = jh - nyr
       else if (jh < 1) then
         jh = jh + nyr
       end if  


       !...00
       wf = 1.d0/sqrt(dirl2 + djrl2 )
       if (wf > minwf) then
       ugrid(il,jl) = ugrid(il,jl) + wf*ua2
       vgrid(il,jl) = vgrid(il,jl) + wf*va2
       Swf(il,jl) = Swf(il,jl) + wf
       end if
    
       !...01
       wf = 1.d0/sqrt(dirl2 + djrh2 )
       if (wf > minwf) then
       ugrid(il,jh) = ugrid(il,jh) + wf*ua2
       vgrid(il,jh) = vgrid(il,jh) + wf*va2
       Swf(il,jh) = Swf(il,jh) + wf
       end if

       !...10
       wf = 1.d0/sqrt(dirh2 + djrl2 )
       if (wf > minwf) then
       ugrid(ih,jl) = ugrid(ih,jl) + wf*ua2
       vgrid(ih,jl) = vgrid(ih,jl) + wf*va2
       Swf(ih,jl) = Swf(ih,jl) + wf
       end if
    
       !...11
       wf = 1.d0/sqrt(dirh2 + djrh2 )
       if (wf > minwf) then   
       ugrid(ih,jh) = ugrid(ih,jh) + wf*ua2
       vgrid(ih,jh) = vgrid(ih,jh) + wf*va2
       Swf(ih,jh) = Swf(ih,jh) + wf
       end if
    end do
    end do
    

    voids = 0
    do j = 1,nyr
    do i = 1,nxr
      if (Swf(i,j) /= 0.d0) then
         ugrid(i,j) = ugrid(i,j)/Swf(i,j)
         vgrid(i,j) = vgrid(i,j)/Swf(i,j)
      else 
         voids = voids + 1
      end if
    end do
    end do
  
    !#########  second search (only for isolated nodes)   #########
    ! isolated nodes mean the nodes with no particles within a distance of R.
    
    maxd = max(CFL,1.d0)     !......because hx = hy = hz

    DO while (voids > 0) ! By YL
      
      write(*,"(' % voids : ',ES10.3)") real(voids)/(nxr*nyr) * 100.
      write(911,"(' % voids : ',ES10.3)") real(voids)/(nxr*nyr) * 100.d0
      R = R + 1.0d0 ! the searching domain increasing
      R2 = R*R
      Ra = R + maxd 
      
      voids = 0
      
      !...calculate limits fpr searching in the mesh 
      if (Ra < 0.5d0*nyr) then
        do i = 1,nxr
          is(i) = ceiling(dble(i) - Ra)
          ie(i) = floor(dble(i) + Ra)  !!! interval [is, ie] is in the sphere with i as the center
        end do
        js(1:nyr) = is(1:nyr) 
        je(1:nyr) = ie(1:nyr)  ! make it a square
      else
        is = 1
        ie = nxr
        js = 1
        je = nyr
        print *,'  Searching in the whole mesh...'
      end if
      
      write(*,"(/)")     
       
      do j = 1,nyr
      do i = 1,nxr
        if (Swf(i,j) /= 0.d0) cycle  !!!! Skip the non-void nodes
    
        Swfs = 0.d0
        Swfu = 0.d0
        Swfv = 0.d0
    
        !!! For each void node
        do jsearch = js(j),je(j)
          jn = jsearch
          if (jn < 1) then 
            jn = jn + nyr
          else if (jn > nyr) then
            jn = jn - nyr
          end if   
    
          do isearch = is(i),ie(i)
            in = isearch
            if (in < 1) then 
              in = in + nxr
            else if (in > nxr) then
              in = in - nxr
            end if   
         
            inh = in/2
            if (inh*2 == in) then  !!! when in is even
              ua1 = aimag(vxred(inh,jn))
              va1 = aimag(vyred(inh,jn))
            else
              ini = (in+1)/2
              ua1 = real(vxred(ini,jn))
              va1 = real(vyred(ini,jn))
            end if
             
            ip = dble(in) + ua1 * dt
            jp = dble(jn) + va1 * dt
         
            if (ip < 1.d0) then
            ip = ip + nxr8                 !!! nxr8 = dble(ldr)
            else if (ip > nxr8) then
            ip = ip - nxr8
            end if
         
            if (jp < 1.d0) then
            jp = jp + nyr8 
            else if (jp > nyr8) then
            jp = jp - nyr8
            end if   
         
         
            dx = dabs(ip - dble(i))
            dx = dmin1(dx,nxr8 - dx)
            dy = dabs(jp - dble(j))
            dy = dmin1(dy,nyr8 - dy)
            dist2 = dx*dx + dy*dy 
         
            if (dist2 <= R2) then   
              ! ????TODO: why is it here R2???
              ! The argument perhaps is: the searching area is Ra = R + CFL, which include all the points
              ! which are possible to enter the area R in the step. 
              wf = 1.d0/sqrt(dist2)
              Swfu = Swfu + wf * ua1
              Swfv = Swfv + wf * va1
              Swfs = Swfs + wf
            end if
          end do
        end do

        if (Swfs /= 0.d0) then  
           !!! The searching method is different from the general case, where the contributions from a
           !!! particle to neighboring points are calculated. Now the contributions from all particle
           !!! around a point are summed up first.
           ugrid(i,j) = Swfu/Swfs
           vgrid(i,j) = Swfv/Swfs
        else
           voids = voids + 1
        end if   
      end do
      end do
    
    END DO
    
    
    !...rescaling of the Urms
    ur = 0.0d0
    do j = 1,nyr
    do i = 1,nxr
       u1 =  ugrid(i,j)
       v1 =  vgrid(i,j)
       ur = ur + u1*u1 + v1*v1
    end do
    end do
    ur = sqrt(ur/(2.0d0*nxr*nyr))
    
    rfac = urmsb/ur
    write(911,"('urmsb = ',ES14.6)") urmsb
    write(911,"('   ur = ',ES14.6)") ur
    write(911,"(' rfac = ',ES14.6)") rfac
    

    do j = 1,nyr
    do i = 1,nxr/2
       i2 = 2*i
       i1 = i2 - 1
       vxred(i,j) = cmplx(rfac*ugrid(i1,j),rfac*ugrid(i2,j),kind=8)
       vyred(i,j) = cmplx(rfac*vgrid(i1,j),rfac*vgrid(i2,j),kind=8)
    end do
    end do
    
END SUBROUTINE MTLM



!============================================================================


real function ran2_true(idum)
    implicit none
    integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    real AM,EPS,RNMX
    parameter(IM1 = 2147483563, IM2 = 2147483399, AM = 1./IM1, IMM1 = IM1 - 1, &
              IA1 = 40014, IA2 = 40692, IQ1 = 53668, IQ2 = 52774, IR1 = 12211, &
              IR2 = 3791, NTAB = 32, NDIV = 1+IMM1/NTAB, EPS=1.2e-7, RNMX=1. - EPS) 
    integer idum2,j,k,iv(NTAB),iy
    save iv,iy,idum2
    data idum2/123456789/, iv/NTAB*0/, iy/0/
    
    if (idum <= 0) then
       idum = max(-idum,1)
       idum2 = idum
       do j = NTAB+8,1,-1
          k = idum/IQ1
          idum = IA1*(idum - k*IQ1) - k*IR1
          if (idum < 0) idum = idum + IM1
          if (j <= NTAB) iv(j) = idum
       end do
       iy = iv(1)
    end if
    
    k = idum/IQ1
    idum = IA1 * (idum - k*IQ1) - k*IR1
    if (idum < 0) idum = idum + IM1
    k = idum2/IQ2
    idum2 = IA2 * (idum2 - k*IQ2) - k*IR2
    j = 1 + iy/NDIV
    iy = iv(j) - idum2
    iv(j) = idum
    if (iy < 1) iy = iy + IMM1
    ran2_true = min(AM*iy,RNMX)
END function





real function ran2(idum)
    !  Random generator from Numerical Recipes, p.271
    implicit none
    integer idum,ia,im,iq,ir,ntab,ndiv
    real am,eps,rnmx
    parameter(ia = 16807, im = 2147483647, am = 1./im, iq = 127773, ir = 2836,  &
              ntab = 32, ndiv = 1 + (im - 1)/ntab, eps = 1.2e-7, rnmx = 1. - eps)
    integer j,k
    integer, dimension(ntab), save :: iv=(/(0,j=1,ntab)/)
    integer, save :: iy = 0
    
    
    if (idum <= 0 .or. iy == 0) then
      idum = max(-idum,1)
      do j = ntab + 8, 1, -1 
        k = idum/iq
        idum = ia*(idum - k*iq) - ir*k
        if (idum < 0) idum = idum + im
        if (j <= ntab) iv(j) = idum
      end do
      iy = iv(1)
    end if
    
    k = idum/iq
    idum = ia*(idum - k*iq) - ir*k
    if (idum < 0) idum = idum + im
    j = 1 + iy/ndiv
    iy = iv(j)
    iv(j) = idum
    ran2 = min(am*iy,rnmx)
End function
    
