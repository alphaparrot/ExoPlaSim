!
!     ***********************
!     * Planet Simulator 17 *
!     ***********************

!     ********************
!     * Frank Lunkeit    *       University of Hamburg
!     * Edilbert Kirk    *      Meteorological Institute
!     * KLaus Fraedrich  *
!     * Valerio Lucarini *
!     ********************

!     **************
!     * Name rules *
!     **************

!     i - local integer
!     j - loop index
!     k - integer dummy parameter
!     l - logical
!     n - global integer

!     g - real gridpoint arrays
!     p - real dummy parameter
!     s - real spectral arrays
!     y - character
!     z - local real

      program plasim_main
      use pumamod

!     based on PUMA (Portable University Model of the Atmosphere)
!     which is based on SGCM (Simple Global Circulation Model)
!     by Ian James and James Dodd (April 1993)
!     Dept of Meteorology, University of Reading, UK.

!     There is a linear drag which can again vary with level,
!     the time scale is entered for each level in days in tfrc(NLEV).
!     the diffusion is del-ndel with a time scale for diffusion of
!     tdiss days on every variable at the truncation wavenumber.
!
!     UPDATE VERSION IDENTIFIER AFTER EACH CODE CHANGE!

plasimversion = "https://github.com/Edilbert/PLASIM/ : 15-Dec-2015"

      call mpstart(-1)       ! -1: Start MPI   >=0 arg = MPI_COMM_WORLD
      call setfilenames
      call opendiag
      if (mrnum == 2) then
         call mrdimensions
      endif
      call allocate_arrays
      call prolog
      call master
      call epilog
      call mpstop

      stop
      end

      ! ***************************
      ! * SUBROUTINE SETFILENAMES *
      ! ***************************
      
      subroutine setfilenames
      use pumamod
      
      character (3) :: mrext
      
      if (mrpid <  0) return ! no multirun
      
      write(mrext,'("_",i2.2)') mrpid
      
      plasim_namelist     = trim(plasim_namelist    ) // mrext
      radmod_namelist     = trim(radmod_namelist    ) // mrext
      miscmod_namelist    = trim(miscmod_namelist   ) // mrext
      fluxmod_namelist    = trim(fluxmod_namelist   ) // mrext
      rainmod_namelist    = trim(rainmod_namelist   ) // mrext
      surfmod_namelist    = trim(surfmod_namelist   ) // mrext
      plasim_output       = trim(plasim_output      ) // mrext
      plasim_diag         = trim(plasim_diag        ) // mrext
      plasim_restart      = trim(plasim_restart     ) // mrext
      plasim_status       = trim(plasim_status      ) // mrext
      planet_namelist     = trim(planet_namelist    ) // mrext
      efficiency_dat      = trim(efficiency_dat     ) // mrext
      icemod_namelist     = trim(icemod_namelist    ) // mrext
      ice_output          = trim(ice_output         ) // mrext
      oceanmod_namelist   = trim(oceanmod_namelist  ) // mrext
      ocean_output        = trim(ocean_output       ) // mrext
      landmod_namelist    = trim(landmod_namelist   ) // mrext
      vegmod_namelist     = trim(vegmod_namelist    ) // mrext
      seamod_namelist     = trim(seamod_namelist    ) // mrext
      aero_namelist       = trim(aero_namelist      ) // mrext
      
      return
      end

!     ***********************
!     * SUBROUTINE OPENDIAG *
!     ***********************

      subroutine opendiag
      use pumamod
      
      if (mypid == NROOT) then
         open(nud,file=plasim_diag)
      endif
      
      return
      end


!     ******************************
!     * SUBROUTINE ALLOCATE_ARRAYS *
!     ******************************

      subroutine allocate_arrays
      use pumamod
      
      if (mrnum == 2) then
         allocate(sdd(nesp,nlev))   ; sdd(:,:)  = 0.0
         allocate(std(nesp,nlev))   ; std(:,:)  = 0.0
         allocate(szd(nesp,nlev))   ; szd(:,:)  = 0.0
         allocate(spd(nesp     ))   ; spd(:  )  = 0.0
      endif
      
      return
      end subroutine allocate_arrays


!     =================
!     SUBROUTINE PROLOG
!     =================

      subroutine prolog
      use pumamod

      logical :: lrestart

!     ************************************************************
!     * Initializations that cannot be run on parallel processes *
!     ************************************************************

      if (mypid == NROOT) then
         call cpu_time(tmstart)
         write(nud,'(54("*"))')
         write(nud,'("* ",17X,"PLANET SIMULATOR",17X," *")')
         write(nud,'("* ",a50," *")') plasimversion
         write(nud,'(54("*"))')
         if (mrnum > 1) then
            write(nud,'("* Instance ",i3," of ",i3,32x,"*")') &
                  mrpid+1, mrnum
            write(nud,'("* My    truncation  :",i5,27x,"*")') NTRU
            write(nud,'("* Other truncation  :",i5,27x,"*")') mrtru(2-mrpid)
            write(nud,'(54("*"))')
         endif
         write(nud,'("* Truncation   NTRU :",i5,27x,"*")') NTRU
         write(nud,'("* Levels       NLEV :",i5,27x,"*")') NLEV
         write(nud,'("* Latitudes    NLAT :",i5,27x,"*")') NLAT
         write(nud,'("* Longitudes   NLON :",i5,27x,"*")') NLON
         write(nud,'(54("*"))')
         if (NPRO > 1) then
            write(nud,'(54("*"))')
            do jpro = 1 , NPRO
              write(nud,'("* CPU",i4,1x,a42," *")') jpro-1,ympname(jpro)
            enddo
            write(nud,'(54("*"))')
         endif
      endif ! (mypid == NROOT)

      call planet_ini               ! Define planet
      call mpbcr(day_24hr)
      
      if (mypid == NROOT) then
         call print_planet
         call restart_ini(lrestart,plasim_restart)
         if (lrestart) then
            nrestart = 1
            nkits    = 0
            ndivdamp = 0
         endif
         call surface_ini           ! Read boundary and other data
         call inigau(NLAT,sid,gwd)  ! Gaussian abscissas and weights
         call inilat                ! Set latitudinal arrays
         call readnl                ! Open and read <plasim_namelist>
         call initpm                ! Several initializations
         call initsi                ! Initialize semi implicit scheme
         call guistart              ! Initialize GUI
         if (nsela > 0) call tracer_ini0 ! initialize tracer data
         if (nsela > 0 .and. l_aero > 0) call aero_ini 
      endif ! (mypid == NROOT)

!     ***********************
!     * broadcast & scatter *
!     ***********************

      call mpscdn(sid ,NLPP)  ! sine of latitude (kind=8)
      call mpscdn(gwd ,NLPP)  ! gaussian weights (kind=8)
      call mpscrn(csq ,NLPP)  ! cosine squared of latitude
      call mpscrn(rcs ,NLPP)  ! 1.0 / cos(lat)

      do jlat = 1 , NLPP
         deglat(jlat) = 180.0 / PI * asin(sid(jlat))
         cola(jlat) = sqrt(csq(jlat))
         rcsq(1+(jlat-1)*NLON:jlat*NLON) = 1.0 / csq(jlat)
      enddo

      call mpbci(nfixorb ) ! Global switch to fix orbit
      call mpbci(ngenkeplerian) ! Global switch for general Keplerian orbits
      call mpbci(ntpal   ) ! color pallet for T

      call mpbci(kick    ) ! add noise for kick > 0
      call mpbci(naqua   ) ! aqua planet switch
      call mpbci(ndesert ) ! desert planet switch
      call mpbci(nveg    ) ! vegetation switch
      call mpbci(noutput ) ! write data switch
      call mpbci(nafter  ) ! write data interval
      call mpbci(nwpd    ) ! number of writes per day
      call mpbci(nsnapshot) ! Switch for writing snapshots
      call mpbci(nstps   ) ! number of steps per snapshot
      call mpbci(nhcadence) ! Switch for a burst of high-cadence output
      call mpbci(hcstartstep) ! Timestep to start high-cadence output
      call mpbci(hcendstep)   ! Timestep to end high-cadence output (exclusive)
      call mpbci(hcinterval)  ! Number of timesteps per high-cadence write
      call mpbci(ncoeff  ) ! number of modes to print
      call mpbci(ndiag   ) ! write diagnostics interval
      call mpbci(ndivdamp) ! divergence damping countdown
      call mpbci(ngui    ) ! GUI on (1) or off (0)
      call mpbci(sellon )  ! index of longitude for column mode
      call mpbci(nkits   ) ! number of initial timesteps
      call mpbci(nrestart) ! 1: read restart file 0: initial run
      call mpbci(nqspec  ) ! 1: spectral q 0: grodpoint q
      call mpbci(nsela   ) ! 1: semi lagrangian advection enabled
      call mpbci(l_aero  ) ! 1: aerosols enabled

      call mpbci(nstep   ) ! current timestep
      call mpbci(mstep   ) ! current timestep in month
      call mpbci(ntspd   ) ! number of timesteps per day
      call mpbci(mtspd   ) ! number of timesteps per standard day

      call mpbcr(mpstep)   ! minutes per timestep
      call mpbci(n_days_per_month)
      call mpbci(n_days_per_year)
      call mpbci(m_days_per_month)
      call mpbci(m_days_per_year)
      call mpbci(mcal_days_per_year)
      call mpbci(n_steps_per_year)
      call mpbci(n_run_steps)
      call mpbci(n_run_days)
      call mpbci(n_run_months)
      call mpbci(n_run_years)
      call mpbci(n_start_step)
      call mpbci(n_start_year)
      call mpbci(n_start_month)

      call mpbci(nflux   ) !
      call mpbci(nadv    ) !
      call mpbci(nhordif ) !
      call mpbci(nrad    ) !
      call mpbci(neqsig  ) ! switch for equidistant sigma levels
      call mpbci(nhdiff  ) ! critical wavenumber for horizonal diffusion

      call mpbci(nprint    ) ! switch for extensive diagnostic pintout (dbug)
      call mpbci(nprhor    ) ! grid point to be printed (dbug)
      call mpbci(ndiaggp   ) ! switch for franks grid point diagnostics
      call mpbci(ndiagsp   ) ! switch for franks spectral diagnostics
      call mpbci(ndiagcf   ) ! switch for cloud forcing diagnostics
      call mpbci(nentropy  ) ! switch for entropy diagnostics
      call mpbci(nentro3d  ) ! switch for 3d entropy diagnostics
      call mpbci(nenergy  )  ! switch for energy diagnostics
      call mpbci(nener3d  )  ! switch for 3d energy diagnostics
      call mpbci(ndiaggp3d ) ! no of 3d gp diagnostic arrays
      call mpbci(ndiaggp2d ) ! no of 2d gp diagnostic arrays
      call mpbci(ndiagsp3d ) ! no of 3d sp diagnostic arrays
      call mpbci(ndiagsp2d ) ! no of 2d sp diagnostic arrays
      call mpbci(ntime     ) ! switch to activate time consuming estimate
      call mpbci(nperpetual) ! day of perpetual integration
      call mpbci(ndheat)     ! switch for heating due to momentum dissipation
      call mpbci(nsponge)    ! switch for top sponge layer
      call mpbci(nstratosponge) ! Switch for stratospheric newtonian cooling

      call mpbcr(acpd    )   ! Specific heat for dry air
      call mpbcr(adv     )
      call mpbcr(akap    )
      call mpbcr(alr     )
      call mpbcr(cv      )
      call mpbcr(ct      )
      call mpbcr(dtep    )
      call mpbcr(dtns    )
      call mpbcr(dtrop   )
      call mpbcr(dttrp   )
      call mpbcr(ga      )
      call mpbcr(gascon  )
      call mpbcr(tgr     )
      call mpbcr(plarad  )
      call mpbcr(pnu     )
      call mpbcr(pnu21   )
      call mpbcr(psurf   )
      call mpbcr(ptop    )
      call mpbcr(ptop2   )
      call mpbcr(ra1     )
      call mpbcr(ra2     )
      call mpbcr(ra4     )
      call mpbcr(rdbrv   )
      call mpbcr(ww      )
      call mpbcr(solar_day)
      call mpbcr(sidereal_day)
      call mpbcr(tropical_year)
      call mpbcr(sidereal_year)
      call mpbcr(rotspd)
      call mpbcr(eccen)
      call mpbcr(obliq)
      call mpbcr(mvelp)
      call mpbcr(meananom0)
      call mpbcr(plavor)
      call mpbcr(dampsp)
      call mpbcr(taucool)
      
      call mpbcr(fixedlon)
      call mpbcr(dttl)
      

      call mpbcin(ndel  ,NLEV) ! ndel
      call mpbcin(ndl   ,NLEV)

      call mpbcrn(tdissd,NLEV)
      call mpbcrn(tdissz,NLEV)
      call mpbcrn(tdisst,NLEV)
      call mpbcrn(tdissq,NLEV)
      call mpbcrn(damp  ,NLEV)
      call mpbcrn(dsigma,NLEV)
      call mpbcrn(rdsig ,NLEV)
      call mpbcrn(restim,NLEV)
      call mpbcrn(sigma ,NLEV)
      call mpbcrn(sigmah,NLEV)
      call mpbcrn(t0    ,NLEV)
      call mpbcrn(t01s2 ,NLEV)
      call mpbcrn(tfrc  ,NLEV)
      call mpbcrn(tkp   ,NLEV)

      call mpbcrn(c     ,NLSQ)
      call mpbcrn(g     ,NLSQ)
      call mpbcrn(tau   ,NLSQ)

      call mpscin(nindex,NSPP)
      call mpscsp(sak,sakpp,NLEV)
      call mpbcrn(sigh  ,NLEV)
      
      call mpbci(nfilter)
      call mpbci(ngptfilter)
      call mpbci(nspvfilter)
      call mpbcr(landhoskn0)
      call mpbci(nfilterexp)
      call mpbcr(filterkappa)
      

!     Copy some calendar variables to calmod

      call calini(n_days_per_month,n_days_per_year,n_start_step,ntspd &
                  ,solar_day,-1,mcal_days_per_year)
!                  ,day_24hr,-1)

      call mpbci(mcal_days_per_year)
      if (nrestart == 0) nstep = n_start_step ! timestep since 01-01-0001
      call updatim(nstep)  ! set date & time array ndatim
!
!     allocate additional diagnostic arrays, if switched on
!

      if(ndiaggp2d > 0) then
       allocate(dgp2d(NHOR,ndiaggp2d))
       dgp2d(:,:)=0.
      end if
      if(ndiaggp3d > 0) then
       allocate(dgp3d(NHOR,NLEV,ndiaggp3d))
       dgp3d(:,:,:)=0.
      end if
      if(ndiagsp2d > 0) then
       allocate(dsp2d(NESP,ndiagsp2d))
       dsp2d(:,:)=0.
      end if
      if(ndiagsp3d > 0) then
       allocate(dsp3d(NESP,NLEV,ndiagsp3d))
       dsp3d(:,:,:)=0.
      end if
      if(ndiagcf > 0) then
       allocate(dclforc(NHOR,7))
       dclforc(:,:)=0.
      end if
      if(nentropy > 0) then
       allocate(dentropy(NHOR,36))
       allocate(dentrot(NHOR,NLEV))
       allocate(dentroq(NHOR,NLEV))
       allocate(dentrop(NHOR))
       allocate(dentro(NHOR))
       dentropy(:,:)=0.
      end if
      if(nentro3d > 0) then
       allocate(dentro3d(NHOR,NLEV,23))
       dentro3d(:,:,:)=0.
      end if
      if(nenergy > 0) then
       allocate(denergy(NHOR,28))
       denergy(:,:)=0.
      end if
      if(nener3d > 0) then
       allocate(dener3d(NHOR,NLEV,28))
       dener3d(:,:,:)=0.
      end if

      call legini

      if (nrestart > 0) then
         call read_atmos_restart
      else
         call initfd
      endif

      if (mypid == NROOT) then
         if (noutput > 0) call outini    ! Open output file <plasim_output>
      endif
!
!*    initialize miscellaneous additional parameterization
!     which are included in *miscmod*
!

      call miscini

!
!*    initialize surface fluxes and vertical diffusion
!

      call fluxini

!
!*    initialize radiation
!

      call radini

!
!*    initialize carbon weathering
!

      call carbonini
      
!
!*    initialize convective and large scale rain and clouds
!

      call rainini

!
!*    initialize surface parameterizations or models
!

      call surfini
      

      if (mypid==NROOT) then
         write(nud,*) "==========Finalized Albedos=========="
         write(nud,*) "-----For lambda < 0.75 microns------ "
         write(nud,*) "Ground:",dgroundalb(1)
         write(nud,*) "Ocean<;",doceanalb(1) 
         write(nud,*) "Snow:",dsnowalb(1)
         write(nud,*) "Snow max:",dsnowalbmx(1)
         write(nud,*) "Snow min:",dsnowalbmn(1)
         write(nud,*) "Sea ice max:",dicealbmx(1) 
         write(nud,*) "Sea ice min:",dicealbmn(1) 
         write(nud,*) "Glacier min:",dglacalbmn(1)
         write(nud,*) "-----For lambda > 0.75 microns------ "
         write(nud,*) "Ground:",dgroundalb(2)
         write(nud,*) "Ocean<<<;",doceanalb(2) 
         write(nud,*) "Snow:",dsnowalb(2)
         write(nud,*) "Snow max:",dsnowalbmx(2)
         write(nud,*) "Snow min:",dsnowalbmn(2)
         write(nud,*) "Sea ice max:",dicealbmx(2) 
         write(nud,*) "Sea ice min:",dicealbmn(2) 
         write(nud,*) "Glacier min:",dglacalbmn(2)
      endif      
  
!
!*    initialize hurricane/storm diagnostics
!
     
      call hurricaneini(gascon)

!
!*    reset psurf according to orography
!

      if (mypid == NROOT) then
         write(nud,'(" Surface pressure with no topography = ",f10.2," [hPa]")') psurf * 0.01
         zmeanoro = so(1) * cv * cv / sqrt(2.0)
         zdlnp    = zmeanoro / gascon / tgr
         psurf    = EXP(LOG(psurf)-zdlnp)
         write(nud,'(" Mean of topographic height          = ",f10.2," [m]")') zmeanoro / ga
         write(nud,'(" Mean of surface pressure            = ",f10.2," [hPa]")') psurf * 0.01
         write(nud,'(" Initial step number : ",i6," ")') nstep
      end if
      call mpbcr(psurf)

!
!*    broadcast and scatter
!

      call mpbcrn(sp,NESP)
      call mpbcrn(sd,NESP*NLEV)
      call mpbcrn(st,NESP*NLEV)
      call mpbcrn(sz,NESP*NLEV)
      call mpbcrn(sq,NESP*NLEV)

      call mpscsp(sd,sdp,NLEV)
      call mpscsp(st,stp,NLEV)
      call mpscsp(sz,szp,NLEV)
      call mpscsp(sq,sqp,NLEV)
      call mpscsp(sr,srp,NLEV)
      call mpscsp(sp,spp,1)
      call mpscsp(so,sop,1)

!
!*    close the namelist file
!

      if(mypid==NROOT) close(11)

!
!*    open efficiency diagnostic file
!

      if(ndheat > 1 .and. mypid == NROOT) then
       open(9,file=efficiency_dat,form='formatted')
      endif

!
!*    start time consuming calculations
!

      if(ntime==1) then
       call mksecond(zsec,0.)
       time0=zsec
      endif

      return
      end


!     =================
!     SUBROUTINE MASTER
!     =================

      subroutine master
      use pumamod

!     ***************************
!     * short initial timesteps *
!     ***************************

      if (mypid == NROOT) write(nud,*) "Preparing to do initial timesteps: ",nkits

      ikits = nkits
      do jkits=1,ikits
         deltsec  = (day_24hr / mtspd) / (2**nkits) !Timestep
         deltsec2 = deltsec + deltsec
         delt     = (TWOPI     / ntspd) / (2**nkits) !Fraction of rotation
         delt2    = delt + delt
         if (mypid == NROOT) then
            write(nud,*) 'Initial timestep ',jkits,'   deltsec = ',deltsec
            write(nud,*) 'Initial timestep ',jkits,'   delt    = ',delt
         endif
         call gridpointa
         call makebm
         call spectrala
         call gridpointd
         call spectrald
         nkits = nkits - 1
      enddo
      
      if (mypid == NROOT) write(nud,*) "Dynamical core is now spun-up."

!     ****************************************************************
!     * The scaling factor "ww" is derived from the rotation "omega" *
!     * with 1 planetary rotation per sidereal day (2 Pi) of Earth   *
!     ****************************************************************

      deltsec  = day_24hr / mtspd   ! timestep in seconds
      deltsec2 = deltsec + deltsec   ! timestep in seconds * 2
      delt     = TWOPI     / ntspd   ! timestep scaled
      delt2    = delt + delt
      call makebm

      if (mypid == NROOT .and. nsela > 0) then
         call tracer_ini
      endif

!     Use either month countdown (n_run_years * 12 + n_run_months)
!     or step-countdown (for debugging purposes)
      
      mocd = n_run_months                     ! month countdown
      nscd = n_run_days * mtspd + n_run_steps ! step countdown

      if (nrestart == 0) nstep = n_start_step ! timestep since 01-01-0001
      call updatim(nstep)  ! set date & time array ndatim

      nstep1 = nstep ! Remember start step for timing stats
      
      nhcstp = 1
      
!       if (mypid == NROOT) write(nud,'("*  nstep1: ",i6,"   *")') nstep1

      do while (mocd > 0 .or. nscd > 0)  ! main loop

 
!        ************************************************************
!        * calculation of non-linear quantities in grid point space *
!        ************************************************************

         call gridpointa

!        ******************************
!        * adiabatic part of timestep *
!        ******************************

         call spectrala

!        *****************************
!        * diabatic part of timestep *
!        *****************************

         call gridpointd

!          call outaccu
!          naccuout = naccuout + 1
         if (mypid == NROOT) then
            if (mod(nhcstp,nafter) == 0 .and. noutput > 0 ) then
               call outsp
            endif
            if (mod(nhcstp,nstps) == 0 .and. nsnapshot > 0) call snapshotsp
            if (nhcadence>0 .and. hcstartstep>0 .and. hcendstep>0) then
            !We use this triple condition to be really sure that this doesn't get
            !turned on by accident--this has the potential to not only create huge
            !amounts of output if misused, but also to actually damage computing infrastructure.
              if ((nhcstp .ge. hcstartstep) .and. (nhcstp<hcendstep)) then
                if (mod(nhcstp-hcstartstep,hcinterval)==0) call hcadencesp(141)
              endif
            endif
            if (nwritehurricane>0 .and. mod(nhcstp,hcinterval)==0) call hcadencesp(142)
            
            if (mod(nstep,ndiag) == 0 ) then
               call diag
            elseif (ngui > 0) then
               if (mod(nstep,ngui) == 0 ) call diag
            endif
         endif

         if (ngui > 0) call guistep_plasim
         call spectrald
         call outaccu

         if (mod(nhcstp,nstps) == 0 .and. nsnapshot > 0) then
           call snapshotsc
           call snapshotgp
           koutdiag=ndiaggp3d+ndiaggp2d+ndiagsp3d+ndiagsp2d+ndiagcf     &
     &             +nentropy+nenergy
           if(koutdiag > 0) call snapshotdiag
         endif
         if (nhcadence>0 .and. hcstartstep>0 .and. hcendstep>0) then
            !We use this triple condition to be really sure that this doesn't get
            !turned on by accident--this has the potential to not only create huge
            !amounts of output if misused, but also to actually damage computing infrastructure.
           if ((nhcstp .ge. hcstartstep) .and. (nhcstp<hcendstep)) then
             write(nud,*) "HC OUTPUT step",nhcstp
             if (mod(nhcstp-hcstartstep,hcinterval)==0) then
                call hcadencegp(141)
                call hcadencegp(141)
!                 koutdiag=ndiaggp3d+ndiaggp2d+ndiagsp3d+ndiagsp2d+ndiagcf     &
!      &                   +nentropy+nenergy
!                 if(koutdiag > 0) call hcadencediag
             endif
           endif
         endif
         if (nwritehurricane>0 .and. mod(nhcstp,hcinterval)==0) then
            write(nud,*) "HC STORM CAPTURE step",nhcstp
            call hcadencesp(142)
            call hcadencegp(142)
         endif
         if (mod(nhcstp,nafter) == 0) then
          if(noutput > 0) then
!            write(nud,*) "High-cadence step",nhcstp
           call outsc
           call outgp
           koutdiag=ndiaggp3d+ndiaggp2d+ndiagsp3d+ndiagsp2d+ndiagcf     &
     &             +nentropy+nenergy
           if(koutdiag > 0) call outdiag
          endif
          call outreset
!           write(nud,*)"Called outreset"
         endif

         iyea  = ndatim(1)    ! current year
         imon  = ndatim(2)    ! current month
         nstep = nstep + 1
         nhcstp = nhcstp + 1
         mstep = mstep + 1
         call updatim(nstep)  ! set date & time array ndatim
         if (imon /= ndatim(2)) then
            mocd = mocd - 1 ! next month
            if (mypid == NROOT) then
               write(nud,"('Completed month ',I2.2,'-',I4.4)") imon,iyea  
            endif
            mstep = 0
         endif
         call stability_check ! aborts if model tends to explode
         if (nshutdown   > 0) return
         if (nscd  > 0) then
            nscd = nscd - 1
            if (nscd == 0) return
         endif
      enddo ! month countdown

      return
      end


!     =================
!     SUBROUTINE EPILOG
!     =================

      subroutine epilog
      use pumamod
      real    (kind=8) :: zut,zst
      integer (kind=8) :: imem,ipr,ipf,isw,idr,idw
!
!     deallocate additional diagnostic arrays, if switched on
!
      if(ndiaggp2d > 0) deallocate(dgp2d)
      if(ndiagsp2d > 0) deallocate(dsp2d)
      if(ndiaggp3d > 0) deallocate(dgp3d)
      if(ndiagsp3d > 0) deallocate(dsp3d)
      if(ndiagcf   > 0) deallocate(dclforc)
      if(nentropy  > 0) then 
       deallocate(dentropy,dentrop,dentrot,dentroq,dentro)
      endif
      if(nenergy   > 0) deallocate(denergy)
      if(nener3d   > 0) deallocate(dener3d)
      if(nentro3d  > 0) deallocate(dentro3d)
!
!     close output file
!
      if (mypid == NROOT) close(40)
      
!
!     close snapshot file
!
      if (nsnapshot > 0 .and. mypid == NROOT) close(140)
      if (nhcadence > 0 .and. mypid == NROOT) close(141)
!
!     close efficiency diagnostic file
!
      if(ndheat > 1 .and. mypid == NROOT) close(9)
!
!     write restart file
!
      if (mypid == NROOT) then
         call restart_prepare(plasim_status)
         
         if (mod(naccuout,nstpw)==0) naccuout = 0
         
         call put_restart_integer('nstep'   ,nstep   )
         call put_restart_integer('naccuout',naccuout)
         call put_restart_integer('nlat'    ,NLAT    )
         call put_restart_integer('nlon'    ,NLON    )
         call put_restart_integer('nlev'    ,NLEV    )
         call put_restart_integer('nrsp'    ,NRSP    )

!        Save current random number generator seed

         call random_seed(get=meed)
         call put_restart_seed('seed',meed,nseedlen)

         call put_restart_array('sz',sz,NRSP,NESP,NLEV)
         call put_restart_array('sd',sd,NRSP,NESP,NLEV)
         call put_restart_array('st',st,NRSP,NESP,NLEV)
         if (nqspec == 1) call put_restart_array('sq',sq,NRSP,NESP,NLEV)
         call put_restart_array('sr',sr,NRSP,NESP,NLEV)
         call put_restart_array('sp',sp,NRSP,NESP,   1)
         call put_restart_array('so',so,NRSP,NESP,   1)
      endif

      call mpputsp('szm',szm,NSPP,NLEV)
      call mpputsp('sdm',sdm,NSPP,NLEV)
      call mpputsp('stm',stm,NSPP,NLEV)
      if (nqspec == 1) call mpputsp('sqm',sqm,NSPP,NLEV)
      call mpputsp('spm',spm,NSPP,   1)
!
!     gridpoint restart
!
      call mpputgp('dls'    ,dls    ,NHOR,1)
      call mpputgp('drhs'   ,drhs   ,NHOR,1)
      call mpputgp('dalb'   ,dalb   ,NHOR,1)
      call mpputgp('dsalb1',dsalb(1,:),NHOR,1)
      call mpputgp('dsalb2',dsalb(2,:),NHOR,1)
      call mpputgp('dz0'    ,dz0    ,NHOR,1)
      call mpputgp('dicec'  ,dicec  ,NHOR,1)
      call mpputgp('diced'  ,diced  ,NHOR,1)
      call mpputgp('dwatc'  ,dwatc  ,NHOR,1)
      call mpputgp('drunoff',drunoff,NHOR,1)
      call mpputgp('dforest',dforest,NHOR,1)
      call mpputgp('dust3'  ,dust3  ,NHOR,1)
      call mpputgp('dcc'    ,dcc    ,NHOR,NLEV)
      call mpputgp('dql'    ,dql    ,NHOR,NLEV)
      call mpputgp('dqsat'  ,dqsat  ,NHOR,NLEV)
      call mpputgp('dt'     ,dt(1,NLEP),NHOR,1)
      if (nqspec == 1) then ! spectral: save only soil humidity
         call mpputgp('dq'  ,dq(1,NLEP),NHOR,1)
      else                  ! semi-langrange: save complete humidity array
         call mpputgp('dq'  ,dq,NHOR,NLEP)
      if (l_aero > 0) then
         call mpputgp('mmr' ,mmr,NHOR,NLEP)
         call mpputgp('nrho',nrho,NHOR,NLEP)
      endif
      endif
!
!     accumulated diagnostics
!
      call mpputgp('aprl'  ,aprl  ,NHOR,1)
      call mpputgp('aprc'  ,aprc  ,NHOR,1)
      call mpputgp('aprs'  ,aprs  ,NHOR,1)
      call mpputgp('aevap' ,aevap ,NHOR,1)
      call mpputgp('ashfl' ,ashfl ,NHOR,1)
      call mpputgp('alhfl' ,alhfl ,NHOR,1)
      call mpputgp('aroff' ,aroff ,NHOR,1)
      call mpputgp('asmelt',asmelt,NHOR,1)
      call mpputgp('asndch',asndch,NHOR,1)
      call mpputgp('acc'   ,acc   ,NHOR,1)
      call mpputgp('assol' ,assol ,NHOR,1)
      call mpputgp('asthr' ,asthr ,NHOR,1)
      call mpputgp('atsol' ,atsol ,NHOR,1)
      call mpputgp('atthr' ,atthr ,NHOR,1)
      call mpputgp('ataux' ,ataux ,NHOR,1)
      call mpputgp('atauy' ,atauy ,NHOR,1)
      call mpputgp('atsolu',atsolu,NHOR,1)
      call mpputgp('assolu',assolu,NHOR,1)
      call mpputgp('asthru',asthru,NHOR,1)
      call mpputgp('aqvi'  ,aqvi  ,NHOR,1)
      call mpputgp('atsa'  ,atsa  ,NHOR,1)
      call mpputgp('ats0'  ,ats0  ,NHOR,1)
      call mpputgp('atsama',atsama,NHOR,1)
      call mpputgp('atsami',atsami,NHOR,1)
      
      call mpputgp('azmuz'       ,azmuz   ,NHOR,1) 
                                           
      call mpputgp('asigrain'    ,asigrain,NHOR,1)
      call mpputgp('tempmax'     ,tempmax ,NHOR,1) 
      call mpputgp('tempmin'     ,tempmin ,NHOR,1) 
                                           
      call mpputgp('aaso'        ,aaso    ,NESP,1)        
      call mpputgp('aasp'        ,aasp    ,NESP,1)        
      call mpputgp('aast'        ,aast    ,NESP,NLEV)   
      call mpputgp('aasqout'     ,aasqout ,NESP,NLEV)
      call mpputgp('aasd'        ,aasd    ,NESP,NLEV)   
      call mpputgp('aasz'        ,aasz    ,NESP,NLEV)   
      call mpputgp('aadq'        ,aadq    ,NHOR,NLEP) 
      call mpputgp('aammr'       ,aammr   ,NHOR,NLEP)
      call mpputgp('aanrho'      ,aanrho  ,NHOR,NLEP)
      call mpputgp('aadmld'      ,aadmld  ,NHOR,1)      
      call mpputgp('aadt'        ,aadt    ,NHOR,NLEP)   
      call mpputgp('aadwatc'     ,aadwatc ,NHOR,1)     
      call mpputgp('aadsnow'     ,aadsnow ,NHOR,1)     
      call mpputgp('aadql'       ,aadql   ,NHOR,NLEP)  
      call mpputgp('aadust3'     ,aadust3 ,NHOR,1)     
      call mpputgp('aadcc'       ,aadcc   ,NHOR,NLEP)  
      call mpputgp('aadtd5'      ,aadtd5  ,NHOR,1)      
      call mpputgp('aadls'       ,aadls   ,NHOR,1)       
      call mpputgp('aadz0'       ,aadz0   ,NHOR,1)       
      call mpputgp('aadalb'      ,aadalb  ,NHOR,1)      
      call mpputgp('aadsalb1'    ,aadsalb1,NHOR,1)    
      call mpputgp('aadsalb2'    ,aadsalb2,NHOR,1)    
      call mpputgp('aadtsoil'    ,aadtsoil,NHOR,1)    
      call mpputgp('aadtd2'      ,aadtd2  ,NHOR,1)      
      call mpputgp('aadtd3'      ,aadtd3  ,NHOR,1)      
      call mpputgp('aadtd4'      ,aadtd4  ,NHOR,1)      
      call mpputgp('aadicec'     ,aadicec ,NHOR,1)     
      call mpputgp('aadiced'     ,aadiced ,NHOR,1)     
      call mpputgp('aadforest'   ,aadforest   ,NHOR,1)   
      call mpputgp('aadwmax'     ,aadwmax     ,NHOR,1)     
      call mpputgp('aadglac'     ,aadglac     ,NHOR,1)     
      call mpputgp('aadqo3'      ,aadqo3      ,NHOR,NLEV) 
      call mpputgp('aagroundoro' ,aagroundoro ,NHOR,1) 
      call mpputgp('aaglacieroro',aaglacieroro,NHOR,1)

!
!*    finish graphical user interface
!

      call guistop

!
!*    finish miscellaneous additional parameterizations
!

      call miscstop

!
!*    finish surface fluxes and vertical diffusion
!

      call fluxstop

!
!*    finish radiation
!

      call radstop

!
!*    finish large scale and convective rain and clouds
!

      call rainstop

!
!*    finish surface parameterizations or models
!

      call surfstop

!
!*    finish carbon-silicate weathering parameterizations
!

      call carbonstop
      
      if (mypid == NROOT) then
         call restart_stop
      endif

      
      call hurricanestop
!
!     time consumption
!
      if (mypid == NROOT) then 
!        Get resource stats from function resources in file pumax.c
         ires = nresources(zut,zst,imem,ipr,ipf,isw,idr,idw)
         call cpu_time(tmstop)
         tmrun = tmstop - tmstart
         if (nstep > nstep1) then 
            zspy = tmrun * m_days_per_year * real(mtspd) / (nstep-nstep1)
            zypd = (24.0 * 3600.0 / zspy)                         ! siy / day
            write(nud,'(/,"****************************************")')
            if (zut > 0.0) &
            write(nud,  '("* User   time         : ", f10.3," sec *")') zut
            if (zst > 0.0) &
            write(nud,  '("* System time         : ", f10.3," sec *")') zst
            if (zut + zst > 0.0) tmrun = zut + zst
            write(nud,  '("* Total CPU time      : ", f10.3," sec *")') tmrun
            if (imem > 0) then
               zmem = imem * 0.000001
               if (zmem < 1.0) zmem = 1000.0 * zmem ! Could be KB or MB
            write(nud,  '("* Memory usage        : ", f10.3," MB  *")') zmem
            endif
            if (ipr > 0) &
            write(nud,  '("* Page reclaims       :", i7," pages   *")') ipr
            if (ipf > 0) &
            write(nud,  '("* Page faults         :", i7," pages   *")') ipf
            if (isw > 0) &
            write(nud,  '("* Page swaps          :", i7," pages   *")') isw
            if (idr > 0) &
            write(nud,  '("* Disk read           :", i7," blocks  *")') idr
            if (idw > 0) &
            write(nud,  '("* Disk write          :", i7," blocks  *")') idw
            write(nud,'("****************************************")')
            if (zspy < 600.0) then
               write(nud,'("* Seconds per sim year: ",i6,9x,"*")') nint(zspy)
            else if (zspy < 900000.0) then
               write(nud,'("* Minutes per sim year  ",i6,9x,"*")') nint(zspy/60.0)
            else    
               write(nud,'("* Days per sim year:    ",i6,5x,"*")') nint(zspy/86400.0)
            endif
               write(nud,'("* Sim years per day   :",i7,9x,"*")') nint(zypd)
            write(nud,'("****************************************")')
         endif   
      endif      
      return
      end


!     =============================
!     SUBROUTINE READ_ATMOS_RESTART
!     =============================

      subroutine read_atmos_restart
      use pumamod

!     read scalars and full spectral arrays

      if (mypid == NROOT) then
         call random_seed(size=nseedlen)
         allocate(meed(nseedlen))
         call get_restart_integer('nstep'   ,nstep)
         call get_restart_integer('naccuout',naccuout)
         call get_restart_seed('seed',meed,nseedlen)
         call get_restart_array('sz',sz,NRSP,NESP,NLEV)
         call get_restart_array('sd',sd,NRSP,NESP,NLEV)
         call get_restart_array('st',st,NRSP,NESP,NLEV)
         if (nqspec == 1) call get_restart_array('sq',sq,NRSP,NESP,NLEV)
         call get_restart_array('sr',sr,NRSP,NESP,NLEV)
         call get_restart_array('sp',sp,NRSP,NESP,   1)
         call get_restart_array('so',so,NRSP,NESP,   1)
         call random_seed(put=meed)
      endif

      call mpbci(nstep)     ! broadcast current timestep
      call mpbci(naccuout)  ! broadcast accumulation timestep for diagnostics

!     read and scatter spectral arrays

      call mpgetsp('szm',szm,NSPP,NLEV)
      call mpgetsp('sdm',sdm,NSPP,NLEV)
      call mpgetsp('stm',stm,NSPP,NLEV)
      if (nqspec == 1) call mpgetsp('sqm',sqm,NSPP,NLEV)
      call mpgetsp('spm',spm,NSPP,   1)

!     read and scatter surface grids

      call mpgetgp('dls'    ,dls    ,NHOR,   1)
      call mpgetgp('drhs'   ,drhs   ,NHOR,   1)
      call mpgetgp('dalb'   ,dalb   ,NHOR,   1)
      call mpgetgp('dsalb1' ,dsalb(1,:),NHOR,1)
      call mpgetgp('dsalb2' ,dsalb(2,:),NHOR,1)
      call mpgetgp('dz0'    ,dz0    ,NHOR,   1)
      call mpgetgp('dicec'  ,dicec  ,NHOR,   1)
      call mpgetgp('diced'  ,diced  ,NHOR,   1)
      call mpgetgp('dwatc'  ,dwatc  ,NHOR,   1)
      call mpgetgp('drunoff',drunoff,NHOR,   1)
      call mpgetgp('dforest',dforest,NHOR,   1)
      call mpgetgp('dust3'  ,dust3  ,NHOR,   1)
      call mpgetgp('dcc'    ,dcc    ,NHOR,NLEV)
      call mpgetgp('dql'    ,dql    ,NHOR,NLEV)
      call mpgetgp('dqsat'  ,dqsat  ,NHOR,NLEV)
      call mpgetgp('dt'     ,dt(1,NLEP),NHOR,1)
      if (nqspec == 1) then ! spectral: read only soil humidity
         call mpgetgp('dq',dq(1,NLEP),NHOR,1)
      else                  ! semi-langrange: read complete humidity array
         call mpgetgp('dq',dq,NHOR,NLEP)
      if (l_aero > 0) then
         call mpgetgp('mmr',mmr,NHOR,NLEP)
         call mpgetgp('nrho',nrho,NHOR,NLEP)
      endif
      endif

!     read and scatter accumulated diagnostics

      call mpgetgp('aprl'  ,aprl  ,NHOR,1)
      call mpgetgp('aprc'  ,aprc  ,NHOR,1)
      call mpgetgp('aprs'  ,aprs  ,NHOR,1)
      call mpgetgp('aevap' ,aevap ,NHOR,1)
      call mpgetgp('ashfl' ,ashfl ,NHOR,1)
      call mpgetgp('alhfl' ,alhfl ,NHOR,1)
      call mpgetgp('aroff' ,aroff ,NHOR,1)
      call mpgetgp('asmelt',asmelt,NHOR,1)
      call mpgetgp('asndch',asndch,NHOR,1)
      call mpgetgp('acc'   ,acc   ,NHOR,1)
      call mpgetgp('assol' ,assol ,NHOR,1)
      call mpgetgp('asthr' ,asthr ,NHOR,1)
      call mpgetgp('atsol' ,atsol ,NHOR,1)
      call mpgetgp('atthr' ,atthr ,NHOR,1)
      call mpgetgp('ataux' ,ataux ,NHOR,1)
      call mpgetgp('atauy' ,atauy ,NHOR,1)
      call mpgetgp('atsolu',atsolu,NHOR,1)
      call mpgetgp('assolu',assolu,NHOR,1)
      call mpgetgp('asthru',asthru,NHOR,1)
      call mpgetgp('aqvi'  ,aqvi  ,NHOR,1)
      call mpgetgp('atsa'  ,atsa  ,NHOR,1)
      call mpgetgp('ats0'  ,ats0  ,NHOR,1)
      call mpgetgp('atsama',atsama,NHOR,1)
      call mpgetgp('atsami',atsami,NHOR,1)
       
      call mpgetgp('azmuz'       ,azmuz   ,NHOR,1) 
                                           
      call mpgetgp('asigrain'    ,asigrain,NHOR,1)
      call mpgetgp('tempmax'     ,tempmax ,NHOR,1) 
      call mpgetgp('tempmin'     ,tempmin ,NHOR,1) 
                                           
      call mpgetgp('aaso'        ,aaso    ,NESP,1)        
      call mpgetgp('aasp'        ,aasp    ,NESP,1)        
      call mpgetgp('aast'        ,aast    ,NESP,NLEV)   
      call mpgetgp('aasqout'     ,aasqout ,NESP,NLEV)
      call mpgetgp('aasd'        ,aasd    ,NESP,NLEV)   
      call mpgetgp('aasz'        ,aasz    ,NESP,NLEV)   
      call mpgetgp('aadq'        ,aadq    ,NHOR,NLEP)
      call mpgetgp('aammr'       ,aammr   ,NHOR,NLEP)
      call mpgetgp('aanrho'      ,aanrho  ,NHOR,NLEP)
      call mpgetgp('aadmld'      ,aadmld  ,NHOR,1)      
      call mpgetgp('aadt'        ,aadt    ,NHOR,NLEP)   
      call mpgetgp('aadwatc'     ,aadwatc ,NHOR,1)     
      call mpgetgp('aadsnow'     ,aadsnow ,NHOR,1)     
      call mpgetgp('aadql'       ,aadql   ,NHOR,NLEP)  
      call mpgetgp('aadust3'     ,aadust3 ,NHOR,1)     
      call mpgetgp('aadcc'       ,aadcc   ,NHOR,NLEP)  
      call mpgetgp('aadtd5'      ,aadtd5  ,NHOR,1)      
      call mpgetgp('aadls'       ,aadls   ,NHOR,1)       
      call mpgetgp('aadz0'       ,aadz0   ,NHOR,1)       
      call mpgetgp('aadalb'      ,aadalb  ,NHOR,1)      
      call mpgetgp('aadsalb1'    ,aadsalb1,NHOR,1)    
      call mpgetgp('aadsalb2'    ,aadsalb2,NHOR,1)    
      call mpgetgp('aadtsoil'    ,aadtsoil,NHOR,1)    
      call mpgetgp('aadtd2'      ,aadtd2  ,NHOR,1)      
      call mpgetgp('aadtd3'      ,aadtd3  ,NHOR,1)      
      call mpgetgp('aadtd4'      ,aadtd4  ,NHOR,1)      
      call mpgetgp('aadicec'     ,aadicec ,NHOR,1)     
      call mpgetgp('aadiced'     ,aadiced ,NHOR,1)     
      call mpgetgp('aadforest'   ,aadforest   ,NHOR,1)   
      call mpgetgp('aadwmax'     ,aadwmax     ,NHOR,1)     
      call mpgetgp('aadglac'     ,aadglac     ,NHOR,1)     
      call mpgetgp('aadqo3'      ,aadqo3      ,NHOR,NLEV) 
      call mpgetgp('aagroundoro' ,aagroundoro ,NHOR,1) 
      call mpgetgp('aaglacieroro',aaglacieroro,NHOR,1)
      

      return
      end subroutine read_atmos_restart


!     ==========================
!     SUBROUTINE STABILITY_CHECK
!     ==========================

      subroutine stability_check
      use pumamod

!     Some operating systems ignore floating point exceptions
!     and continue to run the program after an explosion,
!     e.g. after wrong settings for timestep or other conditions.
!     This subroutine checks some variables for valid ranges
!     and aborts the program if it's obviously broken.

      if (mypid == NROOT) then
         if (gd(1,1) > 1000.0 .or. gd(1,1) < -1000.0 .or. &
             gz(1,1) > 1000.0 .or. gz(1,1) < -1000.0 .or. &
             (gt(1,1)*ct) > 1000.0 .or. gt(1,1) < -t0(1)) then
            open(44,file='Abort_Message')
            write(44,*) 'Planet Simulator aborted'
            write(44,*) 'timestep = ',nhcstp
            write(44,*) 'gd(1,1) = ',gd(1,1)
            write(44,*) 'gz(1,1) = ',gz(1,1)
            write(44,*) 'gt(1,1) = ',gt(1,1)
            close(44)
   
            write(nud,*) 'Planet Simulator aborted'
            write(nud,*) 'timestep = ',nhcstp
            write(nud,*) 'gd(1,1) = ',gd(1,1)
            write(nud,*) 'gz(1,1) = ',gz(1,1)
            write(nud,*) 'gt(1,1) = ',gt(1,1)
   
            stop
         endif
      endif
      end subroutine stability_check
            

!     =================
!     SUBROUTINE INITFD
!     =================

      subroutine initfd
      use pumamod

      if (nkits < 1) nkits = 1
!     ============================================================
!     next subroutine call to set inital temperature.
!     model started from rest with stratification given by
!     trs calculated in setzt(ex). a perturbation is added to
!     constant log(sp) by subroutine noise, called from setzt(ex).
!     ============================================================
!
       
      if (mypid == NROOT) then
      
       call setzt
      endif
      call mpscsp(sp,spm,1)
      if (abs(dttl)>0.0) then    !setzt created a stellar-antistellar dipole
        call mpscsp(sr,srm,NLEV)
        if (mypid == NROOT) then
           st(:,:) = sr(:,:)
          stm(:,:) =srm(:,:)
           sz(3,:) = plavor
          szm(3,:) = plavor
        endif
!         call mpscsp(sr,stm,NLEV)
      else
        if (mypid == NROOT) then
           st(1,:) = sr(1,:)
          stm(1,:) = sr(1,:)
           sz(3,:) = plavor
          szm(3,:) = plavor
        endif
      endif
      return
      end

!     =================
!     SUBROUTINE READNL
!     =================

      subroutine readnl
      use pumamod

      namelist /plasim_nl/ &
                     kick    , mpstep  , nadv    , naqua   , ncoeff     &
                   , ndel    , ndheat  , ndiag   , ndiagcf , ndiaggp    &
                   , ndiaggp2d , ndiaggp3d , ndesert                    &
                   , ndiagsp   , ndiagsp2d , ndiagsp3d, dttl            &
                   , ndl     , nentropy, nentro3d, neqsig  , nflux      &
                   , ngui    , nguidbg , nhdiff  , nhordif , nkits      &
                   , noutput , nlowio  , nstpw   , nsnapshot, nstps     &
                   , npackgp , npacksp , nperpetual        , nprhor     &
                   , nprint  , nqspec  , nrad    , nsela   , nsync      &
                   , ntime   , ntspd   , nveg    , nwpd    &
                   , n_start_year , n_start_month, n_run_steps          &
                   , n_run_years , n_run_months  , n_run_days           &
                   , n_days_per_month, n_days_per_year, fixedlon        &
                   , nhcadence, hcstartstep, hcendstep, hcinterval      &
                   , seed    , sellon  , nfilter , ngptfilter, nspvfilter   &
                   , landhoskn0, nfilterexp, filterkappa                &
                   , syncstr , synctime, nrdrag  , frcmod               &
                   , dtep    , dtns    , dtrop   , dttrp                &
                   , tdissd  , tdissz  , tdisst  , tdissq  , tgr        &
                   , psurf   , ptop    , ptop2   , taucool              &
                   , restim  , t0      , tfrc    , nstratosponge        &
                   , sigh    , nenergy , nener3d , nsponge , dampsp     &
                   , l_aero
!
!     preset namelist parameter according to model set up
!
      if (NLEV==10) then
         tfrc(1)      =  20.0 * day_24hr * frcmod !day_24hr
         tfrc(2)      = 100.0 * day_24hr * frcmod !day_24hr
         tfrc(3:NLEV) =   0.0 * day_24hr * frcmod !day_24hr
      endif
!

      if(NTRU==42) then
       nhdiff=16
       ndel(:)=4
       tdissq(:)=0.1  * day_24hr
       tdisst(:)=0.76 * day_24hr
       tdissz(:)=0.3  * day_24hr
       tdissd(:)=0.06 * day_24hr
      endif

!
!     read namelist
!

      open(11,file=plasim_namelist,form='formatted')
      read (11,plasim_nl)
      
      if ((NLEV==20) .and. (nrdrag==1) .and. (neqsig.ne.5)) then
         tfrc(1)      =  20.0 * day_24hr !day_24hr
         tfrc(2)      =  30.0 * day_24hr
         tfrc(3)      = 100.0 * day_24hr !day_24hr
         tfrc(4)      = 100.0 * day_24hr !day_24hr
         tfrc(5:NLEV) =   0.0 * day_24hr !day_24hr
      else if ((NLEV==20) .and. (nrdrag==1) .and. (neqsig==5)) then
         tfrc(1:NLEV-10) = 0.0 * day_24hr
         tfrc(NLEV-9)  =  20.0 * day_24hr
         tfrc(NLEV-8)  = 100.0 * day_24hr
         tfrc(NLEV-7:NLEV) =   0.0 * day_24hr
      endif

      if ((ndesert == 1) .and. (naqua == 1)) then !If both toggled, turn off both
        naqua = 0
        ndesert = 0
      endif


!     aqua planet settings
      if (naqua == 1) then
         nveg = 0 ! switch off vegetation
      endif
      
!     We won't turn off vegetation for the desert planet, figuring that interested parties
!     might appreciate having the vegetation naturally respond to low moisture.
      solar_day = day_24hr
      sidereal_day = solar_day * (n_days_per_year-1) / n_days_per_year

!     set rotation dependent variables

      if (rotspd /= 1.0) then
!          if (n_days_per_year == 365) n_days_per_year = 360
!          solar_day = solar_day / rotspd
!          solar_day = (n_days_per_year-1.0)*sidereal_day/n_days_per_year
         sidereal_day = day_24hr / rotspd
!          n_days_per_year = n_days_per_year * rotspd
!          sidereal_day =(n_days_per_year*day_24hr)/(n_days_per_year+1.0)
         if (n_days_per_year /= 1) then
           solar_day = sidereal_day * n_days_per_year/(n_days_per_year-1)
         else
           solar_day = sidereal_day !In this case the solar day is infinite, so we set it to 1 year
         endif
!          day_24hr = 86400.0 ! WHY DOESN'T THIS WORK???
!        If there's one day per year, then solar_day is zero. But really, it's infinite
      endif
      
      m_days_per_year = nint(n_days_per_year * sidereal_day / day_24hr) !24-hour days per year
      n_days_per_month = m_days_per_year / 12 !24-hour days per month
      
      
!       It would appear that day_24hr and sidereal_day cannot be inconsistent with rotspd.
!       For some reason, keeping day_24hr to 24 hours produces NaNs. Why is this?? As near
!       as I can tell, it's purely used for time unit conversions.

      ww    = TWOPI / sidereal_day ! Omega (scaling)
      acpd  = gascon / akap        ! Specific heat for dry air
      adv   = ACPV / acpd -1.0     ! Often used
      cv    = plarad * ww          ! cv
      ct    = cv * cv / gascon     ! ct
      pnu21 = 1.0 - 2.0 * pnu      ! Time filter 2
      rdbrv = gascon / RV          ! rd / rv
!
!     calendar and time control
!     set simulation length by using the following parameters:
!     "n_run_years" and "n_run_months"

!     mpstep <= 0 and ntspd <= 0 triggers automatic

      if (mpstep <= 0 .and. ntspd <= 0) then
         if (nlat <= 32) then       ! T21
            mpstep = 45
         else if (nlat <= 48) then  ! T31
            mpstep = 36
         else                       ! T42 and more
            mpstep = (30 * 64) / nlat
         endif
      endif

!     Make sure that (mpstep * 60) * ntspd = day_24hr

      if (mpstep > 0) then             ! timestep given in [min]
         mtspd = nint(day_24hr) / nint(mpstep * 60)
         mtspd = mtspd + mod(mtspd,2)  ! make even
      endif
      
      mpstep = day_24hr  / real(mtspd * 60)
      ntspd = nint(solar_day) / nint(mpstep * 60)
      nafter = mtspd
      if (nwpd > 0 .and. nwpd <= mtspd) then
         nafter = mtspd / nwpd
      endif
      
      n_steps_per_year = nint(sidereal_year / (mpstep*60.0))
      
      if (nlowio > 0 .and. nstpw > 0) nafter = nstpw
      
      if (ndiag < 1) ndiag = 10 * mtspd

      if (nstps == 0) nstps = mtspd
      
!
!     for column runs set horizontal diffusion coefficients to 0
!
      if (nhordif==0) then
       restim = 0.0
       tfrc   = 0.0
       tdissd = 0.0
       tdissz = 0.0
       tdisst = 0.0
       tdissq = 0.0
      endif

      if (synctime > 0.0) syncstr = 1.0 / (TWOPI * synctime)

      write(nud,'(/,"****************************************")')
      write(nud,'("* plasim_nl from    <",a16,"> *")') plasim_namelist
      write(nud,'("****************************************")')
      write(nud,plasim_nl)

!     Convert start date to timesteps since 1-Jan-0000

      call calini(n_days_per_month,n_days_per_year,n_start_step,ntspd &
                   ,solar_day,0,mcal_days_per_year)
!                  ,day_24hr,0)
      
      call cal2step(n_start_step,mtspd,n_start_year,n_start_month,1,0,0)

!     Compute simulation time in [months]

      n_run_months = n_run_years * 12 + n_run_months
      iyea = n_run_months / 12
      imon = mod(n_run_months,12)

!     Print some values

      write(nud,'(/,"*************************************")')
      write(nud,'("* Solar    day      :",f10.1," [s] *")') solar_day
      write(nud,'("* Sidereal day      :",f10.1," [s] *")') sidereal_day
      write(nud,'("* Omega             :",f8.2," [s-6] *")') ww * 1.0e6
      write(nud,'("* Rotation Speed    :",f10.8,"     *")') rotspd
      write(nud,'("* Days / Year       :",i8,"       *")') n_days_per_year
      write(nud,'("* Days / Month      :",i8,"       *")') n_days_per_month
      write(nud,'("* Timestep          :",f8.3," [min] *")') mpstep
      write(nud,'("* Timesteps / write :",i8,"       *")') nafter
      write(nud,'("* Timesteps / day   :",i8,"       *")') ntspd
      if (iyea  > 1 .and. imon == 0) then
         write(nud,'("* Simulation time:  ",i7,"  years *")') iyea
      else if (iyea == 1 .and. imon == 0) then
         write(nud,'("* Simulation time:     one year     *")')
      else if (n_run_months > 1) then
         write(nud,'("* Simulation time:  ",i7," months *")') n_run_months
      else if (n_run_months == 1) then
         write(nud,'("* Simulation time:    one month     *")')
      else if (n_run_days  > 1) then
         write(nud,'("* Simulation time:  ",i7,"   days *")') n_run_days
      else if (n_run_days == 1) then
         write(nud,'("* Simulation time:      one day     *")')
      else if (n_run_steps  > 1) then
         write(nud,'("* Simulation time:  ",i7,"  steps *")') n_run_steps
      else if (n_run_steps == 1) then
         write(nud,'("* Simulation time:   single step    *")')
      endif
      write(nud,'("*************************************")')

!     set sponge layer time scale

      if(dampsp > 0.) then
       if(dampsp < (day_24hr/mtspd)) dampsp=dampsp*day_24hr
       dampsp=day_24hr/(TWOPI*dampsp)
      endif

!     set franks diagnostics

      if(ndiaggp==1) then
       ndiaggp3d=21+ndiaggp3d
      end if
      if(ndiagsp==1) then
       ndiagsp3d=3+ndiagsp3d
      end if

      return
      end

      subroutine dayseccheck(pf,yn)
      use pumamod
      real :: pf(NLEV)
      character (len=*) :: yn

      zmax = maxval(pf(:))
      if (zmax < (day_24hr / mtspd) .and. zmax > 0.0) then
         write(nud,*) 'old maxval(',trim(yn),') = ',zmax
         write(nud,*) 'assuming [days] - converting to [sec]'
         pf(:) = pf(:) * day_24hr
         write(nud,*) 'new maxval(',trim(yn),') = ',maxval(pf(:))
      endif   
      return
      end 
         
!     =================
!     SUBROUTINE INITPM
!     =================

      subroutine initpm
      use pumamod

      real (kind=8) radea,zakk

!     *************************************************************
!     * carries out all initialisation of model prior to running. *
!     * major sections identified with comments.                  *
!     * this s/r sets the model parameters and all resolution     *
!     * dependent quantities.                                     *
!     *************************************************************

      radea = plarad

!     *********************
!     * set vertical grid *
!     *********************

      if(neqsig==-1) then
       sigmah(:)=sigh(:)
       
      elseif(neqsig==1) then
       do jlev = 1 , NLEV
        sigmah(jlev) = real(jlev) / NLEV
       enddo
       
      elseif(neqsig==2) then
      
       zsk0 = log10(ptop/psurf) + 1
       zskf = 0.99
       dzsk = (zskf - zsk0) / REAL(NLEV-1)
       do jlev = 2 , NLEV
        zsk = zsk0 + (jlev-1)*dzsk - 1
        sigma(jlev) = 10**zsk
       enddo
       sigma(1) = 0.5*sigma(2)
       sigmah(NLEV) = 1.0
       sigmah(1:NLEV-1) = 0.5*(sigma(1:NLEV-1)+sigma(2:NLEV))
       
      elseif(neqsig==3) then
      
       zsk0 = log10(ptop/psurf) + 1
       zskf = 0.99
       dzsk = (zskf - zsk0) / REAL(NLEV-1)
       do jlev = 2 , NLEV
        zfk = real(jlev-1)/real(NLEV-1)
        zffk = (real(jlev-1)/real(NLEV))**0.25
        zsk = zsk0 + (jlev-1)*dzsk - 1
        sigma(jlev) = (1.0-exp(-((1-zfk)**2)/0.05))*10**zsk + exp(-((1-zfk)**2)/0.05)*zffk
       enddo
       sigma(1) = 0.5*sigma(2)
       sigmah(NLEV) = 1.0
       sigmah(1:NLEV-1) = 0.5*(sigma(1:NLEV-1)+sigma(2:NLEV))
       
      elseif(neqsig==4) then
      
       do jlev=1,NLEV
        zsk=REAL(jlev)/REAL(NLEV)
        sigmah(jlev)=0.75*zsk+1.75*zsk**3-1.5*zsk**4
       enddo
       zsk = ptop/psurf
       sigmah(:) = sigmah(:) - sigmah(1)
       sigmah(:) = sigmah(:)/sigmah(NLEV)
       sigmah(:) = sigmah(:)*(1.0-zsk) + zsk
       
      elseif ((neqsig==5) .and. (NLEV .gt. 10)) then
       
       !Bottom atmosphere (PlaSim's normal domain)
       do jlev=NLEV-10,NLEV
         zsk=REAL(jlev-(NLEV-9+1))/10.0  !As if it was a 10-layer atmosphere
         sigmah(jlev)=0.75*zsk+1.75*zsk**3-1.5*zsk**4
       enddo
       zsk = ptop/psurf
       sigmah(NLEV-9:) = sigmah(NLEV-9:) - sigmah(NLEV-9)
       sigmah(NLEV-9:) = sigmah(NLEV-9:)/sigmah(NLEV)
       sigmah(NLEV-9:) = sigmah(NLEV-9:)*(1.0-zsk) + zsk
       
       !Upper atmosphere (stratosphere)
       zsk0 = log10(ptop2/psurf) + 1 !Top at 1 hPa
       zskf = sigmah(max(NLEV-10,1)) !Have to add this max statement so it'll compile with debug flag
       dzsk = (zskf - zsk0) / REAL(NLEV-8)
       do jlev = 2 , NLEV-8
          zsk = zsk0 + (jlev-1)*dzsk - 1
          sigma(jlev) = 10**zsk
       enddo
       sigma(NLEV-8) = 0.5*(sigmah(NLEV-9)+sigmah(NLEV-8))
       sigma(1) = 0.5*sigma(2)
       sigmah(1:NLEV-9) = 0.5*(sigma(1:NLEV-9)+sigma(2:NLEV-8))
         
      else
       do jlev=1,NLEV
        zsk=REAL(jlev)/REAL(NLEV)
        sigmah(jlev)=0.75*zsk+1.75*zsk**3-1.5*zsk**4
       enddo
       
      end if

      dsigma(1     ) = sigmah(1)
      dsigma(2:NLEV) = sigmah(2:NLEV) - sigmah(1:NLEV-1)

      rdsig = 0.5 / dsigma
 
      sigma(1     ) = 0.5 * sigmah(1)
      sigma(2:NLEV) = 0.5 * (sigmah(1:NLEV-1) + sigmah(2:NLEV))
        
!     dimensionless coefficient for newtonian cooling
!     friction and timestep. of course a day is 2*PI in non dimensional
!     units using omega as the unit of frquency.
!     
!     dayseccheck assumes units [days] if values < timestep
!     and converts values to [sec] (compatibilty routine)

      call dayseccheck(restim,"restim")
      call dayseccheck(tfrc  ,"tfrc"  )
      call dayseccheck(tdissd,"tdissd")
      call dayseccheck(tdissz,"tdissz")
      call dayseccheck(tdisst,"tdisst")
      call dayseccheck(tdissq,"tdissq")

      where (restim > 0.0)
         damp = day_24hr / (TWOPI * restim)
      elsewhere
         damp = 0.0
      endwhere
         
      where (tfrc > 0.0)
          tfrc = day_24hr / (TWOPI * tfrc)
      elsewhere
          tfrc = 0.0
      endwhere

!     compute internal diffusion parameter (LAUERSON)

      do jlev=1,NLEV
       jdel = ndel(jlev)
       if (tdissd(jlev) > 0.0) then
        tdissd(jlev) = day_24hr/(TWOPI*tdissd(jlev))
       else
        tdissd(jlev)=0.
       endif
       if (tdissz(jlev) > 0.0) then
        tdissz(jlev) = day_24hr/(TWOPI*tdissz(jlev))
       else
        tdissz(jlev)=0.
       endif
       if (tdisst(jlev) > 0.0) then
        tdisst(jlev) = day_24hr/(TWOPI*tdisst(jlev))
       else
        tdisst(jlev) = 0.
       endif
       if (tdissq(jlev) > 0.0) then
        tdissq(jlev) = day_24hr/(TWOPI*tdissq(jlev))
       else
        tdissq(jlev)=0.
       endif
       zakk=1./(real(NTRU-nhdiff)**jdel)
       jr=-1
       do jm=0,NTRU
         do jn=jm,NTRU
            jr=jr+2
            ji=jr+1
            zsq = (jn - nhdiff)
            if(jn >= nhdiff) then
             sak(jr,jlev) = zakk*zsq**jdel
            else
             sak(jr,jlev) = 0.
            endif
            sak(ji,jlev) = sak(jr,jlev)
         enddo
       enddo
      enddo

!     set coefficients which depend on wavenumber

      zrsq2 = 1.0 / sqrt(2.0)

      jr=-1
      do jm=0,NTRU
         do jn=jm,NTRU
            jr=jr+2
            ji=jr+1
            nindex(jr)=jn
            nindex(ji)=jn
            spnorm(jr)=zrsq2
            spnorm(ji)=zrsq2
         enddo
         zrsq2=-zrsq2
      enddo

! finally make temperatures dimensionless

      dtns  = dtns  / ct
      dtep  = dtep  / ct
      dttrp = dttrp / ct
      t0    = t0    / ct ! (NLEV)

!     print out

      zakk=tdisst(NLEV)/(real(NTRU-nhdiff)**ndel(NLEV))
      write(nud,'(/," *************************************************")')
      if (zakk == 0.0) then
      write (nud,'(" * No lateral dissipation *")')
      else
      write(nud,'(" * Lateral dissipation",13x,"NDEL(",i3,") =",i2," *")')&
         NLEV,ndel(NLEV)
      write(nud,'(" * Diffusion coefficient = ",g14.4," [m**",i1,"] *")')&
         zakk*ww*radea**ndel(NLEV),ndel(NLEV)
      write(nud,'(" * e-folding time for smallest scale =",f5.1," days *")')&
         1.0/(TWOPI*tdisst(NLEV))
      endif
      write(nud,'(" *************************************************")')
      return
      end


!     =================
!     SUBROUTINE MAKEBM
!     =================

      subroutine makebm
      use pumamod

      zdeltsq = delt * delt

      do jlev1 = 1 , NLEV
         do jlev2 = 1 , NLEV
            zaq = zdeltsq * (t0(jlev1) * dsigma(jlev2)                  &
     &          + dot_product(g(:,jlev1),tau(jlev2,:)))
            bm1(jlev2,jlev1,:) = zaq
         enddo
      enddo

      do jn=1,NTRU
         do jlev = 1 , NLEV
            bm1(jlev,jlev,jn) = bm1(jlev,jlev,jn) + 1.0 / (jn*(jn+1))
         enddo
         call minvers(bm1(1,1,jn),NLEV)
      enddo
      return
      end

!     =================
!     SUBROUTINE INITSI
!     =================

      subroutine initsi
      use pumamod
!===========================================================
! carries out all initialisation of model prior to running.
! major sections identified with comments.
! this s/r sets the variables and arrays associated with
! the semi-implicit scheme.
!===========================================================
!
      dimension zalp(NLEV),zh(NLEV)
!
! this value, used in setting alpha(1), is irrelevant in the
! angular momentum conserving ecmwf scheme

      tkp = akap * t0
      t01s2(1:NLEM) = t0(2:NLEV) - t0(1:NLEM)
      t01s2(  NLEV) = 0.0

      zalp(2:NLEV) = log(sigmah(2:NLEV)) - log(sigmah(1:NLEM))

      g      = 0.0
      g(1,1) = 1.0
      do jlev = 2 , NLEV
         g(jlev,jlev) = 1.0 - zalp(jlev)*sigmah(jlev-1)/dsigma(jlev)
         g(jlev,1:jlev-1) = zalp(jlev)
      enddo

      do jlev = 1 , NLEV
         c(jlev,:) = g(:,jlev) * (dsigma(jlev) / dsigma(:))
      enddo

      zt01s2   = t01s2(1)
      zsig     = sigmah(1)
      tau(1,1) = 0.5 * zt01s2 * (zsig - 1.0) + tkp(1)
      tau(2:NLEV,1) = 0.5 * zt01s2 * dsigma(2:NLEV)

      do 1410 jlev=2,NLEV
        zttm=zt01s2
        zsigm=zsig
        zt01s2=t01s2(jlev)
        zsig=sigmah(jlev)
        do 1420 jlev2=1,NLEV
          ztm=0.
          ztmm=0.
          if(jlev2.le.jlev) ztm=1
          if(jlev2.lt.jlev) ztmm=1
          ztau=zttm*(zsigm-ztmm)
          if(jlev.lt.NLEV) ztau=ztau+zt01s2*(zsig-ztm)
          ztau=ztau*rdsig(jlev)*dsigma(jlev2)
          if(jlev2.le.jlev) ztau=ztau+tkp(jlev)*c(jlev2,jlev)
          tau(jlev2,jlev)=ztau
 1420   continue
 1410 continue
!
      zfctr=0.001*CV*CV/ga
      do 1500 jlev = 1 , NLEV
         zh(jlev) = dot_product(g(:,jlev),t0) * zfctr
 1500 continue
!
!     **********************************
!     * write out vertical information *
!     **********************************
!
      write(nud,9001)
      write(nud,9003)
      write(nud,9002)
      do jlev = 1 , NLEV
        write(nud,9004) jlev,sigma(jlev),t0(jlev),zh(jlev)
      enddo
      write(nud,9000)

      if (nprint > 0) then
         write(nud,9012)
         write(nud,9013) (jlev,jlev = 1 , 5)
         write(nud,9012)
         do jlev = 1 , NLEV
           write(nud,9014) jlev,(c(i,jlev),i=1,5)
         enddo
         write(nud,9012)
      endif
      return
 9000 format(1x,33('*'),/)
 9001 format(/,1x,33('*'))
 9002 format(1x,33('*'))
 9003 format(' * Lv *    Sigma Basic-T  Height *')
 9004 format(' *',i3,' * ',3f16.8,' *')
 9012 format(1x,69('*'))
 9013 format(' * Lv * C',i11,4i12,' *')
 9014 format(' *',i3,' * ',5f12.8,' *')
      end

!     ==================
!     SUBROUTINE MINVERS
!     ==================

      subroutine minvers(a,n)
      dimension a(n,n),b(n,n),indx(n)

      b = 0.0
      do j = 1 , n
         b(j,j) = 1.0
      enddo
      call ludcmp(a,n,indx)
      do j = 1 , n
         call lubksb(a,n,indx,b(1,j))
      enddo
      a = b
      return
      end

!     =================
!     SUBROUTINE LUBKSB
!     =================

      subroutine lubksb(a,n,indx,b)
      dimension a(n,n),b(n),indx(n)
      k = 0
      do i = 1 , n
         l    = indx(i)
         sum  = b(l)
         b(l) = b(i)
         if (k > 0) then
            do j = k , i-1
               sum = sum - a(i,j) * b(j)
            enddo
         else if (sum /= 0.0) then
            k = i
         endif
         b(i) = sum
      enddo

      do i = n , 1 , -1
         sum = b(i)
         do j = i+1 , n
            sum = sum - a(i,j) * b(j)
         enddo
         b(i) = sum / a(i,i)
      enddo
      return
      end

!     =================
!     SUBROUTINE LUDCMP
!     =================

      subroutine ludcmp(a,n,indx)
      dimension a(n,n),indx(n),vv(n)

      d = 1.0
      vv = 1.0 / maxval(abs(a),2)

      do 19 j = 1 , n
         do i = 2 , j-1
            a(i,j) = a(i,j) - dot_product(a(i,1:i-1),a(1:i-1,j))
         enddo
         aamax = 0.0
         do i = j , n
            if (j > 1)                                                  &
     &      a(i,j) = a(i,j) - dot_product(a(i,1:j-1),a(1:j-1,j))
            dum = vv(i) * abs(a(i,j))
            if (dum .ge. aamax) then
               imax = i
               aamax = dum
            endif
         enddo
         if (j .ne. imax) then
            do 17 k = 1 , n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
   17       continue
            d = -d
            vv(imax) = vv(j)
         endif
         indx(j) = imax
         if (a(j,j) == 0.0) a(j,j) = tiny(a(j,j))
         if (j < n) a(j+1:n,j) = a(j+1:n,j) / a(j,j)
   19 continue
      return
      end

!     =================
!     SUBROUTINE INILAT
!     =================

      subroutine inilat
      use pumamod
      do jlat = 1 , NLAT
         csq(jlat)  = 1.0 - sid(jlat) * sid(jlat)
         rcs(jlat)  = 1.0 / sqrt(csq(jlat))
      enddo
      do jlat = 1 , NLAT/2
         ideg = nint(180.0/PI * asin(sid(jlat)))
         write(chlat(jlat),'(i2,a1)') ideg,'N'
         write(chlat(NLAT+1-jlat),'(i2,a1)') ideg,'S'
      enddo
      return
      end


!     =====================
!     SUBROUTINE INITRANDOM
!     =====================

      subroutine initrandom
      use pumamod
      integer :: i, clock

!     Set random number generator seed

      call random_seed(size=nseedlen)
      allocate(meed(nseedlen))

!     Take seed from namelist parameter 'SEED' ?

      if (seed(1) /= 0) then
         meed(:) = 0
         i = nseedlen
         if (i > 8) i = 8
         meed(1:i) = seed(1:i)
      else
         call system_clock(count=clock)
         meed(:) = clock + 37 * (/(i,i=1,nseedlen)/)
      endif
      call random_seed(put=meed)
      return
      end

!     ================
!     SUBROUTINE NOISE
!     ================

      subroutine noise
      use pumamod

!     if kick is set to 1 or 2
!     adds white noise perturbation to ln(surface pressure)
!     balanced initial state at t=0.
!     for kick=2, the white noise pertubation is
!     symmetric to the equator
!     eps sets magnitude of the noise

      itp1 = NTP1 ! Suppress compiler warnings for T1
      if (itp1 <= 2 .and. kick > 2) kick = 0 ! for T1
      zeps=1.e-4
      zscale=zeps/sqrt(2.0)

      write(nud,'(/," *****************************************")')
      if (kick == 1) then
         jsp1=2*NTP1+1
         do jsp=jsp1,NRSP
            call random_number(zrand)
            if (mrpid > 0) zrand = zrand + mrpid * 0.01
            sp(jsp)=sp(jsp)+zscale*(zrand-0.5)
         enddo
         write(nud,'(" *     White noise added (KICK = 1)      *")')
      elseif (kick == 2) then
         jr=2*NTP1-1
         do jm=1,NTRU
            do jn=jm,NTRU
               jr=jr+2
               ji=jr+1
               if (mod(jn+jm,2) == 0) then
                  call random_number(zrand)
                  if (mrpid > 0) zrand = zrand + mrpid * 0.01
                  sp(jr)=sp(jr)+zscale*(zrand-0.5)
                  sp(ji)=sp(ji)+zscale*(zrand-0.5)
               endif
            enddo
         enddo
         write(nud,'(" * Symmetric white noise added (KICK=2) *")')
      elseif (kick == 3) then
         sp(2*itp1+3) = zscale
         sp(2*itp1+4) = zscale * 0.5
         write(nud,'(" *  Mode sp(1,1) disturbed (KICK = 3)    *")')
      endif
      write(nud,'(" *****************************************")')
      return
      end

!     ================
!     SUBROUTINE SETZT
!     ================
      subroutine setzt
      use pumamod
!
      dimension ztrs(NLEV)
      dimension ztdprs(NUGP,NLEV)
      dimension zfac(NLEV)
      dimension zdpfac(NUGP)
      dimension srlev(NESP)
!
!*********************************************************************
!  this s/r sets up restoration temp field.
! the temperature at sigma = 1 is tgr, entered in kelvin.
! a lapse rate of ALR k/m is assumed under the tropopause and zero
! above. the actual profile tends to this away from then
! tropopause, with smooth interpolation depending on dttrp
! at the model tropopause.the height of
! the tropopause is given as dtrop m.
!*********************************************************************
!
      sr(:,:) = 0.0 ! NESP,NLEV
!
      zdttrp=dttrp*ct
!
      zsigprev=1.
      ztprev=tgr
      zzprev=0.
      do 1100 jlev=NLEV,1,-1
        zzp=zzprev+(gascon*ztprev/ga)*log(zsigprev/sigma(jlev))
        ztp=tgr-dtrop*ALR
        ztp=ztp+sqrt((.5*ALR*(zzp-dtrop))**2+zdttrp**2)
        ztp=ztp-.5*ALR*(zzp-dtrop)
        ztpm=.5*(ztprev+ztp)
        zzpp=zzprev+(gascon*ztpm/ga)*log(zsigprev/sigma(jlev))
        ztpp=tgr-dtrop*ALR
        ztpp=ztpp+sqrt((.5*ALR*(zzpp-dtrop))**2+zdttrp**2)
        ztpp=ztpp-.5*ALR*(zzpp-dtrop)
        ztrs(jlev)=ztpp
        zzprev=zzprev                                                   &
     &        +(.5*(ztpp+ztprev)*gascon/ga)*log(zsigprev/sigma(jlev))
        ztprev=ztpp
        zsigprev=sigma(jlev)
1100  continue
!
!     **********************************
!     * write out vertical information *
!     **********************************
!
      if (nprint > 0) then
      write(nud,9001)
      write(nud,9003)
      write(nud,9002)
      endif
!     
      do 1200 jlev = 1 , NLEV
         if (nprint > 0) write(nud,9004) jlev,sigma(jlev),ztrs(jlev)
         ztrs(jlev)=ztrs(jlev)/ct
 1200 continue
!     
      if (nprint > 0) write(nud,9002)
!
!******************************************************************
! loop to set array zfac - this controls temperature gradients as a
! function of sigma in tres. it is a sine wave from one at
! sigma = 1 to zero at stps (sigma at the tropopause) .
!******************************************************************
! first find sigma at dtrop
!
      zttrop=tgr-dtrop*ALR
      ztps=(zttrop/tgr)**(ga/(ALR*gascon))
!
! now the latitudinal variation in tres is set up ( this being in terms
! of a deviation from t0 which is usually constant with height)
!
      zsqrt2=sqrt(2.)
      zsqrt04=sqrt(0.4)
      zsqrt6=sqrt(6.)
      
      zsubst1 = - 1.0/3.463843223261725/ct * dttl * cos(fixedlon*PI/180.0)
      zsubst2 =   1.0/3.463843223261725/ct * dttl * sin(fixedlon*PI/180.0)
      ksas = NTP1+1
      ksas1 = 2*ksas-1
      ksas2 = 2*ksas
      
      
      do 2100 jlev = 1 , NLEV
        zfac(jlev)=sin(0.5*PI*(sigma(jlev)-ztps)/(1.-ztps))
        if (zfac(jlev).lt.0.0) zfac(jlev)=0.0
        sr(1,jlev)=zsqrt2*(ztrs(jlev)-t0(jlev))
        sr(3,jlev)=(1./zsqrt6)*dtns*zfac(jlev)
        sr(5,jlev)=-2./3.*zsqrt04*dtep*zfac(jlev)
        sr(ksas1,jlev) = zsubst1*zfac(jlev)
        sr(ksas2,jlev) = zsubst2*zfac(jlev)
 2100 continue
!
      call initrandom
      call printseed
      call noise
!
      return
 9001 format(/,1x,26('*'))
 9002 format(1x,26('*'))
 9003 format(' * Lv *    Sigma Restor-T *')
 9004 format(' *',i3,' * ',f8.3,f9.3,' *')
      end


!     ====================
!     SUBROUTINE PRINTSEED
!     ====================

      subroutine printseed
      use pumamod
      integer :: i

      write (nud,9020)
      write (nud,9010)
      do i = 1 , nseedlen
         write (nud,9000) i,meed(i)
      enddo
      write (nud,9010)
      write (nud,9020)
      return
 9000 format('* seed(',i1,') = ',i10,' *')
 9010 format('************************')
 9020 format(/)
      end


!     ===============
!     SUBROUTINE DIAG
!     ===============

      subroutine diag
      use pumamod
      if (mod(nstep,ndiag) == 0) then
         if (ncoeff .gt. 0) call prisp
         call xsect
      endif
      call energy
      return
      end

!     ================
!     SUBROUTINE PRISP
!     ================

      subroutine prisp
      use pumamod

      character(len=30) title

      scale = 100.0
      title = 'Vorticity [10-2]'
      do 100 jlev = 1 , NLEV
         if (ndl(jlev).ne.0) call wrspam(sz(1,jlev),jlev,title,scale)
  100 continue

      title = 'Divergence [10-2]'
      do 200 jlev = 1 , NLEV
         if (ndl(jlev).ne.0) call wrspam(sd(1,jlev),jlev,title,scale)
  200 continue

      scale = 1000.0
      title = 'Temperature [10-3]'
      do 300 jlev = 1 , NLEV
         if (ndl(jlev).ne.0) call wrspam(st(1,jlev),jlev,title,scale)
  300 continue

      if (nqspec == 1) then
         scale = 1000.0
         title = 'Specific Humidity [10-3]'
         do jlev = 1 , NLEV
            if (ndl(jlev).ne.0) call wrspam(sq(1,jlev),jlev,title,scale)
         enddo
      endif

      title = 'Pressure [10-3]'
      call wrspam(sp,0,title,scale)

      return
      end

!     ==============
!     FUNCTION RMSSP
!     ==============

      function rmssp(pf)
      use pumamod
      real pf(NESP,NLEV)

      zsum = 0.0
      do jlev = 1 , NLEV
         zsum = zsum + dsigma(jlev)                                     &
     &        * (dot_product(pf(1:NZOM,jlev),pf(1:NZOM,jlev)) * 0.5     &
     &        +  dot_product(pf(NZOM+1:NRSP,jlev),pf(NZOM+1:NRSP,jlev)))
      enddo
      rmssp = zsum
      return
      end

!     =================
!     SUBROUTINE ENERGY
!     =================

      subroutine energy
      use pumamod

      parameter (idim=6) ! Number of scalars for GUI timeseries
      real (kind=4) ziso(idim)

      ziso(1) = umax   ! maximum value of csu (zonal mean cross section)
      ziso(2) = t2mean - TMELT  ! mean of 2m temperature
      ziso(3) = precip * 1.0e9  ! mean precipitation
      ziso(4) = evap   * 1.0e9  ! mean evaporation
      ziso(5) = olr             ! OLR
      ziso(6) = minval(dt(:,NLEP)) ! Minimum of surface temperature
!     ziso(6) = 1000.0 * sum(dq(:,NLEV))/2048.0 ! Mean of surface wetness

      call guiput("SCALAR" // char(0) ,ziso,idim,1,1)

      return
      end

!     ==================
!     SUBROUTINE UPDATIM
!     ==================

      subroutine updatim(kstep)
      use pumamod
       
      if (n_days_per_year == 365) then
         call step2cal(kstep,ntspd,ndatim)
      else
         call step2cal30(kstep,ndatim)
      endif
      return
      end

!     =================
!     SUBROUTINE WRSPAM
!     =================

      subroutine wrspam(ps,klev,title,scale)
      use pumamod
!
      dimension ps(NRSP)
      character(len=30) title
      character(len=18) datch


      call ntodat(nstep,datch)
      write(nud,'(1x)')
      write(nud,20000)
      write(nud,20030) datch,title,klev
      write(nud,20000)
      write(nud,20020) (i,i=0,9)
      write(nud,20000)
      write(nud,20100) (cab(i),i=1,10)
      write(nud,20200) (cab(i),i=NTRU+2,NTRU+10)
      write(nud,20300) (cab(i),i=2*NTRU+2,2*NTRU+9)
      write(nud,20400) (cab(i),i=3*NTRU+1,3*NTRU+7)
      write(nud,20000)
      write(nud,'(1x)')

      return

20000 format(1x,78('*'))
20020 format(' * n * ',10i7,' *')
20030 format(' *   * ',a18,2x,a30,'  Level ',i2,11x,'*')
20100 format(' * 0 *',f8.2,9f7.2,' *')
20200 format(' * 1 *',8x,9f7.2,' *')
20300 format(' * 2 *',15x,8f7.2,' *')
20400 format(' * 3 *',22x,7f7.2,' *')
      contains
      function cab(i)
      cab=real(scale*sqrt(ps(i+i-1)*ps(i+i-1)+ps(i+i)*ps(i+i)))
      end function cab
      end

!     ===============
!     SUBROUTINE WRZS
!     ===============

      subroutine wrzs(zs,title,scale)
      use pumamod
!
      dimension zs(NLAT,NLEV)
      character(len=30) title
      character(len=18) datch

      ip = NLAT / 16
      ia = ip/2
      ib = ia + 7 * ip
      id = NLAT + 1 - ia
      ic = id - 7 * ip

      call ntodat(nstep,datch)
      write(nud,'(1x)')
      write(nud,20000)
      write(nud,20030) datch,title
      write(nud,20000)
      write(nud,20020) (chlat(i),i=ia,ib,ip),(chlat(j),j=ic,id,ip)
      write(nud,20000)
      do 200 jlev = 1 , NLEV
         write(nud,20100) jlev,((int(zs(i,jlev)*scale)),i=ia,ib,ip),      &
     &                       ((int(zs(j,jlev)*scale)),j=ic,id,ip),jlev
  200 continue
      write(nud,20000)
      write(nud,'(1x)')

20000 format(1x,78('*'))
20020 format(' * Lv * ',16(1x,a3),' * Lv *')
20030 format(' *    * ',a18,2x,a30,20x,'*')
20100 format(' * ',i2,' * ',16i4,' * ',i2,' *')
      end

!     ================
!     SUBROUTINE WRORB
!     ================

      subroutine wrorb(zf,title)
      use pumamod
      
      real zf
      character(len=30) title
      character(len=18) datch
      
      call ntodat(nstep,datch)
      write(nud,20030) datch,title,zf
      
20030 format('>>>   * ',a18,2x,a30,' = ',f7.3,10x,'*') 
      end
      
!     ================
!     SUBROUTINE XSECT
!     ================

      subroutine xsect
      use pumamod
      use radmod
      character(len=30) title

      title = 'Fractional Day'
      call wrorb(zcdayf,title)
      title = 'True Anomaly [deg]'
      call wrorb(orbnu*180./PI,title)
      title = 'Solar Declination [deg]'
      call wrorb(zdeclf*180./PI,title)
      title = 'Ecliptic Longitude [deg]'
      call wrorb(lambm*180./PI,title)
      title = 'Right Ascension'
      call wrorb(rasc*180./PI,title)
      title = 'Distance Modulus'
      call wrorb(eccf,title)
      scale = 10.0
      title = 'Zonal Wind [0.1 m/s]'
      call wrzs(csu,title,scale)
      title = 'Meridional Wind [0.1 m/s]'
      call wrzs(csv,title,scale)
      scale = 1.0
      title = 'Temperature [C]'
      call wrzs(cst,title,scale)
      scale = 10000.0
      title = 'specific humidity [0.1g/Kg]'
      call wrzs(csm,title,scale)
      scale = 100.0
      title = 'cloud cover [%]'
      call wrzs(ccc,title,scale)
      return
      end


!     =====================
!     SUBROUTINE GRIDPOINTA
!     =====================

      subroutine gridpointa
      use pumamod
!
!*    Adiabatic Gridpoint Calculations
!
      real gtn(NHOR,NLEV)
      real gut(NHOR,NLEV)
      real gvt(NHOR,NLEV)
      real gqn(NHOR,NLEV)
      real guq(NHOR,NLEV)
      real gvq(NHOR,NLEV)
      real guz(NHOR,NLEV)
      real gvz(NHOR,NLEV)
      real gke(NHOR,NLEV)
      real gphi(NHOR,NLEV)
      real gvpp(NHOR)
      real gpmt(NLON,NLPP)
      real sdf(NESP,NLEV)
      real stf(NESP,NLEV)
      real sqf(NESP,NLEV)
      real szf(NESP,NLEV)
      real spf(NESP)

      real zgp(NLON,NLAT)

      real (kind=4) zcs(NLAT,NLEV)
      real (kind=4) zsp(NESP)

!
!     inverse Legendre transformation (spectral to fourier domain)
!     st -> gt  sq -> gq  sd -> gd  sz -> gz
!     (sd,sz) -> (gu,gv)
!     sp -> gp  sp -> gpj (dlnps/dphi)
!

      call invlega

      if (ngui > 0 .or. mod(nstep,ndiag) == 0) then
        do jlev = 1 , NLEV
          do jlat = 1 , NLPP
            sec = CV / sqrt(csq(jlat))
            csu(jlat,jlev) =  gu(1+(jlat-1)*NLON,jlev) * sec
            csv(jlat,jlev) =  gv(1+(jlat-1)*NLON,jlev) * sec
            cst(jlat,jlev) =(gt(1+(jlat-1)*NLON,jlev) + t0(jlev))*ct-TMELT
            j1=(jlat-1)*NLON+1
            j2=jlat*NLON
            ccc(jlat,jlev) = SUM(dcc(j1:j2,jlev))/real(NLON)
            if (nqspec == 1) then
               csm(jlat,jlev) = (gq(1+(jlat-1)*NLON,jlev))
            else
               csm(jlat,jlev) = sum(dq(j1:j2,jlev))/real(NLON)
            endif
          enddo
        enddo
        umax = maxval(csu)
      endif

      do jlat = 1 , NLPP
         do jlon = 1 , NLON , 2
           gpmt(jlon  ,jlat) = -gp(jlon+1+(jlat-1)*NLON) * ((jlon-1)/2)
           gpmt(jlon+1,jlat) =  gp(jlon  +(jlat-1)*NLON) * ((jlon-1)/2)
         enddo
      enddo

      call fc2gp(gu  ,NLON,NLPP*NLEV)
      call fc2gp(gv  ,NLON,NLPP*NLEV)
      call fc2gp(gt  ,NLON,NLPP*NLEV)
      call fc2gp(gd  ,NLON,NLPP*NLEV)
      call fc2gp(gz  ,NLON,NLPP*NLEV)
      call fc2gp(gpj ,NLON,NLPP)
      call fc2gp(gpmt,NLON,NLPP)
      call fc2gp(gp  ,NLON,NLPP)
      if (nqspec == 1) call fc2gp(gq  ,NLON,NLPP*NLEV)
      gp = exp(gp)

      call calcgp(gtn,gqn,guz,gvz,gpmt,gvpp,gphi)

      gut = gu * gt
      gvt = gv * gt
      gke = gu * gu + gv * gv
      if (nqspec == 1) then
         guq = gu * gq
         gvq = gv * gq
      endif
!
!     add non linear geopotential terms
!
      gke(:,:) = gke(:,:) + gphi(:,:)

!
!     fft
!

      call gp2fc(gtn ,NLON,NLPP*NLEV)
      call gp2fc(gqn ,NLON,NLPP*NLEV)

      call gp2fc(gut ,NLON,NLPP*NLEV)
      call gp2fc(gvt ,NLON,NLPP*NLEV)
      call gp2fc(guz ,NLON,NLPP*NLEV)
      call gp2fc(gvz ,NLON,NLPP*NLEV)
      call gp2fc(gke ,NLON,NLPP*NLEV)
      call gp2fc(gvpp,NLON,NLPP     )
      if (nqspec == 1) then
         call gp2fc(guq ,NLON,NLPP*NLEV)
         call gp2fc(gvq ,NLON,NLPP*NLEV)
      endif
!
!     direct Legendre transformation (fourier domain to spectral domain)
!
      call fc2sp(gvpp,spf)
      call mktend(sdf,stf,szf,gtn,gvz,guz,gke,gut,gvt)
      call mpsumsc(spf,spt,1)
      call mpsumsc(stf,stt,NLEV)
      call mpsumsc(sdf,sdt,NLEV)
      call mpsumsc(szf,szt,NLEV)
      if (nqspec == 1) then
         call qtend(sqf,gqn,guq,gvq)
         call mpsumsc(sqf,sqt,NLEV)
      endif
!
!     compute entropy
!
      if(nentropy > 0) then
       zp00=100000.
       dentrop(:)=psurf*gp(:)
       dentropy(:,1)=0.
       do jlev=1,NLEV
        dentrot(:,jlev)=ct*(gt(:,jlev)+t0(jlev))
        dentroq(:,jlev)=gq(:,jlev)*psurf/dentrop(:)
        dentro(:)=acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)*dsigma(jlev) &
     &    *log(dentrot(:,jlev)*(zp00/(dentrop(:)*sigma(jlev)))**akap)/ga
        dentropy(:,1)=dentropy(:,1)+dentro(:)
        if(nentro3d > 0) dentro3d(:,jlev,1)=dentro(:)
       enddo
      endif
!
!     save u, v, and ps (at time t) for tracer transport
!
       dp0(:)=psurf*gp(:)
       do jlev=1,NLEV
        du0(:,jlev)=cv*gu(:,jlev)*SQRT(rcsq(:))
        dv0(:,jlev)=cv*gv(:,jlev)*SQRT(rcsq(:))
       enddo
!
!     save u v t q and ps (at time t) for further diagnostics
!
      if(nenergy > 0) then
       dp(:)  =dp0(:)
       du(:,:)=du0(:,:)
       dv(:,:)=dv0(:,:)
       do jlev=1,NLEV
        dt(:,jlev)=ct*(gt(:,jlev)+t0(jlev))
        if (nqspec == 1) dq(:,jlev)=gq(:,jlev)*psurf/dp(:)
       enddo
      endif

!     Diagnostic output and GUI calls

      if (ngui>0 .or. mod(nstep,ndiag)==0 .or. mod(nstep,nafter)==0) then
         call mpgagp(zgp,gp,1);
         call guips(zgp)                ! Send Ps for Hovmoeller
         call guihor("DQVI"//char(0),dqvi,1,1.0,0.0)! Vertically integrated q
         call guigv("GU"  // char(0),gu)            ! Send u to GUI
         call guigv("GV"  // char(0),gv)            ! Send v to GUI
         call mpgagp(zgp,dtsa,1)
         if (mypid == NROOT) t2mean = sum(zgp) / (NLON * NLAT) ! Mean of 2m temperature
         call mpgagp(zgp,dprc,1)        ! Convective precip
         if (mypid == NROOT) precip = sum(zgp)
         call mpgagp(zgp,dprl,1)        ! Large scale precip
         if (mypid == NROOT) precip = (precip + sum(zgp)) / (NLON * NLAT)
         call mpgagp(zgp,devap,1)       ! Evaporation
         if (mypid == NROOT) evap = sum(zgp) / (NLON * NLAT)
         call mpgagp(zgp,dftu,1)        ! OLR = dftu level 1
         if (mypid == NROOT) olr = -sum(zgp) / (NLON * NLAT)! make positive upwards
         gp(:) = gp(:) - 1.0
         call gp2fc(gp,NLON,NLPP)
         call fc2sp(gp,span)
         call mpsum(span,1)
 
         call mpgacs(csu)
         call mpgacs(csv)
         call mpgacs(cst)
         call mpgacs(csm)
         call mpgacs(ccc)
 
         if (mypid == NROOT) then
            zcs(:,:) = csu(:,:)
            call guiput("CSU"  // char(0) ,zcs ,NLAT, NLEV,1)
            zcs(:,:) = csv(:,:)
            call guiput("CSV"  // char(0) ,zcs ,NLAT, NLEV,1)
            zcs(:,:) = cst(:,:)
            call guiput("CST"  // char(0) ,zcs ,NLAT, NLEV,1)
            zsp(:) = span(:)
            call guiput("SPAN" // char(0) ,zsp ,NCSP,-NTP1,1)
         endif

         ! send fields to GUI for column visualization
         call guigvcol("GUCOL"  // char(0),gu,sellon) ! u
         call guigvcol("GVCOL"  // char(0),gv,sellon) ! v 
         call guigtcol(dt,sellon) ! t 
         call guid3dcol("DCCCOL" // char(0),dcc,sellon,NLEP,100.0,0.0) !cl-cov 
         call guid3dcol("DQCOL" // char(0),dq,sellon,NLEP,1000.0,0.0)  !dq
         call guid3dcol("DTDTCOL" // char(0),dtdt,sellon,NLEP,         &
                        day_24hr,0.0)            ! t-tendency
         call guid3dcol("DQDTCOL" // char(0),dqdt,sellon,NLEP,         &
                        1000.0*day_24hr,0.0)     ! q-tendency
      endif
      if (nqspec == 0) then
         do jlev = 1 , NLEV
           ! dq(:,jlev) = dq(:,jlev) + gqn(:,jlev) * psurf / dp(:)
         enddo
      endif
      return
      end

!     =================
!     SUBROUTINE CALCGP
!     =================

      subroutine calcgp(gtn,gqn,guz,gvz,gpm,gvp,gphi)

!     *****************************************************
!     * computes nonlinear tendencies in grid point space *
!     *****************************************************

      use pumamod

      real gtn(NHOR,NLEV)
      real gqn(NHOR,NLEV)                                                       !NEU
      real gphi(NHOR,NLEV)
      real guz(NHOR,NLEV), gvz(NHOR,NLEV)
      real gpm(NHOR)     , gvp(NHOR)
      real zsdotp(NHOR,NLEM),zsumd(NHOR)
      real ztpta(NHOR),ztptb(NHOR)
      real zvgpg(NHOR,NLEV)
      real gtd(NHOR,NLEM)
      real gqm(NHOR,NLEM)                                                       !NEU

      real gud(NHOR,NLEM)
      real gvd(NHOR,NLEM)

      real ztv1(NHOR,NLEV),ztv2(NHOR,NLEV),zq(NHOR)

      do jlev = 1 , NLEV
         zvgpg(:,jlev) = rcsq  * (gu(:,jlev) * gpm + gv(:,jlev) * gpj)
!
!     set pseudo temperatures to use virtual temperatures
!     (note: gq=ps*q)
!
         if (nqspec == 1) then
            zq(:)=AMAX1(gq(:,jlev)/gp(:),0.)
         else
            zq(:)=max(dq(:,jlev),0.0)
         endif
         ztv1(:,jlev)=(gt(:,jlev)+t0(jlev))*(1.+(1./rdbrv-1.)*zq(:))    &
     &               -t0(jlev)
         ztv2(:,jlev)=(gt(:,jlev)+t0(jlev))*(1.+(1./rdbrv-1.)*zq(:))    &
     &               /(1.+ADV*zq(:))                                    &
     &               -t0(jlev)

      enddo

!     *******
!     * gvp *
!     *******

      zsumd = dsigma(1) * gd(:,1)
      gvp   = dsigma(1) * zvgpg(:,1)
      zsdotp(:,1) = zsumd + gvp

      do jlev = 2 , NLEV-1
         zsumd = zsumd + dsigma(jlev) * gd(:,jlev)
         gvp   = gvp   + dsigma(jlev) * zvgpg(:,jlev)
         zsdotp(:,jlev) = zsumd + gvp
      enddo

      zsumd = zsumd + dsigma(NLEV) * gd(:,NLEV)
      gvp   = gvp   + dsigma(NLEV) * zvgpg(:,NLEV)

!     **************
!     * loop  400: *
!     **************

      do jlev = 1 , NLEM
         zsdotp(:,jlev) = (sigmah(jlev) * (zsumd+gvp) - zsdotp(:,jlev))
         gtd(:,jlev) = zsdotp(:,jlev) * (gt(:,jlev+1) - gt(:,jlev))
         gqm(:,jlev) = zsdotp(:,jlev) * (gq(:,jlev+1) + gq(:,jlev))
         gud(:,jlev) = zsdotp(:,jlev) * (gu(:,jlev+1) - gu(:,jlev))
         gvd(:,jlev) = zsdotp(:,jlev) * (gv(:,jlev+1) - gv(:,jlev))
      enddo

!     *************
!     * top level *
!     *************

      zsumd = zvgpg(:,1) * dsigma(1)

      gtn(:,1) = gt(:,1) * gd(:,1) - akap * ztv2(:,1) * gd(:,1)         &
     &   - rdsig(1)*(gtd(:,1) + t01s2(1) * (sigmah(1)*gvp-zsumd))

      gqn(:,1) = - rdsig(1) * gqm(:,1)

      guz(:,1) =-gu(:,1) * gz(:,1) - gpj * ztv1(:,1) - rdsig(1)*gvd(:,1)
      gvz(:,1) = gv(:,1) * gz(:,1) - gpm * ztv1(:,1) - rdsig(1)*gud(:,1)

      dw(:,1)=gd(:,1)*gp(:)*psurf*ww

!     ****************
!     * inner levels *
!     ****************

      do jlev = 2 , NLEV-1
         ztpta = c(1,jlev) *  zvgpg(:,1)
         ztptb = c(1,jlev) * (zvgpg(:,1) + gd(:,1))

         do jlej = 2 , jlev
            ztpta = ztpta + c(jlej,jlev) *  zvgpg(:,jlej)
            ztptb = ztptb + c(jlej,jlev) * (zvgpg(:,jlej) + gd(:,jlej))
         enddo

         zsumd = zsumd + zvgpg(:,jlev) * dsigma(jlev)

         gtn(:,jlev) = gt(:,jlev) * gd(:,jlev)                          &
     &       + akap * ztv2(:,jlev) * (zvgpg(:,jlev) - ztptb)            &
     &       + tkp(jlev) * (zvgpg(:,jlev) - ztpta)                      &
     &       - rdsig(jlev) * (gtd(:,jlev) + gtd(:,jlev-1)               &
     &                       +gvp*(t01s2(jlev)*sigmah(jlev)             &
     &                            +t01s2(jlev-1)*sigmah(jlev-1))        &
     &                       -zsumd*(t01s2(jlev-1)+t01s2(jlev))         &
     &                       +zvgpg(:,jlev)*dsigma(jlev)*t01s2(jlev-1))

         gqn(:,jlev) = - rdsig(jlev)*(gqm(:,jlev) - gqm(:,jlev-1))

         guz(:,jlev) = - gu(:,jlev) * gz(:,jlev) - gpj * ztv1(:,jlev)   &
     &       - rdsig(jlev)*(gvd(:,jlev) + gvd(:,jlev-1))

         gvz(:,jlev) =   gv(:,jlev) * gz(:,jlev) - gpm * ztv1(:,jlev)   &
     &       - rdsig(jlev)*(gud(:,jlev) + gud(:,jlev-1))

         dw(:,jlev)=(zvgpg(:,jlev)-ztptb)*gp(:)*psurf*ww

      enddo

!     ****************
!     * bottom level *
!     ****************

      ztpta = c(1,NLEV) *  zvgpg(:,1)
      ztptb = c(1,NLEV) * (zvgpg(:,1) + gd(:,1))

      do jlej = 2 , NLEV
         ztpta = ztpta + c(jlej,NLEV) *  zvgpg(:,jlej)
         ztptb = ztptb + c(jlej,NLEV) * (zvgpg(:,jlej) + gd(:,jlej))
      enddo

      gtn(:,NLEV) = gt(:,NLEV) * gd(:,NLEV)                             &
     &   + akap*ztv2(:,NLEV)*(zvgpg(:,NLEV)-ztptb)                      &
     &   + tkp(NLEV)*(zvgpg(:,NLEV)-ztpta)                              &
     &   - rdsig(NLEV)*(gtd(:,NLEM)                                     &
     &                 +t01s2(NLEV-1)*(sigmah(NLEV-1)*gvp-zsumd))

      gqn(:,NLEV) = rdsig(NLEV) * gqm(:,NLEM)

      guz(:,NLEV) = -gu(:,NLEV)*gz(:,NLEV) - gpj*ztv1(:,NLEV)           &
     &    - rdsig(NLEV) * gvd(:,NLEM)
      gvz(:,NLEV) =  gv(:,NLEV)*gz(:,NLEV) - gpm*ztv1(:,NLEV)           &
     &    - rdsig(NLEV) * gud(:,NLEM)

      dw(:,NLEV)=(zvgpg(:,NLEV)-ztptb)*gp(:)*psurf*ww

!
!     compute non linear geopotential terms
!

      ztv1(:,:)=ztv1(:,:)-gt(:,:)
      do jlev=1,NLEV
      do jhor=1,NHOR
       gphi(jhor,jlev)=dot_product(g(:,jlev),ztv1(jhor,:))              &
     &                *2./rcsq(jhor)
      enddo
      enddo

      return
      end

!     ====================
!     SUBROUTINE SPECTRALA
!     ====================

      subroutine spectrala
      use pumamod
!
!*    Add adiabatic and diabatic tendencies
!
!     the adiabatic tendencies are added using the semi implicit scheme
!     described in
!     Hoskins and Simmons 1975 (Q.J.R.Meteorol.Soc.,101,637-655) (HS75)
!     To compare the code directly with HS75 the following notes might be
!     helpful (in addition to the comments below):
!
!     - *x* means model variable x
!     - eq.x referres to equation x of HS75
!     - the script T of HS75 (e.g. 1st term rhs of eq.9) is stored in *stt*
!     - the script D of HS75 (e.g. 1st term rhs of eq.8) is stored in *sdt*
!     - the script P of HS75 (e.g. 1st term rhs of eq.10) is stored in -*spt*
!     - the surface geopotential is stored in *so*
!
!     - the vertical scheme has being changed to the ECMWF scheme
!       (see e.g. Simmons and Burridge 1981, Mon.Wea.Rev.,109,758-766).
!       in this scheme,  matrix g differs from that in HS75.
!
!     - in addition to the dry HS75 version also the tendencys of specific
!       humidity are processed
!
      real azm(NSPP,NLEV)
      real atm(NSPP,NLEV)
      real aqm(NSPP,NLEV)
      real adt(NSPP,NLEV)
      real adm(NSPP,NLEV)
      real zgt(NSPP,NLEV)
      real zgm(NSPP,NLEV)
      real apm(NSPP)
!
!     franks diagnostics
!
      real, allocatable :: ztt(:,:)
      real, allocatable :: zttgp(:,:)
      real, allocatable :: zsd(:,:),zsz(:,:),zsq(:,:),zsp(:),zst(:,:)
      real, allocatable :: zugp(:,:),zvgp(:,:),zekin(:,:),zepot(:,:)
      real, allocatable :: zqgp(:,:),zpgp(:)
      real, allocatable :: zqmgp(:,:),zpmgp(:),ztgp(:,:)
!
!*    0. save prognostic variables at (t-dt)
!        and the non-linear divergence tendency terms
!
      apm(:)   = spm(:)   ! log surface pressure
      adm(:,:) = sdm(:,:) ! divergence
      azm(:,:) = szm(:,:) ! (absolut) vorticity
      atm(:,:) = stm(:,:) ! temperature
      adt(:,:) = sdt(:,:) ! divergence tendency
      if (nqspec == 1) aqm(:,:) = sqm(:,:) ! spec.  humidity
!
!*    do the advective time step
!
      if(nadv > 0) then
!
!*    1. calculate divergence on timelevel t (solving eq.17),
!        which will replace the divergence tendency sdt
!
!     1.a precompute (g*script T; *zgt*) and (phi-phi*=g*T, eq.11; *zgm*)
!         needed for the rhs of eq.17.
!         (note that phi is needed in eq.17 and, therefor,
!         the surface geopotential phi* is added later (loop 1.b))
!
       do jlev = 1 , NLEV
         do jsp=1,NSPP
            zgt(jsp,jlev) = dot_product(g(:,jlev),stt(jsp,:))
            zgm(jsp,jlev) = dot_product(g(:,jlev),atm(jsp,:))
         enddo
       enddo
!
!     1.b compute divergence at time t (i.e. solve eq.17)
!         and overwrite *sdt* with the result
!         for comparison with HS75 note:
!
!         - *spt* contains -(!)script P
!         - *spm* contains log(ps)(t-dt)
!         - surface geopotential (*so*) needs to be added (see 1.a)
!         - bm1 is the invers of matrix (1/cn I+B dt**2) (lhs eq.17)
!         - zz is set to the rhs of eq.17
!         - zsum holds the result of the dot product of rhs(eq17) and bm1
!           (therefor jlev2-loop)
!
       do jlev = 1 , NLEV
        do jsp = 1 , NSPP
          jn = nindex(jsp)
          zsum = 0.0
          if (jn > 0) then
            zq = 1.0 / (jn * (jn+1))
            do jlev2 = 1 , NLEV
              z0 = t0(jlev2)
              zt = zgt(jsp,jlev2) - z0 * spt(jsp)
              zm = zgm(jsp,jlev2) + z0 * spm(jsp)
              za = adt(jsp,jlev2) * zq + sop(jsp)
              zb = adm(jsp,jlev2) * zq
              zz = zb + delt * (zm + za + delt * zt)
              zsum = zsum + zz * bm1(jlev2,jlev,jn)
            enddo
          endif
          sdt(jsp,jlev) = zsum
        enddo
       enddo
       if (mypid == NROOT) then
        sdt(1:2,:) = 0.0 ! first mode should be zero (errors due to numerics)
       endif
!
!*    2. calculate (-log) surface pressure tendency from eq.15 (pi=*dsigma*)
!
       do jlev = 1 , NLEV
         spt = spt + dsigma(jlev) * sdt(:,jlev)
       enddo
!
!*    3. calculate temperature tendency from eq.14
!
       do jlev = 1 , NLEV
        do jsp = 1 , NSPP
         stt(jsp,jlev)=stt(jsp,jlev)-dot_product(tau(:,jlev),sdt(jsp,:))
        enddo
       enddo
!
      endif

!     3a. Coupling for synchronization runs

      if (mrnum == 2 .and. nsync > 0) then
         call mrdiff(stp,std,NESP,NLEV)
         call mrdiff(sdp,sdd,NESP,NLEV)
         call mrdiff(szp,szd,NESP,NLEV)
         call mrdiff(spp,spd,NESP,   1)
         stp(:,:) = stp(:,:) + syncstr * std(:,:)
         sdp(:,:) = sdp(:,:) + syncstr * sdd(:,:)
         szp(:,:) = szp(:,:) + syncstr * szd(:,:)
         spp(:  ) = spp(:  ) + syncstr * spd(:  )
      endif
!
!*    4. time stepping and time filtering (1st part)
!        the time filtering is splitted into two parts:
!
!        1st part: xf(prel)=x+eps*(xf(-)-2x)
!        2nd part: xf=xf(prel)+eps*x(+)
!
!        together the complete filter is xf=x+eps(xf(-)-2x+x(+))
!
!        with: x        = the prognostic variable to be filtered
!              xf       = the filtered variable (current time step)
!              xf(-)    = the filtered variable (old time step)
!              xf(prel) = a preliminary filterd variable
!              x(+)     = x at time (t+2dt)
!              eps      = filter constant = *pnu*
!
!     4.a 1st part of time filter
!
      if (nkits == 0) then
         spm = pnu21 * spp + pnu * apm
         sdm = pnu21 * sdp + pnu * adm
         szm = pnu21 * szp + pnu * azm
         stm = pnu21 * stp + pnu * atm
         if (nqspec == 1) sqm = pnu21 * sqp + pnu * aqm
      endif
!
!     4.b add tendencies
!
      if(nadv > 0 ) then
       spp = apm - delt2 * spt   ! log surface pressure (negative tendencies)
       sdp =   2.0 * sdt - adm   ! note that sdt=adm+delt*tendencies (see 1.)
       stp = delt2 * stt + atm   ! temperature
       szp = delt2 * szt + azm   ! vorticity
       if (nqspec == 1) sqp = delt2 * sqt + aqm   ! spec. humidity
      else
       spp = apm   ! log surface pressure (negative tendencies)
       sdp = adm   ! note that sdt=adm+delt*tendencies (see 1.)
       stp = atm   ! temperature
       szp = azm   ! vorticity
       if (nqspec == 1) sqp = aqm   ! spec. humidity
      endif
!
!     conserve p
!
      if (mypid == NROOT) then
       spp(1) = 0.0
       spp(2) = 0.0
      endif

!*    now the adiabatic time step is finished beside the time filtering.
!     2nd part of time filtering is done in subroutine spectrald
!
!     finaly the partial arrays are gathered from all processors (mpi)
!
      call mpgallsp(sd,sdp,NLEV)
      call mpgallsp(sz,szp,NLEV)
      call mpgallsp(st,stp,NLEV)
      call mpgallsp(sp,spp,   1)
      if (nqspec == 1) call mpgallsp(sq,sqp,NLEV)
!
!     franks diagnostic
!
      if(ndiagsp==1) then
       allocate(ztt(NESP,NLEV))
       call mpgallsp(ztt,stt,NLEV)
       dsp3d(:,1:NLEV,3)=ztt(:,1:NLEV)*ct*ww
       deallocate(ztt)
      endif
!
!     entropy/energy diagnostics
!
      if(nenergy > 0 .or. nentropy > 0) then
       allocate(ztt(NESP,NLEV))
       allocate(zttgp(NHOR,NLEV))
       allocate(zst(NESP,NLEV))
       if (nqspec == 1) allocate(zsq(NESP,NLEV))
       allocate(zsp(NESP))
       allocate(ztgp(NHOR,NLEV))
       allocate(zqgp(NHOR,NLEV))
       allocate(zqmgp(NHOR,NLEV))
       allocate(zpgp(NHOR))
       allocate(zpmgp(NHOR))
       call mpgallsp(ztt,stt,NLEV)
       call mpgallsp(zst,atm,NLEV)
       if (nqspec == 1) call mpgallsp(zsq,aqm,NLEV)
       call mpgallsp(zsp,apm,1)
       ztt(:,:)=ztt(:,:)*ct*ww
       call sp2fl(ztt,zttgp,NLEV)
       if (nqspec == 1) call sp2fl(sq,zqgp,NLEV)
       if (nqspec == 1) call sp2fl(zsq,zqmgp,NLEV)
       call sp2fl(zst,ztgp,NLEV)
       call sp2fl(zsp,zpmgp,1)
       call sp2fl(sp,zpgp,1)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
       if (nqspec == 1) then
          call fc2gp(zqgp,NLON,NLPP*NLEV)
          call fc2gp(zqmgp,NLON,NLPP*NLEV)
       else
          zqgp(:,:) = 0.0
          zqmgp(:,:) = 0.0
       endif
       call fc2gp(ztgp,NLON,NLPP*NLEV)
       call fc2gp(zpgp,NLON,NLPP)
       call fc2gp(zpmgp,NLON,NLPP)
       zpmgp(:)=psurf*exp(zpmgp(:))
       zpgp(:)=psurf*exp(zpgp(:))
       do jlev=1,NLEV
        ztgp(:,jlev)=ct*(ztgp(:,jlev)+t0(jlev))
        if (nqspec == 1) zqmgp(:,jlev)=zqmgp(:,jlev)*psurf/zpmgp(:)
        if (nqspec == 1) zqgp(:,jlev)=zqgp(:,jlev)*psurf/zpgp(:)
       enddo
       deallocate(ztt)
       deallocate(zst)
       if (nqspec == 1) deallocate(zsq)
       deallocate(zsp)
      endif       
      if(nentropy > 0) then
       dentropy(:,2)=0.
       dentropy(:,35)=0.
       dentropy(:,36)=0.
       do jlev=1,NLEV
        dentro(:)=zttgp(:,jlev)/dentrot(:,jlev)                         &
     &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)&
     &         -gascon*(1.+(1./rdbrv-1.)*dentroq(:,jlev))               &
     &         *(zpgp(:)-zpmgp(:))/deltsec2/dentrop(:)                  &
     &         *dentrop(:)*dsigma(jlev)/ga
        dentropy(:,2)=dentropy(:,2)+dentro(:)
        if(nentro3d > 0) dentro3d(:,jlev,2)=dentro(:)
        dentro(:)=zttgp(:,jlev)/dentrot(:,jlev)                         &
     &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)
        dentropy(:,35)=dentropy(:,35)+dentro(:)
        dentropy(:,36)=dentropy(:,36)                                   &
     &                +((acpd*(1.+adv*zqgp(:,jlev))                     & 
     &                  *log(ztgp(:,jlev)+zttgp(:,jlev)*deltsec2)       &
     &                  -gascon*(1.+(1./rdbrv-1.)*zqgp(:,jlev))         &
     &                  *log(zpgp(:)))*zpgp(:)                          &
     &                 -(acpd*(1.+adv*zqmgp(:,jlev))                    &
     &                  *log(ztgp(:,jlev))                              &
     &                  -gascon*(1.+(1./rdbrv-1.)*zqmgp(:,jlev))        &
     &                  *log(zpmgp(:)))*zpmgp(:))                       &
     &                /deltsec2/ga*dsigma(jlev) 
       enddo
      endif
      if(nenergy > 0) then
       allocate(zsd(NESP,NLEV))
       allocate(zsz(NESP,NLEV))
       allocate(zugp(NHOR,NLEV))
       allocate(zvgp(NHOR,NLEV))
       allocate(zekin(NHOR,NLEV))
       allocate(zepot(NHOR,NLEV))
       call mpgallsp(zsd,adm,NLEV)
       call mpgallsp(zsz,azm,NLEV)
       call dv2uv(zsd,zsz,zugp,zvgp)
       call fc2gp(zugp,NLON,NLPP*NLEV)
       call fc2gp(zvgp,NLON,NLPP*NLEV)
       do jlev=1,NLEV
        zekin(:,jlev)=0.5*(zugp(:,jlev)*zugp(:,jlev)                    &
     &                    +zvgp(:,jlev)*zvgp(:,jlev))*cv*cv*rcsq(:)     &
     &               *zpmgp(:)   
        zepot(:,jlev)=ztgp(:,jlev)*acpd*(1.+adv*zqmgp(:,jlev))*zpmgp(:) 
       enddo
       call dv2uv(sd,sz,zugp,zvgp)
       call fc2gp(zugp,NLON,NLPP*NLEV)
       call fc2gp(zvgp,NLON,NLPP*NLEV)
       denergy(:,27)=0.
       denergy(:,26)=0.
       denergy(:,2)=0.
       denergy(:,1)=0.
       do jlev=1,NLEV
        zekin(:,jlev)=0.5*(zugp(:,jlev)*zugp(:,jlev)                    &
     &                    +zvgp(:,jlev)*zvgp(:,jlev))*cv*cv*rcsq(:)     &
     &               *zpgp(:)                                           &
     &               -zekin(:,jlev)
        denergy(:,27)=denergy(:,27)                                     &
     &               -zekin(:,jlev)/deltsec2/ga*dsigma(jlev)
        denergy(:,1)=denergy(:,1)                                       &
     &              +(ztgp(:,jlev)+zttgp(:,jlev)*deltsec2)              &
     &              *acpd*(1.+adv*zqgp(:,jlev))*zpgp(:)*dsigma(jlev)/ga
        denergy(:,2)=denergy(:,2)                                       &
     &              +zttgp(:,jlev)                                      &
     &              *acpd*(1.+adv*zqgp(:,jlev))*zpgp(:)/ga*dsigma(jlev)
        denergy(:,26)=denergy(:,26)                                     &
     &               +((ztgp(:,jlev)+zttgp(:,jlev)*deltsec2)            &
     &                 *acpd*(1.+adv*zqgp(:,jlev))*zpgp(:)              &
     &                -zepot(:,jlev))/deltsec2/ga*dsigma(jlev)
        if(nener3d > 0) then 
         dener3d(:,jlev,27)=-zekin(:,jlev)/deltsec2/ga*dsigma(jlev)
         dener3d(:,jlev,1)=(ztgp(:,jlev)+zttgp(:,jlev)*deltsec2)        &
     &               *acpd*(1.+adv*zqgp(:,jlev))*zpgp(:)*dsigma(jlev)/ga
         dener3d(:,jlev,2)=zttgp(:,jlev)                                &
     &              *acpd*(1.+adv*zqgp(:,jlev))*zpgp(:)/ga*dsigma(jlev)
         dener3d(:,jlev,26)=((ztgp(:,jlev)+zttgp(:,jlev)*deltsec2)      &
     &                     *acpd*(1.+adv*zqgp(:,jlev))*zpgp(:)          &
     &                     -zepot(:,jlev))/deltsec2/ga*dsigma(jlev)
        endif  
       enddo
       deallocate(zsd)
       deallocate(zsz)
       deallocate(zugp)
       deallocate(zvgp)
       deallocate(zekin)
       deallocate(zepot)
      endif
      if(nenergy > 0 .or. nentropy > 0) then
       deallocate(zttgp)
       deallocate(zqgp)
       deallocate(zqmgp)
       deallocate(zpgp)
       deallocate(zpmgp)
       deallocate(ztgp)
      endif
!
      return
      end

!-Diabatic-subroutines  (Frank (Larry) 26-Nov-99)

!     =====================
!     SUBROUTINE GRIDPOINTD
!     =====================

      subroutine gridpointd
      use pumamod
!
      real sdf(NESP,NLEV)
      real stf(NESP,NLEV)
      real sqf(NESP,NLEV)
      real szf(NESP,NLEV)

      real zqout(NHOR,NLEV)
      real zgq(NLON,NLAT,NLEV)
      real zmmr(NLON,NLAT,NLEV)
      real znrho(NLON,NLAT,NLEV)
!
!*    Diabatic Gridpoint Calculations
!
      gudt(:,:)=0.
      gvdt(:,:)=0.
      gtdt(:,:)=0.
      gqdt(:,:)=0.
      dudt(:,:)=0.
      dvdt(:,:)=0.
      dtdt(:,:)=0.
      dqdt(:,:)=0.
      mmrt(:,:)=0.

!
!     transform to gridpoint domain
!

      call invlegd

      call fc2gp(gu  ,NLON,NLPP*NLEV)
      call fc2gp(gv  ,NLON,NLPP*NLEV)
      call fc2gp(gt  ,NLON,NLPP*NLEV)
      call fc2gp(gp  ,NLON,NLPP)
      if (nqspec == 1) then
         call fc2gp(gq  ,NLON,NLPP*NLEV)
      endif

!
!     dimensionalize
!

      dp(:)=psurf*exp(gp(:))
      do jlev = 1 , NLEV
       du(:,jlev)=cv*gu(:,jlev)*SQRT(rcsq(:))
       dv(:,jlev)=cv*gv(:,jlev)*SQRT(rcsq(:))
       dt(:,jlev)=ct*(gt(:,jlev)+t0(jlev))
       if (nqspec == 1) dq(:,jlev)=gq(:,jlev)*psurf/dp(:)
      enddo
!
!     tracer transport
!
      if (nsela > 0 .and. nkits == 0 .and. nqspec == 0) then	  
       if (nqspec == 0) then
        write(nud,*) 'Semi-Lagrangian q running'
        flush(nud)
        dqt(:,:) = dq(:,:) ! Save old value of q
        call mpgagp(zgq,dq,NLEV)
        if (mypid == NROOT) then
         do jlat = 1 , NLAT
          dtrace(:,NLAT+1-jlat,:,1) = zgq(:,jlat,:)
         enddo ! jlat
        endif ! mypid
       endif ! nqspec
       call tracer_main
       if (nqspec == 0) then
        if (mypid == NROOT) then
         do jlat = 1 , NLAT
          zgq(:,jlat,:) = dtrace(:,NLAT+1-jlat,:,1)
         enddo
        endif ! mypid
        call mpscgp(zgq,dq,NLEV)
        dqt(:,:) = (dq(:,:) - dqt(:,:)) / deltsec !  q advection term
       endif ! nqspec	  
      endif ! nkits
      call guihor("DQ" // char(0),dq,NLEV,1000.0,0.0)

!     aerosol transport
  
      if (nsela == 1 .and. nkits == 0 .and. l_aero > 0) then
        mmrt(:,:) = mmr(:,:) ! Save old value of mmr
        call mpgagp(zmmr,mmr,NLEV)
        call mpgagp(znrho,nrho,NLEV)
        if (mypid == NROOT) then
         do jlat = 1 , NLAT
          daeros(:,NLAT+1-jlat,:,1) = zmmr(:,jlat,:)
          numrhos(:,NLAT+1-jlat,:,1) = znrho(:,jlat,:)
         enddo ! jlat
        endif ! mypid
        call aero_main
        if (mypid == NROOT) then
         do jlat = 1 , NLAT
          zmmr(:,jlat,:) = daeros(:,NLAT+1-jlat,:,1)
          znrho(:,jlat,:) = numrhos(:,NLAT+1-jlat,:,1)
         enddo
        endif ! mypid
        call mpscgp(zmmr,mmr,NLEV)
        call mpscgp(znrho,nrho,NLEV)
        mmrt(:,:) = (mmr(:,:) - mmrt(:,:)) / deltsec !  q advection term       endif !  
      endif ! nkits
!
!     compute output specific humidity
!

      if (mod(nstep,nafter)==0 .and. nqspec == 1) then
       do jlev=1,NLEV
        zqout(:,jlev)=gq(:,jlev)/exp(gp(:))
        call gp2fc(zqout(1,jlev),NLON,NLPP)
        call fc2sp(zqout(1,jlev),sqout(1,jlev))
       enddo
       call mpsum(sqout,NLEV)
      endif

!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
          write(nud,*)'------------ SUBROUTINE GRIDPOINTD -------------'
          write(nud,*)'NSTEP= ',nstep
          write(nud,*)'YY,MM,DD,HH,MIN= ',ndatim(1:5)
          write(nud,*) 'gp= ',MAXVAL(gp(1:NHOR)),MINVAL(gp(1:NHOR))
          write(nud,*) 'dp= ',MAXVAL(dp(1:NHOR)),MINVAL(dp(1:NHOR))
       endif
       call prdbug1
      endif

!
!     PARAMETERIZATION ROUTINES TO BE INCLUDED BELOW
!
!     a) miscellaneous parameterizations (like q-fixer)
!

      call miscstep

!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        write(nud,*) 'After MISC:'
       endif
       call prdbug2
      endif

!
!     b) surface fluxes and vertical diffusion
!

      call fluxstep
     
!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        write(nud,*)'After FLUX:'
       endif
       call prdbug2
      endif

!
!     c) radiation
!

      if(nrad > 0) call radstep
     
!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        write(nud,*)'After RAD:'
       endif
       call prdbug2
      endif

!
!     d) convective and large scale rain and clouds
!

      call rainstep

!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        write(nud,*)'After RAIN:'
       endif
       call prdbug2
      endif

!
!     e) surface parmeterizations or models
!

      call surfstep

!
!     f) carbon-silicate weathering parameterization
!

      call carbonstep
      
      
!
!     g) newtonian cooling of stratosphere
!
      
      if (nstratosponge > 0) then
         mint(:) = minval(dt(:,:NLEV-12),2)
         do k=1,NLEV-12
            dtdt(:,k) = dtdt(:,k) + (mint(:) - (dt(:,k)+dtdt(:,k)))/(taucool*day_24hr)
         enddo
      endif
      
!
!     h) hurricane/storm diagnostics
!
      
      call hurricanestep   

!
!     END OF PARAMETERISATION ROUTINES
!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        write(nud,*)'FINAL:'
       endif
       call prdbug2
      endif

!
!     franks diagnostic
!

      if(ndiaggp==1) then
       dgp3d(:,1:NLEV,1)=dtdt(:,1:NLEV)
       do jlev=1,NLEV
        dgp3d(:,jlev,11)=dqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000
       enddo
      endif
!
!     entropy/energy diagnostics
!
      if(nentropy > 0) then
       dentropy(:,4)=0.
       dentropy(:,33)=0.
       do jlev=1,NLEV
        dentro(:)=dtdt(:,jlev)/dentrot(:,jlev)                          &
     &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)
        dentropy(:,4)=dentropy(:,4)+dentro(:)
        if(nentro3d > 0) dentro3d(:,jlev,4)=dentro(:)
        zfac=ww*tfrc(jlev)/(1.+tfrc(jlev)*delt2)
        dentro(:)=((du(:,jlev)+dudt(:,jlev)*deltsec2)**2                &
     &             +(dv(:,jlev)+dvdt(:,jlev)*deltsec2)**2)              &
     &             *zfac*dentrop(:)/ga*dsigma(jlev)/dentrot(:,jlev)
        dentropy(:,33)=dentropy(:,33)+dentro(:)
        if(nentro3d > 0) dentro3d(:,jlev,23)=dentro(:)
       enddo
      endif
      if(nenergy > 0) then
       denergy(:,4)=0.
       do jlev=1,NLEV
        denergy(:,4)=denergy(:,4)                                       &
     &              +dtdt(:,jlev)                                       &
     &              *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,4)=dtdt(:,jlev)                                 &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif
       enddo
      endif
!
!     de-dimensionalize
!

      do jlev = 1 , NLEV
       gudt(:,jlev)=dudt(:,jlev)/SQRT(rcsq(:))/cv/ww                    &
     &             +gudt(:,jlev)
       gvdt(:,jlev)=dvdt(:,jlev)/SQRT(rcsq(:))/cv/ww                    &
     &             +gvdt(:,jlev)
       gtdt(:,jlev)=dtdt(:,jlev)/ct/ww                                  &
     &             +gtdt(:,jlev)
       gqdt(:,jlev)=dqdt(:,jlev)*dp(:)/ww/psurf                         &
     &             +gqdt(:,jlev)
      enddo

!
!     transform to spectral space
!

      call gp2fc(gtdt,NLON,NLPP*NLEV)
      call gp2fc(gudt,NLON,NLPP*NLEV)
      call gp2fc(gvdt,NLON,NLPP*NLEV)
      if (nqspec == 1) call gp2fc(gqdt,NLON,NLPP*NLEV)

!
!     direct Legendre transformation (fourier domain to spectral domain)
!
      do jlev = 1 , NLEV
         call fc2sp(gtdt(1,jlev),stf(1,jlev))
         if (nqspec == 1) call fc2sp(gqdt(1,jlev),sqf(1,jlev))
      enddo

      call uv2dv(gudt,gvdt,sdf,szf)

      call mpsumsc(stf,stt,NLEV)
      call mpsumsc(sdf,sdt,NLEV)
      call mpsumsc(szf,szt,NLEV)
      if (nqspec == 1) call mpsumsc(sqf,sqt,NLEV)
      if (nqspec == 0) dq(:,:) = dq(:,:) + dqdt(:,:) * deltsec
      if (nsela == 1 .and. l_aero > 0) mmr(:,:) = mmr(:,:) + mmrt(:,:) * deltsec

      return
      end

!     ==================
!     SUBROUTINE PRDBUG1
!     ==================

      subroutine prdbug1
      use pumamod

      integer kmaxl(1),kminl(1)

      real zfpr1(NLON*NLAT,NLEV)
      real zfpr2(NLON*NLAT,NLEV)
      real zfpr3(NLON*NLAT,NLEV)
      real zfpr4(NLON*NLAT,NLEV)
      real zfpr5(NLON*NLAT,NLEV)
      real zfpr6(NLON*NLAT,NLEV)

      call mpgagp(zfpr1,dt,NLEV)
      call mpgagp(zfpr2,dq,NLEV)
      call mpgagp(zfpr3,du,NLEV)
      call mpgagp(zfpr4,dv,NLEV)
      call mpgagp(zfpr5,dcc,NLEV)
      call mpgagp(zfpr6,dql,NLEV)

      if(mypid==NROOT) then
       if(nprint==1) then
        write(nud,*)'Global diagnostics: '
        do jlev=1,NLEV
         zmaxv=MAXVAL(zfpr1(:,jlev))
         kmaxl=MAXLOC(zfpr1(:,jlev))
         zminv=MINVAL(zfpr1(:,jlev))
         kminl=MINLOC(zfpr1(:,jlev))
         write(nud,*)'L= ',jlev,' MAX T= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN T= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr2(:,jlev))
         kmaxl=MAXLOC(zfpr2(:,jlev))
         zminv=MINVAL(zfpr2(:,jlev))
         kminl=MINLOC(zfpr2(:,jlev))
         write(nud,*)'L= ',jlev,' MAX Q= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN Q= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr3(:,jlev))
         kmaxl=MAXLOC(zfpr3(:,jlev))
         zminv=MINVAL(zfpr3(:,jlev))
         kminl=MINLOC(zfpr3(:,jlev))
         write(nud,*)'L= ',jlev,' MAX U= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN U= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr4(:,jlev))
         kmaxl=MAXLOC(zfpr4(:,jlev))
         zminv=MINVAL(zfpr4(:,jlev))
         kminl=MINLOC(zfpr4(:,jlev))
         write(nud,*)'L= ',jlev,' MAX V= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN V= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr5(:,jlev))
         kmaxl=MAXLOC(zfpr5(:,jlev))
         zminv=MINVAL(zfpr5(:,jlev))
         kminl=MINLOC(zfpr5(:,jlev))
         write(nud,*)'L= ',jlev,' MAX C= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN C= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr6(:,jlev))
         kmaxl=MAXLOC(zfpr6(:,jlev))
         zminv=MINVAL(zfpr6(:,jlev))
         kminl=MINLOC(zfpr6(:,jlev))
         write(nud,*)'L= ',jlev,' MAX L= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN L= ',zminv,' NHOR= ',kminl(1)
        enddo
       elseif(nprint==2) then
        write(nud,*)'Local diagnostics at nhor= ',nprhor,': '
        do jlev=1,NLEV
         write(nud,*)'L= ',jlev,' T= ',zfpr1(nprhor,jlev)                    &
     &                    ,' Q= ',zfpr2(nprhor,jlev)
         write(nud,*)'L= ',jlev,' U= ',zfpr3(nprhor,jlev)                    &
     &                    ,' V= ',zfpr4(nprhor,jlev)
         write(nud,*)'L= ',jlev,' C= ',zfpr5(nprhor,jlev)                    &
     &                    ,' L= ',zfpr6(nprhor,jlev)
        enddo
       endif
      endif

      call mpgagp(zfpr1(1,1),dls,1)
      call mpgagp(zfpr2(1,1),dp,1)
      call mpgagp(zfpr3(1,1),drhs,1)
      call mpgagp(zfpr4(1,1),dalb,1)
      call mpgagp(zfpr5(1,1),dt(1,NLEP),1)
      call mpgagp(zfpr6(1,1),dq(1,NLEP),1)

      if(mypid==NROOT) then
       if(nprint==1) then
        write(nud,*)'Global diagnostic:'
        zmaxv=MAXVAL(zfpr5(:,1))
        kmaxl=MAXLOC(zfpr5(:,1))
        zminv=MINVAL(zfpr5(:,1))
        kminl=MINLOC(zfpr5(:,1))
        write(nud,*)'SURF MAX T= ',zmaxv,' NHOR= ',kmaxl(1)
        write(nud,*)'SURF MIN T= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr6(:,1))
        kmaxl=MAXLOC(zfpr6(:,1))
        zminv=MINVAL(zfpr6(:,1))
        kminl=MINLOC(zfpr6(:,1))
        write(nud,*)'SURF MAX Q= ',zmaxv,' NHOR= ',kmaxl(1)
        write(nud,*)'SURF MIN Q= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr2(:,1))
        kmaxl=MAXLOC(zfpr2(:,1))
        zminv=MINVAL(zfpr2(:,1))
        kminl=MINLOC(zfpr2(:,1))
        write(nud,*)'SURF MAX P= ',zmaxv,' NHOR= ',kmaxl(1)
        write(nud,*)'SURF MIN P= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr3(:,1))
        kmaxl=MAXLOC(zfpr3(:,1))
        zminv=MINVAL(zfpr3(:,1))
        kminl=MINLOC(zfpr3(:,1))
        write(nud,*)'SURF MAX RHS= ',zmaxv,' NHOR= ',kmaxl(1)
        write(nud,*)'SURF MIN RHS= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr4(:,1))
        kmaxl=MAXLOC(zfpr4(:,1))
        zminv=MINVAL(zfpr4(:,1))
        kminl=MINLOC(zfpr4(:,1))
        write(nud,*)'SURF MAX ALB= ',zmaxv,' NHOR= ',kmaxl(1)
        write(nud,*)'SURF MIN ALB= ',zminv,' NHOR= ',kminl(1)
       elseif(nprint==2) then
        write(nud,*)'Local diagnostics at nhor= ',nprhor,': '
        write(nud,*)'SURF: LS= ',zfpr1(nprhor,1)                             &
     &             ,' PS= ',zfpr2(nprhor,1)
        write(nud,*)'SURF: RHS= ',zfpr3(nprhor,1)                            &
     &             ,' ALB= ',zfpr4(nprhor,1)
        write(nud,*)'SURF: TS= ',zfpr5(nprhor,1)                             &
     &             ,' QS= ',zfpr6(nprhor,1)
       endif
      endif

      if(nprint==2) then
       call mpgagp(zfpr1(1,1),dicec,1)

       if(mypid==NROOT) then
        write(nud,*)'SURF: ICEC= ',zfpr1(nprhor,1)
       endif
      endif

      return
      end subroutine prdbug1

!     ==================
!     SUBROUTINE PRDBUG2
!     ==================

      subroutine prdbug2
      use pumamod

      integer kmaxl(1),kminl(1)

      real zfpr1(NLON*NLAT,NLEV)
      real zfpr2(NLON*NLAT,NLEV)
      real zfpr3(NLON*NLAT,NLEV)
      real zfpr4(NLON*NLAT,NLEV)

      call mpgagp(zfpr1,dtdt,NLEV)
      call mpgagp(zfpr2,dqdt,NLEV)
      call mpgagp(zfpr3,dudt,NLEV)
      call mpgagp(zfpr4,dvdt,NLEV)

      if(mypid==NROOT) then
       if(nprint==1) then
        write(nud,*)'Global diagnostics: '
        do jlev=1,NLEV
         zmaxv=MAXVAL(zfpr1(:,jlev))
         kmaxl=MAXLOC(zfpr1(:,jlev))
         zminv=MINVAL(zfpr1(:,jlev))
         kminl=MINLOC(zfpr1(:,jlev))
         write(nud,*)'L= ',jlev,' MAX DT= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN DT= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr2(:,jlev))
         kmaxl=MAXLOC(zfpr2(:,jlev))
         zminv=MINVAL(zfpr2(:,jlev))
         kminl=MINLOC(zfpr2(:,jlev))
         write(nud,*)'L= ',jlev,' MAX DQ= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN DQ= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr3(:,jlev))
         kmaxl=MAXLOC(zfpr3(:,jlev))
         zminv=MINVAL(zfpr3(:,jlev))
         kminl=MINLOC(zfpr3(:,jlev))
         write(nud,*)'L= ',jlev,' MAX DU= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN DU= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr4(:,jlev))
         kmaxl=MAXLOC(zfpr4(:,jlev))
         zminv=MINVAL(zfpr4(:,jlev))
         kminl=MINLOC(zfpr4(:,jlev))
         write(nud,*)'L= ',jlev,' MAX DV= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN DV= ',zminv,' NHOR= ',kminl(1)
        enddo
       elseif(nprint==2) then
        write(nud,*)'Local diagnostics at nhor= ',nprhor,': '
        do jlev=1,NLEV
         write(nud,*)'L= ',jlev,' DT= ',zfpr1(nprhor,jlev)                   &
     &                    ,' DQ= ',zfpr2(nprhor,jlev)
         write(nud,*)'L= ',jlev,' DU= ',zfpr3(nprhor,jlev)                   &
     &                    ,' DV= ',zfpr4(nprhor,jlev)
        enddo
       endif
      endif

      return
      end subroutine prdbug2

!     ====================
!     SUBROUTINE SPECTRALD
!     ====================

      subroutine spectrald
      use pumamod
!
      real :: zsdt1(NSPP,NLEV),zsdt2(NSPP,NLEV)
      real :: zszt1(NSPP,NLEV),zszt2(NSPP,NLEV)
      real :: zsum3(1) ! must be an array because of the NAG compiler
!
!     franks diagnostics
!
      real, allocatable :: ztt(:,:)
      real, allocatable :: zttgp(:,:)
      real, allocatable :: zst(:,:),zsq(:,:),zstt(:,:),zstt2(:,:)
      real, allocatable :: zsttd(:,:)
      real, allocatable :: ztgp(:,:),zdtgp(:,:),zqgp(:,:)
      real, allocatable :: zgw(:),zsum1(:)
!
!     prepare diagnostics of efficiency
!
      if(ndheat > 1) then
       allocate(zst(NESP,NLEV))      
       if (nqspec == 1) allocate(zsq(NESP,NLEV))      
       allocate(zstt(NESP,NLEV)) 
       allocate(zstt2(NESP,NLEV)) 
       allocate(ztgp(NHOR,NLEV)) 
       allocate(zqgp(NHOR,NLEV)) 
       allocate(zdtgp(NHOR,NLEV)) 
       allocate(zsum1(4))
       allocate(zgw(NHOR))
       call mpgallsp(zst,stp,NLEV)
       call mpgallsp(zstt,stt,NLEV)
       if (nqspec == 1) call mpgallsp(zsq,sqp,NLEV)
       jhor=0
       do jlat=1,NLPP
        do jlon=1,NLON
         jhor=jhor+1
         zgw(jhor)=gwd(jlat)
        enddo
       enddo
      endif

!
!     add tendencies from diabatic parameterizations
!

      szp = szp + delt2 * szt
      stp = stp + delt2 * stt
      sdp = sdp + delt2 * sdt
      if (nqspec == 1) sqp = sqp + delt2 * sqt

!
!     franks diagnostic
!

      if(ndiagsp==1) then
       allocate(ztt(NESP,NLEV))
       call mpgallsp(ztt,stt,NLEV)
       dsp3d(:,1:NLEV,1)=ztt(:,1:NLEV)*ct*ww
       deallocate(ztt)
      endif
      if(ndiaggp==1) then
       call mkdqtgp
       do jlev=1,NLEV
        dgp3d(:,jlev,20)=dqt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      endif
!
!     entropy/energy diagnostics
!
      if(nentropy > 0) then
       allocate(ztt(NESP,NLEV))
       allocate(zttgp(NHOR,NLEV))
       call mpgallsp(ztt,stt,NLEV)
       ztt(:,:)=ztt(:,:)*ct*ww
       call sp2fl(ztt,zttgp,NLEV)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
       deallocate(ztt)
       dentropy(:,3)=0.
       do jlev=1,NLEV
        dentro(:)=zttgp(:,jlev)/dentrot(:,jlev)                         &
     &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)
        dentropy(:,3)=dentropy(:,3)+dentro(:)
        if(nentro3d > 0) dentro3d(:,jlev,3)=dentro(:) 
       enddo
       deallocate(zttgp)
      endif
      if(nenergy > 0) then
       allocate(ztt(NESP,NLEV))
       allocate(zttgp(NHOR,NLEV))
       call mpgallsp(ztt,stt,NLEV)
       ztt(:,:)=ztt(:,:)*ct*ww
       call sp2fl(ztt,zttgp,NLEV)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
       deallocate(ztt)
       denergy(:,3)=0.
       do jlev=1,NLEV
        denergy(:,3)=denergy(:,3)                                       &
     &              +zttgp(:,jlev)                                      &
     &              *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,3)=zttgp(:,jlev)                                &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif 
       enddo
       deallocate(zttgp)
      endif

!     calculates spectral tendencies from restoration (if included)
!     and biharmonic diffusion (both implicitely).
!     add newtonian cooling and drag.

      do jlev = 1 , NLEV
       stt(:,jlev)=((-damp(jlev)-tdisst(jlev)*sakpp(1:NSPP,jlev))       &
     &              *stp(:,jlev)+damp(jlev)*srp(:,jlev))                &
     &          /(1.+delt2*(damp(jlev)+tdisst(jlev)*sakpp(1:NSPP,jlev)))
       if (nqspec == 1)                                                 &
     & sqt(:,jlev)=-tdissq(jlev)*sakpp(1:NSPP,jlev)*sqp(:,jlev)         &
     &            /(1.+delt2*tdissq(jlev)*sakpp(1:NSPP,jlev))
       zsdt1(:,jlev)=-tfrc(jlev)*sdp(:,jlev)                            &
     &              /(1.+delt2*tfrc(jlev))
       zsdt2(:,jlev)=-tdissd(jlev)*sakpp(1:NSPP,jlev)*sdp(:,jlev)       &
     &              /(1.+delt2*tdissd(jlev)*sakpp(1:NSPP,jlev))
       sdt(:,jlev)=zsdt1(:,jlev)+zsdt2(:,jlev)
       zszt1(:,jlev)=-tfrc(jlev)*szp(:,jlev)                            &
     &              /(1.+delt2*tfrc(jlev))
       zszt2(:,jlev)=-tdissz(jlev)*sakpp(1:NSPP,jlev)*szp(:,jlev)       &
     &              /(1.+delt2*tdissz(jlev)*sakpp(1:NSPP,jlev))
       szt(:,jlev)=zszt1(:,jlev)+zszt2(:,jlev)
      enddo
!
!     energy/entropy diagnostics
!
      if(nenergy > 0 .or. nentropy > 0) then
       allocate(zsttd(NSPP,NLEV))
       do jlev=1,NLEV
        zsttd(:,jlev)=-tdisst(jlev)*sakpp(1:NSPP,jlev)*stp(:,jlev)      &
     &               /(1.+delt2*tdisst(jlev)*sakpp(1:NSPP,jlev))
       enddo
      endif

!     no friction on planetary vorticity
!     no diffusion on plavor (planetary vorticity)

      if (mypid == NROOT) then
         spp(1) = 0.0
         spp(2) = 0.0
         zszt1(3,:)=zszt1(3,:)+plavor*tfrc(:)/(1.+delt2*tfrc(:))
         zszt2(3,:)=zszt2(3,:)+plavor*tdissz(:)*sakpp(3,:)              &
     &                               /(1.+delt2*tdissz(:)*sakpp(3,:))
         szt(3,:)=szt(3,:)+plavor*(tfrc(:)/(1.+delt2*tfrc(:))           &
     &                            +tdissz(:)*sakpp(3,:)                 &
     &                            /(1.+delt2*tdissz(:)*sakpp(3,:)))
      endif
!
!     sponge layer ad the top
!     (relax t(1) to t(2))
!
      if(nsponge > 0) then
       stt(:,1)=dampsp*(stp(:,2)-stp(:,1))+stt(:,1)
      endif
!
!     add temp-tendencies due to momentum diffusion/dissipation
!     (if switched on)
!
      if(ndheat > 0) call mkdheat(zszt1,zszt2,zsdt1,zsdt2)
!
!     diagnostics of efficiency
!
      if(ndheat > 1) then
       call mpgallsp(zstt2,stt,NLEV)
       zstt(:,:)=zstt(:,:)+zstt2(:,:)
       if (nqspec == 1) call sp2fl(zsq,zqgp,NLEV)
       call sp2fl(zst,ztgp,NLEV)
       call sp2fl(zstt,zdtgp,NLEV)
       if (nqspec == 1) then
          call fc2gp(zqgp,NLON,NLPP*NLEV)
       else
          zqgp(:,:) = 0.0
       endif
       call fc2gp(ztgp,NLON,NLPP*NLEV)
       call fc2gp(zdtgp,NLON,NLPP*NLEV)
       zsum1(:)=0.
       do jlev=1,NLEV
        if (nqspec == 1) zqgp(:,jlev)=zqgp(:,jlev)*psurf/dp(:)
        ztgp(:,jlev)=ct*(ztgp(:,jlev)+t0(jlev))
        zdtgp(:,jlev)=ct*ww*zdtgp(:,jlev)
        zsum1(1)=zsum1(1)+SUM(zdtgp(:,jlev)*zgw(:)                      &
     &                *acpd*(1.+adv*zqgp(:,jlev))*dp(:)/ga*dsigma(jlev) &
     &                ,mask=(zdtgp(:,jlev) >= 0.))
        zsum1(2)=zsum1(2)+SUM(zdtgp(:,jlev)*zgw(:)                      &
     &                *acpd*(1.+adv*zqgp(:,jlev))*dp(:)/ga*dsigma(jlev) &
     &               ,mask=(zdtgp(:,jlev) < 0.))
        zsum1(3)=zsum1(3)+SUM(zdtgp(:,jlev)/ztgp(:,jlev)*zgw(:)         &
     &                *acpd*(1.+adv*zqgp(:,jlev))*dp(:)/ga*dsigma(jlev) &
     &               ,mask=(zdtgp(:,jlev) >= 0.))
        zsum1(4)=zsum1(4)+SUM(zdtgp(:,jlev)/ztgp(:,jlev)*zgw(:)         & 
     &                *acpd*(1.+adv*zqgp(:,jlev))*dp(:)/ga*dsigma(jlev) &
     &               ,mask=(zdtgp(:,jlev) < 0.))
       enddo
       zsum3(1)=SUM(zgw(:))
       call mpsumbcr(zsum1,4)
       call mpsumbcr(zsum3,1)
       zsum1(:)=zsum1(:)/zsum3(1)
       if(mypid == NROOT) then
        ztp=zsum1(1)/zsum1(3)
        ztm=zsum1(2)/zsum1(4)
        write(9,*) zsum1(:),zsum1(1)/zsum1(3),zsum1(2)/zsum1(4)         &
     &            ,(ztp-ztm)/ztp
       endif
       deallocate(zst)
       if (nqspec == 1) deallocate(zsq)
       deallocate(zstt)
       deallocate(zstt2)
       deallocate(ztgp)
       deallocate(zqgp)
       deallocate(zdtgp)
       deallocate(zsum1)
       deallocate(zgw)
      endif
!
!     add diabatic tendencies
!
      szp = szp + delt2 * szt
      stp = stp + delt2 * stt
      sdp = sdp + delt2 * sdt
      if (nqspec == 1) sqp = sqp + delt2 * sqt

!     initial divergence damping for a smooth start of the model in case
!     of steep orography (Mars) or unusual initial conditions

      if (ndivdamp > 0) then
         zdd = 1.0 / (1.0 + ndivdamp * 0.02)
         sdp(:,:) = sdp(:,:) * zdd
         write(nud,*) '### Damping with ',zdd
         ndivdamp = ndivdamp - 1        
      endif

      if (nkits == 0) then
         szm = szm + pnu * szp
         stm = stm + pnu * stp
         sdm = sdm + pnu * sdp
         spm = spm + pnu * spp
         if (nqspec == 1) sqm = sqm + pnu * sqp
      endif

      call mpgallsp(sd,sdp,NLEV)
      call mpgallsp(sz,szp,NLEV)
      call mpgallsp(st,stp,NLEV)
      call mpgallsp(sp,spp,   1)
      if (nqspec == 1) call mpgallsp(sq,sqp,NLEV)
!
!     franks diagnostic
!

      if(ndiagsp==1) then
       allocate(ztt(NESP,NLEV))
       call mpgallsp(ztt,stt,NLEV)
       dsp3d(:,1:NLEV,2)=ztt(:,1:NLEV)*ct*ww
       deallocate(ztt)
      endif
      if(ndiaggp==1) then
       call mkdqtgp
       do jlev=1,NLEV
        dgp3d(:,jlev,21)=dqt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      endif
!
!     entropy/energy diagnostics
!
      if(nentropy > 0) then
       allocate(ztt(NESP,NLEV))
       allocate(zttgp(NHOR,NLEV))
       call mpgallsp(ztt,stt,NLEV)
       ztt(:,:)=ztt(:,:)*ct*ww
       call sp2fl(ztt,zttgp,NLEV)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
       deallocate(ztt)
       dentropy(:,5)=0.
       do jlev=1,NLEV
        dentro(:)=zttgp(:,jlev)/dentrot(:,jlev)                         &
     &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)
        dentropy(:,5)=dentropy(:,5)+dentro(:)
        if(nentro3d > 0) dentro3d(:,jlev,5)=dentro(:)
       enddo
       deallocate(zttgp)
      endif
      if(nenergy > 0) then
       allocate(ztt(NESP,NLEV))
       allocate(zttgp(NHOR,NLEV))
       call mpgallsp(ztt,stt,NLEV)
       ztt(:,:)=ztt(:,:)*ct*ww
       call sp2fl(ztt,zttgp,NLEV)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
       denergy(:,5)=0.
       do jlev=1,NLEV
        denergy(:,5)=denergy(:,5)                                       &
     &              +zttgp(:,jlev)                                      &
     &              *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,5)=zttgp(:,jlev)                                &
     &                  *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif 
       enddo
       call mpgallsp(ztt,zsttd,NLEV)
       ztt(:,:)=ztt(:,:)*ct*ww
       call sp2fl(ztt,zttgp,NLEV)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
       denergy(:,24)=0.
       do jlev=1,NLEV
        denergy(:,24)=denergy(:,24)                                     &
     &              +zttgp(:,jlev)                                      &
     &              *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,24)=zttgp(:,jlev)                               &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif 
       enddo
       deallocate(ztt)
       deallocate(zttgp)
      endif
      if(nenergy > 0 .or. nentropy > 0) then
       deallocate(zsttd)
      endif
!
      return
      end subroutine spectrald

!     ===============
!     SUBROUTINE MKDQ
!     ===============

!      subroutine mkdq(psq,pgq,KLEV)
!      use pumamod
!
!      real psq(NESP,klev),pgq(NHOR,klev)
!
!      real zps(NHOR)
!
!      call sp2fl(psq,pgq,KLEV)
!      call fc2gp(pgq,NLON,NLPP*KLEV)
!      call sp2fl(sp,zps,1)
!      call fc2gp(zps,NLON,NLPP)
!
!      zps(:)=psurf*EXP(zps(:))
!      do jlev=1,KLEV
!       pgq(:,jlev)=pgq(:,jlev)*psurf/zps(:)
!      enddo
!
!      return
!      end subroutine mkdq


!     ===================
!     SUBROUTINE MKSECOND
!     ===================

      subroutine mksecond(ps,poff)
!
      real ps,poff
      integer ivalues(8)
!
      call date_and_time(values=ivalues)
!
      zdays = ivalues(3) - 1
      zhour = 24.0 * zdays + ivalues(5)
      zmins = 60.0 * zhour + ivalues(6)
      zsecs = 60.0 * zmins + ivalues(7) + 0.001 * ivalues(8)
      ps = zsecs - poff
!
      return
      end

!     ==================
!     SUBROUTINE MKDHEAT
!     ==================

      subroutine mkdheat(zszt1,zszt2,zsdt1,zsdt2)
      use pumamod
!
      real zszt1(NSPP,NLEV),zszt2(NSPP,NLEV)
      real zsdt1(NSPP,NLEV),zsdt2(NSPP,NLEV)
!
      real zsd(NESP,NLEV),zsz(NESP,NLEV),zsq(NESP,NLEV)
      real zsdp(NSPP,NLEV),zszp(NSPP,NLEV),zsqp(NSPP,NLEV)
      real zu(NHOR,NLEV),zun(NHOR,NLEV)
      real zv(NHOR,NLEV),zvn(NHOR,NLEV)
      real zq(NHOR,NLEV)
!
      real zdtdt(NHOR,NLEV),zdekin(NHOR,NLEV)
      real zsde(NSPP,NLEV),zsdef(NESP,NLEV)
      real zstt1(NSPP,NLEV),zstf1(NESP,NLEV)
      real zstt2(NSPP,NLEV),zstf2(NESP,NLEV)
!
      zsdp(:,:)=sdp(:,:)
      zszp(:,:)=szp(:,:)
      if (nqspec == 1) then
         zsqp(:,:)=sqp(:,:)
         call mpgallsp(zsq,zsqp,NLEV)
         call sp2fl(zsq,zq,NLEV)
         call fc2gp(zq,NLON,NLPP*NLEV)
      else
         zq(:,:) = 0.0
      endif
      call mpgallsp(zsd,zsdp,NLEV)
      call mpgallsp(zsz,zszp,NLEV)
      call dv2uv(zsd,zsz,zu,zv)
      call fc2gp(zu,NLON,NLPP*NLEV)
      call fc2gp(zv,NLON,NLPP*NLEV)
!
      zsdp(:,:)=sdp(:,:)+zsdt1(:,:)*delt2
      zszp(:,:)=szp(:,:)+zszt1(:,:)*delt2
      call mpgallsp(zsd,zsdp,NLEV)
      call mpgallsp(zsz,zszp,NLEV)
      call dv2uv(zsd,zsz,zun,zvn)
      call fc2gp(zun,NLON,NLPP*NLEV)
      call fc2gp(zvn,NLON,NLPP*NLEV)
!
      do jlev = 1 , NLEV
       zu(:,jlev)=cv*zu(:,jlev)*SQRT(rcsq(:))
       zv(:,jlev)=cv*zv(:,jlev)*SQRT(rcsq(:))
       zun(:,jlev)=cv*zun(:,jlev)*SQRT(rcsq(:))
       zvn(:,jlev)=cv*zvn(:,jlev)*SQRT(rcsq(:))
       if (nqspec == 1) zq(:,jlev)=zq(:,jlev)*psurf/dp(:)

       zdtdt(:,jlev)=-(zun(:,jlev)*zun(:,jlev)                          &
     &                -zu(:,jlev)*zu(:,jlev)                            &
     &                +zvn(:,jlev)*zvn(:,jlev)                          &
     &                -zv(:,jlev)*zv(:,jlev))/deltsec2                  &
     &               *0.5/acpd/(1.+adv*dq(:,jlev))
      enddo
!
      zdtdt(:,:)=zdtdt(:,:)/ct/ww
      call gp2fc(zdtdt,NLON,NLPP*NLEV)
      do jlev=1,NLEV
       call fc2sp(zdtdt(1,jlev),zstf1(1,jlev))
      enddo
      call mpsumsc(zstf1,zstt1,NLEV)
!
      zsdp(:,:)=sdp(:,:)+zsdt2(:,:)*delt2
      zszp(:,:)=szp(:,:)+zszt2(:,:)*delt2
      call mpgallsp(zsd,zsdp,NLEV)
      call mpgallsp(zsz,zszp,NLEV)
      call dv2uv(zsd,zsz,zun,zvn)
      call fc2gp(zun,NLON,NLPP*NLEV)
      call fc2gp(zvn,NLON,NLPP*NLEV)
!
      do jlev = 1 , NLEV
       zun(:,jlev)=cv*zun(:,jlev)*SQRT(rcsq(:))
       zvn(:,jlev)=cv*zvn(:,jlev)*SQRT(rcsq(:))
       zdekin(:,jlev)=(zun(:,jlev)*zun(:,jlev)                          &
     &                -zu(:,jlev)*zu(:,jlev)                            &
     &                +zvn(:,jlev)*zvn(:,jlev)                          &
     &                -zv(:,jlev)*zv(:,jlev))/deltsec2                  &
     &               *dp(:)/ga*dsigma(jlev)   
      enddo
      call gp2fc(zdekin,NLON,NLPP*NLEV)
      do jlev=1,NLEV
       call fc2sp(zdekin(1,jlev),zsdef(1,jlev))
      enddo
      call mpsumsc(zsdef,zsde,NLEV)
      call mpgallsp(zsdef,zsde,NLEV)
      zsdef(2:NESP,:)=0.
      call sp2fl(zsdef,zdekin,NLEV)
      call fc2gp(zdekin,NLON,NLPP*NLEV)
      do jlev=1,NLEV
       zdtdt(:,jlev)=-zdekin(:,jlev)                                    &
     &           *0.5/acpd/(1.+adv*dq(:,jlev))/dp(:)*ga/dsigma(jlev) 
      enddo
      zdtdt(:,:)=zdtdt(:,:)/ct/ww
      call gp2fc(zdtdt,NLON,NLPP*NLEV)
      do jlev=1,NLEV
       call fc2sp(zdtdt(1,jlev),zstf2(1,jlev))
      enddo
      call mpsumsc(zstf2,zstt2,NLEV)
!
!     energy diagnostics
!
      if(nenergy > 0) then
       call mpgallsp(zstf1,zstt1,NLEV)       
       zstf1(:,:)=zstf1(:,:)*ct*ww
       call sp2fl(zstf1,zdtdt,NLEV)
       call fc2gp(zdtdt,NLON,NLPP*NLEV)
       denergy(:,23)=0.
       do jlev=1,NLEV
        denergy(:,23)=denergy(:,23)                                     &
     &               +zdtdt(:,jlev)                                     &
     &               *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,23)=zdtdt(:,jlev)                               &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif
       enddo
       call mpgallsp(zstf2,zstt2,NLEV)
       zstf2(:,:)=zstf2(:,:)*ct*ww
       call sp2fl(zstf2,zdtdt,NLEV)
       call fc2gp(zdtdt,NLON,NLPP*NLEV)
       denergy(:,25)=0.
       do jlev=1,NLEV
        denergy(:,25)=denergy(:,25)                                     &
     &               +zdtdt(:,jlev)                                     &
     &               *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,25)=zdtdt(:,jlev)                               &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif
       enddo

      endif
!
      stt(:,:)=stt(:,:)+zstt1(:,:)+zstt2(:,:)
!
      return
      end subroutine mkdheat
