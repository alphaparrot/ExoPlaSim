!
!     **********************************************************
!     *  Parameters and subroutines for aerosol transport       *
!     **********************************************************

      module aeromod

!     **********************************************************
!     * This module contains the parameters, arrays and        *
!     * subroutines that are needed for transporting aerosols   *
!     * using the Flux-Form Semi-Lagrangian (FFSL) algorithm   *
!     * developed by S.-J. Lin (now at GFDL).                  *
!     **********************************************************
!     * The original transport code (in F77) was written by    *
!     * S.-J. Lin. Adaptation for the Planet Simulator was     *
!     * done by Hui Wan (MPI-M).
!     **********************************************************

      use pumamod

      logical,parameter :: aero_debug  = .FALSE.
      logical,parameter :: aero_zcross = .TRUE.
      logical,parameter :: aero_deform = .FALSE.

      logical,parameter :: aero_fill = .FALSE.
      logical,parameter :: aero_mfct = .FALSE.

      integer,parameter :: aero_iord = 2
      integer,parameter :: aero_jord = 2
      integer,parameter :: aero_kord = 3
      
      integer :: l_source = 1 ! 1 = photochemical haze (source at top level)
                                      ! 2 = dust (source at bottom level)
      integer :: l_bulk = 1 ! 1 = N2 atmosphere, 
                            ! 2 = H2 atmosphere

      integer,parameter :: aero_cnst = 1   ! 1 = constant preserving
                                           ! 2 = mass conserving

!      integer,parameter :: aero_j1  = 2  ! 1st lat. outside polar cap
!      integer,parameter :: aero_j2  = NLAT + 1 - aero_j1 
                                         ! last lat. outside polar cap 
      real :: apart = 50e-9 ! Radius of aerosol particle in m - DECLARED IN RADMOD AS WELL
      real :: rhop = 1000 ! Density of aerosol particle in kg/m3
      real :: fcoeff = 10e-13 ! Haze mass mixing ratio in kg/kg

!     Removal. Both terms are OFF by default, so one executable can run both
!     arms of an A/B and differ only by a namelist key. That is WORKFLOW.md
!     A3 and it is not a style preference: the low-I/O patch changes the
!     restart layout, so two arms built either side of a rebuild cannot share
!     a restart at all, and a term that cannot be switched off cannot be
!     tested.
!
!     ldepvel  0 = the legacy bottom-level sink, mmr = mmr*0.01 every step.
!                  It is a stability hack rather than a parameterisation: the
!                  implied deposition velocity is two orders of magnitude
!                  above any measured one, and because it does not scale with
!                  the timestep the implied velocity is inversely
!                  proportional to it. A rate that depends on the integration
!                  step is not a rate.
!              1 = a deposition velocity. The bottom layer keeps
!                  exp(-vdaero*dt/dz) per step, with dz built from the
!                  layer's own pressure thickness and gas density, so the
!                  rate is a property of the atmosphere and not of the step.
!     vdaero   the NON-gravitational dry deposition velocity, m/s: turbulent
!              transfer and impaction. Sedimentation to the surface is
!              already carried by the settling flux leaving the bottom layer,
!              so putting a settling term here as well would count it twice.
!              There is no default: aero_ini refuses ldepvel = 1 with vdaero
!              at or below zero rather than invent one, because the value
!              belongs with the rest of the removal budget in
!              the run's own configuration.
!     lwetdep  0 = no wet removal, which is what ExoPlaSim has and what the
!                  LMD Generic PCM has as well.
!              1 = below-cloud scavenging at Lambda = scava * p**scavb, with
!                  p the total precipitation rate in mm/h. That is the
!                  Sportisse (2007) form an offline chain
!                  already uses, so the two chains agree in form and not only
!                  in magnitude.
!     scava    that A, in s-1 per (mm/h)**B. No default, same reasoning as
!              vdaero; the run's own configuration carries it and its bracket.
!     scavb    that B, dimensionless. No default.
      integer :: ldepvel = 0
      integer :: lwetdep = 0
      real :: vdaero = 0.0
      real :: scava  = 0.0
      real :: scavb  = 0.0

      end module aeromod

!     ==================
!     SUBROUTINE AERO_INI
!     ==================

      subroutine aero_ini
      use aeromod
      use radmod, only: l_aerorad, aerofile
      
      namelist/aero_nl/l_source,l_bulk,apart,rhop,fcoeff,l_aerorad,aerofile  &
     &                ,ldepvel,vdaero,lwetdep,scava,scavb

      if (mypid==NROOT) then
         open(11,file=aero_namelist)
         read(11,aero_nl)
         close(11)
         write(nud,'(/," *********************************************")')
         write(nud,'(" * AEROMOD ",a34)')
         write(nud,'(" *********************************************")')
         write(nud,'(" * Namelist AERO_NL from <aero_namelist> *")')
         write(nud,'(" *********************************************")')
         write(nud,aero_nl)
!
!        Refuse a switch that is on without the coefficient it needs, rather
!        than carry a plausible-looking default that nothing in the run
!        decided. Both values live in the run's own configuration.
!
         if (ldepvel == 1 .and. vdaero <= 0.0) then
            write(nud,*) '* ldepvel = 1 needs a positive vdaero in m/s.'
            write(nud,*) '* It is the non-gravitational dry deposition'
            write(nud,*) '* velocity; settling is already carried by the'
            write(nud,*) '* sedimentation flux. There is no defensible'
            write(nud,*) '* default, so it must be set explicitly.'
            call mpabort('aero_nl: ldepvel = 1 without vdaero')
         endif
         if (lwetdep == 1 .and. (scava <= 0.0 .or. scavb <= 0.0)) then
            write(nud,*) '* lwetdep = 1 needs positive scava and scavb.'
            write(nud,*) '* Lambda = scava * p**scavb with p in mm/h.'
            write(nud,*) '* Both are aerosol- and rain-specific and have no'
            write(nud,*) '* defensible default, so both must be set.'
            call mpabort('aero_nl: lwetdep = 1 without scava and scavb')
         endif
      endif
      
      return
      end subroutine aero_ini
  
  
!     ======================
!     SUBROUTINE AERO_MAIN
!     ======================

      subroutine aero_main

      use pumamod, only: du,dv,dp,du0,dv0,dp0,daeros,numrhos, &
                         NLON,NLAT,NLEV,NAERO,NHOR,   &
                         mypid,NROOT,sigmah,dt,dls,dswfl,dprl,dprc
      use tracermod
      use aeromod
      use radmod, only: gmu0, l_aerorad ! Use cosine of solar zenith angle from radmod;

      implicit none

      real :: zu   (NLON,NLAT,NLEV)
      real :: zv   (NLON,NLAT,NLEV)
      real :: zps0 (NLON,NLAT)
      real :: zps1 (NLON,NLAT)

      real :: x (NLON+1,NLAT,NLEV,NAERO)  ! for GUI output
      real :: y (NLON+1,NLAT,NLEV)         ! for GUI output
      real ::   angle(NLON,NLAT) ! Array for cosine of solar zenith angle
      real ::   aerosw(NLON,NLAT,NLEV) ! Array for SW flux 
      real ::   land(NLON,NLAT) ! Array for binary land mask

!     Gather buffer. mpgagp returns a field in the MODEL's latitude order,
!     while daeros and numrhos are handed to aerocore flipped south-to-north
!     (see plasim.f90, `daeros(:,NLAT+1-jlat,:,1) = zmmr(:,jlat,:)`). Every
!     field gathered here therefore has to be flipped the same way before it
!     is used against them. Upstream did not, so the land mask and the solar
!     zenith angle drove the aerosol source in the wrong hemisphere.
      real ::   zgath(NLON,NLAT,NLEV)

      real ::   prec(NLON,NLAT) ! Total precipitation rate (m/s), for lwetdep
      real ::   zprec(NHOR)     ! the same before gathering

      integer :: j,jc

      character(len=9) :: aero_name

!     --- 

      call prepare_uvps( zu,zv,zps0,zps1,      & ! output
                         du0,dv0,dp0,du,dv,dp)   ! input

      if (l_aerorad == 0) then ! No radiative transfer
       select case (l_source) ! Choose your aerosol source
       case(1) ! Case 1: photochemical haze
         call solang ! Use subroutine from radmod to calculate solar zenith angle
         call mpgagp(zgath,gmu0,1) ! Gather from nodes
         do j=1,NLAT
            angle(:,NLAT+1-j) = zgath(:,j,1)
         end do
       case(2) ! Case 2: dust
         call mpgagp(zgath,dls,1) ! Import land-sea mask from landmod and reshape to match grid size
         do j=1,NLAT
            land(:,NLAT+1-j) = zgath(:,j,1)
         end do
       end select
      end if
      
      if (l_aerorad == 1) then ! Include radiative transfer
       select case (l_source) ! Choose aerosol source
       case(1) ! Case 1: photochemical haze     
        call mpgagp(zgath,dswfl,NLEV) ! Gather SW flux from nodes
        do j=1,NLAT
           aerosw(:,NLAT+1-j,:) = zgath(:,j,:)
        end do
       case(2) ! Case 2: dust
        call mpgagp(zgath,dls,1) ! Import land-sea mask from landmod and reshape to match grid size
        do j=1,NLAT
           land(:,NLAT+1-j) = zgath(:,j,1)
        end do
       end select
      end if 

!     Precipitation for the wet-scavenging term, large scale plus convective,
!     in m/s. Gathered only when the term is on, so that lwetdep = 0 costs
!     nothing and reproduces the unpatched model exactly. Flipped in latitude
!     like every other field gathered here.

      prec(:,:) = 0.0
      if (lwetdep == 1) then
         zprec(:) = dprl(:) + dprc(:)
         call mpgagp(zgath,zprec,1)
         do j=1,NLAT
            prec(:,NLAT+1-j) = zgath(:,j,1)
         end do
      end if

      if (mypid == NROOT .and. aero_debug) then
         write(nud,'(a,f11.2)') '* max aero u   =',maxval(abs(zu))
         write(nud,'(a,f11.2)') '* max v   =',maxval(abs(zv))
         write(nud,'(a,f11.2)') '* max ps0 =',maxval(zps0)
         write(nud,'(a,f11.2)') '* max ps1 =',maxval(zps1)
       ! write(nud,*)
       ! write(nud,*) 'ps0 NP'
       ! write(nud,*)  dp0(1:NLON)
       ! write(nud,*) 'zps0 NP'
       ! write(nud,*) zps0(1:NLON,NLAT)
       ! write(nud,*) 'ps0 SP'
       ! write(nud,*)  dp0(NLON*(NLAT-1)+1:)
       ! write(nud,*) 'zps0 SP'
       ! write(nud,*) zps0(1:NLON,1)
      end if
    
      if (mypid == NROOT) then

         call aerocore(daeros,numrhos,l_source,sigmah,dt,         &
                      zps0,zps1,zu,zv,                    &
                      dtoa,dtdx,apart,rhop,fcoeff,        &
                      aero_iord,aero_jord,aero_kord,      & 
                      NAERO,NLON,NLAT,NLEV,dap,dbk,       &
                      iml,ffsl_j1,ffsl_j2,js0,jn0,        &
                      colae,colad,rcolad,dlat,rcap,       &
                      aero_cnst,aero_deform,aero_zcross,  &
                      aero_fill,aero_mfct,aero_debug,nud, &
                      angle,land,aerosw,l_aerorad,prec)

!        preparation for the GUI output: 
!        invert the meridional direction and add the 360 deg. longitude

         do j=1,NLAT
            x(1:NLON,j,:,:) = daeros(:,NLAT+1-j,:,:)
            x(NLON+1,j,:,:) = daeros(1,NLAT+1-j,:,:)
         end do

!        send all tracer fields to output

         do jc=1,NAERO
            write(aero_name,'(a,i2.2)') 'DAEROS',jc
            call guiput(aero_name // char(0), x(1,1,1,jc), NLON+1,NLAT,NLEV)
         enddo

!        check the correlation between tracers 3 and 4

!        y(:,:,:) = x(:,:,:,3)+x(:,:,:,4)
!        call guiput('TRC03+04' // char(0), y, NLON+1,NLAT,NLEV)

      end if

      return
      end subroutine aero_main
