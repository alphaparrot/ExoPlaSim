
#define limitwater 1

      module landmod
      use pumamod
!
!     version identifier (date)
!
      character(len=80) :: lversion = '19.09.2019 by Adiv'
!
!     parameters
!
      parameter(NLSOIL=5)
      parameter(WSMAX_EARTH = 0.5) ! Initial value vor Earth
      parameter(WSMAX_MARS  = 0.0) ! Initial value for Mars
!
!     namelist parameters
!
      integer :: nlandt   = 1     ! switch for land model (1/0 : prog./clim)
      integer :: nlandw   = 1     ! switch for soil model (1/0 : prog./clim)
      integer :: newsurf  = 0     ! (dtcl,dwcl) 1: update from file, 2:reset 
      integer :: nwatcini = 0     ! (0/1) initialize water content of soil
      integer :: nwetsoil = 0     ! (0/1) Soil albedo responds to water content
      real    :: alblandnl  = 0.2   ! albedo for land
      real    :: albland  = 0.22
      
      real  :: dgroundalbnl(2)  = 0.22
                   ! Non-spectral versions
      real    :: albsmin  = 0.4   ! min. albedo for snow
      real    :: albsmax  = 0.8   ! max. albedo for snow
      real    :: albsminf = 0.3   ! min. albedo for snow (with forest)
      real    :: albsmaxf = 0.4   ! max. albedo for snow (with forest)
      real    :: albgmin  = 0.6   ! min. albedo for glaciers
      real    :: albgmax  = 0.8   ! max. albedo for glaciers
      real    :: alblandmax = 0.40  ! max land albedo; i.e. desert sand
                   ! lambda < 0.75 microns
      real    :: albsmin1  = 0.4   ! min. albedo for snow
      real    :: albsmax1  = 0.8   ! max. albedo for snow
      real    :: albsminf1 = 0.3   ! min. albedo for snow (with forest)
      real    :: albsmaxf1 = 0.4   ! max. albedo for snow (with forest)
      real    :: albgmin1  = 0.6   ! min. albedo for glaciers
      real    :: albgmax1  = 0.8   ! max. albedo for glaciers
                   ! lambda > 0.75 microns
      real    :: albsmin2  = 0.4   ! min. albedo for snow
      real    :: albsmax2  = 0.8   ! max. albedo for snow
      real    :: albsminf2 = 0.3   ! min. albedo for snow (with forest)
      real    :: albsmaxf2 = 0.4   ! max. albedo for snow (with forest)
      real    :: albgmin2  = 0.6   ! min. albedo for glaciers
      real    :: albgmax2  = 0.8   ! max. albedo for glaciers
      
      real    :: dz0land  = 2.0   ! roughness length land
      real    :: drhsland = 0.25  ! wetness factor land
      real    :: drhsfull = 0.4   ! threshold above which drhs=1 [frac. of wsmax]
      real    :: dzglac   = -1.   ! threshold of orography to be glacier (-1=none)
      real    :: dztop    = 0.20  ! thickness of the uppermost soil layer (m)
      real    :: dsmax    = 5.00  ! maximum snow depth (m-h20; -1 = no limit)

      
      real    :: wsmax    = WSMAX_EARTH ! max field capacity of soil water (m)
      real    :: dwatcini = 0           ! water content of soil (m)

!     SIMBA - fixed parameters

      real    :: rlue     =  3.4E-10 ! Recommended by Pablo Paiewonsky
      real    :: co2conv  =  8.3E-4  ! Recommended by Pablo Paiewonsky
      real    :: tau_veg  = 10.0  ! [years] - in landini scaled to seconds
      real    :: tau_soil = 42.0  ! [years] - in landini scaled to seconds
      real    :: rinifor  =  0.5
      real    :: rnbiocats=  0.0
!
!     global scalars (snow similar to the sea ice module)
!
      real :: rhosnow  = 330.    ! snow density (kg/m**3)
      real :: soildiff = 1.8     ! heat diffusivity of the soil (W/m/K)
      real :: sicediff = 2.03    ! heat diffusivity of ice      (W/m/K)
      real :: snowdiff = 0.31    ! heat diffusivity of snow     (W/m/K)
      real :: soilcap  = 2.4E6   ! heat capacity of the soil  (J/m**3/K)
      real :: sicecap  = 2.07E6  ! heat capacity of ice       (J/m**3/K)
      real :: snowcap  = 0.6897E6! heat capacity of snow      (J/m**3/K)
!
!     global arrays
!
!     a) surface definitions
!
      real :: doro(NHOR) = 0.0     ! orography (m2/s2)
      real :: dts(NHOR)  = 1.0e20  ! surface temperature (K)
      real :: dtsm(NHOR) = 1.0e20  ! surface temperature (K)
      real :: dqs(NHOR)  = 1.0e20  ! surface humidity    (kg/kg)
!
!     b) runoff
!
      real :: driver(NHOR) = 0.0   ! surface water (river runoff) (m)
      real :: duroff(NHOR) = 0.0   ! zonal runoff-velocity  (1/s)
      real :: dvroff(NHOR) = 0.0   ! meridional runoff-velocity  (1/s)
      real :: darea(NHOR)  = 1.0   ! area weights
!
!     c) soil
!
      real :: dsoilz(NLSOIL)=(/0.4,0.8,1.6,3.2,6.4/)  ! soil layer thickness (m)
      real :: dsoilt(NHOR,NLSOIL) = 0.0    ! soil temperatur (K)
      real :: dsnowt(NHOR)        = 0.0    ! snow temperatur (K)
      real :: dtclsoil(NHOR)      = 0.0    ! clim soil temp. (initilization) (K)
      real :: dsnowz(NHOR)        = 0.0    ! snow depth (m water equivalent)
      real :: dwater(NHOR)        = 0.0    ! surface water for soil (m/s)
!
!     e) climatological surface
!
      real :: dtcl(NHOR,0:13)   =  0.0  ! climatological surface temperature
      real :: dwcl(NHOR,0:13)   = -1.0  ! climatological soil wetness
      real :: dalbcl(NHOR,0:13) =  0.22  ! climatological background albedo
      real :: dalbcl1(NHOR,0:13) = 0.22  ! climatological background albedo (<.75 um)
      real :: dalbcl2(NHOR,0:13) = 0.22  ! climatological background albedo (>.75 um)
      real :: dtclim(NHOR)      =  0.0  ! climatological surface temperature
      real :: dwclim(NHOR)      =  0.0  ! climatological soil wetness
      real :: dz0clim(NHOR)     =  2.0  ! climatological z0  (total)
      real :: dz0climo(NHOR)    =  0.0  ! climatological z0  (from topograhpy only)
      real :: dalbclim(NHOR)    =  2.0  ! climatological background albedo
      real :: dalbclim1(NHOR)   =  2.0  ! climatological background albedo (<.75 um)
      real :: dalbclim2(NHOR)   =  2.0  ! climatological background albedo (>.75 um)
!
      end module landmod

!     ==================
!     FUNCTION NCOUNTSEA
!     ==================

      function ncountsea(plsm)
      use pumamod
      real :: plsm(NHOR)

      call mpsumval(plsm,NHOR,1,zsum)
      ilpo    = nint(zsum)
      ispo    = NUGP - ilpo
      ilperc  = nint((100.0 * zsum) / NUGP)
      isperc  = 100 - ilperc
  
      if (mypid == NROOT) then
         write(nud,'(a,i6,a,i6,a,i3,a)') &
              ' Land:',ilpo,' from',NUGP,' = ',ilperc,'%'
         write(nud,'(a,i6,a,i6,a,i3,a)') &
              ' Sea: ',ispo,' from',NUGP,' = ',isperc,'%'
      endif
      ncountsea = ispo

      return
      end

!     ==================
!     SUBROUTINE LANDINI
!     ==================

      subroutine landini
      use landmod
      use radmod
!
!     initialize land surface
!
      namelist/landmod_nl/nlandt,nlandw,albland,dz0land,drhsland        &
     &                ,dsnowalbmn,dsnowalbmx,dglacalbmn,dsnowalb        &
     &                ,dsmax,wsmax,drhsfull,dzglac,dztop,dsoilz         &
     &                ,rlue,co2conv,tau_veg,tau_soil                    &
     &                ,rnbiocats,nwetsoil,soilcap                       &
     &                ,newsurf,rinifor,nwatcini,dwatcini,dgroundalb
!
      dtclsoil(:) = tmelt
      dsoilt(:,:) = tmelt
      dsnowt(:)   = tmelt
      dtcl(:,:)   = tmelt
      dtclim(:)   = tmelt
      

      if (mars == 1) then
         wsmax = WSMAX_MARS
!        nlandt = 0
!        nlandw = 0
!        dtclsoil(:) = tmelt_CO2
!        dsoilt(:,:) = tmelt_CO2
!        dsnowt(:)   = tmelt_CO2
      endif

      if (ndesert == 1 .and. nrestart == 0) then
         dwatcini = 0.0
         dwater(:) = 0.0
      endif
      
!       if (mypid==NROOT) then
!         
!        albsmax1 = dsnowalbmx(1)
!        albgmax1 = dsnowalbmx(1)
!        albsmax2 = dsnowalbmx(2)
!        albgmax2 = dsnowalbmx(2)
!        albsmaxf1 = 0.5*albsmax1
!        albsmaxf2 = 0.5*albsmax2
!        albsmin1 = dsnowalbmn(1)
!        albsmin2 = dsnowalbmn(2)
!        albgmin1 = dglacalbmn(1)
!        albgmin2 = dglacalbmn(2)
!        albsminf1 = 0.75*albsmaxf1
!        albsminf2 = 0.75*albsmaxf2
!        
!       endif
      
      if (mypid == NROOT) then
      open(12,file=landmod_namelist)
      read(12,landmod_nl)
      write(nud,'(/,"***********************************************")')
      write(nud,'("* LANDMOD ",a35," *")') trim(lversion)
      write(nud,'("***********************************************")')
      write(nud,'("* Namelist LANDMOD_NL from <",a17,"> *")') &
            landmod_namelist
      write(nud,'("***********************************************")')
      write(nud,landmod_nl)
      close(12)
      
      albsmax1 = dsnowalbmx(1)
      albgmax1 = dsnowalbmx(1)
      albsmax2 = dsnowalbmx(2)
      albgmax2 = dsnowalbmx(2)
      albsmaxf1 = 0.5*albsmax1
      albsmaxf2 = 0.5*albsmax2
      albsmin1 = dsnowalbmn(1)
      albsmin2 = dsnowalbmn(2)
      albgmin1 = dglacalbmn(1)
      albgmin2 = dglacalbmn(2)
      albsminf1 = 0.75*albsmaxf1
      albsminf2 = 0.75*albsmaxf2
      
      endif

      if (wsmax < 0.0) wsmax = 0.0 ! Catch user error

      
      call mpbci(nlandt)
      call mpbci(nlandw)
      call mpbci(newsurf)
      call mpbci(nwetsoil)
      call mpbci(nwatcini)
      call mpbcr(albland)
      
      call mpbcr(albsmin)
      call mpbcr(albsmax)
      call mpbcr(albsminf)
      call mpbcr(albsmaxf)
      call mpbcr(albgmin)
      call mpbcr(albgmax)
      
      call mpbcr(albsmin1)
      call mpbcr(albsmax1)
      call mpbcr(albsminf1)
      call mpbcr(albsmaxf1)
      call mpbcr(albgmin1)
      call mpbcr(albgmax1)
      
      call mpbcr(albsmin2)
      call mpbcr(albsmax2)
      call mpbcr(albsminf2)
      call mpbcr(albsmaxf2)
      call mpbcr(albgmin2)
      call mpbcr(albgmax2)
      
      call mpbcr(dz0land)
      call mpbcr(drhsland)
      call mpbcr(drhsfull)
      call mpbcr(wsmax)
      call mpbcr(dwatcini)
      call mpbcr(dzglac)
      call mpbcr(dztop)
      call mpbcr(dsmax)
      call mpbcr(rlue)
      call mpbcr(co2conv)
      call mpbcr(tau_veg)
      call mpbcr(tau_soil)
      call mpbcr(rinifor)
      call mpbcr(rnbiocats)
      call mpbcrn(dsoilz,NLSOIL)
      
      call mpbcrn(dsnowalbmn,2)
      call mpbcrn(dsnowalbmx,2)
      call mpbcrn(dglacalbmn,2)
      call mpbcrn(dsnowalb,2)
      call mpbcrn(dgroundalb,2)
      

!     scale taus from years to seconds

      tau_veg  = tau_veg  * m_days_per_year * day_24hr
      tau_soil = tau_soil * m_days_per_year * day_24hr

      if (tau_veg < 1.0 .or. tau_soil < 1.0) then
         write(nud,*)' *** error: tau_veg = ',tau_veg,'  tau_soil = ',tau_soil
         stop 1
       endif

!     preset

      if (nrestart == 0) then
         dforest(:)= rinifor
      endif

      if (nrestart == 0) then
!
!*       preset some fields
!
         dz0clim(:)  = dz0land
         dwmax(:)    = wsmax
         dalbcl(:,:) = albland
         dalbcl1(:,:) = dgroundalb(1)
         dalbcl2(:,:) = dgroundalb(2)
!
!*       read surface parameters
!
         call mpsurfgp('doro'    ,doro    ,NHOR,1)
         call mpsurfgp('dls'     ,dls     ,NHOR,1)
         call mpsurfgp('dz0clim' ,dz0clim ,NHOR,1)
         call mpsurfgp('dz0climo',dz0climo,NHOR,1)
         call mpsurfgp('dglac'   ,dglac   ,NHOR,1)
         call mpsurfgp('dforest' ,dforest ,NHOR,1)
         call mpsurfgp('dwmax'   ,dwmax   ,NHOR,1)
         call mpsurfgp('dtclsoil',dtclsoil,NHOR,1)
   
         call mpsurfgp('dtcl',dtcl,NHOR,14)
         call mpsurfgp('dwcl',dwcl,NHOR,14)
         call mpsurfgp('dalbcl',dalbcl,NHOR,14)
         call mpsurfgp('dalbcl1',dalbcl1,NHOR,14)
         call mpsurfgp('dalbcl2',dalbcl2,NHOR,14)

!        make sure, that dwmax is positive

         where (dwmax(:) < 0.0) dwmax(:) = 0.0

         if (dwcl(1,1) < 0.0) then ! not in file
            do jm = 0 , 13
               dwcl(:,jm) = dwmax(:) * drhsfull * drhsland
            enddo
         endif

!        make sure, that land sea mask values are 0 or 1

         where (dls(:) > 0.5)
            dls(:) = 1.0
         elsewhere
            dls(:) = 0.0
         endwhere
         n_sea_points = ncountsea(dls)
!
!*       modify glacier mask according to dzglac
!

         if (dzglac > 0.0) then
            where (doro(:) > dzglac * ga) dglac(:)=1.0
         endif
!
!*       convert fractional glacier mask to binary mask
!
         where (dglac(:) > 0.5)
            dglac(:) = 1.0
         elsewhere
            dglac(:) = 0.0
         endwhere
!
!*    initialize soil
!
       call soilini
!
!*    initialize runoff (except for land only simulations)
!
       if (n_sea_points > 0) call roffini
!
!*    get new background albedo 
!
       call getalb
!
!*    set other surface variables
!
       do jhor=1,NHOR
        if(dls(jhor) > 0.0) then
         dtsm(jhor)=dts(jhor)
         dqs(jhor)=rdbrv*ra1*EXP(ra2*(dts(jhor)-tmelt)/(dts(jhor)-ra4)) &
     &            /psurf
         dqs(jhor)=dqs(jhor)/(1.-(1./rdbrv-1.)*dqs(jhor))
         dsnow(jhor)=dsnowz(jhor)
         if(dsnow(jhor) > 0.) then
          zalbmax=dforest(jhor)*albsmaxf+(1.-dforest(jhor))*albsmax
          zalbmin=dforest(jhor)*albsminf+(1.-dforest(jhor))*albsmin
          zdalb=(zalbmax-zalbmin)*(dts(jhor)-263.16)/(tmelt-263.16)
          zalbsnow=MAX(zalbmin,MIN(zalbmax,zalbmax-zdalb))
          zalbmax1=dforest(jhor)*albsmaxf1+(1.-dforest(jhor))*albsmax1
          zalbmin1=dforest(jhor)*albsminf1+(1.-dforest(jhor))*albsmin1
          zdalb1=(zalbmax1-zalbmin1)*(dts(jhor)-263.16)/(tmelt-263.16)
          zalbsnow1=MAX(zalbmin1,MIN(zalbmax1,zalbmax1-zdalb1))
          zalbmax2=dforest(jhor)*albsmaxf2+(1.-dforest(jhor))*albsmax2
          zalbmin2=dforest(jhor)*albsminf2+(1.-dforest(jhor))*albsmin2
          zdalb2=(zalbmax2-zalbmin2)*(dts(jhor)-263.16)/(tmelt-263.16)
          zalbsnow2=MAX(zalbmin2,MIN(zalbmax2,zalbmax2-zdalb2))
          dalb(jhor)=dalbclim(jhor)                                     &
     &        +(zalbsnow-dalbclim(jhor))*dsnow(jhor)/(dsnow(jhor)+0.01)
          dsalb(1,jhor) = dalbclim1(jhor)                               &
     &        +(zalbsnow1-dalbclim1(jhor))*dsnow(jhor)/(dsnow(jhor)+0.01)
          dsalb(2,jhor) = dalbclim2(jhor)                               &
     &        +(zalbsnow2-dalbclim2(jhor))*dsnow(jhor)/(dsnow(jhor)+0.01)
          drhs(jhor)=1.
         else
          dalb(jhor)=dalbclim(jhor)
          dsalb(1,jhor)=dalbclim1(jhor)
          dsalb(2,jhor)=dalbclim2(jhor)
          if (dwmax(jhor) > 0.0)                                        &
          drhs(jhor)=AMIN1(1.,dwatc(jhor)/(drhsfull*dwmax(jhor)))
         endif
         dz0(jhor)=dz0clim(jhor)
!
!        diagnostics: soil temps.
!
         dtsoil(jhor)=dsoilt(jhor,1)
         if(NLSOIL > 2) dtd2(jhor)=dsoilt(jhor,2)
         if(NLSOIL > 3) dtd3(jhor)=dsoilt(jhor,3)
         if(NLSOIL > 4) dtd4(jhor)=dsoilt(jhor,4)
         dtd5(jhor)=dsoilt(jhor,NLSOIL)

!
!*       modifications according to glacier mask
!

         if(dglac(jhor) > 0.5) then
          dsnowz(jhor)=AMAX1(dsmax,0.)
          dsnow(jhor)=dsnowz(jhor)
          zdalb=(albgmax-albgmin)*(dts(jhor)-263.16)/(tmelt-263.16)
          zdalb1=(albgmax1-albgmin1)*(dts(jhor)-263.16)/(tmelt-263.16)
          zdalb2=(albgmax2-albgmin2)*(dts(jhor)-263.16)/(tmelt-263.16)
          dalb(jhor)=MAX(albgmin,MIN(albgmax,albgmax-zdalb))
          dsalb(1,jhor)=MAX(albgmin1,MIN(albgmax1,albgmax1-zdalb1))
          dsalb(2,jhor)=MAX(albgmin2,MIN(albgmax2,albgmax2-zdalb2))
          drhs(jhor)=1.0
         end if

!
!*    set PUMA temperature and moisture:
!

         dt(jhor,NLEP)=dts(jhor)
         dq(jhor,NLEP)=dqs(jhor)

        endif
       enddo

      else ! restart > 0

       call mpgetgp('dtsl'    ,dts     ,NHOR,     1)
       call mpgetgp('dtsm'    ,dtsm    ,NHOR,     1)
       call mpgetgp('dqs'     ,dqs     ,NHOR,     1)
       call mpgetgp('driver'  ,driver  ,NHOR,     1)
       call mpgetgp('duroff'  ,duroff  ,NHOR,     1)
       call mpgetgp('dvroff'  ,dvroff  ,NHOR,     1)
       call mpgetgp('darea'   ,darea   ,NHOR,     1)
       call mpgetgp('dwmax'   ,dwmax   ,NHOR,     1)
       call mpgetgp('dtcl'    ,dtcl    ,NHOR,    14)
       call mpgetgp('dwcl'    ,dwcl    ,NHOR,    14)
       call mpgetgp('dsnowt'  ,dsnowt  ,NHOR,     1)
       call mpgetgp('dsnowz'  ,dsnowz  ,NHOR,     1)
       call mpgetgp('dsoilt'  ,dsoilt  ,NHOR,NLSOIL)
       call mpgetgp('dz0clim' ,dz0clim ,NHOR,     1)
       call mpgetgp('dz0climo',dz0climo,NHOR,     1)
       call mpgetgp('dalbcl'  ,dalbcl  ,NHOR,    14)
       call mpgetgp('dalbcl1' ,dalbcl1 ,NHOR,    14)
       call mpgetgp('dalbcl2' ,dalbcl2 ,NHOR,    14)

       n_sea_points = ncountsea(dls)
       
#if limitwater == 1
!	If you had created a restart file, but then need to change the soil water
!	capacity later.
       dwcl(:,:) = wsmax * drhsfull * drhsland
       do jhor=1,NHOR
          dwmax(jhor) = AMAX1(wsmax,0.0)
          if (dsnow(jhor) .le. 0.) then
             drhs(jhor) = 0.
             if (dwmax(jhor) > 0.) &
     &           drhs(jhor) = dwatc(jhor)/(drhsfull*dwmax(jhor))
          endif
       enddo
#endif

      endif ! if (restart == 0)

      if (newsurf == 1) then
         call mpsurfgp('dtcl',dtcl,NHOR,14)
         call mpsurfgp('dwcl',dwcl,NHOR,14)
      endif

      if (newsurf == 2) then ! preset some fields
         dwmax(:)    = wsmax
         dz0clim(:)  = dz0land
         dalbcl(:,:) = albland
         dalbcl1(:,:) = dgroundalb(1)
         dalbcl2(:,:) = dgroundalb(2)
         dwcl(:,:)   = wsmax * drhsfull * drhsland
      endif

      return
      end subroutine landini


!     ===================
!     SUBROUTINE LANDSTEP
!     ===================

      subroutine landstep
      use landmod

!
!     get climatological values if t and/or w are non interactive
!

      if(nlandt==0 .or. nlandw==0 ) call gettcll

!
!     soil
!

      call soilstep

!
!     runoff
!

      call roffstep
!
!     get new background albedo
!
      call getalb
!
!     set surface variables
!
      do jhor=1,NHOR
       if(dls(jhor) > 0.0) then
        dtsm(jhor)=dts(jhor)
        dqs(jhor)=rdbrv*ra1*EXP(ra2*(dts(jhor)-tmelt)/(dts(jhor)-ra4))  &
     &           /dp(jhor)
        dqs(jhor)=dqs(jhor)/(1.-(1./rdbrv-1.)*dqs(jhor))
        dsnow(jhor)=dsnowz(jhor)
        if(dsnow(jhor) > 0.) then
         zalbmax=dforest(jhor)*albsmaxf+(1.-dforest(jhor))*albsmax
         zalbmin=dforest(jhor)*albsminf+(1.-dforest(jhor))*albsmin
         zdalb=(zalbmax-zalbmin)*(dts(jhor)-263.16)/(tmelt-263.16)
         zalbsnow=MAX(zalbmin,MIN(zalbmax,zalbmax-zdalb))
         zalbmax1=dforest(jhor)*albsmaxf1+(1.-dforest(jhor))*albsmax1
         zalbmin1=dforest(jhor)*albsminf1+(1.-dforest(jhor))*albsmin1
         zdalb1=(zalbmax1-zalbmin1)*(dts(jhor)-263.16)/(tmelt-263.16)
         zalbsnow1=MAX(zalbmin1,MIN(zalbmax1,zalbmax1-zdalb1))
         zalbmax2=dforest(jhor)*albsmaxf2+(1.-dforest(jhor))*albsmax2
         zalbmin2=dforest(jhor)*albsminf2+(1.-dforest(jhor))*albsmin2
         zdalb2=(zalbmax2-zalbmin2)*(dts(jhor)-263.16)/(tmelt-263.16)
         zalbsnow2=MAX(zalbmin2,MIN(zalbmax2,zalbmax2-zdalb2))
         dalb(jhor)=dalbclim(jhor)                                     &
     &       +(zalbsnow-dalbclim(jhor))*dsnow(jhor)/(dsnow(jhor)+0.01)
         dsalb(1,jhor) = dalbclim1(jhor)                               &
     &       +(zalbsnow1-dalbclim1(jhor))*dsnow(jhor)/(dsnow(jhor)+0.01)
         dsalb(2,jhor) = dalbclim2(jhor)                               &
     &       +(zalbsnow2-dalbclim2(jhor))*dsnow(jhor)/(dsnow(jhor)+0.01)
         drhs(jhor)=1.
        else
         dalb(jhor)=dalbclim(jhor)
         dsalb(1,jhor)=dalbclim1(jhor)
         dsalb(2,jhor)=dalbclim2(jhor)
         if (dwmax(jhor) > 0.0)                                        &
         drhs(jhor)=AMIN1(1., dwatc(jhor)/(drhsfull *dwmax(jhor)))
        endif

!
!     diagnostics: soil temps.
!

        dtsoil(jhor)=dsoilt(jhor,1)
        if(NLSOIL > 2) dtd2(jhor)=dsoilt(jhor,2)
        if(NLSOIL > 3) dtd3(jhor)=dsoilt(jhor,3)
        if(NLSOIL > 4) dtd4(jhor)=dsoilt(jhor,4)
        dtd5(jhor)=dsoilt(jhor,NLSOIL)

!
!*    set PUMA temperature and moisture:
!

        dt(jhor,NLEP)=dts(jhor)
        dq(jhor,NLEP)=dqs(jhor)

       end if
      enddo

!
!*    vegetation
!

      if (nveg > 0) call vegstep

!
!*     modifications according to glacier mask
!

      where(dglac(:) > 0.5 .and. dls(:) > 0.0)
       dalb(:)=MAX(albgmin,MIN(albgmax                                  &
     &  ,albgmax-(albgmax-albgmin)*(dts(:)-263.16)/(tmelt-263.16)))
       dsalb(1,:)=MAX(albgmin1,MIN(albgmax1                                  &
     &  ,albgmax1-(albgmax1-albgmin1)*(dts(:)-263.16)/(tmelt-263.16)))
       dsalb(2,:)=MAX(albgmin2,MIN(albgmax2                                  &
     &  ,albgmax2-(albgmax2-albgmin2)*(dts(:)-263.16)/(tmelt-263.16)))
       drhs(:)=1.0
      end where

      return
      end subroutine landstep

!     ===================
!     SUBROUTINE LANDSTOP
!     ===================

      subroutine landstop
      use landmod

      if (mypid == NROOT) then
         call put_restart_integer('nlsoil',NLSOIL)
      endif

      call mpputgp('dtsl'    ,dts     ,NHOR, 1)
      call mpputgp('dtsm'    ,dtsm    ,NHOR, 1)
      call mpputgp('dqs'     ,dqs     ,NHOR, 1)
      call mpputgp('driver'  ,driver  ,NHOR, 1)
      call mpputgp('duroff'  ,duroff  ,NHOR, 1)
      call mpputgp('dvroff'  ,dvroff  ,NHOR, 1)
      call mpputgp('darea'   ,darea   ,NHOR, 1)
      call mpputgp('dwmax'   ,dwmax   ,NHOR, 1)
      call mpputgp('dtcl'    ,dtcl    ,NHOR,14)
      call mpputgp('dwcl'    ,dwcl    ,NHOR,14)
      call mpputgp('dsnowt'  ,dsnowt  ,NHOR, 1)
      call mpputgp('dsnowz'  ,dsnowz  ,NHOR, 1)
      call mpputgp('dsoilt'  ,dsoilt  ,NHOR,NLSOIL)
      call mpputgp('dz0clim' ,dz0clim ,NHOR, 1)
      call mpputgp('dz0climo',dz0climo,NHOR, 1)
      call mpputgp('dalbcl'  ,dalbcl  ,NHOR,14)
      call mpputgp('dalbcl1' ,dalbcl1 ,NHOR,14)
      call mpputgp('dalbcl2' ,dalbcl2 ,NHOR,14)
      return
      end subroutine landstop

!     ================
!     SUBROUTINE TANDS
!     ================

      subroutine tands
      use landmod
!
      parameter(zsnowmax=1.)
      parameter(ztop=0.1)
!
      real zsnowz(NHOR)       ! new snow depth
      real zhfls(NHOR)        ! heatflux from soil
      real zhfla(NHOR)        ! heatflux from atmosphere
      real zhflm(NHOR)        ! heatflux used by snowmelt
      real zdsnowz(NHOR)      ! snow depth tendency
      real zcap(NHOR,NLSOIL)  ! heat capacity  of soil layers
      real zdiff(NHOR,NLSOIL) ! thermal conductivity of soil layers
      real zsoilz(NHOR,NLSOIL)! soil layer thicknesses
      real zcap1(NHOR)        ! heat capacity (upper soil layer)
      real zdiff1(NHOR)       ! thermal conductivity (upper soil layer)
      real zsoilz1(NHOR)      ! layer thickness (upper soil layer)
      real zctop(NHOR)        ! heat capacity (top soil/snow layer)
      real zsntop(NHOR)       ! snow depth (top soil/snow layer)
      real zztop(NHOR)        ! depth of the uppermost layer
!
!     debug
!
      integer kmaxl(1),kminl(1)
      real,allocatable :: zfpr1(:)
      real,allocatable :: zfpr2(:)
      real,allocatable :: zfpr3(:)
      real,allocatable :: zfpr4(:)
      real,allocatable :: zfpr5(:)
      real,allocatable :: zfpr6(:)
      real,allocatable :: zfpr7(:)
      real,allocatable :: zfpr8(:)
      real,allocatable :: zfpr9(:)
      real,allocatable :: zfprl(:,:)
!
!     preset
!
      dwater(:)=0.
      dsmelt(:)=0.
      dsndch(:)=dsnowz(:)
      zdsnowz(:)=0.
      zhfls(:)=0.
      zhflm(:)=0.
      zztop(:)=dztop
!
      do jlev=1,NLSOIL
       where(dls(:) > 0.0)
        zcap(:,jlev)=sicecap*dglac(:)+soilcap*(1.-dglac(:))
        zdiff(:,jlev)=sicediff*dglac(:)+soildiff*(1.-dglac(:))
        zsoilz(:,jlev)=dsoilz(jlev)
       endwhere
      enddo
!
      where(dsnowz(:) == 0.0) dsnowt(:)=tmelt
!
!     copy heatflux
!
      zhfla(:)=dshfl(:)+dlhfl(:)+dflux(:,NLEP)
!
!     precip destribution
!
      do jhor=1,NHOR

!
!     surface water
!

       if(dls(jhor) > 0.0) then
        if(dprs(jhor) > 0.) zdsnowz(jhor)=dprs(jhor)
        if(dsnowz(jhor) > 0.) zdsnowz(jhor)=zdsnowz(jhor)+devap(jhor)
        zsnowz(jhor)=AMAX1(0.,dsnowz(jhor)+zdsnowz(jhor)*deltsec)
        zdsnowz(jhor)=(zsnowz(jhor)-dsnowz(jhor))/deltsec
        dwater(jhor)=devap(jhor)+dprl(jhor)+dprc(jhor)-zdsnowz(jhor)
        dsnowz(jhor)=zsnowz(jhor)
       end if
      enddo
!
!     surface temperature and heat flux into the soil
!
      where(dls(:) > 0.0)
!
!     new properties of the top dztop meters (mixed soil/snow layer):
!     1: snow (m of ztop) 2: heat capacity
!
       zsntop(:)=AMIN1(dsnowz(:)*1000./rhosnow,zztop(:))
       zctop(:)=zcap(:,1)*snowcap*zztop(:)                              &
     &         /(zcap(:,1)*zsntop(:)+snowcap*(zztop(:)-zsntop(:)))
!
!     new properties of the uppermost soil layer (mixed soil/snow)
!     1: snow depth (m) 2: total thickness 3: conductivity
!
       zsnowz(:)=AMIN1(dsnowz(:)*1000./rhosnow,zsnowmax)-zsntop(:)
       zsoilz1(:)=zsoilz(:,1)+zsnowz(:)
       zdiff1(:)=zdiff(:,1)*snowdiff*zsoilz1(:)                         &
     &          /(snowdiff*zsoilz(:,1)+zdiff(:,1)*zsnowz(:))
       zcap1(:)=zcap(:,1)*snowcap*zsoilz1(:)                            &
     &         /(zcap(:,1)*zsnowz(:)+snowcap*zsoilz(:,1))
!
!     new surface temp. implicit w.r.t. conductive heat flux
!
       dts(:)=(zctop(:)*zztop(:)*dtsm(:)/deltsec+zhfla(:)               &
              +2.*zdiff1(:)*dsoilt(:,1)/zsoilz1(:))                     &
             /(zctop(:)*zztop(:)/deltsec+2.*zdiff1(:)/zsoilz1(:)) 
!
!     heat flux into the soil
!
       zhfls(:)=2.*zdiff1(:)/zsoilz1(:)*(dts(:)-dsoilt(:,1))
!
      endwhere
!
!     snow melt
!
      where(dls(:) > 0.0 .and. dsnowz(:) > 0. .and. dts(:) > tmelt )
       dts(:)=tmelt
       zhflm(:)=AMAX1(0.,zhfla(:)                                       &
     &                  -zctop(:)*zztop(:)*(dts(:)-dtsm(:))/deltsec     &
     &                  +2.*zdiff1(:)*(dsoilt(:,1)-dts(:))/zsoilz1(:))  
       zsnowz(:)=AMAX1(0.,dsnowz(:)-zhflm(:)*deltsec/((ALS-ALV)*1000.))
       zdsnowz(:)=(zsnowz(:)-dsnowz(:))/deltsec
       zhflm(:)=-zdsnowz(:)*1000.*(ALS-ALV)
!
!     new snow depth (h2o equiv. and snow melt diagnostic)
!
       dsnowz(:)=zsnowz(:)
       dsmelt(:)=-zdsnowz(:)
!
!     heat flux and water flux into the soil
!
       zhfls(:)=2.*zdiff1(:)/zsoilz1(:)*(dts(:)-dsoilt(:,1))
       dwater(:)=dwater(:)-zdsnowz(:)
!
!     new properties of the top dztop meters (mixed soil/snow layer):
!     1: snow (m of dztop) 2: heat capacity
!
       zsntop(:)=AMIN1(dsnowz(:)*1000./rhosnow,zztop(:))
       zctop(:)=zcap(:,1)*snowcap*zztop(:)                              &
     &         /(zcap(:,1)*zsntop(:)+snowcap*(zztop(:)-zsntop(:)))
!
!     new properties of the uppermost soil layer (mixed soil/snow)
!     1: snow depth (m) 2: total thickness 3: conductivity 4: heat capacity
!
       zsnowz(:)=AMIN1(dsnowz(:)*1000./rhosnow,zsnowmax)-zsntop(:)
       zsoilz1(:)=zsoilz(:,1)+zsnowz(:)
       zdiff1(:)=zdiff(:,1)*snowdiff*zsoilz1(:)                         &
     &          /(snowdiff*zsoilz(:,1)+zdiff(:,1)*zsnowz(:))
       zcap1(:)=zcap(:,1)*snowcap*zsoilz1(:)                            &
     &         /(zcap(:,1)*zsnowz(:)+snowcap*zsoilz(:,1))
!
      endwhere
!     write(nud,*) 'landmod 881 zsnowz = ',sum(zsnowz(:))
!
!     correct surface temp and flux where snow is totally melted
!
      where(dls(:) > 0.0 .and. dsnowz(:) == 0. .and. zhflm(:) > 0.)
!
!     new surface (as above)
!
       dts(:)=(zctop(:)*zztop(:)*dtsm(:)/deltsec+zhfla(:)-zhflm(:)      &
     &        +2.*zdiff1(:)*dsoilt(:,1)/zsoilz1(:))                     &
     &       /(zctop(:)*zztop(:)/deltsec+2.*zdiff1(:)/zsoilz1(:)) 
!
!     heat flux into the soil
!
       zhfls(:)=2.*zdiff1(:)/zsoilz1(:)*(dts(:)-dsoilt(:,1))
!
      endwhere
!
!     set new soil properties for the uppermost blended soil/snow layer
!

      where(dls(:) > 0.0)
       zsoilz(:,1)=zsoilz1(:)
       zcap(:,1)=zcap1(:)
       zdiff(:,1)=zdiff1(:)
      endwhere

!
!     deep soil temperatures
!

      call mktsoil(zhfls,zsoilz,zcap,zdiff)

!
!     snow temperature (at the moment set to ts)
!

      where(dls(:) > 0.0 .and. dsnowz(:) > 0.) dsnowt(:)=dts(:)
!
!     limit snow depth to maximum (if switched on)
!     (heat s conserved by artifical cooling of ts (dsnowt)
!      water is conserved)
!
      if(dsmax > 0.) then
       where(dls(:) > 0.0 .and. dsnowz(:) > dsmax)
        zdsnowz(:)=(dsmax-dsnowz(:))/deltsec
!
!     new snow depth (h2o equiv. and snow melt diagnostic)
!
        dsnowz(:)=dsmax
!
!     water flux into the soil
!
        dwater(:)=dwater(:)-zdsnowz(:)
!
!     cool ts (dsnowt)
!
        dts(:)=dts(:)+zdsnowz(:)*1000.*(ALS-ALV)*deltsec/zctop(:)/zztop(:)
!
!       diagnose the lost snow as snow melt
!
        dsmelt(:)=dsmelt(:)-zdsnowz(:)
       endwhere
      endif
!
!     diagnostic: snow depth change
!
      where(dls(:) > 0.) dsndch(:)=(dsnowz(:)-dsndch(:))!/deltsec
!
!     dbug output
!

      if(nprint == 1) then
       allocate(zfpr1(NLON*NLAT))
       allocate(zfpr2(NLON*NLAT))
       allocate(zfpr3(NLON*NLAT))
       allocate(zfpr4(NLON*NLAT))
       allocate(zfprl(NLON*NLAT,NLSOIL))
       call mpgagp(zfpr1,dprs,1)
       call mpgagp(zfpr2,dts,1)
       call mpgagp(zfpr3,dsnowz,1)
       call mpgagp(zfpr4,dls,1)
       call mpgagp(zfprl,dsoilt,NLSOIL)
       if(mypid == NROOT) then
        write(nud,*)'Land Surface Global Diagnostic:'
        zzmax=MAXVAL(zfpr2(:),MASK=(zfpr4(:) > 0.5))
        kmaxl=MAXLOC(zfpr2(:),MASK=(zfpr4(:) > 0.5))
        zzmin=MINVAL(zfpr2(:),MASK=(zfpr4(:) > 0.5))
        kminl=MINLOC(zfpr2(:),MASK=(zfpr4(:) > 0.5))
        write(nud,*)'MAX TS = ',zzmax,' NHOR= ',kmaxl(1)
        write(nud,*)'MIN TS = ',zzmin,' NHOR= ',kminl(1)
        zzmax=MAXVAL(zfpr3(:),MASK=(zfpr4(:) > 0.5))
        kmaxl=MAXLOC(zfpr3(:),MASK=(zfpr4(:) > 0.5))
        zzmin=MINVAL(zfpr3(:),MASK=(zfpr4(:) > 0.5))
        kminl=MINLOC(zfpr3(:),MASK=(zfpr4(:) > 0.5))
        write(nud,*)'MAX ZSNOW = ',zzmax,' NHOR= ',kmaxl(1)
        write(nud,*)'MIN ZSNOW = ',zzmin,' NHOR= ',kminl(1)
        do jlev=1,NLSOIL
         zzmax=MAXVAL(zfprl(:,jlev),MASK=(zfpr4(:) > 0.5))
         kmaxl=MAXLOC(zfprl(:,jlev),MASK=(zfpr4(:) > 0.5))
         zzmin=MINVAL(zfprl(:,jlev),MASK=(zfpr4(:) > 0.5))
         kminl=MINLOC(zfprl(:,jlev),MASK=(zfpr4(:) > 0.5))
         write(nud,*)'MAX TSOIL L= ',jlev,' = ',zzmax,' NHOR= ',kmaxl(1)
         write(nud,*)'MIN TSOIL L= ',jlev,' = ',zzmin,' NHOR= ',kminl(1)
        enddo
        zzmax=MAXVAL(zfpr1(:),MASK=(zfpr4(:) > 0.5))
        kmaxl=MAXLOC(zfpr1(:),MASK=(zfpr4(:) > 0.5))
        write(nud,*)'MAX PRS = ',zzmax,' NHOR= ',kmaxl(1)
       endif
       deallocate(zfpr1)
       deallocate(zfpr2)
       deallocate(zfpr3)
       deallocate(zfpr4)
       deallocate(zfprl)
      endif
      if(nprint == 2) then
       allocate(zfpr1(NLON*NLAT))
       allocate(zfpr2(NLON*NLAT))
       allocate(zfpr3(NLON*NLAT))
       allocate(zfpr4(NLON*NLAT))
       allocate(zfpr5(NLON*NLAT))
       allocate(zfpr6(NLON*NLAT))
       allocate(zfpr7(NLON*NLAT))
       allocate(zfpr8(NLON*NLAT))
       allocate(zfpr9(NLON*NLAT))
       allocate(zfprl(NLON*NLAT,NLSOIL))
       call mpgagp(zfpr1,dprs,1)
       call mpgagp(zfpr2,dts,1)
       call mpgagp(zfpr3,dsnowz,1)
       call mpgagp(zfpr4,zhfla,1)
       call mpgagp(zfpr5,zhflm,1)
       call mpgagp(zfpr6,zhfls,1)
       call mpgagp(zfpr8,dwater,1)
       call mpgagp(zfpr9,dsmelt,1)
       call mpgagp(zfprl,dsoilt,NLSOIL)
       if(mypid == NROOT) then
        write(nud,*)'Land Surface Local Diagnostic:'
        do jlev=1,NLSOIL
         write(nud,*)'TSOIL L= ',jlev,' = ',zfprl(NPRHOR,jlev)
        enddo
        write(nud,*)'TS             = ',zfpr2(NPRHOR)
        write(nud,*)'ZSNOW          = ',zfpr3(NPRHOR)
        write(nud,*)'HFLA           = ',zfpr4(NPRHOR)
        write(nud,*)'HFLX MELTING   = ',zfpr5(NPRHOR)
        write(nud,*)'HFLX TO SOIL   = ',zfpr6(NPRHOR)
        write(nud,*)'PRS (mm/d)     = ',zfpr1(NPRHOR)*1000.*deltsec*mtspd
        write(nud,*)'DWATER (mm/d)  = ',zfpr8(NPRHOR)*1000.*deltsec*mtspd
        write(nud,*)'SNOWMELT (mm/d)= ',zfpr9(NPRHOR)*1000.*deltsec*mtspd
       endif
       deallocate(zfpr1)
       deallocate(zfpr2)
       deallocate(zfpr3)
       deallocate(zfpr4)
       deallocate(zfpr5)
       deallocate(zfpr6)
       deallocate(zfpr7)
       deallocate(zfpr8)
       deallocate(zfpr9)
       deallocate(zfprl)
      endif
!
!     entropy diagnostics
!
      if(nentropy > 0) then
       where(dls(:) > 0.)
        dentropy(:,18)=zctop(:)*zztop(:)*(dts(:)-dtsm(:))/(dtsm(:)*deltsec)
       endwhere
      endif
!
      return
      end subroutine tands

!     ==================
!     SUBROUTINE MKTSOIL
!     ==================

      subroutine mktsoil(pftop,psoilz,pcap,pdiff)
      use landmod
!
!     input
!
      real pftop(NHOR)         ! surface heat flux
      real psoilz(NHOR,NLSOIL) ! layer thickness
      real pcap(NHOR,NLSOIL)   ! heat capacities
      real pdiff(NHOR,NLSOIL)  ! heat conductivities
!
!     local
!
      real ztn(NHOR,NLSOIL)    ! new temperatures
      real zebs(NHOR,NLSOIL)
      real zcap(NHOR,NLSOIL)
      real zdiff(NHOR,NLSOIL-1)
      real ztold(NHOR,NLSOIL)
!
      if(nentropy > 0) ztold(:,:)=dsoilt(:,:)
!
!     implicit scheme for soiltemp
!
!     a) deep layer elimination (zero flux at bottom)
!

      jlev=NLSOIL
      jlem=NLSOIL-1
      where(dls(:) > 0.0)
       zdiff(:,jlem)=2.*pdiff(:,jlev)*pdiff(:,jlem)                     &
     &              /(pdiff(:,jlev)*psoilz(:,jlem)                      &
     &               +pdiff(:,jlem)*psoilz(:,jlev))
       zcap(:,jlev)=pcap(:,jlev)*psoilz(:,jlev)/deltsec
       zebs(:,jlev)=1./(zcap(:,jlev)+zdiff(:,jlem))
       ztn(:,jlev)=zcap(:,jlev)*dsoilt(:,jlev)*zebs(:,jlev)
      endwhere

!
!     middle layer elimination
!

      do jlev=NLSOIL-1,2,-1
       jlep=jlev+1
       jlem=jlev-1
       where(dls(:) > 0.0)
        zdiff(:,jlem)=2.*pdiff(:,jlev)*pdiff(:,jlem)                    &
     &                /(pdiff(:,jlev)*psoilz(:,jlem)                    &
     &                 +pdiff(:,jlem)*psoilz(:,jlev))
        zcap(:,jlev)=pcap(:,jlev)*psoilz(:,jlev)/deltsec
        zebs(:,jlev)=1./(zcap(:,jlev)+zdiff(:,jlem)                     &
     &                  +zdiff(:,jlev)*(1.-zdiff(:,jlev)*zebs(:,jlep)))
        ztn(:,jlev)=zebs(:,jlev)*(zcap(:,jlev)*dsoilt(:,jlev)           &
     &                           +zdiff(:,jlev)*ztn(:,jlep))
       endwhere
      enddo

!
!     top layer elimination
!

       jlev=1
       jlep=2
       where(dls(:) > 0.0)
        zcap(:,jlev)=pcap(:,jlev)*psoilz(:,jlev)/deltsec
        zebs(:,jlev)=1./(zcap(:,jlev)+zdiff(:,jlev)                     &
     &                               *(1.-zdiff(:,jlev)*zebs(:,jlep)))
        ztn(:,jlev)=zebs(:,jlev)*(zcap(:,jlev)*dsoilt(:,jlev)           &
     &                           +zdiff(:,jlev)*ztn(:,jlep)+pftop(:))
       endwhere

!
!     back-substitution
!

      do jlev=2,NLSOIL
       jlem=jlev-1
       where(dls(:) > 0.0)
        ztn(:,jlev)=ztn(:,jlev)+zebs(:,jlev)*zdiff(:,jlem)*ztn(:,jlem)
       endwhere
      enddo

!
!     new temperatures
!

      do jlev=1,NLSOIL
       where(dls(:) > 0.0)
        dsoilt(:,jlev)=ztn(:,jlev)
       endwhere

!
!      do not allow warming of permanet glaciers
!

       where(dls(:) > 0.0 .and. dglac(:) > 0.5)
        dsoilt(:,jlev)=AMIN1(dsoilt(:,jlev),tmelt)
       endwhere
      enddo
!
!     entropy diagnostics
!
      if(nentropy > 0) then
       dentropy(:,19)=0.
       do jlev=1,NLSOIL
        where(dls(:) > 0.)
         dentropy(:,19)=dentropy(:,19)                                  &
     &                 +zcap(:,jlev)*(dsoilt(:,jlev)-ztold(:,jlev))     &
     &                 /dsoilt(:,jlev)
        endwhere
       enddo
      endif
!
      return
      end subroutine mktsoil

!     ================
!     SUBROUTINE WANDR
!     ================

      subroutine wandr
      use landmod
!
!     calculate soil water and runoff
!
      drunoff(:)=0.
      where (dls(:) > 0.0)
       dwatc(:)=dwatc(:)+deltsec*dwater(:)
       drunoff(:)=AMAX1(0.,dwatc(:)-dwmax(:))/deltsec
       dwatc(:)=AMAX1(AMIN1(dwmax(:),dwatc(:)),0.)
      end where
!
      return
      end subroutine wandr

!     ==================
!     SUBROUTINE SOILINI
!     ==================

      subroutine soilini
      use landmod
!
      if(nrestart == 0) then
!
!     initialize snow and temperatures
!
       where(dls(:) > 0.0)
        dsnowz(:)=0.
        dsnowt(:)=tmelt
        dts(:)=dtclsoil(:)
       endwhere
       do jlev=1,NLSOIL
        where(dls(:) > 0.0) dsoilt(:,jlev)=dts(:)
       enddo
!
!     initialize soil water and runoff
!
       if (nwatcini > 0) then
        where(dls(:) > 0.0)
         dwatc(:)=dwatcini
        endwhere
       else
        where(dls(:) > 0.0)
         dwatc(:)=drhsland*drhsfull*dwmax(:)
        endwhere
       endif
       drunoff(:)=0.
!
      endif
!
      return
      end subroutine soilini

!     ===================
!     SUBROUTINE SOILSTEP
!     ===================

      subroutine soilstep
      use landmod
!
!     calculate snow and temperatures
!
      if(nlandt==1) then
       call tands
      else
       where(dls(:) > 0.0)
        dwater(:)=devap(:)+dprl(:)+dprc(:)
        dts(:)=dtclim(:)
       endwhere
      end if
!
!     calculate soil water and runoff
!
      if(nlandw==1) then
       call wandr
      else
       drunoff(:)=0.
       where(dls(:) > 0.)
        drunoff(:)=AMAX1(0.,dwater(:))
        dwatc(:)=dwclim(:)
       endwhere
      endif
!
      return
      end subroutine soilstep

!     ==================
!     SUBROUTINE ROFFINI
!     ==================

      subroutine roffini
      use landmod
      parameter(zcvel=4.2)
      parameter(zcexp=0.18)
!
      real zuroff(NLON,NLAT)
      real zvroff(NLON,NLAT)
      real zoro(NLON,NLAT)
      real zoron(NLON,NLAT)
      real zlsm(NLON,NLAT)
      real zsi(NLON,NLAT)
      real zsir(NLON,NLPP)

!
      ilat = NLAT ! using ilat suppresses compiler warnings for T1
!
      if(NRESTART==0) then

       do jlat=1,NLPP
        do jlon=1,NLON
         jhor=(jlat-1)*NLON+jlon
         darea(jhor)=gwd(jlat)
         zsir(jlon,jlat)=sid(jlat)
        enddo
       enddo

       call mpgagp(zsi,zsir,1)
       call mpgagp(zoro,doro,1)
       call mpgagp(zlsm,dls,1)

       if(mypid==NROOT) then

        zoro(:,:)=MAX(zoro(:,:),0.)
        where(zlsm(:,:) < 1.) zoro(:,:)=zlsm(:,:)-1.

!
!     iterate to remove local minima
!

 1000   continue
        jconv=0
        do jlat=2,ilat-1
         do jlon=2,NLON-1
          if(zlsm(jlon,jlat) > 0.                                       &
     &      .and. zoro(jlon,jlat) <= zoro(jlon,jlat-1)                  &
     &      .and. zoro(jlon,jlat) <= zoro(jlon,jlat+1)                  &
     &      .and. zoro(jlon,jlat) <= zoro(jlon+1,jlat)                  &
     &      .and. zoro(jlon,jlat) <= zoro(jlon-1,jlat)) then
           zoron(jlon,jlat)=1.+MIN(zoro(jlon+1,jlat),zoro(jlon-1,jlat)  &
     &                            ,zoro(jlon,jlat+1),zoro(jlon,jlat-1))
           jconv=jconv+1
          else
           zoron(jlon,jlat)=zoro(jlon,jlat)
          endif
         enddo
         if(zlsm(1,jlat) > 0.                                           &
     &     .and. zoro(1,jlat) <= zoro(1,jlat-1)                         &
     &     .and. zoro(1,jlat) <= zoro(1,jlat+1)                         &
     &     .and. zoro(1,jlat) <= zoro(2,jlat)                           &
     &     .and. zoro(1,jlat) <= zoro(NLON,jlat)) then
          zoron(1,jlat)=1.+MIN(zoro(2,jlat),zoro(NLON,jlat)             &
     &                        ,zoro(1,jlat+1),zoro(1,jlat-1))
          jconv=jconv+1
         else
          zoron(1,jlat)=zoro(1,jlat)
         endif
         if(zlsm(NLON,jlat) > 0.                                        &
     &     .and. zoro(NLON,jlat) <= zoro(NLON,jlat-1)                   &
     &     .and. zoro(NLON,jlat) <= zoro(NLON,jlat+1)                   &
     &     .and. zoro(NLON,jlat) <= zoro(1,jlat)                        &
     &     .and. zoro(NLON,jlat) <= zoro(NLON-1,jlat)) then
          zoron(NLON,jlat)=1.+MIN(zoro(1,jlat),zoro(NLON-1,jlat)        &
     &                           ,zoro(NLON,jlat+1),zoro(NLON,jlat-1))
          jconv=jconv+1
         else
          zoron(NLON,jlat)=zoro(NLON,jlat)
         endif
        enddo
        do jlon=2,NLON-1
         if(zlsm(jlon,1) > 0.                                           &
     &    .and. zoro(jlon,1) <= zoro(jlon,2)                            &
     &    .and. zoro(jlon,1) <= zoro(jlon+1,1)                          &
     &    .and. zoro(jlon,1) <= zoro(jlon-1,1)) then
          zoron(jlon,1)=1.+MIN(zoro(jlon+1,1),zoro(jlon-1,1)            &
     &                        ,zoro(jlon,2))
          jconv=jconv+1
         else
          zoron(jlon,1)=zoro(jlon,1)
         endif
         if(zlsm(jlon,NLAT) > 0.                                        &
     &    .and. zoro(jlon,NLAT) <= zoro(jlon,NLAT-1)                    &
     &    .and. zoro(jlon,NLAT) <= zoro(jlon+1,NLAT)                    &
     &    .and. zoro(jlon,NLAT) <= zoro(jlon-1,NLAT)) then
          zoron(jlon,NLAT)=1.+MIN(zoro(jlon+1,NLAT),zoro(jlon-1,NLAT)   &
     &                           ,zoro(jlon,NLAT-1))
          jconv=jconv+1
         else
          zoron(jlon,NLAT)=zoro(jlon,NLAT)
         endif
        enddo
        if(zlsm(1,1) > 0.                                               &
     &    .and. zoro(1,1) <= zoro(1,2)                                  &
     &    .and. zoro(1,1) <= zoro(2,1)                                  &
     &    .and. zoro(1,1) <= zoro(NLON,1)) then
         zoron(1,1)=1.+MIN(zoro(2,1),zoro(NLON,1)                       &
     &                    ,zoro(1,2))
         jconv=jconv+1
        else
         zoron(1,1)=zoro(1,1)
        endif
        if(zlsm(NLON,NLAT) > 0.                                         &
     &    .and. zoro(NLON,NLAT) <= zoro(NLON,NLAT-1)                    &
     &    .and. zoro(NLON,NLAT) <= zoro(1,NLAT)                         &
     &    .and. zoro(NLON,NLAT) <= zoro(NLON-1,NLAT)) then
         zoron(NLON,NLAT)=1.+MIN(zoro(1,NLAT),zoro(NLON-1,NLAT)         &
     &                          ,zoro(NLON,NLAT-1))
         jconv=jconv+1
        else
         zoron(NLON,NLAT)=zoro(NLON,NLAT)
        endif
        if(zlsm(NLON,1) > 0.                                            &
     &    .and. zoro(NLON,1) <= zoro(NLON,2)                            &
     &    .and. zoro(NLON,1) <= zoro(1,1)                               &
     &    .and. zoro(NLON,1) <= zoro(NLON-1,1)) then
         zoron(NLON,1)=1.+MIN(zoro(1,1),zoro(NLON-1,1)                  &
     &                       ,zoro(NLON,2))
         jconv=jconv+1
        else
         zoron(NLON,1)=zoro(NLON,1)
        endif
        if(zlsm(1,NLAT) > 0.                                            &
     &    .and. zoro(1,NLAT) <= zoro(1,NLAT-1)                          &
     &    .and. zoro(1,NLAT) <= zoro(2,NLAT)                            &
     &    .and. zoro(1,NLAT) <= zoro(NLON,NLAT)) then
         zoron(1,NLAT)=1.+MIN(zoro(NLON,NLAT),zoro(2,NLAT)              &
     &                       ,zoro(1,NLAT-1))
         jconv=jconv+1
        else
         zoron(1,NLAT)=zoro(1,NLAT)
        endif
        zoro(:,:)=zoron(:,:)
        if(jconv > 0 ) goto 1000

        do jlat=1,NLAT
         do jlon=1,NLON-1
          zdx=TWOPI*cos(ASIN(zsi(jlon,jlat)))*plarad/real(NLON)
          zdh=(zoro(jlon+1,jlat)-zoro(jlon,jlat))/zdx
          zfac=1.
          if(zdh > 0.) zfac=-1.
          zuroff(jlon,jlat)=zfac*zcvel/zdx*ABS(zdh)**zcexp
         enddo
         zdx=TWOPI*cos(ASIN(zsi(NLON,jlat)))*plarad/real(NLON)
         zdh=(zoro(1,jlat)-zoro(NLON,jlat))/zdx
         zfac=1.
         if(zdh > 0.) zfac=-1.
         zuroff(NLON,jlat)=zfac*zcvel/zdx*ABS(zdh)**zcexp
        enddo

        do jlat=1,NLAT-1
         do jlon=1,NLON
          zdy=(ASIN(zsi(jlon,jlat+1))-ASIN(zsi(jlon,jlat)))*plarad
          zdh=(zoro(jlon,jlat+1)-zoro(jlon,jlat))/zdy
          zfac=1.
          if(zdh < 0.) zfac=-1.
          zvroff(jlon,jlat)=zfac*zcvel/zdy*ABS(zdh)**zcexp
         enddo
        enddo
        zvroff(1:NLON,NLAT)=0.
!!FL        
!        open(62,file='roffuv.srv',form='unformatted')
!        write(62) 131,0,0,0,NLON,NLAT,0,0
!        write(62) zuroff
!        write(62) 132,0,0,0,NLON,NLAT,0,0
!        write(62) zvroff
!        write(62) 129,0,0,0,NLON,NLAT,0,0
!        write(62) zoro
!        write(62) 172,0,0,0,NLON,NLAT,0,0
!        write(62) zlsm
!        close(62)
!!
       endif

       call mpscgp(zuroff,duroff,1)
       call mpscgp(zvroff,dvroff,1)

       driver(:)=0.
       drunoff(:)=0.

      endif

      return
      end subroutine roffini

!     ===================
!     SUBROUTINE ROFFSTEP
!     ===================

      subroutine roffstep
      use landmod

!
!     make advection
!

      call mkradv(driver,drunoff,duroff,dvroff,darea,dls)

!
!     new water
!

      where(dls(:) > 0.) driver(:)=driver(:)+drunoff(:)*deltsec

      return
      end subroutine roffstep

!     =================
!     SUBROUTINE MKRADV
!     =================

      subroutine mkradv(zriver,zroff,zuroff,zvroff,zarea,zls)
      use landmod

      real zarea(NLON,NLPP)
      real zuroff(NLON,NLPP)
      real zvroff(NLON,NLPP)
      real zriver(NLON,NLPP)
      real zrivin(NLON,NLPP)
      real zroff(NLON,NLPP)
      real zls(NLON,NLPP)

      real zrip(NLON,NLPP)
      real zrop(NLON,NLPP)

      real zrigl(NLON,NLAT),zriglp(NLON,NLAT)
      real zrogl(NLON,NLAT),zroglp(NLON,NLAT)
!

      zrivin(:,:)=zriver(:,:)*zarea(:,:)*zls(:,:)

!
!     make advection in zonal direction (no partitioning in mpp)
!

      do jlon=1,NLON-1
       where(zuroff(jlon,:) > 0.)
        zroff(jlon,:)=zroff(jlon,:)                                     &
     &               -zuroff(jlon,:)*zrivin(jlon,:)/zarea(jlon,:)
        zroff(jlon+1,:)=zroff(jlon+1,:)                                 &
     &                 +zuroff(jlon,:)*zrivin(jlon,:)/zarea(jlon+1,:)
       elsewhere
        zroff(jlon,:)=zroff(jlon,:)                                     &
     &               -zuroff(jlon,:)*zrivin(jlon+1,:)/zarea(jlon,:)
        zroff(jlon+1,:)=zroff(jlon+1,:)                                 &
     &                 +zuroff(jlon,:)*zrivin(jlon+1,:)/zarea(jlon+1,:)
       endwhere
      enddo
      where(zuroff(NLON,:) > 0.)
       zroff(NLON,:)=zroff(NLON,:)                                      &
     &              -zuroff(NLON,:)*zrivin(NLON,:)/zarea(NLON,:)
       zroff(1,:)=zroff(1,:)                                            &
     &           +zuroff(NLON,:)*zrivin(NLON,:)/zarea(1,:)
      elsewhere
       zroff(NLON,:)=zroff(NLON,:)                                      &
     &              -zuroff(NLON,:)*zrivin(1,:)/zarea(NLON,:)
       zroff(1,:)=zroff(1,:)                                            &
     &           +zuroff(NLON,:)*zrivin(1,:)/zarea(1,:)
      endwhere

!
!     make advection in meridional direction (partitioning in mpp)
!     simple schema (exchange all gps)
!

      call mpgagp(zrigl,zrivin,1)
      if(mypid==NROOT) then
       zriglp(:,1:NLAT-1)=zrigl(:,2:NLAT)
       zriglp(:,NLAT)=0.
      endif
      call mpscgp(zriglp,zrip,1)
      where(zvroff(:,:) < 0.)
       zroff(:,:)=zroff(:,:)+zvroff(:,:)*zrivin(:,:)/zarea(:,:)
       zrop(:,:)=-zvroff(:,:)*zrivin(:,:)
      elsewhere
       zroff(:,:)=zroff(:,:)+zvroff(:,:)*zrip(:,:)/zarea(:,:)
       zrop(:,:)=-zvroff(:,:)*zrip(:,:)
      endwhere
      call mpgagp(zroglp,zrop,1)
      if(mypid==NROOT) then
       zrogl(:,2:NLAT)=zroglp(:,1:NLAT-1)
       zrogl(:,1)=0.
      endif
      call mpscgp(zrogl,zrop,1)

      zroff(:,:)=zroff(:,:)+zrop(:,:)/zarea(:,:)

      return
      end subroutine mkradv

!     =================
!     SUBROUTINE GETTCL
!     =================

      subroutine gettcll
      use landmod
!
!     get surface temperature annual cycle
!
      call momint(nperpetual,nstep+1,jm1,jm2,zgw2)
      zgw1 = 1.0 - zgw2
      dtclim(:)=zgw1*dtcl(:,jm1)+zgw2*dtcl(:,jm2)
      dwclim(:)=zgw1*dwcl(:,jm1)+zgw2*dwcl(:,jm2)
      return
      end subroutine gettcll

!     =================
!     SUBROUTINE GETALB
!     =================

      subroutine getalb
      use landmod
      
      real :: aa = 5.2
      real :: yy = 4.0
      real :: bf = 0.0
      real :: al = 0.0
!
!     get surface background albedo from  annual cycle
!
      call momint(nperpetual,nstep+1,jm1,jm2,zgw2)
      zgw1 = 1.0 - zgw2
      dalbclim(:)=zgw1*dalbcl(:,jm1)+zgw2*dalbcl(:,jm2)
      dalbclim1(:)=zgw1*dalbcl1(:,jm1)+zgw2*dalbcl1(:,jm2)
      dalbclim2(:)=zgw1*dalbcl2(:,jm1)+zgw2*dalbcl2(:,jm2)
      
!      
!     Modify the surface background albedo according to soil water capacity
!
      if (nwetsoil > 0.5) then
        do jhor = 1,NHOR
          if (dls(jhor)>0.5) then
            if (dts(jhor)>273.15) then !We have wet soil. The wetter, the darker.
              bf = dwater(jhor)/wsmax
              al = 1.0/aa*(max(bf,1.0e-3)**((1.0-yy)/yy) - 1)**(1.0/yy)
              dalbclim(jhor) = alblandmax * exp(-(bf*25)**6)+ &
&                              al * (1 - exp(-(bf*25)**6) - exp(-((1-bf)*30)**9))+ &
&                              doceanalb(1) * exp(-((1-bf)*30)**9)
              dalbclim1(jhor) = dalbclim(jhor)
              dalbclim2(jhor) = alblandmax * exp(-(bf*25)**6)+ &
&                              al * (1 - exp(-(bf*25)**6) - exp(-((1-bf)*30)**9))+ &
&                              doceanalb(2) * exp(-((1-bf)*30)**9)
            else !Frozen soil.
              dalbclim(jhor) = albland
              dalbclim1(jhor) = dgroundalb(1)
              dalbclim2(jhor) = dgroundalb(2)
            endif
          endif
        enddo
      endif
      
      return
      end subroutine getalb

