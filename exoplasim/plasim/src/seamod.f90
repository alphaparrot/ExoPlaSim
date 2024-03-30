      module seamod
      use pumamod
!
!     version identifier (date)
!
      character(len=80) :: version = '29.07.2004 by Larry'

!
!     namelist parameters
!
      integer :: ncpl_atmos_ice = 1     ! ice-atmosphere coupling timesteps
!
      real    :: albsea   = 0.069 ! albedo for free ocean
      real    :: albice   = 0.7   ! max. albedo for sea ice
      real    :: dz0sea   = 1.5E-5! roughness length sea
      real    :: dz0ice   = 0.001 !  "          "    ice
      real    :: drhssea  = 1.    ! wetness factor sea
      real    :: drhsice  = 1.    !  "         "   ice
      real    :: charnock = 0.018 ! albedo for free ocean
      
!
!     global arrays
!
!     surface definitions
!
      real :: dts(NHOR)            ! surface temperature (K)
      real :: dqs(NHOR)            ! surface humidity    (kg/kg)
      real :: dsst(NHOR)           ! sea surface temperature (K)
!
!     variables and arrays for ocean coupling
!
      integer :: naccua      = 0   ! counter for accumulation puma
!
      real :: csst(NHOR)     = 0.  ! sst (K)  (ocean -> puma)
      real :: cmld(NHOR)     = 0.  ! mixed-layer depth (m) (ocean -> puma)
      real :: cts(NHOR)      = 0.  ! surface temp (K)  (ice -> puma)
      real :: cicec(NHOR)    = 0.  ! sea ice compactness (ice ->puma)
      real :: ciced(NHOR)    = 0.  ! sea ice thickness  (ice -> puma)
      real :: csnow(NHOR)    = 0.  ! snow depth  (ice -> puma)
      real :: csmelt(NHOR)   = 0.  ! snow melt  (ice -> puma)
      real :: csndch(NHOR)   = 0.  ! snow depth change (ice -> puma)
!
!     accumulated fluxes  (puma -> ice)
!
      real :: cheata(NHOR)   = 0.  ! total heat flux (w/m2)
      real :: cpmea(NHOR)    = 0.  ! P-E (m/s)
      real :: cprsa(NHOR)    = 0.  ! snow precipitation (m/s)
      real :: croffa(NHOR)   = 0.  ! Runoff (m/s) (puma -> ocean)
      real :: ctauxa(NHOR)   = 0.  ! u-stress (pa)
      real :: ctauya(NHOR)   = 0.  ! v-stress (pa)
      real :: cust3a(NHOR)   = 0.  ! ustar**3 (m3/s3)
      real :: cshfla(NHOR)   = 0.  ! surface sens. heat flx (w/m2) (puma->ice)
      real :: cshdta(NHOR)   = 0.  ! deriv. of dshfl (w/(m2 s))(puma->ice)
      real :: clhfla(NHOR)   = 0.  ! surface latent heat flx (w/m2) (puma->ice)
      real :: clhdta(NHOR)   = 0.  ! deriv. of dlhfl (w/(m2 s))(puma->ice)
      real :: cswfla(NHOR)   = 0.  ! net solar radiation (w/m2) (puma -> ice)
      real :: clwfla(NHOR)   = 0.  ! surface thermal radiation (w/m2) (puma -> ice)
!
      end module seamod

!     =================
!     SUBROUTINE SEAINI
!     =================

      subroutine seaini
      use seamod
!
      namelist/seamod_nl/albsea,albice,dz0sea,dz0ice,drhssea,drhsice       &
     &               ,ncpl_atmos_ice,charnock,doceanalb,dicealbmn,dicealbmx
!
!     read namelist
!
      if(mypid == NROOT) then
         open(12,file=seamod_namelist)
         read(12,seamod_nl)
         write(nud,'(/,"*********************************************")')
         write(nud,'("* SEAMOD ",a34," *")') trim(version)
         write(nud,'("*********************************************")')
         write(nud,'("* Namelist SEAMOD_NL from <",a15,"> *")') &
               trim(seamod_namelist)
         write(nud,'("*********************************************")')
         write(nud,seamod_nl)
         close(12)
      end if
!
      call mpbcr(albsea)
      call mpbcr(albice)
      call mpbcr(dz0sea)
      call mpbcr(dz0ice)
      call mpbcr(drhssea)
      call mpbcr(drhsice)
      call mpbci(ncpl_atmos_ice)
      call mpbcr(charnock)
      call mpbcrn(doceanalb,2*NLIGHTS)
      call mpbcrn(dicealbmn,2*NLIGHTS)
      call mpbcrn(dicealbmx,2*NLIGHTS)
!
!     initialize ice (and ocean)
!
      call iceini(n_start_step,nrestart,noutput,m_days_per_year         &
     &     ,ngui,cts,csst,cmld,cicec,ciced,csnow,mtspd,day_24hr,deglat &
     &     ,icemod_namelist,oceanmod_namelist,ice_output,ocean_output)
!
!     set puma surface variables
!
      where(dls(:) < 0.5)
       dts(:)=cts(:)
       dicec(:)=cicec(:)
       diced(:)=ciced(:)
       dsnow(:)=csnow(:)
       dsst(:)=csst(:)
       dmld(:)=cmld(:)
       dt(:,NLEP)=dts(:)
       dqs(:)=dq(:,NLEP)
      endwhere
!
!     read coupling parameters if restart
!
      if (nrestart == 1) then
         if (mypid == NROOT) then
            call get_restart_integer('naccua',naccua)
         endif
         call mpbci(naccua)

         call mpgetgp('dts'   ,dts   ,NHOR,1)
         call mpgetgp('cheata',cheata,NHOR,1)
         call mpgetgp('cpmea' ,cpmea ,NHOR,1)
         call mpgetgp('cprsa' ,cprsa ,NHOR,1)
         call mpgetgp('croffa',croffa,NHOR,1)
         call mpgetgp('ctauxa',ctauxa,NHOR,1)
         call mpgetgp('ctauya',ctauya,NHOR,1)
         call mpgetgp('cust3a',cust3a,NHOR,1)
         call mpgetgp('cshfla',cshfla,NHOR,1)
         call mpgetgp('cshdta',cshdta,NHOR,1)
         call mpgetgp('clhfla',clhfla,NHOR,1)
         call mpgetgp('clhdta',clhdta,NHOR,1)
         call mpgetgp('cswfla',cswfla,NHOR,1)
         call mpgetgp('clwfla',clwfla,NHOR,1)
      else
      where(dls(:) < 0.5)
       dqs(:)  = rdbrv*ra1*EXP(ra2*(dt(:,NLEP)-TMELT)   &
     &          /(dt(:,NLEP)-ra4))/psurf
       dqs(:)  = dqs(:)/(1.-(1./rdbrv-1.)*dqs(:))
       dq(:,NLEP) = dqs(:)
       drhs(:)=drhssea*(1.-dicec(:))+drhsice*dicec(:)
       dz0(:)=dz0sea*(1.-dicec(:))+dz0ice*dicec(:)
      endwhere
      do klight=1,NLIGHTS
        where(dls(:) < 0.5)
         dalb(:,klight)=albsea*(1.-dicec(:))   &
     &           +dicec(:)*AMIN1(albice,0.5+0.025*(273.-dts(:)))
         dsalb(:,2*klight-1)=doceanalb(2*klight-1)*(1.-dicec(:))   &
     &           +dicec(:)*AMIN1(dicealbmx(2*klight-1),dicealbmn(2*klight-1)+0.025*(273.-dts(:)))
         dsalb(:,2*klight  )=doceanalb(2*klight)*(1.-dicec(:))   &
     &           +dicec(:)*AMIN1(dicealbmx(2*klight),dicealbmn(2*klight)+0.025*(273.-dts(:)))
        endwhere
      enddo
      endif
!
      return
      end subroutine seaini

!     ===================
!     SUBROUTINE SEASTEP
!     ===================

      subroutine seastep
      use seamod
!
      real :: zz0(NHOR) = 0.
!
!     coupling to sea ice
!
!     a) accumulate fluxes
!
      where(dls(:) < 0.5)
       cheata(:)=dshfl(:)+dswfl(:,NLEP)+dlwfl(:,NLEP)+dlhfl(:)+cheata(:)
       cpmea(:)=dprl(:)+dprc(:)+devap(:)+cpmea(:)
       cprsa(:)=dprs(:)+cprsa(:)
       croffa(:)=drunoff(:)+croffa(:)
       ctauxa(:)=dtaux(:)+ctauxa(:)
       ctauya(:)=dtauy(:)+ctauya(:)
       cust3a(:)=dust3(:)+cust3a(:)
       cshfla(:)=dshfl(:)+cshfla(:)
       cshdta(:)=dshdt(:)+cshdta(:)
       clhfla(:)=dlhfl(:)+clhfla(:)
       clhdta(:)=dlhdt(:)+clhdta(:)
       cswfla(:)=dswfl(:,NLEP)+cswfla(:)
       clwfla(:)=dlwfl(:,NLEP)+clwfla(:)
      end where
      naccua=naccua+1
!
!     call for ice (and ocean)
!
      if (mod(nstep,ncpl_atmos_ice) == 0 .and. nkits==0) then
       where(dls(:) < 0.5)
        cheata(:)=cheata(:)/real(naccua)
        cpmea(:)=cpmea(:)/real(naccua)
        cprsa(:)=cprsa(:)/real(naccua)
        croffa(:)=croffa(:)/real(naccua)
        ctauxa(:)=ctauxa(:)/real(naccua)
        ctauya(:)=ctauya(:)/real(naccua)
        cust3a(:)=cust3a(:)/real(naccua)
        cshfla(:)=cshfla(:)/real(naccua)
        cshdta(:)=cshdta(:)/real(naccua)
        clhfla(:)=clhfla(:)/real(naccua)
        clhdta(:)=clhdta(:)/real(naccua)
        cswfla(:)=cswfla(:)/real(naccua)
        clwfla(:)=clwfla(:)/real(naccua)
       endwhere
!
       call icestep(cheata,cshfla,cshdta,clhfla,clhdta,clwfla,cswfla    &
     &             ,cpmea,croffa,cprsa,ctauxa,ctauya,cust3a,cts,cicec   &
     &             ,ciced,csnow,csmelt,csndch,csst,cmld)
!
       cheata(:)=0.
       cpmea(:)=0.
       cprsa(:)=0.
       croffa(:)=0.
       ctauxa(:)=0.
       ctauya(:)=0.
       cust3a(:)=0.
       cshfla(:)=0.
       cshdta(:)=0.
       clhfla(:)=0.
       clhdta(:)=0.
       cswfla(:)=0.
       clwfla(:)=0.
       naccua=0
!
!     set puma variable
!
       where (dls(:) < 0.5)
        dsst(:)  = csst(:)
        dts(:)   = cts(:)
        dmld(:)  = cmld(:)
        dicec(:) = cicec(:)
        diced(:) = ciced(:)
        dsnow(:) = csnow(:)
        dsmelt(:) = csmelt(:)
        dsndch(:) = csndch(:)
       endwhere
      endif
!
!     set dependent surface variables
!
      where(dls(:) < 0.5)
       dt(:,NLEP)=dts(:)
       dqs(:)=rdbrv*ra1*EXP(ra2*(dt(:,NLEP)-TMELT)                      &
     &       /(dt(:,NLEP)-ra4))/dp(:)
       dqs(:)=dqs(:)/(1.-(1./rdbrv-1.)*dqs(:))
       dq(:,NLEP)=dqs(:)
       drhs(:)=1.
!
!     z0 after charnock
!
       zz0(:)=charnock*SQRT(dtaux(:)**2+dtauy(:)**2)*gascon*dt(:,NLEP)  &
     &       /(ga*dp(:))
       zz0(:)=AMAX1(zz0(:),dz0sea)
!
       drhs(:)=drhssea*(1.-dicec(:))+drhsice*dicec(:)
       dz0(:)=zz0(:)*(1.-dicec(:))+dz0ice*dicec(:)
      endwhere
      do klight=1,NLIGHTS
       where(dls(:) < 0.5)
          dalb(: ,klight)=albsea*(1.-dicec(:))                                     &
     &              +dicec(:)*AMIN1(albice,0.5+0.025*(273.-dts(:)))
          dsalb(:,2*klight-1)=doceanalb(2*klight-1)*(1.-dicec(:))                                     &
     &              +dicec(:)*AMIN1(dicealbmx(2*klight-1),dicealbmn(2*klight-1)+0.025*(273.-dts(:)))
          dsalb(:,2*klight  )=doceanalb(2*klight  )*(1.-dicec(:))                                     &
     &              +dicec(:)*AMIN1(dicealbmx(2*klight  ),dicealbmn(2*klight  )+0.025*(273.-dts(:)))
       endwhere
      enddo
!
      return
      end subroutine seastep

!     ===================
!     SUBROUTINE SEASTOP
!     ===================

      subroutine seastop
      use seamod
!
!     finalize ice (and ocean)
!
      call icestop
!
!     write restart file
!
      if (mypid == NROOT) then
         call put_restart_integer('naccua',naccua)
      endif

      call mpputgp('dts'   ,dts   ,NHOR,1)
      call mpputgp('cheata',cheata,NHOR,1)
      call mpputgp('cpmea' ,cpmea ,NHOR,1)
      call mpputgp('cprsa' ,cprsa ,NHOR,1)
      call mpputgp('croffa',croffa,NHOR,1)
      call mpputgp('ctauxa',ctauxa,NHOR,1)
      call mpputgp('ctauya',ctauya,NHOR,1)
      call mpputgp('cust3a',cust3a,NHOR,1)
      call mpputgp('cshfla',cshfla,NHOR,1)
      call mpputgp('cshdta',cshdta,NHOR,1)
      call mpputgp('clhfla',clhfla,NHOR,1)
      call mpputgp('clhdta',clhdta,NHOR,1)
      call mpputgp('cswfla',cswfla,NHOR,1)
      call mpputgp('clwfla',clwfla,NHOR,1)

      return
      end subroutine seastop
