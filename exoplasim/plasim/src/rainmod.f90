      module rainmod
      use pumamod
!
!     version identifier (date)
!
      character(len=80) :: version = '11.07.2010 by Larry'
!
!     namelist parameter:
!

      integer :: kbeta    = 1 ! switch for beta in kuo (1/0=yes/no, 2=expl.)
      integer :: nprl     = 1 ! switch for large scale precip (1/0=yes/no)
      integer :: nprc     = 1 ! switch for large convective precip (1/0=yes/no)
      integer :: nclouds  = 1 ! switch for cloud generation (1/0=yes/no)
      integer :: ndca     = 1 ! switch for dry convec. adjustment (1/0=yes/no)
      integer :: ncsurf   = 1 ! switch for conv. starts at surface (1/0=yes/no)
      integer :: nmoment  = 0 ! switch for momentum mixing (1/0=yes/no)
      integer :: nshallow = 1 ! switch for shallow convection (1/0=yes/no)
      integer :: nstorain = 0 ! switch for stochastic rain (1/0=yes/no)
      integer :: nevapprec= 1 ! switch for evaporation of precip.(1/0=yes/no)
      integer :: nbeta    = 3 ! exp. to compute beta
!
      real :: rhbeta   =  0.  ! critical humidity (comp. of beta)
      real :: rbeta    =  0.  ! beta if kbeta=2 (sensitivity tests)
      real :: clwcrit1 = -0.1 ! 1st critical vertical velocity for clouds
      real :: clwcrit2 =  0.0 ! 2nd critical vertical velocity for clouds
                              ! not active if clwcrit1 < clwcrit2
      real :: pdeep   =999999.! pressure to destinguish deep and shallow con. 
      real :: pdeepth = 70000.! pressure threshold for shallow con. 
      real :: rkshallow = 10. ! diffusivity for shallow convection (m*m/s)
      real :: gamma   = 0.01  ! tuning parameter for evaporation of precip.

      real :: rcritmod = 1.0   ! Modifier for cloud critical relative humidity
      real :: rcritslope = 0.0 !By-level modifier for changing cloud height bias
      real :: rcrit(NLEV)    ! critical relative hum. for non conv. clouds

!
!     global scalras
!

      real :: clwfac    = 0. ! smothing for cloud suppression
      real :: time4rain = 0. ! CPU time needed for rainmod
      real :: time4prl  = 0. ! CPU time spent for large scale precip
      real :: time4prc  = 0. ! CPU time spent for convective precip
      real :: time4dca  = 0. ! CPU time spent dry convection
      real :: time4cl   = 0. ! CPU time spent for cloud formation

!
!     global arrays
!

      integer :: icclev(NHOR,NLEV)  ! flag for convective cloud layer
      integer :: icctot(NHOR)       ! total convective cloud layers
!
!     note: the following snow/rain will be only dq/dt*dsigma (units!)
!
      real :: dprcl(NHOR,NLEV)      ! convective rain for each level
      real :: dprll(NHOR,NLEV)      ! large scale rain for each level
      real :: dprscl(NHOR,NLEV)     ! convective snow for each level
      real :: dprsll(NHOR,NLEV)     ! large scale snow for each level
!
      end module rainmod

!     ==================
!     SUBROUTINE RAININI
!     ==================

      subroutine rainini
      use rainmod
!
      namelist/rainmod_nl/kbeta,nprl,nprc,ndca,ncsurf,nmoment,nshallow  &
     &       ,nstorain,rcrit,clwcrit1,clwcrit2,pdeep,rkshallow,gamma    &
     &       ,nclouds,pdeepth,nevapprec,nbeta,rhbeta,rbeta,rcritmod,rcritslope      
!
!     reset defaults (according to general setup... tuning)
!
      if(NTRU==21 .and. NLEV==5) then
       nshallow=0
      endif 
      if(NTRU==42 .and. NLEV==10) then
       gamma=0.007
      endif    
!
      rcrit(:)=MAX(0.85,MAX(sigma(:),1.-sigma(:)))
!
      if(mypid==NROOT) then
         open(11,file=rainmod_namelist)
         read(11,rainmod_nl)
         close(11)
         write(nud,'(/," ***********************************************")')
         write(nud,'(" * RAINMOD ",a35," *")') trim(version)
         write(nud,'(" ***********************************************")')
         write(nud,'(" * Namelist RAINMOD_NL from <rainmod_namelist> *")')
         write(nud,'(" ***********************************************")')
         write(nud,rainmod_nl)
         
         do jlev=1,NLEV
            rcrit(jlev) = rcrit(jlev)*rcritmod*(1-2*rcritslope*(float(NLEV)/2-jlev)/(float(NLEV)))
         enddo
      endif
!
!     broadcast namelist parameter
!
      call mpbci(kbeta)
      call mpbci(nprl)
      call mpbci(nprc)
      call mpbci(nclouds)
      call mpbci(ndca)
      call mpbci(ncsurf)
      call mpbci(nmoment)
      call mpbci(nshallow)
      call mpbci(nstorain)
      call mpbci(nevapprec)
      call mpbci(nbeta)
      call mpbcr(rbeta)
      call mpbcr(rhbeta)
      call mpbcr(clwcrit1)
      call mpbcr(clwcrit2)
      call mpbcr(pdeep)
      call mpbcr(pdeepth)
      call mpbcr(rkshallow)
      call mpbcr(gamma)
      call mpbcrn(rcrit,NLEV)
      call mpbcrn(rcrit,NLEV)
!
!     smoothing factor for cloud suppression
!
      if(clwcrit1 > clwcrit2) then
       clwfac=1./(clwcrit1-clwcrit2)
      elseif(clwcrit1 < clwcrit2) then
       clwfac=-1.
      else
       clwfac=1.
      endif
!
      return
      end subroutine rainini

!     ===================
!     SUBROUTINE RAINSTEP
!     ===================

      subroutine rainstep
      use rainmod
!
!**   performance estimate
!
      if(ntime==1) call mksecond(zsec,0.)
!
!
!**   0) preset precipitation
!

      dprl(:)=0.
      dprc(:)=0.
      dprs(:)=0.

!
!**   1) transfer adiabatic q-tendencies to grid point space
!        (needed for kuo scheme)
!

      if(nprc == 1) call mkdqtgp

!
!**   2) kuo convection
!

      if(ntime ==1) call mksecond(zsec1,0.)

      if(nprc == 1) call kuo
     
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4prc=time4prc+zsec1
      endif

!
!**   3) dry convective adjustment
!

      if(ntime ==1) call mksecond(zsec1,0.)
      if(ndca == 1) call mkdca
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4dca=time4dca+zsec1
      endif

!
!**   4) large scale precipitation
!

      if(ntime ==1) call mksecond(zsec1,0.)
      if(nprl == 1) call mklsp
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4prl=time4prl+zsec1
      endif
  
!
!*    rain and snow
!

      call mkrain

!
!**   5) diagnose cloud cover
!

      if(ntime ==1) call mksecond(zsec1,0.)
      if(nclouds ==1) call mkclouds
      if(ntime ==1) then
       call mksecond(zsec1,zsec1)
       time4cl=time4cl+zsec1
       call mksecond(zsec,zsec)
       time4rain=time4rain+zsec
      endif
!
      return
      end subroutine rainstep

!     ===================
!     SUBROUTINE RAINSTOP
!     ===================

      subroutine rainstop
      use rainmod
!
      if(mypid == NROOT .and. ntime == 1) then
       write(nud,*)'*********************************************'
       write(nud,*)' CPU usage in RAINSTEP (ROOT process only):  '
       write(nud,*)'    All routines   : ',time4rain,' s'
       write(nud,*)'    Large scale p  : ',time4prl,' s'
       write(nud,*)'    Convective p   : ',time4prc,' s'
       write(nud,*)'    Cloud formation: ',time4cl,' s'
       write(nud,*)'    Dry convection : ',time4dca,' s'
       write(nud,*)'*********************************************'
      endif
!
      return
      end subroutine rainstop


!     ==================
!     SUBROUTINE MKDQTGP
!     ==================

      subroutine mkdqtgp
      use rainmod
!
      real zsqt(NESP,NLEV)
!
      if (nqspec == 1) then
         call mpgallsp(zsqt,sqt,NLEV)
         call sp2fl(zsqt,dqt,NLEV)
         call fc2gp(dqt,NLON,NLPP*NLEV)
!
         do jlev=1,NLEV
          dqt(:,jlev)=dqt(:,jlev)*psurf*ww/dp(:)
         enddo
      endif
!
      return
      end subroutine mkdqtgp

!     ==========================
!     SUBROUTINE GET_RANDOM_ROOT
!     ==========================

      subroutine get_random_root(pnoise,n)
      use pumamod
      real :: pnoise(n)
      real :: znoise(n * nproc)

      if (mypid == NROOT) then
         do j = 1 , n * nproc
            call random_number(znoise(j))
         enddo
      endif
      call mpscrn(znoise,n)
      pnoise(:) = znoise(1:n)
      return
      end subroutine get_random_root
         
      
!     ================
!     SUBROUTINE MKLSP
!     ================

      subroutine mklsp
      use rainmod
!
!     local arrays
!

      real zq(NHOR)          ! preliminary q at t+dt
      real zt(NHOR)          ! preliminary t at t+dt
      real zqn(NHOR)         ! q after condensation
      real ztn(NHOR)         ! t after condensation
      real zdqdt(NHOR)       ! local q-tendency
      real zdtdt(NHOR)       ! local t-tendency
      real zqsat(NHOR)       ! saturation humidity
      real zlcpe(NHOR)       ! L/CP
      real zcor(NHOR)        ! correction term for saturation humidity
      real zqsto(NHOR)       ! saturation humidity lowered by stochastic
      real zrand(100)        ! array for random numbers
!
!     allocatable arrys for dbug diagnostic
!
      real, allocatable :: zprf1(:)
      real, allocatable :: zprf2(:)
      real, allocatable :: zprf3(:)
      real, allocatable :: zprf4(:)
      real, allocatable :: zprf5(:)
      real, allocatable :: zprf6(:)
      real, allocatable :: zprf7(:)

      if(nentropy > 0) then
       dentropy(:,11)=0.
      endif
      if(nentro3d > 0) then
       dentro3d(:,:,11)=0.
      endif
      if(nenergy > 0) then
       denergy(:,11)=0.
       denergy(:,15)=0.
      endif
      if(nener3d > 0) then
       dener3d(:,:,11)=0.
       dener3d(:,:,15)=0.
      endif
!
!*    initialize
!
      dprll(:,:)=0.
      dprsll(:,:)=0.

!
!*    make condensation
!

      do jlev=1,NLEV

!
!     preset tendencies
!

       zdqdt(:)=0.
       zdtdt(:)=0.

!
!     compute preliminary t and q at t+dt
!

       zq(:)=dq(:,jlev)+dqdt(:,jlev)*deltsec2
       zt(:)=dt(:,jlev)+dtdt(:,jlev)*deltsec2

!
!*    set l/cp
!

       zlcpe(:)=ALV/(acpd*(1.+ADV*dq(:,jlev)))

!
!     saturation humidity
!

       zqsat(:)=rdbrv*ra1*exp(ra2*(zt(:)-TMELT)/(zt(:)-ra4))            &
     &         /(dp(:)*sigma(jlev))
!
!      avoid negative qsat
!
       zqsat(:)=AMIN1(zqsat(:),rdbrv)
!
       zcor(:)=1./(1.-(1./rdbrv-1.)*zqsat(:))
       zqsat(:)=zqsat(:)*zcor(:)

!      Fix for problem with higher Mars atmosphere

       where (zqsat(:) < 0.0)
          zqsat(:) = 0.0
       endwhere

!      Add noise to precipitation

       zqsto(:) = zqsat(:) ! modifiable zqsat
       if (nstorain > 0) then 
         do jhor = 1, NHOR 
           call get_random_root(zrand,100) ! get 100 random numbers
           if (zq(jhor)<=zqsat(jhor).and.zq(jhor)>0.9*zqsat(jhor)) then
             zqdiffm = 0.0
             do j = 1 , 100
               zqdiff = (zrand(j)-0.5) * 0.2 * zqsat(jhor)
               if (zqdiff > zqsat(jhor)-zq(jhor)) then 
                 zqdiffm = zqdiffm + zqdiff
               endif
             enddo ! j
             zqsto(jhor) = zqsat(jhor) - zqdiffm * 0.01 ! / 100.0
           endif ! zqsat <= zq < 0.9 zqsat
         enddo ! jhor
       endif ! (nstorain > 0)

!     diagnostics: check # of cells with supersaturation

!     ic = 0
!     do jhor = 1 , NHOR
!        if (zq(jhor) > zqsat(jhor)) ic=ic+1
!     enddo
!     write (77,'(i5,2f12.4)') ic,sum(zq),sum(zqsat)

!
!     look for supersaturation
!
      where(zq(:) > zqsto(:))
!
!     make condensation (2 iterations)
!

        zdqdt(:)=(zqsto(:)-zq(:))                                       &
     &          /(1.0+zlcpe(:)*ra2*(TMELT-ra4)                          &
     &           *zqsat(:)*zcor(:)/(zt(:)-ra4)**2)
        ztn(:)=zt(:)-zdqdt(:)*zlcpe(:)
        zqn(:)=zq(:)+zdqdt(:)
        zqsat(:)=rdbrv*ra1*exp(ra2*(ztn(:)-TMELT)/(ztn(:)-ra4))         &
     &          /(dp(:)*sigma(jlev))
!
!      avoid negative qsat
!
        zqsat(:)=AMIN1(zqsat(:),rdbrv)
!
        zcor(:)=1./(1.-(1./rdbrv-1.)*zqsat(:))
        zqsat(:)=zqsat(:)*zcor(:)
        zdqdt(:)=(zqsat(:)-zqn(:))                                      &
     &          /(1.0+zlcpe(:)*ra2*(TMELT-ra4)                          &
     &           *zqsat(:)*zcor(:)/(ztn(:)-ra4)**2)
        ztn(:)=ztn(:)-zdqdt(:)*zlcpe(:)
        zqn(:)=zqn(:)+zdqdt(:)

!
!     final tendencies
!

        zdqdt(:)=(zqn(:)-zq(:))/deltsec2
        zdtdt(:)=(ztn(:)-zt(:))/deltsec2

!
!     add tendencies
!

        dqdt(:,jlev)=dqdt(:,jlev)+zdqdt(:)
        dtdt(:,jlev)=dtdt(:,jlev)+zdtdt(:)

!
!     add precipitation
!

        dprll(:,jlev)=-AMIN1(zdqdt(:),0.0)*dsigma(jlev)

       endwhere ! (zq(:) > zqsto(:))

!
!     dbug print out
!

       if(nprint==2) then
        allocate(zprf1(NLON*NLAT))
        allocate(zprf2(NLON*NLAT))
        allocate(zprf3(NLON*NLAT))
        allocate(zprf4(NLON*NLAT))
        allocate(zprf5(NLON*NLAT))
        allocate(zprf6(NLON*NLAT))
        allocate(zprf7(NLON*NLAT))
        call mpgagp(zprf1,zq,1)
        call mpgagp(zprf2,zt,1)
        call mpgagp(zprf3,zqsat,1)
        call mpgagp(zprf4,zdqdt,1)
        call mpgagp(zprf5,zdtdt,1)
        call mpgagp(zprf6,dprll(1,jlev),1)
        call mpgagp(zprf7,dp,1)
        if(mypid==NROOT) then
         if(jlev==1) write(nud,*)'In mklsp:'
         write(nud,*)'L= ',jlev,' zt= ',zprf2(nprhor)                        &
     &         ,' zq= ',zprf1(nprhor),' zqsat= ',zprf3(nprhor)
         write(nud,*)'L= ',jlev,' dtdt= ',zprf5(nprhor)                      &
     &         ,' dqdt= ',zprf4(nprhor)
         write(nud,*)'L= ',jlev,' prl (mm/day)= ',zprf6(nprhor)*zprf7(nprhor)&
     &                                      *0.5*deltsec2*mtspd/ga
        endif
        deallocate(zprf1)
        deallocate(zprf2)
        deallocate(zprf3)
        deallocate(zprf4)
        deallocate(zprf5)
        deallocate(zprf6)
        deallocate(zprf7)
       endif

!
!     franks dbug
!

       if(ndiaggp==1) then
        dgp3d(:,jlev,7)=zdtdt(:)
        dgp3d(:,jlev,15)=zdqdt(:)*dp(:)*dsigma(jlev)/ga/1000.
       endif
!
!     entropy/energy diagnostics
!
       if(nentropy > 0) then
        dentro(:)=zdtdt(:)/dentrot(:,jlev)                              &
     &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev) 
        dentropy(:,11)=dentropy(:,11)+dentro(:)
        if(nentro3d > 0) dentro3d(:,jlev,11)=dentro(:)
       endif
       if(nenergy > 0) then
        denergy(:,11)=denergy(:,11)                                     &
     &               +zdtdt(:)                                          &
     &               *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then 
         dener3d(:,jlev,11)=zdtdt(:)                                    &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif
        do jhor=1,NHOR
!        Phase choice for the latent heat released by condensation. Below the
!        melting point the vapour goes to ice and sublimation applies; above it
!        the vapour goes to liquid and vaporisation applies.
!
!        This was inverted, assigning als above TMELT and alv below. It sits
!        inside if(nenergy > 0) and feeds only denergy(:,15) and
!        dener3d(:,:,15), so the prognostic physics was never affected and no
!        completed run is wrong. But term 15 was, by als - alv at every
!        condensing gridpoint, which is exactly the quantity anyone enabling
!        these diagnostics is trying to measure.
!
!        The correct convention is used by the prognostic code in this same file
!        (rainmod.f90:731, 850, 944, 1047, all `if(ztnew < TMELT) zlcp=ALS`) and
!        by the surface latent heat flux in fluxmod.f90:687.
         if(zt(jhor) < TMELT) then
          zzal=als
         else
          zzal=alv
         endif
         denergy(jhor,15)=denergy(jhor,15)+zzal*zdqdt(jhor)             &
     &                                    *dsigma(jlev)*dp(jhor)/ga
         if(nener3d > 0) then
          dener3d(jhor,jlev,15)=zzal*zdqdt(jhor)                        &
     &                         *dsigma(jlev)*dp(jhor)/ga
         endif
        enddo
       endif
      enddo

      return
      end subroutine mklsp

!     ==============
!     SUBROUTINE KUO
!     ==============

      subroutine kuo
      use rainmod

!
!     local arrays
!

      integer ilift(NHOR)    ! lifting level (= cloud base + 1)
      integer itop(NHOR)     ! cloud top level
      integer ishallow(NHOR) ! flag for shallow convection
      integer ideep(NHOR)    ! flag for deep convection
      real zqe(NHOR,NLEV)    ! moisture of environmental air
      real zte(NHOR,NLEV)    ! temperatur of environmental air
      real zdqtot(NHOR,NLEV) ! total moisture accession in resp. layer
      real zdqdt(NHOR,NLEV)  ! moisture diff. (parcel-environment)
      real zdtdt(NHOR,NLEV)  ! temperature diff. (parcel-environment)
      real zalcpe(NHOR,NLEV) ! help array for factors
      real zale(NHOR,NLEV)   ! help array for factors
      real zacpe(NHOR,NLEV)  ! help array for factors
      real zqp(NHOR)         ! humidity of lifted parcel
      real ztp(NHOR)         ! temperature of lifted parcel
      real zi(NHOR)          ! moisture supply in the collumn
      real zpt(NHOR)         ! energy surplus of the cloud (t-part)
      real zpq(NHOR)         ! energy surplus of the cloud (q-part)
      real zbeta(NHOR)       ! beta factor (new kuo scheme)
      real zpbeta(NHOR)      ! p-scaling for beta factor
      real zat(NHOR)         ! distribution factor (temperature)
      real zaq(NHOR)         ! distribution factor (moisture)
      real zptop(NHOR)       ! cloud top pressure
      real zdqdts(NHOR,NLEV) ! dq/dt due to shallow convection
      real zdtdts(NHOR,NLEV) ! dt/dt due to shallow convection
      real zrhm(NHOR)        ! mean relative humidity (for beta)

      real zdqdtold(NHOR,NLEV)    ! old q-tendency (for energy diagn.)

!
!     array for momentum mixing
!
      real zue(NHOR,NLEV)    ! u of environmental air
      real zve(NHOR,NLEV)    ! v of environmental air
      real zdudt(NHOR,NLEV)  ! u-momentum difference (cloud-environment)
      real zdvdt(NHOR,NLEV)  ! v-momentum difference (cloud-environment)
      real zau(NHOR)         ! distribution factor (momentum)
      real zup(NHOR)         ! reference (cloud) u
      real zvp(NHOR)         ! reference (cloud) v
      real zsums(NHOR)       ! sum(dsigma) for mixed levels

!
!     allocatable arrays for dbug diagnostic
!

      real, allocatable :: zprf1(:,:)
      real, allocatable :: zprf2(:,:)
      real, allocatable :: zprf3(:,:)
      real, allocatable :: zprf4(:,:)
      real, allocatable :: zprf5(:)
      real, allocatable :: zprf6(:)
      real, allocatable :: zprf7(:,:)
      real, allocatable :: zprf8(:)
      real, allocatable :: zprf9(:)
      real, allocatable :: zprf10(:)
      real, allocatable :: zprf11(:)
      real, allocatable :: zprf12(:)
      real, allocatable :: zprf13(:)
      real, allocatable :: zprf14(:)

!
!*    0) initialize some fields
!

      ishallow(:)=0
      ilift(:)=-1
      itop(:)=-1
      icctot(:)=0
      zi(:)=0.
      zpt(:)=0.
      zpq(:)=0.
      zat(:)=0.
      zaq(:)=0.
      zbeta(:)=0.
      zpbeta(:)=0.
      zdtdt(:,:)=0.
      zdqdt(:,:)=0.
      icclev(:,:)=0
      dprcl(:,:)=0.
      dprscl(:,:)=0.
      zdtdts(:,:)=0.
      zdqdts(:,:)=0.
      zptop(:)=dp(:)

      if(nenergy > 0) then
       zdqdtold(:,:)=dqdt(:,1:NLEV)
      endif

!
!     moisture accession (level and column)
!

      zdqtot(:,:)=(dqt(:,1:NLEV)+dqdt(:,1:NLEV))*deltsec2

!
!     for original kuo (whole column) use the next loop
!
!     do jlev=1,NLEV
!      zi(:)=zi(:)+zdqtot(:,jlev)*dsigma(jlev)
!     enddo

!
!     franks dbug
!

      if(ndiaggp==1) then
       do jlev=1,NLEV
        dgp3d(:,jlev,16)=-dqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
        dgp3d(:,jlev,18)=dqt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      end if

!
!*    1) update q and t
!

      zqe(:,:)=dq(:,1:NLEV)+dqdt(:,1:NLEV)*deltsec2
      zte(:,:)=dt(:,1:NLEV)+dtdt(:,1:NLEV)*deltsec2

!
!     set L/CP (only rain here, phase transition in mkrain)
!

      zalcpe(:,:)=ALV/(acpd*(1.+ADV*dq(:,1:NLEV)))
      zacpe(:,:)=acpd*(1.+ADV*dq(:,1:NLEV))
      zale(:,:)=ALV

!
!*    2) search for lifting level
!
      if(ncsurf == 1) then

!     a) start at surface (lifting level only if dry unstable)
!

      zfac=sigma(NLEV)**akap
      do jhor=1,NHOR
       ztnew=dt(jhor,NLEP)*zfac
       if(zdqtot(jhor,NLEV) > 0. .and. ztnew > zte(jhor,NLEV)) then

!
!     surface air specific humidity
!

        zqsnl=rdbrv*ra1*exp(ra2*(zte(jhor,NLEV)-TMELT)                  &
     &                          /(zte(jhor,NLEV)-ra4))                  &
     &       /(dp(jhor)*sigma(NLEV))
!
!      avoid negative qsat
!
        zqsnl=AMIN1(zqsnl,rdbrv)
!
        zqsnl=zqsnl/(1.-(1./rdbrv-1.)*zqsnl)
        zqnew=dq(jhor,NLEP)*AMIN1(1.,zqe(jhor,NLEV)/zqsnl)

!
!     saturation humidity:
!

        zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))              &
     &       /(dp(jhor)*sigma(NLEV))
!
!      avoid negative qsat
!
        zqsnl=AMIN1(zqsnl,rdbrv)
!
        zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
        zqsnl=zqsnl*zcor

!
!     condensation if supersaturated
!

        if(zqsnl < zqnew) then

!
!     update constants (ice/water phase)
!

         if(ztnew < TMELT) then
          zlcp=ALS/(acpd*(1.+ADV*zqnew))
         else
          zlcp=ALV/(acpd*(1.+ADV*zqnew))
         endif

!
!     temperature correction term
!

         zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))


!
!     update t and q
!

         zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
         ztnew=ztnew+zlcp*zrain
         zqnew=zqnew-zrain

!
!     second iteration:
!

         zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))             &
     &        /(dp(jhor)*sigma(NLEV))
!
!      avoid negative qsat
!
         zqsnl=AMIN1(zqsnl,rdbrv)
!
         zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
         zqsnl=zqsnl*zcor
         zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))
         zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
         ztnew=ztnew+zlcp*zrain
         zqnew=zqnew-zrain

!
!     if convection: set lifting level and parcel q,t (zpq,ztp)
!                    calculate prelim. dq,dt (zdqdt,zdtdt)
!                    update surplus of total energy (zp)
!

         ilift(jhor)=NLEP
         itop(jhor)=NLEV
         ztp(jhor)=ztnew
         zqp(jhor)=zqnew
         zqsa=rdbrv*ra1*exp(ra2*(zte(jhor,NLEV)-TMELT)                  &
     &                          /(zte(jhor,NLEV)-ra4))                  &
     &       /(dp(jhor)*sigma(NLEV))
!
!      avoid negative qsat
!
         zqsa=AMIN1(zqsa,rdbrv)
!
         zqsa=zqsa/(1.-(1./rdbrv-1.)*zqsa)
         zdqdt(jhor,NLEV)=AMAX1(zqsa-zqe(jhor,NLEV),0.)
         zdtdt(jhor,NLEV)=ztp(jhor)-zte(jhor,NLEV)
         zpt(jhor)=zdtdt(jhor,NLEV)*dsigma(NLEV)*zacpe(jhor,NLEV)
         zpq(jhor)=zdqdt(jhor,NLEV)*dsigma(NLEV)*zale(jhor,NLEV)
         if(kbeta > 0 ) then
          zbeta(jhor)=AMIN1(1.,zqe(jhor,NLEV)/zqsa)*dsigma(NLEV)
          zpbeta(jhor)=dsigma(NLEV)
         endif

!
!     moisture accession (delete for original kuo)
!

         zi(jhor)=zdqtot(jhor,NLEV)*dsigma(NLEV)*zale(jhor,NLEV)
        endif
       endif
      enddo

      endif
!
!     b) continue search for free atmosphere
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       zfac=(sigma(jlev)/sigma(jlep))**akap
       do jhor=1,NHOR
        if(ilift(jhor) == -1                                            &
!    &    .and. zi(jhor) > 0.                   ! use for original kuo
     &    .and. zdqtot(jhor,jlep) > 0.) then

!
!     adiabatic uplift of parcel from below
!

         ztnew=zte(jhor,jlep)*zfac
         zqnew=zqe(jhor,jlep)

!
!     saturation humidity:
!

         zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))             &
     &        /(dp(jhor)*sigma(jlev))
!
!      avoid negative qsat
!
         zqsnl=AMIN1(zqsnl,rdbrv)
! 
         zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
         zqsnl=zqsnl*zcor

!
!     condensation if supersaturated
!

         if(zqsnl < zqnew) then

!
!     update constants (ice/water phase)
!

          if(ztnew < TMELT) then
           zlcp=ALS/(acpd*(1.+ADV*zqnew))
          else
           zlcp=ALV/(acpd*(1.+ADV*zqnew))
          endif

!
!     temperature correction term
!

          zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))


!
!     update t and q
!

          zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
          ztnew=ztnew+zlcp*zrain
          zqnew=zqnew-zrain

!
!     second iteration:
!

          zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))            &
     &         /(dp(jhor)*sigma(jlev))
!
!      avoid negative qsat
!
          zqsnl=AMIN1(zqsnl,rdbrv)
!
          zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
          zqsnl=zqsnl*zcor
          zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))
          zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
          ztnew=ztnew+zlcp*zrain
          zqnew=zqnew-zrain

!
!     check for buyoncy
!

          if(ztnew > zte(jhor,jlev)) then

!
!     if convection: set lifting level and parcel q,t (zpq,ztp)
!                    calculate prelim. dq,dt (zdqdt,zdtdt)
!                    update surplus of total energy (zp)
!

           ilift(jhor)=jlep
           itop(jhor)=jlev
           ztp(jhor)=ztnew
           zqp(jhor)=zqnew
           zqsa=rdbrv*ra1*exp(ra2*(zte(jhor,jlev)-TMELT)                &
     &                           /(zte(jhor,jlev)-ra4))                 &
     &         /(dp(jhor)*sigma(jlev))
!
!      avoid negative qsat
!
           zqsa=AMIN1(zqsa,rdbrv)
!
           zqsa=zqsa/(1.-(1./rdbrv-1.)*zqsa)
           zdqdt(jhor,jlev)=AMAX1(zqsa-zqe(jhor,jlev),0.)
           zdtdt(jhor,jlev)=ztp(jhor)-zte(jhor,jlev)
           zpt(jhor)=zdtdt(jhor,jlev)*dsigma(jlev)*zacpe(jhor,jlev)
           zpq(jhor)=zdqdt(jhor,jlev)*dsigma(jlev)*zale(jhor,jlev)
           if(kbeta > 0 ) then
            zbeta(jhor)=AMIN1(1.,zqe(jhor,jlev)/zqsa)*dsigma(jlev)
            zpbeta(jhor)=dsigma(jlev)
           endif

!
!          check if lifting level is supersaturated
!

           ztnew=zte(jhor,jlep)
           zqnew=zqe(jhor,jlep)
           zqsnl=rdbrv*ra1*exp(ra2*(zte(jhor,jlep)-TMELT)               &
     &                            /(zte(jhor,jlep)-ra4))                &
     &         /(dp(jhor)*sigma(jlep))
!
!      avoid negative qsat
!
           zqsnl=AMIN1(zqsnl,rdbrv)
!
           zqsnl=zqsnl/(1.-(1./rdbrv-1.)*zqsnl)
           if(zqe(jhor,jlep) > zqsnl) then

!
!     update constants (ice/water phase)
!

           if(ztnew < TMELT) then
            zlcp=ALS/(acpd*(1.+ADV*zqnew))
           else
            zlcp=ALV/(acpd*(1.+ADV*zqnew))
           endif

!
!     temperature correction term
!

           zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))

!
!     update t and q
!

           zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
           ztnew=ztnew+zlcp*zrain
           zqnew=zqnew-zrain

!
!     second iteration:
!

           zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))            &
     &          /(dp(jhor)*sigma(jlep))
!
!      avoid negative qsat
!
           zqsnl=AMIN1(zqsnl,rdbrv)
!
           zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
           zqsnl=zqsnl*zcor
           zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))
           zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
           ztnew=ztnew+zlcp*zrain
           zqnew=zqnew-zrain
           zdqdt(jhor,jlep)=AMAX1(zqnew-zqe(jhor,jlep),0.)
           zdtdt(jhor,jlep)=ztnew-zte(jhor,jlep)
           zpt(jhor)=zdtdt(jhor,jlep)*dsigma(jlep)*zacpe(jhor,jlep)    &
     &              +zpt(jhor)
           zpq(jhor)=zdqdt(jhor,jlep)*dsigma(jlep)*zale(jhor,jlep)        &
     &              +zpq(jhor)
           if(kbeta > 0 ) then
            zbeta(jhor)=AMIN1(1.,zqe(jhor,jlep)/zqsnl)*dsigma(jlep)    &
     &                  +zbeta(jhor)
            zpbeta(jhor)=dsigma(jlep)+zpbeta(jhor)
           endif
           endif

!
!     moisture accession (delete for original kuo)
!

           zi(jhor)=zdqtot(jhor,jlev)*dsigma(jlev)*zale(jhor,jlev)      &
     &             +zdqtot(jhor,jlep)*dsigma(jlep)*zale(jhor,jlep)
          endif
         endif
        endif
       enddo
      enddo

!
!*    3) search for top level
!

      do jlev=NLEM,1,-1
       jlep=jlev+1
       zfac=(sigma(jlev)/sigma(jlep))**akap
       do jhor=1,NHOR
        if(jlep == itop(jhor)) then

!
!     3a) parcel uplift
!

         ztnew=ztp(jhor)*zfac
         zqnew=zqp(jhor)

!
!     saturation humidity
!

         zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))             &
     &        /(dp(jhor)*sigma(jlev))
!
!      avoid negative qsat
!
         zqsnl=AMIN1(zqsnl,rdbrv)
!
         zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
         zqsnl=zqsnl*zcor

!
!     3b) more condensation
!

         if(zqsnl < zqnew) then

!
!     update constants (ice/water phase)
!

          if(ztnew < TMELT) then
           zlcp=ALS/(acpd*(1.+ADV*zqnew))
          else
           zlcp=ALV/(acpd*(1.+ADV*zqnew))
          endif

!
!     temperature correction term
!

          zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))

!
!     update t,q
!

          zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
          ztnew=ztnew+zlcp*zrain
          zqnew=zqnew-zrain

!
!     second iteration:
!

          zqsnl=rdbrv*ra1*exp(ra2*(ztnew-TMELT)/(ztnew-ra4))            &
     &         /(dp(jhor)*sigma(jlev))
!
!      avoid negative qsat
!
          zqsnl=AMIN1(zqsnl,rdbrv)
!
          zcor=1./(1.-(1./rdbrv-1.)*zqsnl)
          zqsnl=zqsnl*zcor
          zdqsdt=ra2*(TMELT-ra4)*zqsnl*zcor/((ztnew-ra4)*(ztnew-ra4))
          zrain=(zqnew-zqsnl)/(1.+zlcp*zdqsdt)
          ztnew=ztnew+zlcp*zrain
          zqnew=zqnew-zrain

!
!     3c) check for buyoncy
!

          if(ztnew > zte(jhor,jlev)) then

!
!     3d) if further convection: update parcel q,t (zpq,ztp)
!                                calculate prelim. dq,dt (zdqdt,zdtdt)
!                                update surplus of total energy (zp)
!

           itop(jhor)=jlev
           ztp(jhor)=ztnew
           zqp(jhor)=zqnew
           zqsa=rdbrv*ra1*exp(ra2*(zte(jhor,jlev)-TMELT)                &
     &                            /(zte(jhor,jlev)-ra4))                &
     &         /(dp(jhor)*sigma(jlev))
!
!      avoid negative qsat
!
           zqsa=AMIN1(zqsa,rdbrv)
!
           zqsa=zqsa/(1.-(1./rdbrv-1.)*zqsa)
           zdqdt(jhor,jlev)=AMAX1(zqsa-zqe(jhor,jlev),0.)
           zdtdt(jhor,jlev)=ztp(jhor)-zte(jhor,jlev)
           zpt(jhor)=zdtdt(jhor,jlev)*dsigma(jlev)*zacpe(jhor,jlev)    &
     &              +zpt(jhor)
           zpq(jhor)=zdqdt(jhor,jlev)*dsigma(jlev)*zale(jhor,jlev)     &
     &              +zpq(jhor)
           if(kbeta > 0 ) then
            zbeta(jhor)=AMIN1(1.,zqe(jhor,jlev)/zqsa)*dsigma(jlev)     &
     &                  +zbeta(jhor)
            zpbeta(jhor)=zpbeta(jhor)+dsigma(jlev)
           endif

!
!     moisture accession (delete for original kuo)
!

           zi(jhor)=zdqtot(jhor,jlev)*dsigma(jlev)*zale(jhor,jlev)  &
     &             +zi(jhor)
          endif
         endif
        endif
       enddo
      enddo

!
!*    4) calc. moisture distribution factors
!

      zi(:)=AMAX1(0.,zi(:))

      if(kbeta==0) then

!
!     kuo without beta
!

       where(zpt(:)+zpq(:) > 0.)
        zat(:)=zi(:)/(zpt(:)+zpq(:))
        zaq(:)=zat(:)
       end where
      else
!
!     kuo with beta
!
!     a) beta computed from rel. humidity       
!
       zrhm(:)=0.
       where(zpbeta(:) > 0.) zrhm(:)=zbeta(:)/zpbeta(:)
       zbeta(:)=1.
       where(zrhm(:) >= rhbeta) 
        zbeta(:)=((1.-zrhm(:))/(1.-rhbeta))**nbeta
       endwhere
!
!     b) beta give explicitely (for sensitivity studies, etc.)

       if(kbeta==2) zbeta(:)=rbeta        

       where(zpt(:) > 0.) zat(:)=zi(:)*(1.-zbeta(:))/zpt(:)
       where(zpq(:) > 0.) zaq(:)=zi(:)*zbeta(:)/zpq(:)

      endif

!
!*    5) update dqdt and dtdt, compute precip.
!
!
!     select points for shallow convection: top pressure needed
!     shallow convection if: 
!     a) not enought moisture convergence but unstable (zi=0)
!     b) convection to shallow (ptop > pdeep)
!
      if(nshallow == 1) then
       where(itop(:) > 0) zptop(:)=dp(:)*sigma(itop(:))
       where(itop(:) > 0 .and. zptop(:) > pdeep) ishallow(:)=1
       where(itop(:) > 0 .and. zi(:) <= 0.) ishallow(:)=1
      endif

      do jlev=1,NLEV
!
!     select points with deep convection
!
       ideep(:)=0
       where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.    &
     &      .and. ishallow(:) == 0.)
        ideep(:)=1
       endwhere
!
!     5a) in cloud: add tendencies due to convection and sub. advection
!         lifting level: sub. advection. (dtdt != 0 only if level supersat)
!
       where(ideep(:) > 0)
         dqdt(:,jlev)=zaq(:)*zdqdt(:,jlev)/deltsec2                     &
     &               -dqt(:,jlev)
         dtdt(:,jlev)=zat(:)*zdtdt(:,jlev)/deltsec2                     &
     &               +dtdt(:,jlev)
!
!     5b) calc. precipitation
!
         dprcl(:,jlev)=zat(:)*zdtdt(:,jlev)*dsigma(jlev)/zalcpe(:,jlev) &
     &                /deltsec2
       endwhere
!
!     5c) original kuo:
!         rest of column: subtrac all q-tendencies (i.e. reset dqdt)
!
!      where(zi(:) > 0 .and. (jlev >= ilift(:) .or. jlev < itop(:)))
!       dqdt(:,jlev)=-1.*dqt(:,jlev)
!      endwhere

!
!     5d) flag and counter for clouds:
!

       where(dprcl(:,jlev) > 0.)
        icclev(:,jlev)=1
        icctot(:)=icctot(:)+1
       endwhere
       if(sigma(jlev) < 0.4) then
        where(icclev(:,jlev)==1 .and. itop(:)==jlev)
         icclev(:,jlev)=2
        endwhere
       endif

      enddo
!
!*    make shallow convection
!
      if(nshallow == 1) then
       call mkshallow(zte,zqe,ishallow,ilift,itop,zdtdts,zdqdts)
       dtdt(:,1:NLEV)=dtdt(:,1:NLEV)+zdtdts(:,1:NLEV)
       dqdt(:,1:NLEV)=dqdt(:,1:NLEV)+zdqdts(:,1:NLEV)
      endif
!
!*    6) momentum mixing
!

      if(nmoment == 1) then

!
!     6a) update u and v
!

       zue(:,:)=du(:,1:NLEV)+dudt(:,1:NLEV)*deltsec2
       zve(:,:)=dv(:,1:NLEV)+dvdt(:,1:NLEV)*deltsec2

!
!     6a) make distribution factor
!

       zau(:)=0.
       where(zpt(:)+zpq(:) > 0.)
        zau(:)=AMIN1(1.,zi(:)/(zpt(:)+zpq(:)))
       end where

!
!     6b) reference u and v
!

       zup(:)=0.
       zvp(:)=0.
       zsums(:)=0.
       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.   &
     &       .and. ishallow(:) == 0)
         zup(:)=zup(:)+zue(:,jlev)*dsigma(jlev)
         zvp(:)=zvp(:)+zve(:,jlev)*dsigma(jlev)
         zsums(:)=zsums(:)+dsigma(jlev)
        endwhere
       enddo
       where(zsums(:) > 0.)
        zup(:)=zup(:)/zsums(:)
        zvp(:)=zvp(:)/zsums(:)
       endwhere

!
!     6c) make and add tendencies
!

       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.   &
     &       .and. ishallow(:) == 0)
          zdudt(:,jlev)=zau(:)*(zup(:)-zue(:,jlev))
          zdvdt(:,jlev)=zau(:)*(zvp(:)-zve(:,jlev))
          dudt(:,jlev)=dudt(:,jlev)+zdudt(:,jlev)/deltsec2
          dvdt(:,jlev)=dvdt(:,jlev)+zdvdt(:,jlev)/deltsec2
        endwhere
       enddo

      endif
!
!*    dbug diagnostic for nprint==2 (see pumamod)
!

      if(nprint==2) then
       allocate(zprf1(NLON*NLAT,NLEV))
       allocate(zprf2(NLON*NLAT,NLEV))
       allocate(zprf3(NLON*NLAT,NLEV))
       allocate(zprf4(NLON*NLAT,NLEV))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       allocate(zprf7(NLON*NLAT,NLEV))
       allocate(zprf8(NLON*NLAT))
       allocate(zprf9(NLON*NLAT))
       allocate(zprf10(NLON*NLAT))
       allocate(zprf13(NLON*NLAT))
       allocate(zprf11(NHOR))
       allocate(zprf12(NHOR))
       allocate(zprf14(NHOR))
       zprf11(:)=real(ilift(:))
       zprf12(:)=real(itop(:))
       zprf14(:)=real(ishallow(:))
       call mpgagp(zprf1,zdqdt,NLEV)
       call mpgagp(zprf2,zdtdt,NLEV)
       call mpgagp(zprf3,zdqtot,NLEV)
       call mpgagp(zprf4,dqt,NLEV)
       call mpgagp(zprf5,zprf11,1)
       call mpgagp(zprf6,zprf12,1)
       call mpgagp(zprf13,zprf14,1)
       call mpgagp(zprf7,dprcl,NLEV)
       call mpgagp(zprf8,zaq,1)
       call mpgagp(zprf9,zat,1)
       call mpgagp(zprf10,zi,1)
       deallocate(zprf11)
       deallocate(zprf12)
       deallocate(zprf14)
       if(mypid==NROOT) then
        write(nud,*)'In Kuo:'
        write(nud,*)'ishallow= ',zprf13(nprhor)
        write(nud,*)'ilift= ',zprf5(nprhor),' itop= ',zprf6(nprhor)
        do jlev=1,NLEV
         write(nud,*)'L= ',jlev,' dqtot= ',zprf3(nprhor,jlev)                &
     &                    ,' dqt= ',zprf4(nprhor,jlev)
         if(zprf5(nprhor) >= jlev .and. zprf6(nprhor) <= jlev           &
     &     .and. zprf10(nprhor) > 0. .and. zprf13(nprhor) == 0) then
          zzdt=zprf9(nprhor)*zprf2(nprhor,jlev)/deltsec2
          zzdq=zprf8(nprhor)*zprf1(nprhor,jlev)/deltsec2
          write(nud,*)'L= ',jlev,' dtdt= ',zzdt,' dqdt= ',zzdq
          write(nud,*)'L= ',jlev,' prc (mm/day)= '                           &
     &                     ,zprf7(nprhor,jlev)*1000.*mtspd*deltsec2*0.5
         endif
        enddo
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
       deallocate(zprf6)
       deallocate(zprf7)
       deallocate(zprf8)
       deallocate(zprf9)
       deallocate(zprf10)
       deallocate(zprf13)
      endif
!
!     franks dbug
!
      if(ndiaggp==1) then
       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.   &
     &       .and. ishallow(:) == 0)
         dgp3d(:,jlev,8)=zat(:)*zdtdt(:,jlev)/deltsec2
         dgp3d(:,jlev,16)=(zaq(:)*zdqdt(:,jlev)-zdqtot(:,jlev))/deltsec2&
     &                   *dp(:)*dsigma(jlev)/ga/1000.
        elsewhere
         dgp3d(:,jlev,8)=0.
         dgp3d(:,jlev,16)=0.
        endwhere
        if(nshallow > 0) then
         dgp3d(:,jlev,8)=dgp3d(:,jlev,8)+zdtdts(:,jlev)
         dgp3d(:,jlev,16)=dgp3d(:,jlev,16)+zdqdts(:,jlev)
        endif
       enddo
      end if
!
!     entropy/energy diagnostics
!
      if(nentropy > 0) then
       dentropy(:,12)=0.
       if(nentro3d > 0) dentro3d(:,:,12)=0.
       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.   &
     &       .and.ishallow(:) == 0)
         dentro(:)=+zat(:)*zdtdt(:,jlev)/deltsec2/dentrot(:,jlev)       &
     &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)
         dentropy(:,12)=dentropy(:,12)+dentro(:)
        endwhere
        if(nentro3d > 0) then
         where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.  &
     &        .and.ishallow(:) == 0)
          dentro3d(:,jlev,12)=dentro(:)
         endwhere
        endif
        if(nshallow == 1) then
         dentro(:)=zdtdts(:,jlev)/dentrot(:,jlev)                       &
     &        *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)
         dentropy(:,12)=dentropy(:,12)+dentro(:)
         if(nentro3d > 0) dentro3d(:,jlev,12)=dentro3d(:,jlev,12)       &
     &                                       +dentro(:)
        endif
       enddo
      endif
      if(nenergy > 0) then
       denergy(:,12)=0.
       denergy(:,16)=0.
       do jlev=1,NLEV
        where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.   &
     &       .and.ishallow(:) == 0)
         denergy(:,12)=denergy(:,12)                                    &
     &                +zat(:)*zdtdt(:,jlev)/deltsec2                    &
     &                *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endwhere
        denergy(:,16)=denergy(:,16)                                     &
     &               +alv*(dqdt(:,jlev)-zdqdtold(:,jlev))               &
     &               *dsigma(jlev)*dp(:)/ga
        if(nshallow == 1) then
         denergy(:,12)=denergy(:,12)                                    &
     &                +zdtdts(:,jlev)/deltsec2                          &
     &                *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
         denergy(:,16)=denergy(:,16)                                    &
     &                +alv*zdqdts(:,jlev)                               &
     &                *dsigma(jlev)*dp(:)/ga
        endif
       enddo
       if(nener3d > 0) then
        dener3d(:,:,12)=0.
        dener3d(:,:,16)=0.
        do jlev=1,NLEV
         where(ilift(:) >= jlev .and. itop(:) <= jlev .and. zi(:) > 0.  &
     &        .and.ishallow(:) == 0)
          dener3d(:,jlev,12)=zat(:)*zdtdt(:,jlev)/deltsec2              &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
         endwhere
         dener3d(:,jlev,16)=alv*(dqdt(:,jlev)-zdqdtold(:,jlev))         &
     &                     *dsigma(jlev)*dp(:)/ga
         if(nshallow == 1) then
          dener3d(:,jlev,12)=dener3d(:,jlev,12)                         &
     &                      +zdtdts(:,jlev)/deltsec2                    &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
          dener3d(:,jlev,16)=dener3d(:,jlev,16)                         &
     &                      +alv*zdqdts(:,jlev)                         &
     &                      *dsigma(jlev)*dp(:)/ga
         endif
        enddo
       endif
      endif
!
      return
      end subroutine kuo

!     ====================
!     SUBROUTINE MKSHALLOW
!     ====================

      subroutine mkshallow(pt,pq,kshallow,klift,ktop,pdtdt,pdqdt)
      use rainmod
!
!     calculate t,q tendencies due to shallow convection
!     (i.e. vertical diffusion using the ECHAM semi-implicit scheme)
!
      integer kshallow(NHOR)
      integer klift(NHOR)
      integer ktop(NHOR)
      integer kktop(NHOR)
!
      real pt(NHOR,NLEV)
      real pq(NHOR,NLEV)
      real pdtdt(NHOR,NLEV)
      real pdqdt(NHOR,NLEV)
!
      real zdtdt(NHOR,NLEV)
      real zdudt(NHOR,NLEV)
      real zdvdt(NHOR,NLEV)
!
      real zkdiff(NHOR,NLEV)
      real ztn(NHOR,NLEV)
      real zqn(NHOR,NLEV)
      real zun(NHOR,NLEV)
      real zvn(NHOR,NLEV)
      real zt(NHOR,NLEV)
      real zq(NHOR,NLEV)
      real zu(NHOR,NLEV)
      real zv(NHOR,NLEV)
      real zebs(NHOR,NLEM)
      real zskap(NLEV),zskaph(NLEV)
      real zke(NHOR,NLEV),zken(NHOR,NLEV)
!
      real, allocatable :: zprf1(:,:)
      real, allocatable :: zprf2(:,:)
      real, allocatable :: zprf3(:,:)
      real, allocatable :: zprf4(:,:)
      real, allocatable :: zprf5(:)
      real, allocatable :: zprf6(:)
      real, allocatable :: zprf7(:)
      real, allocatable :: zprf8(:)
      real, allocatable :: zprf9(:)
      real, allocatable :: zprf10(:)
!
      zkonst1=ga*deltsec2/gascon
      zkonst2=rkshallow*zkonst1*ga/gascon
!
      zskap(:)=sigma(:)**akap
      zskaph(:)=sigmah(:)**akap
!
      zt(:,:)=pt(:,1:NLEV)
      zq(:,:)=pq(:,1:NLEV)
      pdqdt(:,:)=0.
      pdtdt(:,:)=0.
!
!     modified ktop according to pdeepth
!
      kktop(:)=1
      do jlev=2,NLEM
       where(kshallow(:) > 0 .and. dp(:)*sigma(jlev) < pdeepth) 
        kktop(:)=jlev
       endwhere
      enddo
      ktop(:)=MAX(ktop(:),kktop(:))
!
!     0) modified diffusion coefficents
!
      zkdiff(:,:)=0.
      do jlev=1,NLEM
       jlp=jlev+1
       zrdsig=1./(sigma(jlp)-sigma(jlev))
       zfac=zkonst2*sigmah(jlev)*sigmah(jlev)*zrdsig
       do jhor=1,NHOR
!
!     within cloud + cloud base
!
        if(kshallow(jhor) > 0 .and. jlev >= ktop(jhor)                  &
     &    .and. jlev < klift(jhor)) then
         zrthl=(dsigma(jlp)+dsigma(jlev))                               &
     &        /(zt(jhor,jlev)*dsigma(jlp)+zt(jhor,jlp)*dsigma(jlev))
         zkdiff(jhor,jlev)=zfac*zrthl*zrthl
        endif
!
!     top entrainment
!
        if(kshallow(jhor) > 0 .and. jlev+1 == ktop(jhor)) then
         zrthl=(dsigma(jlp)+dsigma(jlev))                               &
     &        /(zt(jhor,jlev)*dsigma(jlp)+zt(jhor,jlp)*dsigma(jlev))
         zzqs=rdbrv*ra1*exp(ra2*(zt(jhor,jlev)-TMELT)                   &
     &                     /(zt(jhor,jlev)-ra4))                        &
     &       /(dp(jhor)*sigma(jlev))
!
!      avoid negative qsat
!
         zzqs=AMIN1(zzqs,rdbrv)
!
         zcor=1./(1.-(1./rdbrv-1.)*zzqs)
         zzqs=zzqs*zcor
         zzqsp=rdbrv*ra1*exp(ra2*(zt(jhor,jlp)-TMELT)                   &
     &                     /(zt(jhor,jlp)-ra4))                         &
     &       /(dp(jhor)*sigma(jlp))
!
!      avoid negative qsat
!
         zzqsp=AMIN1(zzqsp,rdbrv)
!
         zcor=1./(1.-(1./rdbrv-1.)*zzqsp)
         zzqsp=zzqsp*zcor
         zrh=zq(jhor,jlev)/zzqs
         zrhp=zq(jhor,jlp)/zzqsp
         zzf=ABS((zrh-0.8)/0.2*(zrhp-zrh))
         zkdiff(jhor,jlev)=zfac*zrthl*zrthl*zzf
        endif
       enddo
      enddo
!
!     1. semi implicit scheme
!
!     1a moisture diffusion
!
!     top layer elimination
!
      where(kshallow(:) > 0)
       zebs(:,1)=zkdiff(:,1)/(dsigma(1)+zkdiff(:,1))
       zqn(:,1)=dsigma(1)*zq(:,1)/(dsigma(1)+zkdiff(:,1))
      endwhere
!
!     middle layer elimination
!
      do jlev=2,NLEM
       jlem=jlev-1
       where(kshallow(:) > 0)
        zebs(:,jlev)=zkdiff(:,jlev)                                     &
     &              /(dsigma(jlev)+zkdiff(:,jlev)                       &
     &               +zkdiff(:,jlem)*(1.-zebs(:,jlem)))
        zqn(:,jlev)=(zq(:,jlev)*dsigma(jlev)                            & 
     &              +zkdiff(:,jlem)*zqn(:,jlem))                        &
     &             /(dsigma(jlev)+zkdiff(:,jlev)                        &
     &              +zkdiff(:,jlem)*(1.-zebs(:,jlem)))
       endwhere
      enddo
!
!     bottom layer elimination
!
      where(kshallow(:) > 0)
       zqn(:,NLEV)=(zq(:,NLEV)*dsigma(NLEV)+zkdiff(:,NLEM)*zqn(:,NLEM)) &
     &            /(dsigma(NLEV)+zkdiff(:,NLEM)*(1.-zebs(:,NLEM)))
      endwhere
!
!     back-substitution
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       where(kshallow(:) > 0)
        zqn(:,jlev)=zqn(:,jlev)+zebs(:,jlev)*zqn(:,jlep)
       endwhere
      enddo
!
!     tendencies
!
      do jlev=1,NLEV
       where(kshallow(:) > 0)      
        pdqdt(:,jlev)=(zqn(:,jlev)-zq(:,jlev))/deltsec2
       endwhere
      enddo
!
!     2b potential temperature
!
      do jlev=1,NLEM
       where(kshallow(:) > 0)
        zkdiff(:,jlev)=zkdiff(:,jlev)*zskaph(jlev)
       endwhere
      enddo
!
!     semi implicit scheme
!
!     top layer elimination
!
      where(kshallow(:) > 0)
       zebs(:,1)=zkdiff(:,1)/(dsigma(1)+zkdiff(:,1)/zskap(1))
       ztn(:,1)=dsigma(1)*zt(:,1)/(dsigma(1)+zkdiff(:,1)/zskap(1))
      endwhere
!
!     middle layer elimination
!
      do jlev=2,NLEM
       jlem=jlev-1
       where(kshallow(:) > 0)
        zebs(:,jlev)=zkdiff(:,jlev)                                     &
     &              /(dsigma(jlev)+(zkdiff(:,jlev)                      &
     &               +zkdiff(:,jlem)*(1.-zebs(:,jlem)/zskap(jlem)))     &
     &               /zskap(jlev))
        ztn(:,jlev)=(zt(:,jlev)*dsigma(jlev)                            &
     &              +zkdiff(:,jlem)/zskap(jlem)*ztn(:,jlem))            &
     &             /(dsigma(jlev)+(zkdiff(:,jlev)                       &
     &              +zkdiff(:,jlem)*(1.-zebs(:,jlem)/zskap(jlem)))      &
     &              /zskap(jlev))
       endwhere
      enddo
!
!     bottom layer elimination
!
      where(kshallow(:) > 0) 
       ztn(:,NLEV)=(zt(:,NLEV)*dsigma(NLEV)                             &
     &             +zkdiff(:,NLEM)*ztn(:,NLEM)/zskap(NLEM))             &
     &            /(dsigma(NLEV)+zkdiff(:,NLEM)/zskap(NLEV)             &
     &                          *(1.-zebs(:,NLEM)/zskap(NLEM)))
      endwhere
!
!     back-substitution
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       where(kshallow(:) > 0)
        ztn(:,jlev)=ztn(:,jlev)+zebs(:,jlev)*ztn(:,jlep)/zskap(jlep)
       endwhere
      enddo
!
!     tendencies
!
      do jlev=1,NLEV
       where(kshallow(:) > 0)
        pdtdt(:,jlev)=(ztn(:,jlev)-zt(:,jlev))/deltsec2
       endwhere
      enddo

      if(nmoment == 1) then
!
!     3 momentum
!
       zu(:,:)=du(:,1:NLEV)
       zv(:,:)=dv(:,1:NLEV)
!
!     top layer elimination
!
       where(kshallow(:) > 0)
        zebs(:,1)=zkdiff(:,1)/(dsigma(1)+zkdiff(:,1))
        zun(:,1)=dsigma(1)*zu(:,1)/(dsigma(1)+zkdiff(:,1))
        zvn(:,1)=dsigma(1)*zv(:,1)/(dsigma(1)+zkdiff(:,1))
       endwhere
!
!     middle layer elimination
!
       do jlev=2,NLEM
        jlem=jlev-1
        where(kshallow > 0)
         zebs(:,jlev)=zkdiff(:,jlev)                                     &
     &               /(dsigma(jlev)+zkdiff(:,jlev)                       &
     &              +zkdiff(:,jlem)*(1.-zebs(:,jlem)))
         zun(:,jlev)=(zu(:,jlev)*dsigma(jlev)+zkdiff(:,jlem)*zun(:,jlem))&
     &              /(dsigma(jlev)+zkdiff(:,jlev)                        &
     &             +zkdiff(:,jlem)*(1.-zebs(:,jlem)))
         zvn(:,jlev)=(zv(:,jlev)*dsigma(jlev)+zkdiff(:,jlem)*zvn(:,jlem))&
     &              /(dsigma(jlev)+zkdiff(:,jlev)                        &
     &               +zkdiff(:,jlem)*(1.-zebs(:,jlem)))
        endwhere
       enddo
!
!     bottom layer elimination
! 
       where(kshallow(:) > 0)
        zun(:,NLEV)=(zu(:,NLEV)*dsigma(NLEV)+zkdiff(:,NLEM)*zun(:,NLEM))&
     &             /(dsigma(NLEV)+zkdiff(:,NLEM)*(1.-zebs(:,NLEM)))
        zvn(:,NLEV)=(zv(:,NLEV)*dsigma(NLEV)+zkdiff(:,NLEM)*zvn(:,NLEM))&
     &             /(dsigma(NLEV)+zkdiff(:,NLEM)*(1.-zebs(:,NLEM)))
       endwhere
!
!     back-substitution
!
       do jlev=NLEM,1,-1
        jlep=jlev+1
        where(kshallow(:) > 0)
         zun(:,jlev)=zun(:,jlev)+zebs(:,jlev)*zun(:,jlep)
         zvn(:,jlev)=zvn(:,jlev)+zebs(:,jlev)*zvn(:,jlep)
        endwhere
       enddo
!
!     tendencies
!
       do jlev=1,NLEV
        where(kshallow(:) > 0)
         zdudt(:,jlev)=(zun(:,jlev)-zu(:,jlev))/deltsec2
         zdvdt(:,jlev)=(zvn(:,jlev)-zv(:,jlev))/deltsec2
         dudt(:,jlev)=dudt(:,jlev)+zdudt(:,jlev)
         dvdt(:,jlev)=dvdt(:,jlev)+zdvdt(:,jlev)
        endwhere
       enddo
!
!     temp.-tendencies due to the loss of kin.energy (energy-conservation)
!
       if(ndheat > 0) then
!
!     compute diffusion of ekin (needs to be substracted)
!
        do jlev=1,NLEV
         where(kshallow(:) > 0)
          zke(:,jlev)=0.5*(zu(:,jlev)*zu(:,jlev)+zv(:,jlev)*zv(:,jlev))
         endwhere
        enddo
        where(kshallow(:) > 0)
         zebs(:,1)=zkdiff(:,1)/(dsigma(1)+zkdiff(:,1))
         zken(:,1)=dsigma(1)*zke(:,1)/(dsigma(1)+zkdiff(:,1))
        endwhere
        do jlev=2,NLEM
         jlem=jlev-1
         where(kshallow(:) > 0)
          zebs(:,jlev)=zkdiff(:,jlev)                                   &
     &                /(dsigma(jlev)+zkdiff(:,jlev)                     &
     &                 +zkdiff(:,jlem)*(1.-zebs(:,jlem)))
          zken(:,jlev)=(zke(:,jlev)*dsigma(jlev)                        &
     &                 +zkdiff(:,jlem)*zken(:,jlem))                    &
     &                /(dsigma(jlev)+zkdiff(:,jlev)                     &
     &                 +zkdiff(:,jlem)*(1.-zebs(:,jlem)))
         endwhere
        enddo
        where(kshallow(:) > 0)
         zken(:,NLEV)=(zke(:,NLEV)*dsigma(NLEV)                         &
     &                +zkdiff(:,NLEM)*zken(:,NLEM))                     &
     &               /(dsigma(NLEV)+zkdiff(:,NLEM)*(1.-zebs(:,NLEM)))
        endwhere
        do jlev=NLEM,1,-1
         jlep=jlev+1
         where(kshallow(:) > 0)
          zken(:,jlev)=zken(:,jlev)+zebs(:,jlev)*zken(:,jlep)
         endwhere
        enddo
!
!     temperature tendencies du to ekin dissipation
!
        do jlev=1,NLEV
         where(kshallow(:) > 0)
          zdtdt(:,jlev)=-((zun(:,jlev)*zun(:,jlev)                      &
     &                    -zu(:,jlev)*zu(:,jlev)                        &
     &                    +zvn(:,jlev)*zvn(:,jlev)                      &
     &                    -zv(:,jlev)*zv(:,jlev))/deltsec2              &
     &                    -(zken(:,jlev)-zke(:,jlev))/deltsec2)         &
     &                 *0.5/acpd/(1.+adv*dq(:,jlev))
          pdtdt(:,jlev)=pdtdt(:,jlev)+zdtdt(:,jlev)
          zdtdt(:,jlev)=0.
         endwhere
        enddo
       endif
      endif
!
!*    dbug diagnostic for nprint==2 (see pumamod)
!
      if(nprint==2) then
       allocate(zprf1(NLON*NLAT,NLEV))
       allocate(zprf2(NLON*NLAT,NLEV))
       allocate(zprf3(NLON*NLAT,NLEV))
       allocate(zprf4(NLON*NLAT,NLEV))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       allocate(zprf7(NLON*NLAT))
       allocate(zprf8(NHOR))
       allocate(zprf9(NHOR))
       allocate(zprf10(NHOR))
       zprf8(:)=real(klift(:))
       zprf9(:)=real(ktop(:))
       zprf10(:)=real(kshallow(:))
       call mpgagp(zprf1,pt,NLEV)
       call mpgagp(zprf2,pq,NLEV)
       call mpgagp(zprf3,pdtdt,NLEV)
       call mpgagp(zprf4,pdqdt,NLEV)
       call mpgagp(zprf5,zprf8,1)
       call mpgagp(zprf6,zprf9,1)
       call mpgagp(zprf7,zprf10,1)
       deallocate(zprf8)
       deallocate(zprf9)
       deallocate(zprf10)
       if(mypid==NROOT) then
        write(nud,*)'In mkshallow:'
        write(nud,*)'kshallow= ',zprf7(nprhor),' total= ',SUM(nint(zprf7))
        write(nud,*)'ilift= ',zprf5(nprhor),' itop= ',zprf6(nprhor)
        do jlev=1,NLEV
         zztn=zprf1(nprhor,jlev)+zprf3(nprhor,jlev)*deltsec2
         zzqn=zprf2(nprhor,jlev)+zprf4(nprhor,jlev)*deltsec2
         write(nud,*)'L= ',jlev,' t= ',zprf1(nprhor,jlev)                    &
     &                    ,' q= ',zprf2(nprhor,jlev)                    
         write(nud,*)'L= ',jlev,' tnew= ',zztn,' dtdt= ',zprf3(nprhor,jlev)  &
     &                    ,' qnew= ',zzqn,' dqdt= ',zprf4(nprhor,jlev)  
        enddo
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
       deallocate(zprf6)
       deallocate(zprf7)
      endif
!
      return
      end subroutine mkshallow

!     ===================
!     SUBROUTINE MKCLOUDS
!     ===================

      subroutine mkclouds
      use rainmod
!
!     compute cloud properties (cover and liquid water content)
!

      parameter(zclmax=1.)
      parameter(zcca=0.245,zccb=0.125)
      parameter(zccmax=0.8,zccmin=0.05)

      real zq(NHOR),zt(NHOR),zqsat(NHOR),zrh(NHOR)
      real zrcrit(NHOR),zccconv(NHOR),zcctot(NHOR)
      real zwfac(NHOR),zcc(NHOR)
      real zzf(NHOR,NLEV),zzh(NHOR),zdh(NHOR)

!
!*    preset
!

      dql(:,:)=0.
      dcc(:,1:NLEV)=0.
      zcctot(:)=zccmin
      zccconv(:)=zccmin

!
!*    convective cloud cover for random overlab:
!

      zrfac = solar_day * 1000.0 ! convert m/s into mm/day
      where(dprc(:) > 0.)
       zcctot(:)=zcca+zccb*log(dprc(:)*zrfac)
       zcctot(:)=AMIN1(zccmax,AMAX1(zccmin,zcctot(:)))
       zccconv(:)=1.-(1.-zcctot(:))**(1./real(icctot(:)))
      end where

!
!*    cloud cover for the levels
!

      do jlev=1,NLEV
       zt(:)=dt(:,jlev)+dtdt(:,jlev)*deltsec2
       zq(:)=dq(:,jlev)+dqdt(:,jlev)*deltsec2
       zqsat(:)=rdbrv*ra1*exp(ra2*(zt(:)-TMELT)/(zt(:)-ra4))            &
     &         /(dp(:)*sigma(jlev))
!
!      avoid negative qsat
!
       zqsat(:)=AMIN1(zqsat(:),rdbrv)
!
       zqsat(:)=zqsat(:)/(1.-(1./rdbrv-1.)*zqsat(:))
       zrh(:)=zq(:)/zqsat(:)
!
!     convective clouds
!
!
!     a) cu/cb
!
       where(icclev(:,jlev) > 0)
        dcc(:,jlev)=zccconv(:)
       endwhere
!
!     b) ci
!
       where(icclev(:,jlev)==2 .and. zcctot(:) >= 0.3)
        dcc(:,jlev)=AMIN1(dcc(:,jlev)+2.*(zcctot(:)-0.3),zclmax)
       endwhere
!
!     reduce relative humidity
!
        zrh(:)=zrh(:)*(1.-dcc(:,jlev))
!
!     c) non convective clouds
!
       zrcrit(:)=rcrit(jlev)
       zwfac(:)=1.
       if(sigma(jlev) > 0.8 .and. clwfac > 0.)  then
        zwfac(:)=clwfac*AMAX1(0.,clwcrit1-dw(:,jlev))
        where(clwcrit2 > dw(:,jlev)) zwfac(:)=1.
       endif
       zcc(:)=zwfac(:)*AMAX1(0.,(zrh(:)-zrcrit(:))/(1.-zrcrit(:)))**2
       dcc(:,jlev)=AMIN1(dcc(:,jlev)+zcc(:),zclmax)

      enddo

!
!*    precipitable water (and z which is later used to compute dql)
!

      dqvi(:)=0.
      zzh(:)=0.
      do jlev=NLEV,2,-1
       dqvi(:)=dq(:,jlev)*dsigma(jlev)+dqvi(:)
       zdh(:)=-dt(:,jlev)*gascon/ga*ALOG(sigmah(jlev-1)/sigmah(jlev))
       zzf(:,jlev)=zzh(:)+zdh(:)*0.5
       zzh(:)=zzh(:)+zdh(:)
      enddo
      dqvi(:)=dq(:,1)*dsigma(1)+dqvi(:)
      zdh(:)=-dt(:,1)*gascon/ga*ALOG(sigma(1)/sigmah(1))*0.5
      zzf(:,1)=zzh(:)+zdh(:)

!
!     cloud liquid water (see CCM3 description: Kiehl et al. 1996 pp49+50)
!     and total cloud cover for random overlap
!

      dcc(:,NLEP)=1.
      dqvi(:)=dqvi(:)*dp(:)/ga
      zzh(:)=700.*ALOG(1.+dqvi(:))
      do jlev=1,NLEV
       where(zzh(:) > 0. .and. dcc(:,jlev) > 0.)
        dql(:,jlev)=0.00021*EXP(-zzf(:,jlev)/zzh(:))*gascon*dt(:,jlev)  &
     &             /(sigma(jlev)*dp(:))
        dql(:,jlev)=MAX(dql(:,jlev),1.E-9)
       endwhere
       dcc(:,NLEP)=dcc(:,NLEP)*(1.-dcc(:,jlev))
      enddo
      dcc(:,NLEP)=1.-dcc(:,NLEP)
!
      return
      end subroutine mkclouds

!     ==============
!     SUBROUTINE DCA
!     ==============

      subroutine mkdca
      use rainmod
!
!     calculate t and q tendencies due to dry convection
!
      integer itop(NHOR),ibase(NHOR),icon(NHOR)
      real zdtdt(NHOR,NLEV)
      real zdqdt(NHOR,NLEV)
      real zt(NHOR,NLEV),zth(NHOR,NLEV)
      real zq(NHOR,NLEV),zqn(NHOR,NLEV)
      real ztht(NHOR)
      real zqt(NHOR)
      real zsum1(NHOR),zsum2(NHOR)
      real zsumq1(NHOR),zsumq2(NHOR)
      real zskap(NLEV)
!
      do jlev = 1,NLEV
       zskap(jlev)=sigma(jlev)**akap
       zt(:,jlev)=dt(:,jlev)+dtdt(:,jlev)*deltsec2
       zth(:,jlev)=zt(:,jlev)/zskap(jlev)
       zq(:,jlev)=dq(:,jlev)+dqdt(:,jlev)*deltsec2
       zqn(:,jlev)=zq(:,jlev)
      enddo
!
      jiter=0
 1000 continue
       icon(:)=0
       ibase(:)=0
       itop(:)=NLEV
       ztht(:)=0.
       zqt(:)=0.
       do jlev=1,NLEM
        jlep=jlev+1
        where(zth(:,jlev) < zth(:,jlep))
         ibase(:)=jlep
         itop(:)=jlev-1
         zsum1(:)=zskap(jlep)*dsigma(jlep)                              &
     &              +zskap(jlev)*dsigma(jlev)
         zsum2(:)=zskap(jlep)*dsigma(jlep)*zth(:,jlep)                  &
     &           +zskap(jlev)*dsigma(jlev)*zth(:,jlev)
         zsumq1(:)=dsigma(jlep)+dsigma(jlev)
         zsumq2(:)=dsigma(jlep)*zqn(:,jlep)+dsigma(jlev)*zqn(:,jlev)
         ztht(:)=zsum2(:)/zsum1(:)
         zqt(:)=zsumq2(:)/zsumq1(:)
         icon(:)=1
        endwhere
       enddo
       if(SUM(icon(:)) > 0) then
        do jlev=NLEM-1,1,-1
         where(itop(:) == jlev .and. zth(:,jlev) < ztht(:))
          itop(:)=jlev-1
          zsum1(:)=zsum1(:)+zskap(jlev)*dsigma(jlev)
          zsum2(:)=zsum2(:)+zskap(jlev)*dsigma(jlev)*zth(:,jlev)
          zsumq1(:)=zsumq1(:)+dsigma(jlev)
          zsumq2(:)=zsumq2(:)+dsigma(jlev)*zqn(:,jlev)
          ztht(:)=zsum2(:)/zsum1(:)
          zqt(:)=zsumq2(:)/zsumq1(:)
         endwhere
        enddo
        do jlev=1,NLEV
         where(jlev > itop(:) .and. jlev <= ibase(:))
          zth(:,jlev)=ztht(:)
          zqn(:,jlev)=zqt(:)
         endwhere
        enddo
        jiter=jiter+1
        if(jiter > 2*NLEV) goto 1001
        goto 1000
      endif
!
      do jlev=1,NLEV
       zdtdt(:,jlev)=(zth(:,jlev)*zskap(jlev)-zt(:,jlev))/deltsec2
       zdqdt(:,jlev)=(zqn(:,jlev)-zq(:,jlev))/deltsec2
       dtdt(:,jlev)=dtdt(:,jlev)+zdtdt(:,jlev)
       dqdt(:,jlev)=dqdt(:,jlev)+zdqdt(:,jlev)
      enddo
!
!     franks dbug
!
      if(ndiaggp==1) then
       dgp3d(:,1:NLEV,9)=zdtdt(:,1:NLEV)
       do jlev=1,NLEV
        dgp3d(:,jlev,17)=zdqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      endif
!
!     entropy/energy diagnostics
!
      if(nentropy > 0) then
       dentropy(:,13)=0.
       do jlev=1,NLEV
        dentro(:)=zdtdt(:,jlev)/dentrot(:,jlev)                         &
     &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev) 
        dentropy(:,13)=dentropy(:,13)+dentro(:)
        if(nentro3d > 0) dentro3d(:,jlev,13)=dentro(:)
       enddo
      endif
      if(nenergy > 0) then
       denergy(:,13)=0.
       do jlev=1,NLEV
        denergy(:,13)=denergy(:,13)                                     &
     &               +zdtdt(:,jlev)                                     &
     &               *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,13)=zdtdt(:,jlev)                               &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif
       enddo
      endif
!
      return
 1001 continue
      end subroutine mkdca

      subroutine mkrain
      use rainmod
      parameter(zeps=0.)
!
      real zq(NHOR),zt(NHOR),zqsat(NHOR)
      real zlcpe(NHOR)       ! L/CP
      real zcor(NHOR)        ! correction term for saturation humidity
!
!*    convert rain <-> snow during its way down
!
      real :: zprc(NHOR)
      real :: zprl(NHOR)
      real :: zprsc(NHOR)
      real :: zprsl(NHOR)
      real :: zdqdt(NHOR)               ! q-change due to evap of rain
      real :: zdtdt(NHOR,NLEV)          ! t-change
      real :: zdqdts(NHOR,NLEV)         ! sublimation rate
      real :: zdqdtl(NHOR,NLEV)         ! evaporation rate
      real :: zlcp(NHOR)
!
!*    initialize
!
      zprc(:)=0.
      zprl(:)=0.
      zprsc(:)=0.
      zprsl(:)=0.
      zdtdt(:,:)=0.
      zdqdts(:,:)=0.
      zdqdtl(:,:)=0.
!
!*    convert rain <-> snow according to temperature
!
      do jlev=1,NLEV
       zprc(:)=zprc(:)+dprcl(:,jlev)-dprscl(:,jlev)
       zprl(:)=zprl(:)+dprll(:,jlev)-dprsll(:,jlev)
       zprsc(:)=zprsc(:)+dprscl(:,jlev)
       zprsl(:)=zprsl(:)+dprsll(:,jlev)
       where(zprsc(:)+zprsl(:) > 0. .and. dt(:,jlev) > TMELT)
        zlcp(:)=(ALS-ALV)/(acpd*(1.+ADV*dq(:,jlev)))
        zdtdt(:,jlev)=-(zprsc(:)+zprsl(:))/dsigma(jlev)*zlcp(:)
        zprc(:)=zprsc(:)+zprc(:)
        zprl(:)=zprsl(:)+zprl(:)
        zprsc(:)=0.
        zprsl(:)=0.
       elsewhere(zprc(:)+zprl(:) > 0. .and. dt(:,jlev) <= TMELT)
        zlcp(:)=(ALS-ALV)/(acpd*(1.+ADV*dq(:,jlev)))
        zdtdt(:,jlev)=(zprc(:)+zprl(:))/dsigma(jlev)*zlcp(:)
        zprsc(:)=zprc(:)+zprsc(:)
        zprsl(:)=zprl(:)+zprsl(:)
        zprc(:)=0.
        zprl(:)=0.
       endwhere
!
!*     convert also according to surface temp.
!      (since melting of falling snow and, in particular,
!       freezing of rain is mayy not be incl. in surface scheme)
!
       if(jlev==NLEV) then
        where(zprsc(:)+zprsl(:) > 0. .and. dt(:,NLEP) > TMELT)
         zlcp(:)=(ALS-ALV)/(acpd*(1.+ADV*dq(:,jlev)))
         zdtdt(:,jlev)=-(zprsc(:)+zprsl(:))/dsigma(jlev)*zlcp(:)        &
     &                +zdtdt(:,jlev)
         zprc(:)=zprsc(:)+zprc(:)
         zprl(:)=zprsl(:)+zprl(:)
         zprsc(:)=0.
         zprsl(:)=0.
        elsewhere(zprc(:)+zprl(:) > 0. .and. dt(:,NLEP) <= TMELT)
         zlcp(:)=(ALS-ALV)/(acpd*(1.+ADV*dq(:,jlev)))
         zdtdt(:,jlev)=(zprc(:)+zprl(:))/dsigma(jlev)*zlcp(:)           &
     &                +zdtdt(:,jlev) 
         zprsc(:)=zprc(:)+zprsc(:)
         zprsl(:)=zprl(:)+zprsl(:)
         zprc(:)=0.
         zprl(:)=0.
        endwhere
       endif
!
!*     evaporation/sublimation of rain/snow
!
       if(nevapprec.eq.1) then
         zq(:)=dq(:,jlev)+dqdt(:,jlev)*deltsec2
         zt(:)=dt(:,jlev)+dtdt(:,jlev)*deltsec2
!
!        saturation humidity
! 
         zqsat(:)=rdbrv*ra1*exp(ra2*(zt(:)-TMELT)/(zt(:)-ra4))          &
     &         /(dp(:)*sigma(jlev))
!
!        avoid negative qsat
!
         zqsat(:)=AMIN1(zqsat(:),rdbrv)
!
         zcor(:)=1./(1.-(1./rdbrv-1.)*zqsat(:))
         zqsat(:)=zqsat(:)*zcor(:)
!        Fix for problem with higher Mars atmosphere
         where (zqsat(:) < 0.0)
           zqsat(:) = 0.0
         endwhere
!
         where((zq(:).lt.zqsat(:)).AND.zprl(:) > zeps)
           zlcpe(:)=ALV/(acpd*(1.+ADV*dq(:,jlev)))
           zdqdt(:)=AMIN1(gamma*(zqsat(:)-zq(:))*dsigma(jlev)/deltsec2  &
     &                   /(1.0+zlcpe(:)*ra2*(TMELT-ra4)                 &
     &                   *zqsat(:)*zcor(:)/(zt(:)-ra4)**2)              &
     &                   ,zprl(:))
           zprl(:)=zprl(:)-zdqdt(:)
           zdqdt(:)=zdqdt(:)/dsigma(jlev)
           zdtdt(:,jlev)=zdtdt(:,jlev)-zdqdt(:)*zlcpe(:)
           zdqdtl(:,jlev)=zdqdtl(:,jlev)+zdqdt(:)
         endwhere
         where((zq(:).lt.zqsat(:)).AND.zprc(:) > zeps)
           zlcpe(:)=ALV/(acpd*(1.+ADV*dq(:,jlev)))
           zdqdt(:)=AMIN1(gamma*(zqsat(:)-zq(:))*dsigma(jlev)/deltsec2  &
     &                   /(1.0+zlcpe(:)*ra2*(TMELT-ra4)                 &
     &                   *zqsat(:)*zcor(:)/(zt(:)-ra4)**2)              &
     &                   ,zprc(:))
           zprc(:)=zprc(:)-zdqdt(:)
           zdqdt(:)=zdqdt(:)/dsigma(jlev)
           zdtdt(:,jlev)=zdtdt(:,jlev)-zdqdt(:)*zlcpe(:)
           zdqdtl(:,jlev)=zdqdtl(:,jlev)+zdqdt(:)
         endwhere
         where((zq(:).lt.zqsat(:)).AND.zprsl(:) > zeps)
           zlcpe(:)=ALS/(acpd*(1.+ADV*dq(:,jlev)))
           zdqdt(:)=AMIN1(gamma*(zqsat(:)-zq(:))*dsigma(jlev)/deltsec2  &
     &                   /(1.0+zlcpe(:)*ra2*(TMELT-ra4)                 &
     &                   *zqsat(:)*zcor(:)/(zt(:)-ra4)**2)              &
     &                   ,zprsl(:))
           zprsl(:)=zprsl(:)-zdqdt(:)
           zdqdt(:)=zdqdt(:)/dsigma(jlev)
           zdtdt(:,jlev)=zdtdt(:,jlev)-zdqdt(:)*zlcpe(:)
           zdqdts(:,jlev)=zdqdts(:,jlev)+zdqdt(:)
         endwhere
         where((zq(:).lt.zqsat(:)).AND.zprsc(:) > zeps)
           zlcpe(:)=ALS/(acpd*(1.+ADV*dq(:,jlev)))
           zdqdt(:)=AMIN1(gamma*(zqsat(:)-zq(:))*dsigma(jlev)/deltsec2  &
     &                   /(1.0+zlcpe(:)*ra2*(TMELT-ra4)                 &
     &                   *zqsat(:)*zcor(:)/(zt(:)-ra4)**2)              &
     &                   ,zprsc(:))
           zprsc(:)=zprsc(:)-zdqdt(:)
           zdqdt(:)=zdqdt(:)/dsigma(jlev)
           zdtdt(:,jlev)=zdtdt(:,jlev)-zdqdt(:)*zlcpe(:)
           zdqdts(:,jlev)=zdqdts(:,jlev)+zdqdt(:)
         endwhere
       endif
!
      enddo
!
!*    compute final rain and snow (im m/s)
!
      dprc(:)=(zprc(:)+zprsc(:))*dp(:)/ga/1000.
      dprl(:)=(zprl(:)+zprsl(:))*dp(:)/ga/1000.
      dprs(:)=(zprsc(:)+zprsl(:))*dp(:)/ga/1000.
!
!*    add tendencies
!
      dqdt(:,1:NLEV)=dqdt(:,1:NLEV)+zdqdtl(:,1:NLEV)+zdqdts(:,1:NLEV)
      dtdt(:,1:NLEV)=dtdt(:,1:NLEV)+zdtdt(:,1:NLEV)
!
!*    franks diagnostics
!
      if(ndiaggp==1) then
       dgp3d(:,1:NLEV,10)=zdtdt(:,1:NLEV)
       dgp3d(:,1:NLEV,19)=zdqdtl(:,1:NLEV)+zdqdts(:,1:NLEV)
      endif
!
!*    entropy/energy diagnostics
!
      if(nentropy > 0) then
       dentropy(:,20)=0.
       do jlev=1,NLEV
        dentro(:)=zdtdt(:,jlev)/dentrot(:,jlev)                         &
     &         *acpd*(1.+adv*dentroq(:,jlev))*dentrop(:)/ga*dsigma(jlev)
        dentropy(:,20)=dentropy(:,20)+dentro(:)
        if(nentro3d > 0) dentro3d(:,jlev,14)=dentro(:)
       enddo
      endif
      if(nenergy > 0) then
       denergy(:,14)=0.
       do jlev=1,NLEV
         denergy(:,14)=denergy(:,14)                                    &
     &                +zdtdt(:,jlev)                                    &
     &                *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
         if(nener3d > 0) then
          dener3d(:,jlev,14)=zdtdt(:,jlev)                              &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
         endif 
       enddo
      endif

      return
      end subroutine mkrain
