!     =================
!     SUBROUTINE OUTINI
!     =================

      subroutine outini
      use pumamod

!     Output has always 32 bit precision

      integer (kind=4) :: ihead(8)   ! header of first data set
      integer (kind=4) :: ihead2(8)  ! header for source variable
      real    (kind=4) :: zsig(NUGP) ! first block contains settings
      real    (kind=4) :: sources(NLIGHTS)
      
      call ntomin(nstep,nmin,nhour,nday,nmonth,nyear)

      ihead(1) = 333  ! ID for PUMA/PLASIM parameter block
      ihead(2) = 0
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
      ihead(4) = 0
      ihead(5) = NLON
      ihead(6) = NLAT
      ihead(7) = NLEV
      ihead(8) = NTRU

!     The first data block with Code = 333 is used for transferring
!     planet properties to the postprocessor "burn"

      zsig(:)      = 0.0       ! initialize
      zsig(1:NLEV) = sigmah(:) ! vertical coordinate table
      zsig(NLEV+1) = m_days_per_year
      
      ihead2(1) = 49 !ID for light source block
      ihead2(2) = 0
      ihead2(3) = nday + 100 * nmonth + 10000 * nyear
      ihead2(4) = 0
      ihead2(5) = NLIGHTS
      ihead2(6) = 1
      ihead2(7) = 1
      ihead2(8) = 1
      
      do k=1,NLIGHTS
        sources(k) = k
      enddo

      open  (40,file=plasim_output,form='unformatted')
      write (40) ihead(:)
      write (40) zsig(:)
      write (40) ihead2(:)
      write (40) sources(:)
      
      if (nsnapshot > 0) then
        open (140,file=plasim_snapshot,form='unformatted')
        write(140) ihead(:)
        write(140) zsig(:)
        write(140) ihead2(:)
        write(140) sources(:)
      endif
      
      if (nhcadence > 0) then
        open (141,file=plasim_hcadence,form='unformatted')
        write(141) ihead(:)
        write(141) zsig(:)
        write(141) ihead2(:)
        write(141) sources(:)
      endif
    
      return
      end
  
! !
!     ============================
!     SUBROUTINE HCOUTPUTINI
!     ============================
!
!     Set up output file
!    
      
      subroutine hcoutputini
      use pumamod
      use hurricanemod
 
!     Output has always 32 bit precision

      integer (kind=4) :: ihead(8)   ! header of first data set
      integer (kind=4) :: ihead2(8)
      real    (kind=4) :: zsig(NUGP) ! first block contains settings
      real    (kind=4) :: sources(NLIGHTS)
      
      call ntomin(nstep,nmin,nhour,nday,nmonth,nyear)

      ihead(1) = 333  ! ID for PUMA/PLASIM parameter block
      ihead(2) = 0
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
      ihead(4) = 0
      ihead(5) = NLON
      ihead(6) = NLAT
      ihead(7) = NLEV
      ihead(8) = NTRU

!     The first data block with Code = 333 is used for transferring
!     planet properties to the postprocessor "burn"

      zsig(:)      = 0.0       ! initialize
      zsig(1:NLEV) = sigmah(:) ! vertical coordinate table
      zsig(NLEV+1) = m_days_per_year     
      
      ihead2(1) = 49 !ID for light source block
      ihead2(2) = 0
      ihead2(3) = nday + 100 * nmonth + 10000 * nyear
      ihead2(4) = 0
      ihead2(5) = NLIGHTS
      ihead2(6) = 1
      ihead2(7) = 1
      ihead2(8) = 1
      
      do k=1,NLIGHTS
        sources(k) = k
      enddo
      
      write(hc_output,"(A17,I0.3)") "hurricane_output.",kstorms
      
      open(142,file=hc_output,form='unformatted')
      write(142) ihead(:)
      write(142) zsig(:)
      write(142) ihead2(:)
      write(142) sources(:)
      
      return
      end subroutine hcoutputini
!     ==================
!     SUBROUTINE WRITEGP
!     ==================

      subroutine writegp(kunit,pf,kcode,klev)
      use pumamod
      real :: pf(NHOR)
      real :: zf(NUGP)
      integer(kind=4) :: ihead(8)
      integer(kind=4) :: la(NPGP)
      real(kind=4) :: zmin
      real(kind=4) :: zsca
      real(kind=4) :: zzf(NUGP)

! We'll need to change how we write this header. It wants the timestep to translate to
! multiple months per year, and multiple days per year, as well as minutes and hours per
! day, which resets every day. There should be multiple days per year. It might make the
! most sense to rewrite ntomin.
      
      
      
      istep = nstep
      call ntomin(istep,nmin,nhour,nday,nmonth,nyear)

      ihead(1) = kcode
      ihead(2) = klev
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
      ihead(4) = nmin + 100 * nhour
      ihead(5) = NLON
      ihead(6) = NLAT
      ihead(7) = nstep - nstep1
      ihead(8) = m_days_per_year

      call mpgagp(zf,pf,1)

      if (mypid == NROOT) then
         write (kunit) ihead
         zzf(:) = zf(:)
         write (kunit) zzf
      endif

      return
      end

!     ==================
!     SUBROUTINE WRITESP
!     ==================

      subroutine writesp(kunit,pf,kcode,klev,pscale,poff)
      use pumamod
      real :: pf(NRSP)
      integer(kind=4) :: ihead(8)
      integer(kind=4) :: la(NTP1+1:NCSP)
      real(kind=4) :: zf(NRSP)
      real(kind=4) :: za(NTP1+1)

      istep = nstep
      call ntomin(istep,nmin,nhour,nday,nmonth,nyear)

      ihead(1) = kcode
      ihead(2) = klev
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
      ihead(4) = nmin + 100 * nhour
      ihead(5) = NRSP
      ihead(6) = 1
      ihead(7) = nstep - nstep1
      ihead(8) = m_days_per_year

!     normalize ECHAM compatible and scale to physical dimensions

      zf(:) = pf(:) * spnorm(1:NRSP) * pscale
      zf(1) = zf(1) + poff ! Add offset if necessary
      write (kunit) ihead
      write (kunit) zf

      return
      end

!     ======================
!     SUBROUTINE WRITESCALAR
!     ======================

!     This will output a single scalar per timestep

      subroutine writescalar(kunit,pf,kcode)
      use pumamod
      real :: pf
      integer(kind=4) :: ihead(8)

      istep = nstep
      call ntomin(istep,nmin,nhour,nday,nmonth,nyear)

      ihead(1) = kcode
      ihead(2) = 0
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
      ihead(4) = nmin + 100 * nhour
      ihead(5) = 1
      ihead(6) = 1
      ihead(7) = nstep - nstep1
      ihead(8) = m_days_per_year

      write (kunit) ihead
      write (kunit) pf

      return
      end
      
!     ======================
!     SUBROUTINE WRITECOLUMN
!     ======================

!     This will output a single column per timestep. Currently unused, but it could be.

      subroutine writecolumn(kunit,pf,kcode,klev)
      use pumamod
      real :: pf
      integer(kind=4) :: ihead(8)

      istep = nstep
      call ntomin(istep,nmin,nhour,nday,nmonth,nyear)

      ihead(1) = kcode
      ihead(2) = klev
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
      ihead(4) = nmin + 100 * nhour
      ihead(5) = 1
      ihead(6) = 1
      ihead(7) = nstep - nstep1
      ihead(8) = m_days_per_year

      write (kunit) ihead
      write (kunit) pf

      return
      end
      
!     ================
!     SUBROUTINE OUTSC
!     ================
!     !Output scalars of interest

      subroutine outsc
      use pumamod
      use radmod
      
      if (nlowio .eq. 0) then
        
        !True anomaly
        call writescalar(40,orbnu*180./PI,50)
        
        !Solar ecliptic longitude
        call writescalar(40,lambm*180./PI,51)
        
        !Solar declination
        call writescalar(40,zdeclf*180./PI,52)
        
        !Solar distance modulus
        call writescalar(40,1.0/sqrt(eccf),53)
        
        !Solar right ascension
        call writescalar(40,rasc*180./PI,54)
        
      else !low IO mode
        !True anomaly
        aorbnu = aorbnu/real(max(1,naccuout))
        call writescalar(40,aorbnu*180./PI,50)
        
        !Solar ecliptic longitude
        alambm = alambm/real(max(1,naccuout))
        call writescalar(40,alambm*180./PI,51)
        
        !Solar declination
        azdecl = azdecl/real(max(1,naccuout))
        call writescalar(40,azdecl*180./PI,52)
        
        !Solar distance modulus
        ardist = ardist/real(max(1,naccuout))
        call writescalar(40,ardist,53)
        
        !Solar right ascension
        arasc = arasc/real(max(1,naccuout))
        call writescalar(40,rasc*180./PI,54)
      endif
      
      return
      end
      
      
!     ================
!     SUBROUTINE OUTSP
!     ================

      subroutine outsp
      use pumamod

      
      if (nlowio .eq. 0) then
!       ************
!       * orograpy *
!       ************
       
        call writesp(40,so,129,0,CV*CV,0.)
       
!       ************
!       * pressure *
!       ************
       
        call writesp(40,sp,152,0,1.0,log(psurf))
       
!       ***************
!       * temperature *
!       ***************
       
        do jlev = 1 , NLEV
           call writesp(40,st(1,jlev),130,jlev,ct,t0(jlev) * ct)
        enddo
       
!       *********************
!       * specific humidity *
!       *********************
       
        if (nqspec == 1) then
           do jlev = 1 , NLEV
              call writesp(40,sqout(1,jlev),133,jlev,1.0,0.0)
           enddo
        endif
       
!       **************
!       * divergence *
!       **************
       
        do jlev = 1 , NLEV
           call writesp(40,sd(1,jlev),155,jlev,ww,0.0)
        enddo
       
!       *************
!       * vorticity *
!       *************
       
        do jlev = 1 , NLEV
           zsave = sz(3,jlev)
           sz(3,jlev) = sz(3,jlev) - plavor
           call writesp(40,sz(1,jlev),138,jlev,ww,0.0)
           sz(3,jlev) = zsave
        enddo

      else !Low-I/O mode
!       ************
!       * orograpy *
!       ************
        
        
!         write(nud,'("* nstep ",i6,"  *")') nstep
!         write(nud,*) real(max(1,naccuout))
!         write(nud,*) aaso(:)
        aaso(:) = aaso(:) / real(max(1,naccuout))
        call writesp(40,aaso,129,0,CV*CV,0.)
        
!       ************
!       * pressure *
!       ************
        
        aasp(:) = aasp(:) / real(max(1,naccuout))
        call writesp(40,aasp,152,0,1.0,log(psurf))
        
!       ***************
!       * temperature *
!       ***************
        
        do jlev = 1 , NLEV
           aast(:,jlev) = aast(:,jlev) / real(max(1,naccuout))
           call writesp(40,aast(1,jlev),130,jlev,ct,t0(jlev) * ct)
        enddo
        
!       *********************
!       * specific humidity *
!       *********************
        
        if (nqspec == 1) then
           do jlev = 1 , NLEV
              aasqout(:,jlev) = aasqout(:,jlev) / real(max(1,naccuout))
              call writesp(40,aasqout(1,jlev),133,jlev,1.0,0.0)
           enddo
        endif
        
!       **************
!       * divergence *
!       **************
        
        do jlev = 1 , NLEV
           aasd(:,jlev) = aasd(:,jlev) / real(max(1,naccuout))
           call writesp(40,aasd(1,jlev),155,jlev,ww,0.0)
        enddo
        
!       *************
!       * vorticity *
!       *************
        
        do jlev = 1 , NLEV
           aasz(:,jlev) = aasz(:,jlev) / real(max(1,naccuout))
           zsave = aasz(3,jlev)
           aasz(3,jlev) = aasz(3,jlev) - plavor
           call writesp(40,aasz(1,jlev),138,jlev,ww,0.0)
           aasz(3,jlev) = zsave
        enddo      
        
      
      endif
      
      return
      end


!     ================
!     SUBROUTINE OUTGP
!     ================

      subroutine outgp
      use pumamod
      use carbonmod
      use radmod
      use glaciermod

      if (nlowio .eq. 0) then
!       *********************
!       * specific humidity *
!       *********************
        
        if (nqspec == 0) then ! Semi Langrangian advection active
           do jlev = 1 , NLEV
              call writegp(40,dq(1,jlev),133,jlev)
           enddo
        endif
        
!       **********************************
!       * mixed-layer depth (from ocean) *
!       **********************************
        
        call writegp(40,dmld,110,0)
        
!       ***********************
!       * surface temperature *
!       ***********************
        
        call writegp(40,dt(1,NLEP),139,0)
        
!       ****************
!       * soil wetness *
!       ****************
        
        call writegp(40,dwatc,140,0)
        
!       **************
!       * snow depth *
!       **************
        
        call writegp(40,dsnow,141,0)
      
      else ! Low I/O mode
      
!       *********************
!       * specific humidity *
!       *********************
        
        if (nqspec == 0) then ! Semi Langrangian advection active
           do jlev = 1 , NLEV
              aadq(:,jlev) = aadq(:,jlev)/real(naccuout)
              call writegp(40,aadq(1,jlev),133,jlev)
           enddo
        endif
        
!       **********************************
!       * mixed-layer depth (from ocean) *
!       **********************************
        
        aadmld(:) = aadmld(:)/real(naccuout)
        call writegp(40,aadmld,110,0)
        
!       ***********************
!       * surface temperature *
!       ***********************
        
!         thyng = SUM(aadt(:,NLEP))/NHOR/real(naccuout)
!         if (thyng .lt. 15.0) then
!            open(unit=79,file='surftemp_error.log',position='append',status='unknown')
!            write(79,*) nstep, thyng
!            write(79,*) aadt(1,NLEP),naccuout,dt(1,NLEP),mypid
!            close(79)
!         endif
        if (aadt(1,NLEP) .ne. dt(1,NLEP)) then
          do jlev=1,NLEP
            aadt(:,jlev) = aadt(:,jlev)/real(naccuout)
          enddo
        endif
!         thyng = SUM(aadt(:,NLEP))/NHOR
!         if (thyng .lt. 15.0) then
!            open(unit=79,file='surftemp_error.log',position='append',status='unknown')
!            write(79,*) nstep, thyng
!            write(79,*) aadt(1,NLEP),naccuout,dt(1,NLEP),mypid
!            close(79)
!         endif
        call writegp(40,aadt(1,NLEP),139,0)
        
!       ****************
!       * soil wetness *
!       ****************
        
        aadwatc(:) = aadwatc(:)/real(naccuout)
        call writegp(40,aadwatc,140,0)
        
!       **************
!       * snow depth *
!       **************
        
        aadsnow(:) = aadsnow(:)/real(naccuout)
        call writegp(40,aadsnow,141,0)
      
      endif
      
!     **********************
!     * large scale precip *
!     **********************

      aprl(:)=aprl(:)/real(naccuout)
      call writegp(40,aprl,142,0)

!     *********************
!     * convective precip *
!     *********************

      aprc(:)=aprc(:)/real(naccuout)
      call writegp(40,aprc,143,0)

!     *************
!     * snow fall *
!     *************

      aprs(:)=aprs(:)/real(naccuout)
      call writegp(40,aprs,144,0)

!     **********************
!     * sensible heat flux *
!     **********************

      ashfl(:)=ashfl(:)/real(naccuout)
      call writegp(40,ashfl,146,0)

!     ********************
!     * latent heat flux *
!     ********************

      alhfl(:)=alhfl(:)/real(naccuout)
      call writegp(40,alhfl,147,0)

      if (nlowio .eq. 0) then
      
!       ************************
!       * liquid water content *
!       ************************
     
        do jlev = 1 , NLEV
           call writegp(40,dql(1,jlev),161,jlev)
        enddo
     
!       *************
!       * u-star**3 *
!       *************
     
        call writegp(40,dust3,159,0)
     
      else !Low I/O mode
      
!       ************************
!       * liquid water content *
!       ************************
     
        do jlev = 1 , NLEV
           aadql(:,jlev) = aadql(:,jlev)/real(naccuout)
           call writegp(40,aadql(1,jlev),161,jlev)
        enddo
     
!       *************
!       * u-star**3 *
!       *************
     
        aadust3(:) = aadust3(:)/real(naccuout)
        call writegp(40,aadust3,159,0)
     
     endif
      
     
!     **********
!     * runoff *
!     **********

      aroff(:)=aroff(:)/real(naccuout)
      call writegp(40,aroff,160,0)

!     ***************
!     * cloud cover *
!     ***************

      if (nlowio .eq. 0) then
        do jlev = 1 , NLEV
          call writegp(40,dcc(1,jlev),162,jlev)
        enddo
      else
        do jlev = 1 , NLEV
          aadcc(:,jlev) = aadcc(:,jlev)/real(naccuout)
          call writegp(40,aadcc(1,jlev),162,jlev)
        enddo
      endif
      acc(:)=acc(:)/real(naccuout)
      call writegp(40,acc,164,0)

!     ***************************
!     * surface air temperature *
!     ***************************

      atsa(:)=atsa(:)/real(naccuout)
      call writegp(40,atsa,167,0)

!     ******************************
!     * surface temperature (accu) *
!     ******************************

      ats0(:)=ats0(:)/real(naccuout)
      call writegp(40,ats0,169,0)

      if (nlowio .eq. 0) then
      
!       *************************
!       * deep soil temperature *
!       *************************
      
        call writegp(40,dtd5,170,0)
      
!       *****************
!       * land sea mask *
!       *****************
      
        call writegp(40,dls,172,0)
      
!       *********************
!       * surface roughness *
!       *********************
      
        call writegp(40,dz0,173,0)
      
!       **********
!       * albedo *
!       **********
        
        do klight=1,NLIGHTS
           call writegp(40,dalb(:,klight),175,klight)
        enddo
        do klight=1,NLIGHTS
           call writegp(40,dsalb(:,2*klight-1),174,klight)
        enddo
        do klight=1,NLIGHTS
           call writegp(40,dsalb(:,2*klight),184,klight)
        enddo
        
      else
      
!       *************************
!       * deep soil temperature *
!       *************************
        
        aadtd5(:) = aadtd5(:)/real(naccuout)
        call writegp(40,aadtd5,170,0)
      
!       *****************
!       * land sea mask *
!       *****************
      
        aadls(:) = aadls(:)/real(naccuout)
        call writegp(40,aadls,172,0)
      
!       *********************
!       * surface roughness *
!       *********************
      
        aadz0(:) = aadz0(:)/real(naccuout)
        call writegp(40,aadz0,173,0)
      
!       **********
!       * albedo *
!       **********
      
        do klight=1,NLIGHTS
           aadalb(:,klight) = aadalb(:,klight)/real(naccuout)
           call writegp(40,aadalb(:,klight),175,klight)
        enddo
        do klight=1,NLIGHTS
           aadsalb1(:,klight) = aadsalb1(:,klight)/real(naccuout)
           call writegp(40,aadsalb1(:,klight),174,klight)
        enddo
        do klight=1,NLIGHTS
           aadsalb2(:,klight) = aadsalb2(:,klight)/real(naccuout)
           call writegp(40,aadsalb2(:,klight),184,klight)
        enddo
        
      endif

!     ***************************
!     * surface solar radiation *
!     ***************************

      assol(:)=assol(:)/real(naccuout)
      call writegp(40,assol,176,0)

!     *****************************
!     * surface thermal radiation *
!     *****************************

      asthr(:)=asthr(:)/real(naccuout)
      call writegp(40,asthr,177,0)

!     ***********************
!     * top solar radiation *
!     ***********************

      atsol(:)=atsol(:)/real(naccuout)
      call writegp(40,atsol,178,0)

!     *************************
!     * top thermal radiation *
!     *************************

      atthr(:)=atthr(:)/real(naccuout)
      call writegp(40,atthr,179,0)

!     ************
!     * u-stress *
!     ************

      ataux(:)=ataux(:)/real(naccuout)
      call writegp(40,ataux,180,0)

!     *************
!     * v- stress *
!     *************

      atauy(:)=atauy(:)/real(naccuout)
      call writegp(40,atauy,181,0)

!     ***************
!     * evaporation *
!     ***************

      aevap(:)=aevap(:)/real(naccuout)
      call writegp(40,aevap,182,0)

!     *********************
!     * soil temperature *
!     *********************
      
      if (nlowio .eq. 0) then
        call writegp(40,dtsoil,183,0)
      else
        aadtsoil(:) = aadtsoil(:)/real(naccuout)
        call writegp(40,aadtsoil,183,0)
      endif
      
      
!     ***********************************
!     * maximum surface air temperature *
!     ***********************************

      call writegp(40,atsama,201,0)

!     ***********************************
!     * minimum surface air temperature *
!     ***********************************

      call writegp(40,atsami,202,0)

!     ********************
!     * top solar upward *
!     ********************

      atsolu(:)=atsolu(:)/real(naccuout)
      call writegp(40,atsolu,203,0)

!     ************************
!     * surface solar upward *
!     ************************

      assolu(:)=assolu(:)/real(naccuout)
      call writegp(40,assolu,204,0)

!     **************************
!     * surface thermal upward *
!     **************************

      asthru(:)=asthru(:)/real(naccuout)
      call writegp(40,asthru,205,0)

      if (nlowio .eq. 0) then
      
!       *******************************
!       * soil temperatures level 2-4 *
!       *******************************
       
        call writegp(40,dtd2,207,0)
        call writegp(40,dtd3,208,0)
        call writegp(40,dtd4,209,0)
       
!       *****************
!       * sea ice cover *
!       *****************
       
        call writegp(40,dicec,210,0)
       
!       *********************
!       * sea ice thickness *
!       *********************
       
        call writegp(40,diced,211,0)
       
!       ****************
!       * forest cover *
!       ****************
       
        call writegp(40,dforest,298,0)
       
      else
      
!       *******************************
!       * soil temperatures level 2-4 *
!       *******************************
        
        aadtd2(:) = aadtd2(:) / real(naccuout)
        aadtd3(:) = aadtd3(:) / real(naccuout)
        aadtd4(:) = aadtd4(:) / real(naccuout)
        
        call writegp(40,aadtd2,207,0)
        call writegp(40,aadtd3,208,0)
        call writegp(40,aadtd4,209,0)
       
!       *****************
!       * sea ice cover *
!       *****************
       
        aadicec(:) = aadicec(:) / real(naccuout)
        call writegp(40,aadicec,210,0)
       
!       *********************
!       * sea ice thickness *
!       *********************
       
        aadiced(:) = aadiced(:) / real(naccuout)
        call writegp(40,aadiced,211,0)
       
!       ****************
!       * forest cover *
!       ****************
       
        aadforest(:) = aadforest(:) / real(naccuout)
        call writegp(40,aadforest,298,0)
      
      endif
       
!     *************
!     * snow melt *
!     *************

      asmelt(:)=asmelt(:)/real(naccuout)
      call writegp(40,asmelt,218,0)

!     *********************
!     * snow depth change *
!     *********************

!       asndch(:)=asndch(:)/real(naccuout)
      call writegp(40,asndch,221,0)

!     ******************
!     * field capacity *
!     ******************
      if (nlowio .eq. 0) then
        call writegp(40,dwmax,229,0)
      else
        aadwmax(:) = aadwmax(:)/real(naccuout)
        call writegp(40,aadwmax,229,0)
      endif
        
        
!     *****************************************
!     * vertical integrated specific humidity *
!     *****************************************

      aqvi(:)=aqvi(:)/real(naccuout)
      call writegp(40,aqvi,230,0)

!     ****************
!     * glacier mask *
!     ****************

      if (nlowio .eq. 0) then
        call writegp(40,dglac,232,0)
      else
        aadglac(:) = aadglac(:) / real(naccuout)
        call writegp(40,aadglac,232,0)
      endif
        
!     *********************
!     ***   S I M B A   ***
!     *********************

      
      if (nveg > 0) then
        if (nlowio .eq. 0) then
          call vegsnapout(40)
        else
          call vegout
        endif
      endif

!     ***********************
!     * Ozone concentration *
!     ***********************

      if (nlowio .eq. 0) then
        do jlev = 1 , NLEV
           call writegp(40,dqo3(1,jlev),265,jlev)
        enddo
      else
        do jlev = 1 , NLEV
           aadqo3(:,jlev) = aadqo3(:,jlev)/real(naccuout)
           call writegp(40,aadqo3(1,jlev),265,jlev)
        enddo
      endif
      
!     ********************
!     * local weathering *
!     ********************

      aweathering(:) = aweathering(:)/real(naccuout)
      call writegp(40,aweathering,266,0)

      
      if (nlowio .eq. 0) then
!       ********************
!       * ground elevation *
!       ********************
       
        call writegp(40,groundoro,267,0)
       
!       *********************
!       * glacier elevation *
!       *********************
       
        call writegp(40,glacieroro,268,0)
       
        
!       *********************
!       *    net elevation  *
!       *********************
        netoro(:) = groundoro(:) + glacieroro(:)
        call writegp(40,netoro,269,0)
        
      else
!       ********************
!       * ground elevation *
!       ********************
        
        aagroundoro(:) = aagroundoro(:) / real(naccuout)
        call writegp(40,aagroundoro,267,0)
       
!       *********************
!       * glacier elevation *
!       *********************
       
        aaglacieroro(:) = aaglacieroro(:) / real(naccuout)
        call writegp(40,aaglacieroro,268,0)
       
        
!       *********************
!       *    net elevation  *
!       *********************
        netoro(:) = aagroundoro(:) + aaglacieroro(:)
        call writegp(40,netoro,269,0)
      
      endif
       
!     ********************
!     * Cos Solar Zenith *
!     ********************
         
      azmuz(:) = azmuz(:)/real(naccuout)   
      call writegp(40,azmuz,318,0)
      
!     *****************************
!     * Weatherable Precipitation *
!     *****************************
        
      asigrain(:) = asigrain(:)/real(naccuout)
      call writegp(40,asigrain,319,0)

!     ***********************
!     * Minimum Temperature *
!     ***********************
         
      call writegp(40,tempmin,320,0)

!     ***********************
!     * Maximum Temperature *
!     ***********************
         
      call writegp(40,tempmax,321,0)
      
      
      ! Hurricane fields
      if (nlowio .eq. 0) then
         !Convective Available Potential Energy at the surface
         call writegp(40,capen,322,0)
         !Level of neutral buoyancy
         call writegp(40,lnb,323,0)
         ! Tropospheric entropy deficit
         call writegp(40,chim,324,0)
         ! Tropospheric entropy deficit
         call writegp(40,laav,325,0)
         ! Maximum potential intensity 
         call writegp(40,mpoti,326,0)
         ! Ventilation index
         call writegp(40,venti,327,0)
         ! Ventilation-reduced Maximum Potential Intensity
         call writegp(40,vrmpi,328,0)
         ! Genesis Potential Index
         call writegp(40,gpi,329,0)
      else
         !Convective Available Potential Energy at the surface
         acapen(:) = acapen(:) / real(naccuout)
         call writegp(40,acapen,322,0)
         !Level of neutral buoyancy
         alnb(:) = alnb(:) / real(naccuout)
         call writegp(40,alnb,323,0)
         ! Tropospheric entropy deficit
         achim(:) = achim(:) / real(naccuout)
         call writegp(40,achim,324,0)
         ! Tropospheric entropy deficit
         alaav(:) = alaav(:) / real(naccuout)
         call writegp(40,alaav,325,0)
         ! Maximum potential intensity 
         ampoti(:) = ampoti(:) / real(naccuout)
         call writegp(40,ampoti,326,0)
         ! Ventilation index
         aventi(:) = aventi(:) / real(naccuout)
         call writegp(40,aventi,327,0)
         ! Ventilation-reduced Maximum Potential Intensity
         avrmpi(:) = avrmpi(:) / real(naccuout)
         call writegp(40,avrmpi,328,0)
         ! Genesis Potential Index
         agpi(:) = agpi(:) / real(naccuout)
         call writegp(40,agpi,329,0)
       endif
            
      return
      end

!     ==================
!     SUBROUTINE OUTDIAG
!     ==================

      subroutine outdiag
      use pumamod

!     *****************************************
!     * 2-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiagsp2d > 0 .and. mypid == NROOT) then
       do jdiag=1,ndiagsp2d
        jcode=50+jdiag
        call writesp(40,dsp2d(1,jdiag),jcode,0,1.,0.0)
       enddo
      end if

!     *****************************************
!     * 3-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiagsp3d > 0 .and. mypid == NROOT) then
       do jdiag=1,ndiagsp3d
        jcode=60+jdiag
        do jlev=1,NLEV
         call writesp(40,dsp3d(1,jlev,jdiag),jcode,jlev,1.,0.0)
        enddo
       enddo
      end if

!     *****************************************
!     * 2-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiaggp2d > 0) then
       do jdiag=1,ndiaggp2d
        jcode=jdiag
        call writegp(40,dgp2d(1,jdiag),jcode,0)
       enddo
      end if

!     *****************************************
!     * 3-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiaggp3d > 0) then
       do jdiag=1,ndiaggp3d
        jcode=20+jdiag
        do jlev=1,NLEV
         call writegp(40,dgp3d(1,jlev,jdiag),jcode,jlev)
        enddo
       enddo
      end if

!     ************************************************
!     * cloud forcing (clear sky fluxes) diagnostics *
!     ************************************************

      if(ndiagcf > 0) then
       call writegp(40,dclforc(1,1),101,0)
       call writegp(40,dclforc(1,2),102,0)
       call writegp(40,dclforc(1,3),103,0)
       call writegp(40,dclforc(1,4),104,0)
       call writegp(40,dclforc(1,5),105,0)
       call writegp(40,dclforc(1,6),106,0)
       call writegp(40,dclforc(1,7),107,0)
      end if

!     **************************************
!     * entropy diagnostics if switched on *
!     **************************************

      if(nentropy > 0) then
       do jdiag=1,36
        jcode=319+jdiag
        if(jcode == 333) cycle                      !333 is reserved
        call writegp(40,dentropy(1,jdiag),jcode,0)
       enddo
      end if
      if(nentro3d > 0) then
       do jdiag=1,23
        jcode=419+jdiag
        do jlev=1,NLEV
         call writegp(40,dentro3d(1,jlev,jdiag),jcode,jlev)
        enddo
       enddo
      end if

!     *************************************
!     * energy diagnostics if switched on *
!     *************************************

      if(nenergy > 0) then
       do jdiag=1,28
        jcode=359+jdiag
        call writegp(40,denergy(1,jdiag),jcode,0)
       enddo
      end if
      if(nener3d > 0) then
       do jdiag=1,28
        jcode=459+jdiag
        do jlev=1,NLEV
         call writegp(40,dener3d(1,jlev,jdiag),jcode,jlev)
        enddo
       enddo
      end if
!
      return
      end


!     ======================
!     SUBROUTINE SNAPSHOTSP
!     ======================

      subroutine snapshotsp
      use pumamod

      
!       ************
!       * orograpy *
!       ************
       
        call writesp(140,so,129,0,CV*CV,0.)
       
!       ************
!       * pressure *
!       ************
       
        call writesp(140,sp,152,0,1.0,log(psurf))
       
!       ***************
!       * temperature *
!       ***************
       
        do jlev = 1 , NLEV
           call writesp(140,st(1,jlev),130,jlev,ct,t0(jlev) * ct)
        enddo
       
!       *********************
!       * specific humidity *
!       *********************
       
        if (nqspec == 1) then
           do jlev = 1 , NLEV
              call writesp(140,sqout(1,jlev),133,jlev,1.0,0.0)
           enddo
        endif
       
!       **************
!       * divergence *
!       **************
       
        do jlev = 1 , NLEV
           call writesp(140,sd(1,jlev),155,jlev,ww,0.0)
        enddo
       
!       *************
!       * vorticity *
!       *************
       
        do jlev = 1 , NLEV
           zsave = sz(3,jlev)
           sz(3,jlev) = sz(3,jlev) - plavor
           call writesp(140,sz(1,jlev),138,jlev,ww,0.0)
           sz(3,jlev) = zsave
        enddo

      
      return
      end

!     ======================
!     SUBROUTINE HCADENCESP
!     ======================

      subroutine hcadencesp(kunit)
      use pumamod

      integer, intent(in) :: kunit
      
!       ************
!       * orograpy *
!       ************
       
        call writesp(kunit,so,129,0,CV*CV,0.)
       
!       ************
!       * pressure *
!       ************
       
        call writesp(kunit,sp,152,0,1.0,log(psurf))
       
!       ***************
!       * temperature *
!       ***************
       
        do jlev = 1 , NLEV
           call writesp(kunit,st(1,jlev),130,jlev,ct,t0(jlev) * ct)
        enddo
       
!       *********************
!       * specific humidity *
!       *********************
       
        if (nqspec == 1) then
           do jlev = 1 , NLEV
              call writesp(kunit,sqout(1,jlev),133,jlev,1.0,0.0)
           enddo
        endif
       
!       **************
!       * divergence *
!       **************
       
        do jlev = 1 , NLEV
           call writesp(kunit,sd(1,jlev),155,jlev,ww,0.0)
        enddo
       
!       *************
!       * vorticity *
!       *************
       
        do jlev = 1 , NLEV
           zsave = sz(3,jlev)
           sz(3,jlev) = sz(3,jlev) - plavor
           call writesp(kunit,sz(1,jlev),138,jlev,ww,0.0)
           sz(3,jlev) = zsave
        enddo

      
      return
      end

!     =====================
!     SUBROUTINE SNAPSHOTGP
!     =====================

      subroutine snapshotgp
      use pumamod
      use carbonmod
      use radmod
      use glaciermod

!       *********************
!       * specific humidity *
!       *********************
        
        if (nqspec == 0) then ! Semi Langrangian advection active
           do jlev = 1 , NLEV
              call writegp(140,dq(1,jlev),133,jlev)
           enddo
        endif
        
!       **********************************
!       * mixed-layer depth (from ocean) *
!       **********************************
        
        call writegp(140,dmld,110,0)
        
!       ***********************
!       * surface temperature *
!       ***********************
        
        call writegp(140,dt(1,NLEP),139,0)
        
!       ****************
!       * soil wetness *
!       ****************
        
        call writegp(140,dwatc,140,0)
        
!       **************
!       * snow depth *
!       **************
        
        call writegp(140,dsnow,141,0)
     
      
!     **********************
!     * large scale precip *
!     **********************

      call writegp(140,dprl,142,0)

!     *********************
!     * convective precip *
!     *********************

      call writegp(140,dprc,143,0)

!     *************
!     * snow fall *
!     *************

      call writegp(140,dprs,144,0)

!     **********************
!     * sensible heat flux *
!     **********************

      call writegp(140,dshfl,146,0)

!     ********************
!     * latent heat flux *
!     ********************

      call writegp(140,dlhfl,147,0)

!       ************************
!       * liquid water content *
!       ************************
     
        do jlev = 1 , NLEV
           call writegp(140,dql(1,jlev),161,jlev)
        enddo
     
!       *************
!       * u-star**3 *
!       *************
     
        call writegp(140,dust3,159,0)
     
     
!     **********
!     * runoff *
!     **********

      call writegp(140,drunoff,160,0)

!     ***************
!     * cloud cover *
!     ***************

!         cl
        do jlev = 1 , NLEV
          call writegp(140,dcc(1,jlev),162,jlev)
        enddo
        
!         clt        
      call writegp(140,dcc(:,NLEP),164,0)
     

!     ***************************
!     * surface air temperature *
!     ***************************

      call writegp(140,dtsa,167,0)

!     ******************************
!     * surface temperature (accu) *
!     ******************************

      call writegp(140,dt(:,NLEP),169,0)

      
!       *************************
!       * deep soil temperature *
!       *************************
      
        call writegp(140,dtd5,170,0)
      
!       *****************
!       * land sea mask *
!       *****************
      
        call writegp(140,dls,172,0)
      
!       *********************
!       * surface roughness *
!       *********************
      
        call writegp(140,dz0,173,0)
      
!       **********
!       * albedo *
!       **********
        do klight = 1, NLIGHTS
          call writegp(140,dalb(:,klight),175,klight)
        enddo
        do klight = 1, NLIGHTS
          call writegp(140,dsalb(:,2*klight-1),174,klight)
        enddo
        do klight = 1, NLIGHTS
          call writegp(140,dsalb(:,2*klight),184,klight)
        enddo
        

!     ***************************
!     * surface solar radiation *
!     ***************************

      call writegp(140,dswfl(:,NLEP),176,0)

!     *****************************
!     * surface thermal radiation *
!     *****************************

      call writegp(140,dlwfl(:,NLEP),177,0)

!     ***********************
!     * top solar radiation *
!     ***********************

      call writegp(140,dswfl(:,1),178,0)

!     *************************
!     * top thermal radiation *
!     *************************

      call writegp(140,dlwfl(:,1),179,0)

!     ************
!     * u-stress *
!     ************

      call writegp(140,dtaux,180,0)

!     *************
!     * v- stress *
!     *************

      call writegp(140,dtauy,181,0)

!     ***************
!     * evaporation *
!     ***************

      call writegp(140,devap,182,0)

!     *********************
!     * soil temperature *
!     *********************
      
        call writegp(140,dtsoil,183,0)
      
      
!     ***********************************
!     * maximum surface air temperature *
!     ***********************************

      call writegp(140,atsama,201,0)

!     ***********************************
!     * minimum surface air temperature *
!     ***********************************

      call writegp(140,atsami,202,0)

!     ********************
!     * top solar upward *
!     ********************

      call writegp(140,dfu(:,1),203,0)

!     ************************
!     * surface solar upward *
!     ************************

      call writegp(140,dfu(:,NLEP),204,0)
      
!     ************************
!     * radiation level data *
!     ************************

      do jlev=1,NLEP
        call writegp(140,dfu(:,jlev),404,jlev)
        call writegp(140,dfd(:,jlev),405,jlev)
        call writegp(140,dftu(:,jlev),406,jlev)
        call writegp(140,dftd(:,jlev),407,jlev)
      enddo
      do jlev=1,NLEV
        call writegp(140,dtdtlwr(:,jlev)+dtdtswr(:,jlev),408,jlev)
        call writegp(140,dconv(:,jlev),409,jlev)
      enddo

!     **************************
!     * surface thermal upward *
!     **************************

      call writegp(140,dftu(:,NLEP),205,0)

!       *******************************
!       * soil temperatures level 2-4 *
!       *******************************
       
        call writegp(140,dtd2,207,0)
        call writegp(140,dtd3,208,0)
        call writegp(140,dtd4,209,0)
       
!       *****************
!       * sea ice cover *
!       *****************
       
        call writegp(140,dicec,210,0)
       
!       *********************
!       * sea ice thickness *
!       *********************
       
        call writegp(140,diced,211,0)
       
!       ****************
!       * forest cover *
!       ****************
       
        call writegp(140,dforest,298,0)
       
      
       
!     *************
!     * snow melt *
!     *************

      call writegp(140,dsmelt,218,0)

!     *********************
!     * snow depth change *
!     *********************

!       asndch(:)=asndch(:)/real(naccuout)
      call writegp(140,dsndch,221,0)

!     ******************
!     * field capacity *
!     ******************
        call writegp(140,dwmax,229,0)
        
        
!     *****************************************
!     * vertical integrated specific humidity *
!     *****************************************

      call writegp(140,dqvi,230,0)

!     ****************
!     * glacier mask *
!     ****************

      call writegp(140,dglac,232,0)
        
        
!     *********************
!     ***   S I M B A   ***
!     *********************

      
      if (nveg > 0) call vegsnapout(140)

!     ***********************
!     * Ozone concentration *
!     ***********************

        do jlev = 1 , NLEV
           call writegp(140,dqo3(1,jlev),265,jlev)
        enddo
      
!     ********************
!     * local weathering *
!     ********************

      call writegp(140,localweathering,266,0)

      
!       ********************
!       * ground elevation *
!       ********************
       
        call writegp(140,groundoro,267,0)
       
!       *********************
!       * glacier elevation *
!       *********************
       
        call writegp(140,glacieroro,268,0)
       
        
!       *********************
!       *    net elevation  *
!       *********************
        netoro(:) = groundoro(:) + glacieroro(:)
        call writegp(140,netoro,269,0)
        
       
!     ********************
!     * Cos Solar Zenith *
!     ********************
         
      call writegp(140,gmu0,318,0)
      
!     *****************************
!     * Weatherable Precipitation *
!     *****************************
        
      call writegp(140,asigrain,319,0)

!     ***********************
!     * Minimum Temperature *
!     ***********************
         
      call writegp(140,tempmin,320,0)

!     ***********************
!     * Maximum Temperature *
!     ***********************
         
      call writegp(140,tempmax,321,0)
      
      ! Hurricane quantities
      !Convective Available Potential Energy at the surface
      call writegp(140,capen,322,0)
      !Level of neutral buoyancy
      call writegp(140,lnb,323,0)
      ! Tropospheric entropy deficit
      call writegp(140,chim,324,0)
      ! Tropospheric entropy deficit
      call writegp(140,laav,325,0)
      ! Maximum potential intensity 
      call writegp(140,mpoti,326,0)
      ! Ventilation index
      call writegp(140,venti,327,0)      
      ! Ventilation-reduced Maximum Potential Intensityreturn
      call writegp(140,vrmpi,328,0)
      ! Genesis Potential Index
      call writegp(140,gpi,329,0)
      
      return
      end
      
!     =====================
!     SUBROUTINE SNAPSHOTSC
!     =====================
!     !Output scalars of interest

      subroutine snapshotsc
      use pumamod
      use radmod
      
        !True anomaly
        call writescalar(140,orbnu*180./PI,50)
        
        !Solar ecliptic longitude
        call writescalar(140,lambm*180./PI,51)
        
        !Solar declination
        call writescalar(140,zdeclf*180./PI,52)
        
        !Solar distance modulus
        call writescalar(140,1.0/sqrt(eccf),53)
        
        !Solar right ascension
        call writescalar(140,rasc*180./PI,54)
      
      return
      end
      
!     =====================
!     SUBROUTINE HCADENCEGP
!     =====================

      subroutine hcadencegp(kunit)
      use pumamod
      use carbonmod
      use radmod
      use glaciermod
      
      integer, intent(in) :: kunit

!       *********************
!       * specific humidity *
!       *********************
        
        if (nqspec == 0) then ! Semi Langrangian advection active
           do jlev = 1 , NLEV
              call writegp(kunit,dq(1,jlev),133,jlev)
           enddo
        endif
        
!       **********************************
!       * mixed-layer depth (from ocean) *
!       **********************************
        
        call writegp(kunit,dmld,110,0)
        
!       ***********************
!       * surface temperature *
!       ***********************
        
        call writegp(kunit,dt(1,NLEP),139,0)
        
!       ****************
!       * soil wetness *
!       ****************
        
        call writegp(kunit,dwatc,140,0)
        
!       **************
!       * snow depth *
!       **************
        
        call writegp(kunit,dsnow,kunit,0)
     
      
!     **********************
!     * large scale precip *
!     **********************

      call writegp(kunit,dprl,142,0)

!     *********************
!     * convective precip *
!     *********************

      call writegp(kunit,dprc,143,0)

!     *************
!     * snow fall *
!     *************

      call writegp(kunit,dprs,144,0)

!     **********************
!     * sensible heat flux *
!     **********************

      call writegp(kunit,dshfl,146,0)

!     ********************
!     * latent heat flux *
!     ********************

      call writegp(kunit,dlhfl,147,0)

!       ************************
!       * liquid water content *
!       ************************
     
        do jlev = 1 , NLEV
           call writegp(kunit,dql(1,jlev),161,jlev)
        enddo
     
!       *************
!       * u-star**3 *
!       *************
     
        call writegp(kunit,dust3,159,0)
     
     
!     **********
!     * runoff *
!     **********

      call writegp(kunit,drunoff,160,0)

!     ***************
!     * cloud cover *
!     ***************

!         cl
        do jlev = 1 , NLEV
          call writegp(kunit,dcc(1,jlev),162,jlev)
        enddo
        
!         clt        
      call writegp(kunit,dcc(:,NLEP),164,0)
     

!     ***************************
!     * surface air temperature *
!     ***************************

      call writegp(kunit,dtsa,167,0)

!     ******************************
!     * surface temperature (accu) *
!     ******************************

      call writegp(kunit,dt(:,NLEP),169,0)

      
!       *************************
!       * deep soil temperature *
!       *************************
      
        call writegp(kunit,dtd5,170,0)
      
!       *****************
!       * land sea mask *
!       *****************
      
        call writegp(kunit,dls,172,0)
      
!       *********************
!       * surface roughness *
!       *********************
      
        call writegp(kunit,dz0,173,0)
      
!       **********
!       * albedo *
!       **********
        do klight=1,NLIGHTS
           call writegp(kunit,dalb(:,klight),175,klight)
        enddo
        do klight=1,NLIGHTS
           call writegp(kunit,dsalb(:,2*klight-1),174,0)
        enddo
        do klight=1,NLIGHTS
           call writegp(kunit,dsalb(:,2*klight),184,0)
        enddo
        

!     ***************************
!     * surface solar radiation *
!     ***************************

      call writegp(kunit,dswfl(:,NLEP),176,0)

!     *****************************
!     * surface thermal radiation *
!     *****************************

      call writegp(kunit,dlwfl(:,NLEP),177,0)

!     ***********************
!     * top solar radiation *
!     ***********************

      call writegp(kunit,dswfl(:,1),178,0)

!     *************************
!     * top thermal radiation *
!     *************************

      call writegp(kunit,dlwfl(:,1),179,0)

!     ************
!     * u-stress *
!     ************

      call writegp(kunit,dtaux,180,0)

!     *************
!     * v- stress *
!     *************

      call writegp(kunit,dtauy,181,0)

!     ***************
!     * evaporation *
!     ***************

      call writegp(kunit,devap,182,0)

!     *********************
!     * soil temperature *
!     *********************
      
        call writegp(kunit,dtsoil,183,0)
      
      
!     ***********************************
!     * maximum surface air temperature *
!     ***********************************

      call writegp(kunit,atsama,201,0)

!     ***********************************
!     * minimum surface air temperature *
!     ***********************************

      call writegp(kunit,atsami,202,0)

!     ********************
!     * top solar upward *
!     ********************

      call writegp(kunit,dfu(:,1),203,0)

!     ************************
!     * surface solar upward *
!     ************************

      call writegp(kunit,dfu(:,NLEP),204,0)
      
!     ************************
!     * radiation level data *
!     ************************

      do jlev=1,NLEP
        call writegp(kunit,dfu(:,jlev),404,jlev)
        call writegp(kunit,dfd(:,jlev),405,jlev)
        call writegp(kunit,dftu(:,jlev),406,jlev)
        call writegp(kunit,dftd(:,jlev),407,jlev)
      enddo
      do jlev=1,NLEV
        call writegp(kunit,dtdtlwr(:,jlev)+dtdtswr(:,jlev),408,jlev)
        call writegp(kunit,dconv(:,jlev),409,jlev)
      enddo

!     **************************
!     * surface thermal upward *
!     **************************

      call writegp(kunit,dftu(:,NLEP),205,0)

!       *******************************
!       * soil temperatures level 2-4 *
!       *******************************
       
        call writegp(kunit,dtd2,207,0)
        call writegp(kunit,dtd3,208,0)
        call writegp(kunit,dtd4,209,0)
       
!       *****************
!       * sea ice cover *
!       *****************
       
        call writegp(kunit,dicec,210,0)
       
!       *********************
!       * sea ice thickness *
!       *********************
       
        call writegp(kunit,diced,211,0)
       
!       ****************
!       * forest cover *
!       ****************
       
        call writegp(kunit,dforest,298,0)
       
      
       
!     *************
!     * snow melt *
!     *************

      call writegp(kunit,dsmelt,218,0)

!     *********************
!     * snow depth change *
!     *********************

!       asndch(:)=asndch(:)/real(naccuout)
      call writegp(kunit,dsndch,221,0)

!     ******************
!     * field capacity *
!     ******************
        call writegp(kunit,dwmax,229,0)
        
        
!     *****************************************
!     * vertical integrated specific humidity *
!     *****************************************

      call writegp(kunit,dqvi,230,0)

!     ****************
!     * glacier mask *
!     ****************

      call writegp(kunit,dglac,232,0)
        

!     *********************
!     ***   S I M B A   ***
!     *********************

      
      if (nveg > 0) call vegsnapout(kunit)

!     ***********************
!     * Ozone concentration *
!     ***********************

        do jlev = 1 , NLEV
           call writegp(kunit,dqo3(1,jlev),265,jlev)
        enddo
      
!     ********************
!     * local weathering *
!     ********************

      call writegp(kunit,localweathering,266,0)

      
!       ********************
!       * ground elevation *
!       ********************
       
        call writegp(kunit,groundoro,267,0)
       
!       *********************
!       * glacier elevation *
!       *********************
       
        call writegp(kunit,glacieroro,268,0)
       
        
!       *********************
!       *    net elevation  *
!       *********************
        netoro(:) = groundoro(:) + glacieroro(:)
        call writegp(kunit,netoro,269,0)
        
       
!     ********************
!     * Cos Solar Zenith *
!     ********************
         
      call writegp(kunit,gmu0,318,0)
      
!     *****************************
!     * Weatherable Precipitation *
!     *****************************
        
      call writegp(kunit,asigrain,319,0)

!     ***********************
!     * Minimum Temperature *
!     ***********************
         
      call writegp(kunit,tempmin,320,0)

!     ***********************
!     * Maximum Temperature *
!     ***********************
         
      call writegp(kunit,tempmax,321,0)
         
      ! Hurricane quantities
      !Convective Available Potential Energy at the surface
      call writegp(kunit,capen,322,0)
      !Level of neutral buoyancy
      call writegp(kunit,lnb,323,0)
      ! Tropospheric entropy deficit
      call writegp(kunit,chim,324,0)
      ! Tropospheric entropy deficit
      call writegp(kunit,laav,325,0)
      ! Maximum potential intensity 
      call writegp(kunit,mpoti,326,0)
      ! Ventilation index
      call writegp(kunit,venti,327,0)      
      ! Ventilation-reduced Maximum Potential Intensityreturn
      call writegp(kunit,vrmpi,328,0)
      ! Genesis Potential Index
      call writegp(kunit,gpi,329,0)               
            
      return
      end
         
!     =====================
!     SUBROUTINE HCADENCESC
!     =====================
!     !Output scalars of interest

      subroutine hcadencesc(kunit)
      use pumamod
      use radmod
      
      integer, intent(in) :: kunit
      
        !True anomaly
        call writescalar(kunit,orbnu*180./PI,50)
        
        !Solar ecliptic longitude
        call writescalar(kunit,lambm*180./PI,51)
        
        !Solar declination
        call writescalar(kunit,zdeclf*180./PI,52)
        
        !Solar distance modulus
        call writescalar(kunit,1.0/sqrt(eccf),53)
        
        !Solar right ascension
        call writescalar(kunit,rasc*180./PI,54)
        
      
      return
      end
      
      
!     ==================
!     SUBROUTINE SNAPSHOTDIAG
!     ==================

      subroutine snapshotdiag
      use pumamod

!     *****************************************
!     * 2-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiagsp2d > 0 .and. mypid == NROOT) then
       do jdiag=1,ndiagsp2d
        jcode=50+jdiag
        call writesp(140,dsp2d(1,jdiag),jcode,0,1.,0.0)
       enddo
      end if

!     *****************************************
!     * 3-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiagsp3d > 0 .and. mypid == NROOT) then
       do jdiag=1,ndiagsp3d
        jcode=60+jdiag
        do jlev=1,NLEV
         call writesp(140,dsp3d(1,jlev,jdiag),jcode,jlev,1.,0.0)
        enddo
       enddo
      end if

!     *****************************************
!     * 2-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiaggp2d > 0) then
       do jdiag=1,ndiaggp2d
        jcode=jdiag
        call writegp(140,dgp2d(1,jdiag),jcode,0)
       enddo
      end if

!     *****************************************
!     * 3-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiaggp3d > 0) then
       do jdiag=1,ndiaggp3d
        jcode=20+jdiag
        do jlev=1,NLEV
         call writegp(140,dgp3d(1,jlev,jdiag),jcode,jlev)
        enddo
       enddo
      end if

!     ************************************************
!     * cloud forcing (clear sky fluxes) diagnostics *
!     ************************************************

      if(ndiagcf > 0) then
       call writegp(140,dclforc(1,1),101,0)
       call writegp(140,dclforc(1,2),102,0)
       call writegp(140,dclforc(1,3),103,0)
       call writegp(140,dclforc(1,4),104,0)
       call writegp(140,dclforc(1,5),105,0)
       call writegp(140,dclforc(1,6),106,0)
       call writegp(140,dclforc(1,7),107,0)
      end if

!     **************************************
!     * entropy diagnostics if switched on *
!     **************************************

      if(nentropy > 0) then
       do jdiag=1,36
        jcode=319+jdiag
        if(jcode == 333) cycle                      !333 is reserved
        call writegp(140,dentropy(1,jdiag),jcode,0)
       enddo
      end if
      if(nentro3d > 0) then
       do jdiag=1,23
        jcode=419+jdiag
        do jlev=1,NLEV
         call writegp(140,dentro3d(1,jlev,jdiag),jcode,jlev)
        enddo
       enddo
      end if

!     *************************************
!     * energy diagnostics if switched on *
!     *************************************

      if(nenergy > 0) then
       do jdiag=1,28
        jcode=359+jdiag
        call writegp(140,denergy(1,jdiag),jcode,0)
       enddo
      end if
      if(nener3d > 0) then
       do jdiag=1,28
        jcode=459+jdiag
        do jlev=1,NLEV
         call writegp(140,dener3d(1,jlev,jdiag),jcode,jlev)
        enddo
       enddo
      end if
!
      return
      end
      
     
!     ==================
!     SUBROUTINE HCADENCEDIAG
!     ==================

      subroutine hcadencediag(kunit)
      use pumamod

      integer, intent(in) :: kunit
      
!     *****************************************
!     * 2-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiagsp2d > 0 .and. mypid == NROOT) then
       do jdiag=1,ndiagsp2d
        jcode=50+jdiag
        call writesp(kunit,dsp2d(1,jdiag),jcode,0,1.,0.0)
       enddo
      end if

!     *****************************************
!     * 3-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiagsp3d > 0 .and. mypid == NROOT) then
       do jdiag=1,ndiagsp3d
        jcode=60+jdiag
        do jlev=1,NLEV
         call writesp(kunit,dsp3d(1,jlev,jdiag),jcode,jlev,1.,0.0)
        enddo
       enddo
      end if

!     *****************************************
!     * 2-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiaggp2d > 0) then
       do jdiag=1,ndiaggp2d
        jcode=jdiag
        call writegp(kunit,dgp2d(1,jdiag),jcode,0)
       enddo
      end if

!     *****************************************
!     * 3-D diagnostic arrays, if switched on *
!     *****************************************

      if(ndiaggp3d > 0) then
       do jdiag=1,ndiaggp3d
        jcode=20+jdiag
        do jlev=1,NLEV
         call writegp(kunit,dgp3d(1,jlev,jdiag),jcode,jlev)
        enddo
       enddo
      end if

!     ************************************************
!     * cloud forcing (clear sky fluxes) diagnostics *
!     ************************************************

      if(ndiagcf > 0) then
       call writegp(kunit,dclforc(1,1),101,0)
       call writegp(kunit,dclforc(1,2),102,0)
       call writegp(kunit,dclforc(1,3),103,0)
       call writegp(kunit,dclforc(1,4),104,0)
       call writegp(kunit,dclforc(1,5),105,0)
       call writegp(kunit,dclforc(1,6),106,0)
       call writegp(kunit,dclforc(1,7),107,0)
      end if

!     **************************************
!     * entropy diagnostics if switched on *
!     **************************************

      if(nentropy > 0) then
       do jdiag=1,36
        jcode=319+jdiag
        if(jcode == 333) cycle                      !333 is reserved
        call writegp(kunit,dentropy(1,jdiag),jcode,0)
       enddo
      end if
      if(nentro3d > 0) then
       do jdiag=1,23
        jcode=419+jdiag
        do jlev=1,NLEV
         call writegp(kunit,dentro3d(1,jlev,jdiag),jcode,jlev)
        enddo
       enddo
      end if

!     *************************************
!     * energy diagnostics if switched on *
!     *************************************

      if(nenergy > 0) then
       do jdiag=1,28
        jcode=359+jdiag
        call writegp(kunit,denergy(1,jdiag),jcode,0)
       enddo
      end if
      if(nener3d > 0) then
       do jdiag=1,28
        jcode=459+jdiag
        do jlev=1,NLEV
         call writegp(kunit,dener3d(1,jlev,jdiag),jcode,jlev)
        enddo
       enddo
      end if
!
      return
      end
      
      
!     ===================
!     SUBROUTINE OUTRESET
!     ===================

      subroutine outreset
      use pumamod
      use carbonmod
!
!     reset accumulated arrays and counter
!
      azmuz(:)=0.
      aprl(:)=0.
      aprc(:)=0.
      aprs(:)=0.
      aevap(:)=0.
      ashfl(:)=0.
      alhfl(:)=0.
      acc(:)=0.
      assol(:)=0.
      asthr(:)=0.
      atsol(:)=0.
      atthr(:)=0.
      assolu(:)=0.
      asthru(:)=0.
      atsolu(:)=0.
      ataux(:)=0.
      atauy(:)=0.
      aroff(:)=0.
      asmelt(:)=0.
!       asndch(:)=0. !Let the net snow change keep accumulating
      aqvi(:)=0.
      atsa(:)=0.
      ats0(:)=0.
      atsami(:)=1.E10
      atsama(:)=0.
      aweathering(:) = 0.
      
      asigrain(:) = 0.
      tempmax(:) = 0.
      tempmin(:) = 1.0e3
      
      if (nlowio > 0) then
        aorbnu = 0.
        alambm = 0.
        azdecl = 0.
        ardist = 0.
        arasc = 0.
      
        aaso(:)  = 0.
        aasp(:)  = 0.
        do j=1,NLEV
          aast(:,j)    = 0.
          aasqout(:,j) = 0.
          aasd(:,j)    = 0.
          aasz(:,j)    = 0.
          aadqo3(:,j)  = 0.
        enddo
        do j=1,NLEP
          aadq(:,j)  = 0.
          aadt(:,j)  = 0.
          aadql(:,j) = 0.
          aadcc(:,j) = 0.
        enddo
        aadmld(:)       = 0.
        aadwatc(:)      = 0.
        aadsnow(:)      = 0.
        aadust3(:)      = 0.
        aadtd5(:)       = 0.
        aadls(:)        = 0.
        aadz0(:)        = 0.
        do k=1,NLIGHTS
           aadalb(:,k)       = 0.
           aadsalb1(:,k)     = 0.
           aadsalb2(:,k)     = 0.
        enddo
        aadtsoil(:)     = 0.
        aadtd2(:)       = 0.
        aadtd3(:)       = 0.
        aadtd4(:)       = 0.
        aadicec(:)      = 0.
        aadiced(:)      = 0.
        aadforest(:)    = 0.
        aadwmax(:)      = 0.
        aadglac(:)      = 0.
        aagroundoro(:)  = 0.
        aaglacieroro(:) = 0.
        
        agpi(:)   = 0.
        aventi(:) = 0.
        alaav(:)  = 0.
        ampoti(:) = 0.
        avrmpi(:) = 0.
        acapen(:) = 0.
        alnb(:)   = 0.
        achim(:)  = 0.

      endif
      
      naccuout=0

!     ************************************************
!     * cloud forcing (clear sky fluxes) diagnostics *
!     ************************************************

      if(ndiagcf > 0) then
       dclforc(:,:)=0.
      end if
!
      return
      end

!     ==================
!     SUBROUTINE OUTACCU
!     ==================

      subroutine outaccu
      use pumamod
      use carbonmod
      use radmod
      use glaciermod
!
!     accumulate diagnostic arrays
!

      where (dls(:) .gt. 0.5) !weatherable precipitation
         where (dt(:,NLEP) .gt. 273.15)
            asigrain(:)=asigrain(:)+(dprl(:)+dprc(:))*8.64e7 ![mm/day]
         endwhere
      endwhere
      do i=1,NHOR
        tempmax(i) = MAX(tempmax(i),dt(i,NLEP))
        tempmin(i) = MIN(tempmin(i),dt(i,NLEP))
      enddo
      
      aprl(:)=aprl(:)+dprl(:)
      aprc(:)=aprc(:)+dprc(:)
      aprs(:)=aprs(:)+dprs(:)
      aevap(:)=aevap(:)+devap(:)
      ashfl(:)=ashfl(:)+dshfl(:)
      alhfl(:)=alhfl(:)+dlhfl(:)
      acc(:)=acc(:)+dcc(:,NLEP)
      assol(:)=assol(:)+dswfl(:,NLEP)
      asthr(:)=asthr(:)+dlwfl(:,NLEP)
      atsol(:)=atsol(:)+dswfl(:,1)
      atthr(:)=atthr(:)+dlwfl(:,1)
      assolu(:)=assolu(:)+dfu(:,NLEP)
      asthru(:)=asthru(:)+dftu(:,NLEP)
      atsolu(:)=atsolu(:)+dfu(:,1)
      ataux(:)=ataux(:)+dtaux(:)
      atauy(:)=atauy(:)+dtauy(:)
      aroff(:)=aroff(:)+drunoff(:)
      asmelt(:)=asmelt(:)+dsmelt(:)
      asndch(:)=asndch(:)+dsndch(:)
      aqvi(:)=aqvi(:)+dqvi(:)
      atsa(:)=atsa(:)+dtsa(:)
      ats0(:)=ats0(:)+dt(:,NLEP)
      atsami(:)=AMIN1(atsami(:),dtsa(:))
      atsama(:)=AMAX1(atsama(:),dtsa(:))
      aweathering(:)=aweathering(:)+localweathering(:)
      azmuz(:) = azmuz(:)+gmu0(:)

      if (nlowio > 0) then
      
        aorbnu = aorbnu + orbnu
        alambm = alambm + lambm
        azdecl = azdecl + zdeclf
        ardist = ardist + 1.0/sqrt(eccf)
        arasc = arasc + rasc
      
        aaso(:) = aaso(:) + so(:)        
        aasp(:) = aasp(:) + sp(:)
        do j=1,NLEV
          aast(:,j)    = aast(:,j)    + st(:,j)        
          aasqout(:,j) = aasqout(:,j) + sqout(:,j)     
          aasd(:,j)    = aasd(:,j)    + sd(:,j)        
          aasz(:,j)    = aasz(:,j)    + sz(:,j) 
          aadqo3(:,j)  = aadqo3(:,j)  + dqo3(:,j)   
        enddo
        do j=1,NLEP
          aadq(:,j)  = aadq(:,j)  + dq(:,j)   
          aadt(:,j)  = aadt(:,j)  + dt(:,j)     
          aadql(:,j) = aadql(:,j) + dql(:,j)     
          aadcc(:,j) = aadcc(:,j) + dcc(:,j)            
        enddo
        aadmld(:)       = aadmld(:)       + dmld(:)     
        aadwatc(:)      = aadwatc(:)      + dwatc(:)     
        aadsnow(:)      = aadsnow(:)      + dsnow(:)       
        aadust3(:)      = aadust3(:)      + dust3(:)         
        aadtd5(:)       = aadtd5(:)       + dtd5(:)      
        aadls(:)        = aadls(:)        + dls(:)       
        aadz0(:)        = aadz0(:)        + dz0(:)
        do klight=1,NLIGHTS
           aadalb(:  ,klight)       = aadalb(:,klight)        + dalb(:,klight  )
           aadsalb1(:,klight)     = aadsalb1(:,klight)     + dsalb(:,2*klight-1)
           aadsalb2(:,klight)     = aadsalb2(:,klight)     + dsalb(:,2*klight  )
        enddo
        aadtsoil(:)     = aadtsoil(:)     + dtsoil(:)    
        aadtd2(:)       = aadtd2(:)       + dtd2(:)      
        aadtd3(:)       = aadtd3(:)       + dtd3(:)      
        aadtd4(:)       = aadtd4(:)       + dtd4(:)      
        aadicec(:)      = aadicec(:)      + dicec(:)     
        aadiced(:)      = aadiced(:)      + diced(:)     
        aadforest(:)    = aadforest(:)    + dforest(:)   
        aadwmax(:)      = aadwmax(:)      + dwmax(:)     
        aadglac(:)      = aadglac(:)      + dglac(:)       
        aagroundoro(:)  = aagroundoro(:)  + groundoro(:) 
        aaglacieroro(:) = aaglacieroro(:) + glacieroro(:)
        
        agpi(:)   = agpi(:)   + gpi(:)  
        aventi(:) = aventi(:) + venti(:)
        alaav(:)  = alaav(:)  + laav(:) 
        ampoti(:) = ampoti(:) + mpoti(:)
        avrmpi(:) = avrmpi(:) + vrmpi(:)
        acapen(:) = acapen(:) + capen(:)
        alnb(:)   = alnb(:)   + lnb(:)  
        achim(:)  = achim(:)  + chim(:) 
      
      endif
      
      naccuout=naccuout+1
!
      return
      end
