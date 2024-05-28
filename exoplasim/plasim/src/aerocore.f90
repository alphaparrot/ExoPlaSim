!=======================================================================
! Subroutines in this file have been adapted or rewritten 
! for the PlanetSimulator of University Hamburg 
! by Hui Wan (Max Planck Institute for Meteorology, 2009).
!=======================================================================
!
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine aerocore(mmr,numrhos,l_source,sigmah,dt,ps1,ps2, &
                        u,v,dtoa,dtdx,apart,rhop,fcoeff,  &
                        iord,jord,kord,                   &
                        NAERO,im,jm,nl,dap,dbk,           &
                        iml,j1,j2,js0,jn0,                &
                        cose,cosp,acosp,dlat,rcap,        &
                        cnst,deform,zcross,               &
                        fill,mfct,debug,nud,angle,land,   &
                        aerosw,l_aerorad)
!****6***0*********0*********0*********0*********0*********0**********72
!
! The subroutine aerocore is a duplicate of the tracer transport
! subroutine tpcore, with an added gravitational settling term in
! the vertical direction to account for the effect of gravity on
! solid phase particles (haze, dust, etc). The comments below are
! retained from the tpcore subroutine.
!
!
! Purpose: perform the transport of  3-D mixing ratio fields using 
!          externally specified winds on the hybrid Eta-coordinate.
!          One call to aerocore updates the 3-D mixing ratio
!          fields for one time step. [vertical mass flux is computed
!          internally using a center differenced hydrostatic mass
!          continuity equation].
!
! Schemes: Multi-dimensional Flux Form Semi-Lagrangian (FFSL) schemes
!          (Lin and Rood 1996, MWR) with a modified MFCT option
!          (Zalesak 1979).
!
! Multitasking version: 6m
! History: Original version by S.-J. Lin, Oct. 1, 1997
!          Modified for the PlanetSimulator by Hui Wan, May 2009
!
! Author of the algorithm: S.-J. Lin
! Address:
!                 Code 910.3, NASA/GSFC, Greenbelt, MD 20771
!                 Phone: 301-286-9540
!                 E-mail: lin@dao.gsfc.nasa.gov
!
! Affiliation:
!                 Joint Center for Earth Systems Technology
!                 The University of Maryland Baltimore County
!                 and  NASA - Goddard Space Flight Center
!
! The algorithm is based on the following papers:
! 
! 1. Lin, S.-J., and R. B. Rood, 1996: Multidimensional flux form semi-
!    Lagrangian transport schemes. Mon. Wea. Rev., 124, 2046-2070.
!
! 2. Lin, S.-J., W. C. Chao, Y. C. Sud, and G. K. Walker, 1994: A class of
!    the van Leer-type transport schemes and its applications to the moist-
!    ure transport in a General Circulation Model. Mon. Wea. Rev., 122,
!    1575-1593.
!
! ======
! INPUT:
! ======
! MMR(IM,JM,NL,NAERO): mixing ratios of aerosols at current time (t)
! NAERO: total # of constituent aerosols (1 by default)
! IM: first (E-W) dimension; # of Grid intervals in E-W is IM
! JM: 2nd (N-S) dimension;   # of Gaussian latitudes is JM
! NL: 3rd dimension (# of layers); vertical index increases from 1 at
!       the model top to NL near the surface (see fig. below).
!       It is assumed that NL > 5.
!
! PS1(IM,JM): surface pressure at current time (t)
! PS2(IM,JM): surface pressure at next time level (t+dt) 
!             predicted by the host model
! PS2 is replaced by the PS predicted in this routine (at t+dt) on output.
! Note: surface pressure can have any unit or can be multiplied by any
!       const.
!
! The hybrid ETA-coordinate:
!
! pressure at layer edges are defined as follows:
!
!        p(i,j,k) = AP(k)*PT  +  BP(k)*PS(i,j)          (1)
!
! Where PT is a constant having the same unit as PS.
! AP and BP are unitless constants given at layer edges.
! In all cases  BP(1) = 0., BP(NL+1) = 1.
! The pressure at the model top is PTOP = AP(1)*PT
!
! Assume the upper and lower boundaries of vertical layer k are 
! layer edges k and k+1, respectively (see the sketch below). 
! The layer thickness at grid point (i,j,k) reads 
!
!   delp(i,j,k) = dap(k) + dbk(k)*PS(i,j)
!
! where    dap(k) = ( AP(k+1) - AP(k) )*PT
!          dbk(k) =   BP(k+1) - BP(k)
!
! *********************
! For pure sigma system
! *********************
! AP(k) = 1 for all k, PT = PTOP,
! BP(k) = sige(k) (sigma at edges), PS = Psfc - PTOP, where Psfc
! is the true surface pressure.
! In this implementation, we calculate the values of dap and dbk in 
! the initialization phase of a model run, and pass them into this
! subroutine. Since the PlanetSimulator uses a pure sigma coordinate 
! with PT = 0, we have
!
!         dap(:) = 0
!         dbk(:) = dsigma(:) 
!
!
!
!                  /////////////////////////////////
!              / \ ------ Model top P=PTOP ---------  AP(1), BP(1)
!               |
!    delp(1)    |  ........... Q(i,j,1) ............
!               |
!     W(k=1)   \ / ---------------------------------  AP(2), BP(2)
!
!
!
!     W(k-1)   / \ ---------------------------------  AP(k), BP(k)
!               |
!    delp(K)    |  ........... Q(i,j,k) ............
!               |
!      W(k)    \ / ---------------------------------  AP(k+1), BP(k+1)
!
!
!
!              / \ ---------------------------------  AP(NL), BP(NL)
!               |
!    delp(NL)   |  ........... Q(i,j,NL) .........
!               |
!     W(NL)=0  \ / -----Earth's surface P=Psfc ------ AP(NL+1), BP(NL+1)
!                 //////////////////////////////////
!
! U(IM,JM,NL) & V(IM,JM,NL):winds (m/s) at mid-time-level (t+dt/2)
! Note that on return U and V are destroyed.
!
! dtoa(real): time step (in seconds) divided by the planet radius (in meters). 
!      Suggested value for time step: 30 min. for 4x5, 15 min. for 2x2.5
!      (Lat-Lon) resolution. Smaller values maybe needed if the model
!      has a well-resolved stratosphere and Max(V) > 225 m/s
!
! dtdx: real array of shape (JM). time step (in seconds) divided by the 
!       zonal grid size (in meters). The values are initialized in subroutine
!       tracer_ini.
!
! J1, J2 are the starting and ending indices of the Gaussian latitudes 
!        outside the polar caps. Note that in the tracer transport routines
!        we count from south to north.
!        The South polar cap edge is located between the first and second
!        Gaussian latitudes from the south; The North polar cap edge is 
!        located between the last two Gaussian latitudes. 
!
! js0, jn0: [1, js0] and [jn0, JM] are the latitude ranges in which 
!           the semi-Lagrangian method may be needed for transport in 
!           the x-direction. Their values are initialized in subroutine
!           tracer_ini.
!
! cose,cosp,acosp,dlat,rcap: these parameters are related to the 
!        locations of the Gaussiang latitudes. They are initialized in
!        subroutine tracer_ini.
!
! IORD, JORD, and KORD are integers controlling various options in E-W, N-S,
! and vertical transport, respectively. 
!
!
!  _ORD=
!     1: 1st order upstream scheme (too diffusive, not a real option; it
!        can be used for debugging purposes; this is THE only known "linear"
!        monotonic advection scheme.).
!     2: 2nd order van Leer (full monotonicity constraint;
!        see Lin et al 1994, MWR)
!     3: monotonic PPM* (Collela & Woodward 1984)
!     4: semi-monotonic PPM (same as 3, but overshoots are allowed)
!     5: positive-definite PPM (constraint on the subgrid distribution is
!        only strong enough to prevent generation of negative values;
!        both overshoots & undershootes are possible).
!     6: un-constrained PPM (nearly diffusion free; faster but
!        positivity of the subgrid distribution is not quaranteed. Use
!        this option only when the fields and winds are very smooth or
!        when MFCT=.true.)
!
! *PPM: Piece-wise Parabolic Method
!
! Recommended values:
! IORD=JORD=3 for high horizontal resolution.
! KORD=6       if MFCT=.true.
! KORD=3 or 5  if MFCT=.false.
!
! The implicit numerical diffusion decreases as _ORD increases.
! DO not use option 4 or 5 for non-positive definite scalars
! (such as Ertel Potential Vorticity).
!
! If numerical diffusion is a problem (particularly at low horizontal
! resolution) then the following setup is recommended:
! IORD=JORD=KORD=6 and MFCT=.true.
!
! FILL (logical):   flag to do filling for negatives (see note below).
! MFCT (logical):   flag to do a Zalesak-type Multidimensional Flux
!                   correction. It shouldn't be necessary to call the
!                   filling routine when MFCT is true.
!
! ======
! Output
! ======
!
! MMR: the updated mixing ratios at t+NDT (original values are over-written)
! W(;;NL): large-scale vertical mass flux as diagnosed from the hydrostatic
!          relationship. W will have the same unit as PS1 and PS2 (eg, mb).
!          W must be divided by dt to get the correct mass-flux unit.
!          The vertical Courant number C = W/delp_UPWIND, where delp_UPWIND
!          is the pressure thickness in the "upwind" direction. For example,
!          C(k) = W(k)/delp(k)   if W(k) > 0;
!          C(k) = W(k)/delp(k+1) if W(k) < 0.
!              ( W > 0 is downward, ie, toward surface)
! PS2: predicted PS at t+dt (original values are over-written)
!
! =====
! NOTES:
! =====
!
! This forward-in-time upstream-biased transport scheme degenerates to
! the 2nd order center-in-time center-in-space mass continuity eqn.
! if Q = 1 (constant fields will remain constant). This degeneracy ensures
! that the computed vertical velocity to be identical to GEOS-1 GCM
! for on-line transport.
!
! A larger polar cap is used if j1=3 (recommended for C-Grid winds or when
! winds are noisy near poles).
!
! PPM is 4th order accurate when grid spacing is uniform (x & y); 3rd
! order accurate for non-uniform grid (vertical sigma coord.).
!
! Time step is limitted only by transport in the meridional direction.
! (the FFSL scheme is not implemented in the meridional direction).
!
! Since only 1-D limiters are applied, negative values could
! potentially be generated when large time step is used and when the
! initial fields contain discontinuities.
! This does not necessarily imply the integration is unstable.
! These negatives are typically very small. A filling algorithm is
! activated if the user set "fill" to be true.
! Alternatively, one can use the MFCT option to enforce monotonicity.
!
      use pumamod, only: NLAT,NLON,NLEV,ga ! Use planet's gravity from pumamod
      implicit none

! Input-Output variables

      integer,intent(in) :: NAERO,im,jm,nl,iml,j1,j2,js0,jn0 ! Input variables passed from other modules
      integer,intent(in) :: iord,jord,kord,cnst,l_source,l_aerorad

! Input-Output arrays

      real,intent(inout) ::   mmr(im,jm,nl,NAERO) ! Mixing ratio of aerosol
      real,intent(inout) ::   numrhos(im,jm,nl,NAERO) ! Number density of aerosol in particles/m3 (converted from mass mixing ratio)
      real,intent(inout) :: ps1(im,jm) ! Surface pressure at time t
      real,intent(inout) :: ps2(im,jm)! Surface pressure at time t+dt
      real,intent(inout) ::   u(im,jm,nl) ! Zonal wind at t+dt/2
      real,intent(inout) ::   v(im,jm,nl) ! Meridional wind at t+dt/2
      real,intent(in)    ::  dtdx(jm), dtoa ! Ratio of timestep to zonal grid step; ratio of timestep to planetary radius
      real,intent(in)    :: sigmah(nl) ! Sigma at half level
      real,intent(in)    :: dt(im,jm,nl) ! Air temperature (dt in plasimmod)
      real,intent(in)    :: apart ! Particle radius
      real,intent(in)    :: rhop ! Particle density
      real,intent(in)    :: fcoeff ! Haze mass production rate at solar zenith in kg/m2/s

      real,intent(in)    :: cose(jm+1),cosp(jm), acosp(jm), dlat(jm), rcap ! Geometric factors calculated in tracer_ini subroutine in tracermod
      real,intent(in)    :: dap(nl), dbk(nl) ! Coefficients which control layer thickness for different height schemes (pure sigma, hybrid)
      logical,intent(in) :: zcross, fill, mfct, deform, debug ! Logicals set in aeromod

! Local dynamic arrays

      integer :: js(nl),jn(nl) !

      real :: dtdx5(jm) 
      real ::   w(im,jm,nl) ! Vertical mass flux
      real :: crx(im,jm,nl),cry(im,jm,nl) ! Courant number in x and y direction
      real :: delp(im,jm,nl),delp1(im,jm,nl),delp2(im,jm,nl) ! Layer thickness
      real :: delp2dyn(im,jm,nl) ! Layer thickness predicted by dynamical core outside this subroutine
      real :: xmass(im,jm,nl),ymass(im,jm,nl) ! Mass fluxes in x and y direction
      real ::   dg1(im),dg2(im,jm),dpi(im,jm,nl) ! 
      real ::  qlow(im,jm,nl), daero(im,jm,nl)!
      real ::    qz(im,jm,nl),qmax(im,jm,nl),qmin(im,jm,nl) !
      real ::    wk(im,jm,nl),pu(im,jm,nl) !
      real ::    fx(im+1,jm,nl),fy(im,jm,nl),fz(im,jm,nl+1),gz(im,jm,nl+1) ! x, y, z density-weighted mmr fluxes + gravitational settling flux
      real ::    mu(im,jm,nl) ! Viscosity of bulk gas calculated by viscos subroutine
      real ::   beta(im,jm,nl) ! Cunningham factor used in terminal velocity, calculated by chamfac subroutine
      real ::   vels(im,jm,nl) ! Terminal velocity of falling aerosol particles
      real ::    temp(im,jm,nl) ! Local temperature array, remove negative values and replace with 0 to avoid errors
      real ::   rhog(im,jm,nl) ! Array for gas density calculated in routine
      real ::   angle(im,jm) ! Array for cosine of solar zenith angle
      real ::   land(im,jm) ! Array for binary land mask
      real ::   aerosw(im,jm,nl) ! Array for net SW flux

! scalars

      integer :: imh,i,j,k,jt,ic,nud ! nud = I/O unit for diagnostic output, defined in pumamod
      real    :: d5, dtoa5, sum1, sum2 ! Values calculated below

!---------------------------------------------------------------
      
      imh      = im/2 ! Number of longitudes divided by 2
      dtdx5(:) = 0.5*dtdx(:) ! Ratio of timestep to zonal grid step divided by 2
      dtoa5    = 0.5*dtoa ! Ratio of timestep to planetary radius divided by 2

      if (debug) then ! For debugging, output existing max and min values for particle MMR
         write(nud,*) '* Entered routine aerocore'
         do ic=1,naero
            write(nud,*) '* tracer',ic,'* max. & min. mmr =', &
                       maxval(mmr(:,:,:,ic)),minval(mmr(:,:,:,ic))
         end do
      endif

! save the vertical layer thickness (at t+dt) predicted by 
! the dynamical core 

      do k=1,nl
         delp2dyn(:,:,k) = dap(k) + dbk(k)*ps2(:,:) ! Calculate layer thickness for t+dt using surface pressure at t+dt and coefficients
      end do

! estimate the value at t+0.5dt

      ps2(:,:) = 0.5*( ps1(:,:) + ps2(:,:) ) ! Calculate new ps2 (surface pressure at t+0.5dt)

!****6***0*********0*********0*********0*********0*********0**********72
! Compute Courant number
!****6***0*********0*********0*********0*********0*********0**********72
! Convert winds on A-Grid to Courant # on C-Grid.

      do j=2,jm-1 ! Latitude loop, skips polar caps
      do i=2,im ! Longitude loop
         crx(i,j,:) = dtdx5(j)*(u(i,j,:)+u(i-1,j,:)) ! Zonal Courant number for stability
      end do
         crx(1,j,:) = dtdx5(j)*(u(1,j,:)+u(im,j,:))
      end do
 
      do j=2,jm ! Latitude loop, skips polar caps
         cry(:,j,:) = dtoa5*(v(:,j,:)+v(:,j-1,:)) ! Meridional Courant number
      enddo
 
!****6***0*********0*********0*********0*********0*********0**********72
! Find JN and JS
!****6***0*********0*********0*********0*********0*********0**********72
! [2,js(k)] and [jn(k),jm-1] are the latitudes at which semi-Lagrangian
! treatment is needed.
 
!MIC$ do all autoscope shared(JS,JN,CRX,CRY,PS2,U,V,DPI,ymass,delp2,PU)
!MIC$* shared(xmass)
!MIC$* private(i,j,k,sum1,sum2,D5)

      do k=1,nl ! Vertical loop

        js(k) = j1 ! Southernmost lat outside polar cap
        jn(k) = j2 ! Northernmost lat outside polar cap
    
        outer1: do j=js0,j1+1,-1 ! Just for these lats
            do i=1,im ! For all lons
                if(abs(crx(i,j,k)) .gt. 1.) then
                        js(k) = j ! If zonal Courant number is >1, extend limit of semi-Lagrangian treatment at southern polar cap
                        exit outer1
                endif
            enddo
        enddo outer1

        outer2: do j=jn0,j2-1 ! For these lats
            do i=1,im ! For all lons
                if(abs(crx(i,j,k)) .gt. 1.) then 
                        jn(k) = j ! If zonal Courant number is >1, extend limit of semi-Lagrangian treatment at northern polar cap
                        exit outer2
                endif
            enddo
        enddo outer2

      enddo  !end vertical layer loop. We have set the limits of the caps for this timestep to ensure stability.
 
!****6***0*********0*********0*********0*********0*********0**********72
! ***** Compute horizontal mass fluxes *****
!****6***0*********0*********0*********0*********0*********0**********72
! for the polar caps: replace the values at the northest (southest)
! latitude with the zonal mean. 

      sum1 = sum(ps1(:,1 ))/im ! Surface pressure at time t at southernmost lat, summed at all lons, then divided by number of lons (zonal mean)
      sum2 = sum(ps1(:,jm))/im ! Same as above but for northernmost lat
      ps1(:,1)  = sum1 ! Replace southernmost lat with zonal mean in ps1 array
      ps1(:,jm) = sum2 ! Replace northernmost last with zonal mean in ps1 array

      sum1 = sum(ps2(:,1 ))/im ! Same as above but for surface pressure at time t+dt
      sum2 = sum(ps2(:,jm))/im
      ps2(:,1)  = sum1
      ps2(:,jm) = sum2


      do k=1,nl ! Vertical levels loop

! delp = pressure thickness: the psudo-density in a hydrostatic system.

      delp2(:,:,k) = dap(k) + dbk(k)*ps2(:,:) ! Pressure thickness at time t+dt, coefficients set by choice of vertical levels scheme (pure sigma, hybrid, etc.)
 
! N-S componenet
 
      do j=j1,j2+1 ! For latitudes excluding the polar caps
         d5 = 0.5 * cose(j) ! cose defined in trc_routines, cosine at edges
         ymass(:,j,k) = cry(:,j,k)*d5*(delp2(:,j,k)+delp2(:,j-1,k)) ! Meridional mass flux, prognostic, uses delp for t+dt
      enddo
 
      do j=j1,j2 ! For lats at edges of polar caps
         dpi(:,j,k) = (ymass(:,j,k)-ymass(:,j+1,k)) *acosp(j)/dlat(j) ! Difference in mass flux at adjacent lats times scaled geometrically somehow.
      end do
 
! polar caps

      sum1 = -sum(ymass(:,j1  ,k))*rcap ! (negative) Sum of y mass flux at southernmost lat x 1/cosphi/deltaphi (scaled by some geometric derivative)
      sum2 =  sum(ymass(:,j2+1,k))*rcap ! Sum of y mass flux at northernmost lat x 1/cosphi/deltaphi

      dpi(:, 1,k) = sum1 ! Tendency of pi (mass flux) at southernmost lat
      dpi(:,jm,k) = sum2 ! Tendency of pi at northernmost lat
 
! E-W component

      do j=j1,j2 ! Between edges of caps
      do i=2,im ! Lon loop
         pu(i,j,k) = 0.5 * (delp2(i,j,k) + delp2(i-1,j,k)) ! Slope of pressure thickness at t+dt
      enddo
      enddo
 
      do j=j1,j2
         pu(1,j,k) = 0.5 * (delp2(1,j,k) + delp2(im,j,k)) ! Same but close circle of lons
      enddo
 
      do j=j1,j2 ! Between edges of caps
      do i=1,im
         xmass(i,j,k) = pu(i,j,k)*crx(i,j,k) ! Zonal mass flux
      enddo
      enddo
 
      do j=j1,j2 ! Between edges of caps
      do i=1,im-1
         dpi(i,j,k) = dpi(i,j,k) + xmass(i,j,k) - xmass(i+1,j,k) ! Add derivative of zonal mass flux to tendency array
      enddo
      enddo
 
      do j=j1,j2
         dpi(im,j,k) = dpi(im,j,k) + xmass(im,j,k) - xmass(1,j,k) ! Same but close circle of lons
      enddo

      enddo ! vertical layer loop 

!****6***0*********0*********0*********0*********0*********0**********72
! Compute Courant number at cell center (upwinding)
!****6***0*********0*********0*********0*********0*********0**********72

! E-W direction

      do k=1,nl ! Vertical level loop

      do j=j1,j2 ! Between edges of caps
      do i=1,im-1 ! Lons loop
         if(crx(i,j,k)*crx(i+1,j,k) .gt. 0.) then ! If product of Courant numbers of adjacent lons is greater than 0
            if(crx(i,j,k) .gt. 0.) then ! And if Courant number at longitude i is greater than 0
            u(i,j,k) = crx(i,j,k) ! Zonal wind is equal to the Courant number in cell i
            else
            u(i,j,k) = crx(i+1,j,k) ! Otherwise zonal wind is equal to the Courant number in cell i+1
            endif
         else
            u(i,j,k) = 0. ! Otherwise zonal wind is 0
         endif
      enddo
      enddo
 
      i=im ! At final lon
      do j=j1,j2 ! Between edges of caps
      if(crx(i,j,k)*crx(1,j,k) .gt. 0.) then
         if(crx(i,j,k) .gt. 0.) then
         u(i,j,k) = crx(i,j,k) ! Zonal wind is equal to Courant number at im (final lon)
         else
         u(i,j,k) = crx(1,j,k) ! Or zonal wind is equal to Courant number at i=1
         endif
      else
         u(i,j,k) = 0. ! Or zonal wind equals 0
      endif
      enddo

      enddo ! vertical layer loop 

! N-S direction

      do k=1,nl ! Vertical level loop

      do j=j1,j2 ! Between edges of caps
      do i=1,im ! Lon loop
         if(cry(i,j,k)*cry(i,j+1,k) .gt. 0.) then
            if(cry(i,j,k) .gt. 0.) then
            v(i,j,k) = cry(i,j,k)/dlat(j) ! Meridional wind equals y Courant number divided by delta phi at j
            else
            v(i,j,k) = cry(i,j+1,k)/dlat(j) ! Meridional wind equals y Courant number divided by delta phi at j+1
            endif
         else
            v(i,j,k) = 0. ! Meridional wind equals 0
         endif
      enddo
      enddo

!++ to be checked ++ 
!     do 139 i=1,imh
!     v(i,     1,k) = 0.5*(cry(i,2,k)-cry(i+imh,2,k))
!     v(i+imh, 1,k) = -v(i,1,k)
!     v(i,    jm,k) = 0.5*(cry(i,jm,k)-cry(i+imh,jm1,k))
!139   v(i+imh,jm,k) = -v(i,jm,k)
!== to be checked == 

      enddo ! vertical layer loop 
 
!****6***0*********0*********0*********0*********0*********0**********72
! Compute vertical mass flux (same dimensional unit as PS)
!****6***0*********0*********0*********0*********0*********0**********72
 
! compute total column mass CONVERGENCE.
 
!MIC$ do all autoscope shared(im,jm,DPI,PS1,PS2,W,DBK)
!MIC$* shared(DPI,PS1,PS2,W,DBK)
!MIC$* private(i,j,k,DG1)

      do j=1,jm ! For all lats

         dg1(1:im) = 0. ! Initialise dg1 to 0. dg1 has length = number of lons. Reinitalised for each latitude.
         do k=1,nl ! Vertical level loop
            dg1(1:im) = dg1(1:im) + dpi(1:im,j,k) ! Now fill with values for dpi tendency array
         enddo
 
! Compute PS2 (PS at n+1) using the hydrostatic assumption.
! Changes (increases) to surface pressure = total column mass convergence
 
         ps2(1:im,j)  = ps1(1:im,j) + dg1(1:im) ! Prognostic surface pressure at t+dt calculated from surface press at t plus tendency of mass flux (dpi)
 
! compute vertical mass flux from mass conservation principle:
! lower boundary of the first (uppermost) layer 

         w(1:im,j,1) = dpi(1:im,j,1) - dbk(1)*dg1(1:im) ! For level 1 (the top), for each lat, vertical mass flux equals delta mass flux minus delta sigma*delta mass flux
 
         do k=2,nl-1 ! For other levels
         w(1:im,j,k) = w(1:im,j,k-1) + dpi(1:im,j,k) - dbk(k)*dg1(1:im) ! Vertical mass flux calculated using mass flux from layer above plus delta mass flux minus delta sigma*delta mass flux
         enddo

! Earth's surface 

         w(1:im,j,nl) = 0. ! Vertical mass flux is zero at nl, the lowest level

      enddo 
 
!MIC$ do all
!MIC$* shared(deform,NL,im,jm,delp,delp1,delp2,DPI,DAP,DBK,PS1,PS2)
!MIC$* private(i,j,k)

      do k=1,nl ! Vertical level loop

         delp1(:,:,k) = dap(k) + dbk(k)*ps1(:,:) ! Pressure thickness at time t, depends on choice of levels scheme via coefficients dap and dbk (see intro comment)
         delp2(:,:,k) = dap(k) + dbk(k)*ps2(:,:) ! Pressure thickness at time t+dt, depends on choice of levels scheme
         delp (:,:,k) = delp1(:,:,k) + dpi(:,:,k) ! Pressure thickness at time t plus tendency of horizontal mass flux
 
! Check deformation of the flow fields

        if(deform) then
          do j=1,jm
          do i=1,im
          if(delp(i,j,k) .le. 0.) then
             write(nud,'(a)') '* FFSL transport'
             write(nud,'(a,i3,a)') ' Vertical layer',k, &
                         'Noisy wind fields -> delp* is negative!'
             write(nud,*) '* Smooth the wind fields or reduce time step'
             stop
          endif
          enddo
          enddo
        endif

      enddo !vertical layer loop


    temp = dt
    where(temp .le. 1.0) temp = 1.0
      
!****6***0*********0*********0*********0*********0*********0**********72
! Calculate viscosity of bulk gas, store as mu (3D)
!****6***0*********0*********0*********0*********0*********0**********72

    call viscos(temp,im,jm,nl,nud,mu)
    
!****6***0*********0*********0*********0*********0*********0**********72
! Calculate Cunningham factor, store as beta (3D) + call calculation of
! gas density in SI units; both used in calculation of terminal velocity
!****6***0*********0*********0*********0*********0*********0**********72

    call chamfac(temp,im,jm,nl,sigmah,ps2,apart,nud,beta)

    call density(im,jm,nl,temp,ps2,sigmah,nud,rhog)
    
!****6***0*********0*********0*********0*********0*********0**********72
! Calculate terminal velocity, store as vels (3D)
!****6***0*********0*********0*********0*********0*********0**********72

    call vterm(im,jm,nl,beta,apart,rhop,mu,temp,sigmah,ps2,rhog,nud,vels)	

!****6***0*********0*********0*********0*********0*********0**********72
! Do transport one tracer at a time.
!****6***0*********0*********0*********0*********0*********0**********72
 
      DO ic=1,naero ! From aerosol 1 to however many there are
 
!MIC$ do all autoscope
!MIC$* shared(q,DQ,delp1,U,V,j1,j2,JS,JN,im,jm,IML,IC,IORD,JORD)
!MIC$* shared(CRX,CRY,PU,xmass,ymass,fx,fy,acosp,rcap,qz)
!MIC$* private(i,j,k,jt,wk,DG2)

      if (l_aerorad == 0) then
       select case (l_source) ! Choose your aerosol source
       case(1) ! Case 1: photochemical haze
         mmr(:,:,1,ic) = fcoeff*angle ! The coefficient fcoeff sets the haze mass production rate at the solar zenith at k=1
       case(2)
         mmr(:,:,NLEV,ic) = fcoeff*land ! At k=surface, land grid boxes are given the abundance fcoeff (kg/kg) and sea is given 0
       end select
      end if
      
      if (l_aerorad == 1) then
       select case (l_source) ! Choose your aerosol source
       case(1) ! Case 1: photochemical haze
         mmr(:,:,1,ic) = fcoeff*(abs(aerosw(:,:,1))/(maxval(abs(aerosw(:,:,1))) + 1.E-6)) ! At the max SW flux, the source strength is the input source
       case(2) ! Case 2: dust
         mmr(:,:,NLEV,ic) = fcoeff*land ! At k=surface, land grid boxes are given the abundance fcoeff (kg/kg) and sea is given 0
       end select
      end if
      
      do k=1,nl ! Vertical levels loop

! for the polar caps: replace the mixing ratio at the northest (southest)
! latitude with the zonal mean. 

         sum1 = sum(mmr(:,1 ,k,ic))/im ! Zonal mean of mixing ratio at southernmost lat
         sum2 = sum(mmr(:,jm,k,ic))/im ! Zonal mean of mixing ratio at northernmost lat
         mmr(:,1 ,k,ic) = sum1 ! Put these back into mixing ratio array
         mmr(:,jm,k,ic) = sum2
      
! Initialize DAERO
 
         daero(:,:,k) = mmr(:,:,k,ic)*delp1(:,:,k)  ! Mixing ratio multipled by pressure thickness (quasi-density) at time t

! E-W advective cross term
      call xadv(im,jm,j1,j2,mmr(1,1,k,ic),u(1,1,k),js(k),jn(k),iml, &
                wk(1,1,1))

      wk(:,:,1) = mmr(:,:,k,ic) + 0.5*wk(:,:,1) 
 
! N-S advective cross term
      do j=j1,j2
        do i=1,im
          jt = float(j) - v(i,j,k)
          wk(i,j,2) = v(i,j,k) * (mmr(i,jt,k,ic) - mmr(i,jt+1,k,ic))
        enddo
      enddo
 
      do j=j1,j2
        do i=1,im
          wk(i,j,2) = mmr(i,j,k,ic) + 0.5*wk(i,j,2)
        enddo
      enddo

!****6***0*********0*********0*********0*********0*********0**********72
! compute flux in  E-W direction - 2-D fluxes untouched from tpcore
      call xtp(im,jm,iml,j1,j2,jn(k),js(k),pu(1,1,k),daero(1,1,k), &
               wk(1,1,2),crx(1,1,k),fx(1,1,k),xmass(1,1,k),iord)

! compute flux in  N-S direction
      call ytp(im,jm,j1,j2,acosp,dlat,rcap,daero(1,1,k),wk(1,1,1), &
               cry(1,1,k),dg2,ymass(1,1,k),wk(1,1,3),wk(1,1,4), &
               wk(1,1,5),wk(1,1,6),fy(1,1,k),jord,nud)
!****6***0*********0*********0*********0*********0*********0**********72

      if(ZCROSS) then

! qz is the horizontal advection modified value for input to the
! vertical transport operator FZPPM
! Note: DQ contains only first order upwind contribution.

        do j=1,JM
            do i=1,IM
            qz(i,j,k) = daero(i,j,k) / delp(i,j,k)
            enddo
        enddo

      else

        do j=1,JM
            do i=1,IM
            qz(i,j,k) = mmr(i,j,k,IC)
            enddo
        enddo

      endif
       
      enddo     ! end of k-loop

 
!****6***0*********0*********0*********0*********0*********0**********72
! Compute fluxes in the vertical direction
      call FZPPM(qz,fz,IM,JM,NL,daero,W,delp,KORD)
!****6***0*********0*********0*********0*********0*********0**********72
      if (debug) write(nud,*) 'Max and min mmr =',maxval(mmr), minval(mmr)
!****6***0*********0*********0*********0*********0*********0**********72
! Compute vertical flux due to gravitational settling 
      call gsettle(qz,im,jm,nl,rhog,vels,nud,gz)
      
!****6***0*********0*********0*********0*********0*********0**********72
 
      if( MFCT ) then ! Some flux correction if needed

      if (debug) write(nud,*) '* mfct on' 
       
! qlow is the low order "monotonic" solution
 
!MIC$ do all
!MIC$* shared(NL,im,jm,j1,jm1,qlow,DQ,delp2)
!MIC$* private(i,j,k)

      DO k=1,NL ! For all levels

      DO j=1,JM ! For all lats
        DO i=1,IM ! For all lons
          qlow(i,j,k) = daero(i,j,k) / delp2(i,j,k) !
        enddo
      enddo
 
      if(j1.ne.2) then
        DO i=1,IM
            qlow(i,   2,k) = qlow(i, 1,k)
            qlow(i,jm-1,k) = qlow(i,jm,k)
        enddo
      endif

      enddo
 
!****6***0*********0*********0*********0*********0*********0**********72
      call FCT3D(mmr(1,1,1,IC),qlow,fx,fy,fz,IM,JM,NL,j1,j2,delp2, &
                 DPI,qz,wk,Qmax,Qmin,DG2,U,V,acosp,dlat,RCAP)
! Note: Q is destroyed!!!
!****6***0*********0*********0*********0*********0*********0**********72
      ENDIF
 
! Final update

!MIC$ do all autoscope
!MIC$* private(i,j,k,sum1,sum2)

! For k=1, update top level by subtracting downward flux from this level 
      do j=j1,j2 
        do i=1,IM
            daero(i,j,1) = daero(i,j,1) +  fx(i,j,1) - fx(i+1,j,1)                   &
                            + (fy(i,j,1) - fy(i,j+1,1))*acosp(j)/dlat(j) &
                            +  fz(i,j,1) - fz(i,j,2)  &
                            - gz(i,j,1)*ga   
        enddo ! Lon loop
      enddo ! Lat loop
      
      ! poles:
      sum1 = fy(IM,j1  ,1) ! fy at last lon, southernmost lat, level k
      sum2 = fy(IM,J2+1,1) ! fy at last lon, northernmost lat, level k

      do i=1,IM-1 ! For all other lons
        sum1 = sum1 + fy(i,j1  ,1) ! Add fy for last lon at southernmost lat
        sum2 = sum2 + fy(i,J2+1,1) ! Add fy for last lon at northernmost lat
      enddo
 
      daero(1, 1,1) = daero(1, 1,1) - sum1*RCAP + fz(1, 1,1) - fz(1, 1,2) ! First lon, southernmost lat
      daero(1,JM,1) = daero(1,JM,1) + sum2*RCAP + fz(1,JM,1) - fz(1,JM,2) ! First lon, northernmost lat
 
      do i=2,IM ! All other lons except first
        daero(i, 1,1) = daero(1, 1,1) ! At southernmost lat, daero is equal to daero at first lon
        daero(i,JM,1) = daero(1,JM,1) ! At northernmost lat, daero is equal to daero at first lon
      enddo							

! For all levels except the top and bottom, add flux from level above and subtract flux falling out of current level
      do k=2,NL

        do j=j1,j2 ! Between the caps
            do i=1,IM ! For all lons
                daero(i,j,k) = daero(i,j,k) +  fx(i,j,k) - fx(i+1,j,k)             &
                                        + (fy(i,j,k) - fy(i,j+1,k))*acosp(j)/dlat(j) &
                                        +  fz(i,j,k) - fz(i,j,k+1)                   &
                                        + (gz(i,j,k-1) - gz(i,j,k))*ga                
            enddo
        enddo
    
    ! poles:
        sum1 = fy(IM,j1  ,k) ! fy at last lon, southernmost lat, level k
        sum2 = fy(IM,J2+1,k) ! fy at last lon, northernmost lat, level k

        do i=1,IM-1 ! For all other lons
            sum1 = sum1 + fy(i,j1  ,k) ! Add fy for last lon at southernmost lat
            sum2 = sum2 + fy(i,J2+1,k) ! Add fy for last lon at northernmost lat
        enddo
    
        daero(1, 1,k) = daero(1, 1,k) - sum1*RCAP + fz(1, 1,k) - fz(1, 1,k+1) ! First lon, southernmost lat
        daero(1,JM,k) = daero(1,JM,k) + sum2*RCAP + fz(1,JM,k) - fz(1,JM,k+1) ! First lon, northernmost lat
    
        do i=2,IM ! All other lons except first
            daero(i, 1,k) = daero(1, 1,k) ! At southernmost lat, daero is equal to daero at first lon
            daero(i,JM,k) = daero(1,JM,k) ! At northernmost lat, daero is equal to daero at first lon
        enddo

      enddo

! Now for k=nl, only add the flux coming from the level above
      do j=j1,j2 
        do i=1,IM
            daero(i,j,nl) = daero(i,j,nl) +  fx(i,j,nl) - fx(i+1,j,nl)                   &
                            + (fy(i,j,nl) - fy(i,j+1,nl))*acosp(j)/dlat(j) &
                            +  fz(i,j,nl) - fz(i,j,nl+1) &
                            + gz(i,j,nl)*ga
        enddo ! Lon loop
      enddo ! Lat loop
      
      ! poles:
      sum1 = fy(IM,j1  ,nl) ! fy at last lon, southernmost lat, level k
      sum2 = fy(IM,J2+1,nl) ! fy at last lon, northernmost lat, level k

      do i=1,IM-1 ! For all other lons
        sum1 = sum1 + fy(i,j1  ,nl) ! Add fy for last lon at southernmost lat
        sum2 = sum2 + fy(i,J2+1,nl) ! Add fy for last lon at northernmost lat
      enddo
 
      daero(1, 1,nl) = daero(1, 1,nl) - sum1*RCAP + fz(1, 1,nl) - fz(1, 1,nl+1) ! First lon, southernmost lat
      daero(1,JM,nl) = daero(1,JM,nl) + sum2*RCAP + fz(1,JM,nl) - fz(1,JM,nl+1) ! First lon, northernmost lat
 
      do i=2,IM ! All other lons except first
        daero(i, 1,nl) = daero(1, 1,nl) ! At southernmost lat, daero is equal to daero at first lon
        daero(i,JM,nl) = daero(1,JM,nl) ! At northernmost lat, daero is equal to daero at first lon
      enddo							

 
!****6***0*********0*********0*********0*********0*********0**********72
      if(FILL) call qckxyz(daero,DG2,IM,JM,NL,j1,j2,cosp,acosp,IC)
!****6***0*********0*********0*********0*********0*********0**********72
 
!MIC$ do all
!MIC$* shared(q,IC,NL,j1,im,jm,jm1,DQ,delp2)
!MIC$* private(i,j,k)

!******************************************************************
! finally, convert daero to mmr
!******************************************************************
! We have two options here:
! - to conserve the total mass of each tracer (cnst=2), use the 
!   surface pressure predicted by the dynamical core
! - to preserve constant mixing ration, use the surface pressure
!   calculated in this subroutine (i.e., the updated ps2)

      select case (cnst)
      case(1) ! Set to 1 here
        mmr(:,:,:,IC) = daero(:,:,:) / delp2(:,:,:) ! Mixing ratio for this particle is (updated mmr*pressure thickness divided by pressure thickness at time t+dt)
      case(2)
        mmr(:,:,:,IC) = daero(:,:,:) / delp2dyn(:,:,:)
      end select

! Finally, put in a sink term at the bottom level to avoid infinite build-up of haze particles	  
      mmr(:,:,nl,ic) = mmr(:,:,nl,ic)*10e-3
      where(mmr .lt. 0.) mmr = 0.0
      
      call mmr2n(mmr(:,:,:,ic),apart,rhop,rhog,im,jm,nl,numrhos(:,:,:,ic))

      enddo !tracer loop

      if (debug) then
         write(nud,*) '* Leaving routine aerocore'
         do ic=1,naero
            write(nud,*) '* tracer',ic,'* max. & min. mmr =', &
                       maxval(mmr(:,:,:,ic)),minval(mmr(:,:,:,ic))
         end do
      endif

      RETURN
      END ! End of aerocore subroutine
 

 
! The subroutines below calculate the terminal velocity of aerosol particles and the flux due to gravitational settling.
! Sources:
!
! Rosner, D. (2000) Transport Processes in Chemically Reacting Flow Systems
! Parmentier, V. et al. (2013) 3D mixing in hot Jupiter atmospheres. I: Application to the day/night cold trap in HD 209458b
! Steinrueck, M. et al. (2021) 3D simulations of photochemical hazes in the atmosphere of hot Jupiter HD 189733b
!

 !****6***0*********0*********0*********0*********0*********0**********72
    SUBROUTINE viscos(temp,im,jm,nl,nud,    &
                mu)
!   Calculate the viscosity of the bulk gas (molecular nitrogen, N2)
!****6***0*********0*********0*********0*********0*********0**********72
    USE pumamod, ONLY: PI
    USE aeromod, only: l_bulk
    
    IMPLICIT NONE
    
! Dimensions of arrays

    INTEGER,INTENT(IN) :: im,jm,nl ! Length of x, y, and z dimensions
    INTEGER,INTENT(IN) :: nud
    
! Define arrays

    REAL,INTENT(IN) :: temp(im,jm,nl) ! Air temperature
    REAL,INTENT(OUT) :: mu(im,jm,nl) ! Viscosity
    
! Define constants

    REAL :: kb ! Boltzmann constant
    REAL :: mair ! Molecular mass of bulk gas
    REAL :: dair ! Molecular diameter of bulk gas
    REAL :: eps ! Lennard-Jones potential well of N2
    
! Local arrays

    REAL :: rt_temp(im,jm,nl), temp_eps(im,jm,nl)
    REAL ::	rt_mair, rt_pi, rt_kb, sq_dair, coeff
    
    kb = 1.3806E-23 ! m2 kg s-2 K-1
    select case (l_bulk)
    case(1)
     mair = 4.652E-26 ! kg
     dair = 3.64E-10 ! m
     eps = 95.5 ! K
    case(2)
     mair = 3.34E-27
     dair = 2.827E-10
     eps = 59.7
    end select
    
!---------------------------------------------------------------
! Calculation
    
    rt_temp = SQRT(temp)
    rt_mair = SQRT(mair)
    rt_pi = SQRT(PI)
    rt_kb = SQRT(kb)
    sq_dair = dair*dair
    temp_eps = (temp/eps)**(4/25)
    coeff = (5./16.)*(1./1.22)*(1./PI)*rt_pi*rt_mair*rt_kb/sq_dair
    
    mu = coeff*rt_temp*temp_eps
    
    RETURN
    END
    
!****6***0*********0*********0*********0*********0*********0**********72
    SUBROUTINE chamfac(temp,im,jm,nl,sigmah,ps,apart,nud,   &
                        beta)
    
!   Calculate the Cunningham factor used in terminal velocity calculation
!****6***0*********0*********0*********0*********0*********0**********72
    USE pumamod, ONLY: PI
    USE aeromod, only: l_bulk
    
    IMPLICIT NONE

! Dimensions of arrays

    INTEGER,INTENT(IN) :: im,jm,nl,nud ! Length of x, y, and z dimensions
    
! Define inputs and outputs

    REAL,INTENT(IN) :: temp(im,jm,nl) ! Air temperature
    REAL,INTENT(IN) :: sigmah(nl) ! Sigma at half-level
    REAL,INTENT(IN) :: ps(im,jm) ! Surface air pressure
    REAL,INTENT(OUT) :: beta(im,jm,nl) ! Cunningham factor
    
! Local arrays
    REAL :: lambda(im,jm,nl) ! Mean free path
    REAL :: airp(im,jm,nl) ! Air pressure 
    REAL :: kn(im,jm,nl) ! Knudsen number
    REAL :: pwer(im,jm,nl) ! Power in exponent
    INTEGER :: k ! For vertical loop
    REAL :: coeff
    
! Define constants

    REAL :: kb ! Boltzmann constant
    REAL :: dm ! Molecular diameter of bulk gas
    REAL :: apart ! Particle radius
    
    kb = 1.3806E-23 ! J/K
    select case(l_bulk) ! N2
    case(1)
     dm = 3.64E-10 ! m
    case(2) ! H2
     dm = 2.827E-10 ! m
    end select

!---------------------------------------------------------------
! Calculation of mean free path of the bulk gas

    DO k=1,nl
        airp(:,:,k) = ps(:,:)*sigmah(k) ! Air pressure at mid level in Pa
    ENDDO
        
    coeff = (kb/dm)*(1./(SQRT(2.)*PI*dm))
    lambda = coeff*(temp/airp)
    
! Calculation of Cunningham factor

    kn = lambda/apart ! Knudsen number
    pwer = -1.1/kn
    beta = 1. + kn*(1.256 + 0.4*EXP(pwer)) ! Cunningham factor	
    
    RETURN
    END


!****6***0*********0*********0*********0*********0*********0**********72
    SUBROUTINE density(im,jm,nl,temp,ps,sigmah,nud,   &
                        rhog)

!   Calculate the gas density in kg/m3 (input for terminal velocity calculation)
!****6***0*********0*********0*********0*********0*********0**********72

    USE aeromod, only: l_bulk
    
    IMPLICIT NONE
    
! Dimensions of arrays

    INTEGER,INTENT(IN) :: im,jm,nl,nud ! Length of x, y, and z dimensions
    
! Define arrays

    REAL,INTENT(IN) :: temp(im,jm,nl) ! Temperature
    REAL,INTENT(IN) :: sigmah(nl) ! Sigma
    REAL,INTENT(IN) :: ps(im,jm) ! Surface air pressure
    REAL,INTENT(OUT) :: rhog(im,jm,nl) ! Density of surrounding gas

! Local arrays	
    REAL :: airp(im,jm,nl) ! Air pressure 
    INTEGER :: k ! For vertical loop
    
! Constants

    REAL :: R_gas ! Gas constant for ideal gas law
    REAL :: M_mass ! Molar mass of nitrogen
        
    R_gas = 8.314 ! SI units : J/K mol
    
    select case (l_bulk)
    case(1)
     M_mass = 0.0280 ! kg/mol for N2
    case(2)
     M_mass = 0.00201 ! kg/mol for H2
    end select
    
! Calculations
    
    DO k=1,nl
        airp(:,:,k) = ps(:,:)*sigmah(k) ! Air pressure at mid level in Pa
    ENDDO
        
    rhog = (airp*M_mass)/(R_gas*temp)
    
    RETURN
    END


!****6***0*********0*********0*********0*********0*********0**********72
    SUBROUTINE vterm(im,jm,nl,beta,apart,rhop,mu,   &
                    temp,sigmah,ps,rhog,nud,vels)
    
!   Calculate the terminal velocity of aerosol particles in the vertical direction
!****6***0*********0*********0*********0*********0*********0**********72
    USE pumamod, ONLY: ga ! Import planet's gravity from pumamod (plasimmod.f90)
    
    IMPLICIT NONE

! Dimensions of arrays
    
    INTEGER,INTENT(IN) :: im,jm,nl,nud ! Length of x, y, and z dimensions
    
! Define arrays

    REAL,INTENT(INOUT) :: beta(im,jm,nl) ! Cunningham factor
    REAL,INTENT(INOUT) :: mu(im,jm,nl) ! Viscosity
    REAL,INTENT(IN) :: temp(im,jm,nl) ! Temperature
    REAL,INTENT(IN) :: sigmah(nl) ! Sigma
    REAL,INTENT(IN) :: ps(im,jm,nl) ! Surface air pressure
    REAL,INTENT(IN) :: rhog(im,jm,nl) ! Gas density
    REAL,INTENT(OUT) :: vels(im,jm,nl) ! Terminal velocity

! Define constants

    REAL :: apart ! Aerosol particle radius
    REAL :: rhop ! Density of aerosol particle
    
    where(mu .le. 0.) mu = 1E-07
    
    vels = 2*beta*(apart**2)*ga*(rhop - rhog)/(9*mu)
    
    RETURN
    END
    
    
!****6***0*********0*********0*********0*********0*********0**********72
    SUBROUTINE gsettle(mmr,im,jm,nl,rhog,vterm,nud,gz)
    
!   Calculate vertical flux of haze particles due to gravitational settling at each timestep
!****6***0*********0*********0*********0*********0*********0**********72

    IMPLICIT NONE

! Dimensions of arrays
    
    INTEGER,INTENT(IN) :: im,jm,nl,nud ! Length of x and y dimensions
    
! Define arrays
    
    REAL,INTENT(IN) :: rhog(im,jm,nl)
    REAL,INTENT(IN) :: vterm(im,jm,nl)
    REAL,INTENT(INOUT) :: mmr(im,jm,nl)
    REAL,INTENT(OUT)   :: gz(im,jm,nl)
    
    gz = rhog*mmr*vterm
    
    where(gz .lt. 0.) gz = 0.0
    
    
    RETURN
    END
    
!****6***0*********0*********0*********0*********0*********0**********72
    SUBROUTINE mmr2n(mmr,apart,rhop,rhog,im,jm,nl,   &
                      numrho)
! Convert mass mixing ratio (kg/kg) from aerocore into number density (particles/m3)
!****6***0*********0*********0*********0*********0*********0**********72
    USE pumamod, ONLY: PI

    IMPLICIT NONE
    
    integer,intent(in) :: im,jm,nl
    real,intent(in) ::   mmr(im,jm,nl) ! Mixing ratio of aerosol
    real,intent(in)    :: apart ! Particle radius
    real,intent(in)    :: rhop ! Particle density
    real,intent(in)    :: rhog(im,jm,nl) ! Bulk gas density
    real,intent(out) :: numrho(im,jm,nl) ! Number density of aerosol
    
    REAL :: svol
    REAL :: mpart
    
    svol = (4/3)*PI*(apart**3) ! Sphere volume
    mpart = svol*rhop
    numrho = mmr*(1/mpart)*rhog
    RETURN
    END
    

!****6***0*********0*********0*********0*********0*********0**********72
    SUBROUTINE n2mmr(numrho,apart,rhop,rhog,im,jm,nl,   &
                      mmr)
! Convert number density (particles/m3) into mass mixing ratio (kg/kg)
! Not currently in use but just in case...
!****6***0*********0*********0*********0*********0*********0**********72
    USE pumamod, ONLY: PI

    IMPLICIT NONE
    
    integer,intent(in) :: im,jm,nl
    real,intent(in) ::   numrho(im,jm,nl) ! Mixing ratio of aerosol
    real,intent(in)    :: apart ! Particle radius
    real,intent(in)    :: rhop ! Particle density
    real,intent(in)    :: rhog(im,jm,nl) ! Bulk gas density
    real,intent(out) :: mmr(im,jm,nl) ! Number density of aerosol
    
    REAL :: svol
    REAL :: mpart
    
    svol = (4/3)*PI*(apart**3) ! Sphere volume
    mpart = svol*rhop ! Mass of one particle
    
    mmr = numrho*mpart/rhog
    
    RETURN
    END