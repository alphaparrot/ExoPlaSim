!=======================================================================
! Subroutines in this file have been adapted or rewritten 
! for the PlanetSimulator of University Hamburg 
! by Hui Wan (Max Planck Institute for Meteorology, 2009).
!=======================================================================
!
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine tpcore(q,ps1,ps2,u,v,                &
                        dtoa,dtdx,                    &
                        iord,jord,kord,               &
                        nc,im,jm,nl,dap,dbk,          &
                        iml,j1,j2,js0,jn0,            &
                        cose,cosp,acosp,dlat,rcap,    &
                        cnst,deform,zcross,           &
                        fill,mfct,debug,nud)
!****6***0*********0*********0*********0*********0*********0**********72
!
! Purpose: perform the transport of  3-D mixing ratio fields using 
!          externally specified winds on the hybrid Eta-coordinate.
!          One call to tpcore updates the 3-D mixing ratio
!          fields for one time step. [vertical mass flux is computed
!          internally using a center differenced hydrostatic mass
!          continuity equation].
!
! Schemes: Multi-dimensional Flux Form Semi-Lagrangian (FFSL) schemes
!          (Lin and Rood 1996, MWR) with a modified MFCT option
!          (Zalesak 1979).
!
! Multitasking version: 6m
! History: Mriginal vertion by S.-J. Lin, Oct. 1, 1997
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
! Q(IM,JM,NL,NC): mixing ratios at current time (t)
! NC: total # of constituents
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
! Q: the updated mixing ratios at t+NDT (original values are over-written)
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
      implicit none

! Input-Output variables

      integer,intent(in) :: nc,im,jm,nl,iml,j1,j2,js0,jn0
      integer,intent(in) :: iord,jord,kord,cnst

! Input-Output arrays

      real,intent(inout) ::   q(im,jm,nl,nc)
      real,intent(inout) :: ps1(im,jm)
      real,intent(inout) :: ps2(im,jm)
      real,intent(inout) ::   u(im,jm,nl)
      real,intent(inout) ::   v(im,jm,nl)
      real,intent(in)    ::  dtdx(jm), dtoa

      real,intent(in)    :: cose(jm+1),cosp(jm), acosp(jm), dlat(jm), rcap
      real,intent(in)    :: dap(nl), dbk(nl)
      logical,intent(in) :: zcross, fill, mfct, deform, debug

! Local dynamic arrays

      integer :: js(nl),jn(nl)

      real :: dtdx5(jm) 
      real ::   w(im,jm,nl)
      real :: crx(im,jm,nl),cry(im,jm,nl)
      real :: delp(im,jm,nl),delp1(im,jm,nl),delp2(im,jm,nl)
      real :: delp2dyn(im,jm,nl)
      real :: xmass(im,jm,nl),ymass(im,jm,nl)
      real ::   dg1(im),dg2(im,jm),dpi(im,jm,nl)
      real ::  qlow(im,jm,nl), dq(im,jm,nl)
      real ::    qz(im,jm,nl),qmax(im,jm,nl),qmin(im,jm,nl)
      real ::    wk(im,jm,nl),pu(im,jm,nl)
      real ::    fx(im+1,jm,nl),fy(im,jm,nl),fz(im,jm,nl+1)

! scalars

      integer :: imh, i,j,k,jt,ic,nud
      real    :: d5, dtoa5, sum1, sum2 

!---------------------------------------------------------------

      imh      = im/2
      dtdx5(:) = 0.5*dtdx(:)
      dtoa5    = 0.5*dtoa

      if (debug) then
         write(nud,*) '* Entered routine tpcore'
         do ic=1,nc
            write(nud,*) '* tracer',ic,'* max. & min. q =', &
                       maxval(q(:,:,:,ic)),minval(q(:,:,:,ic))
         end do
      endif

! save the vertical layer thickness (at t+dt) predicted by 
! the dynamical core 

      do k=1,nl
         delp2dyn(:,:,k) = dap(k) + dbk(k)*ps2(:,:)
      end do

! estimate the value at t+0.5dt

      ps2(:,:) = 0.5*( ps1(:,:) + ps2(:,:) )

!****6***0*********0*********0*********0*********0*********0**********72
! Compute Courant number
!****6***0*********0*********0*********0*********0*********0**********72
! Convert winds on A-Grid to Courant # on C-Grid.

      do j=2,jm-1
      do i=2,im
         crx(i,j,:) = dtdx5(j)*(u(i,j,:)+u(i-1,j,:))
      end do
         crx(1,j,:) = dtdx5(j)*(u(1,j,:)+u(im,j,:))
      end do
 
      do j=2,jm
         cry(:,j,:) = dtoa5*(v(:,j,:)+v(:,j-1,:))
      enddo
 
!****6***0*********0*********0*********0*********0*********0**********72
! Find JN and JS
!****6***0*********0*********0*********0*********0*********0**********72
! [2,js(k)] and [jn(k),jm-1] are the latitudes at which semi-Lagrangian
! treatment is needed.
 
!MIC$ do all autoscope shared(JS,JN,CRX,CRY,PS2,U,V,DPI,ymass,delp2,PU)
!MIC$* shared(xmass)
!MIC$* private(i,j,k,sum1,sum2,D5)

      do k=1,nl

      js(k) = j1
      jn(k) = j2
 
      outer1: do j=js0,j1+1,-1
      do i=1,im
      if(abs(crx(i,j,k)) .gt. 1.) then
            js(k) = j
            exit outer1
      endif
      enddo
      enddo outer1
 
      outer2: do j=jn0,j2-1
      do i=1,im
      if(abs(crx(i,j,k)) .gt. 1.) then
            jn(k) = j
            exit outer2
      endif
      enddo
      enddo outer2

      enddo  !vertical layer loop
 
!****6***0*********0*********0*********0*********0*********0**********72
! ***** Compute horizontal mass fluxes *****
!****6***0*********0*********0*********0*********0*********0**********72
! for the polar caps: replace the values at the northest (southest)
! latitude with the zonal mean. 

      sum1 = sum(ps1(:,1 ))/im
      sum2 = sum(ps1(:,jm))/im
      ps1(:,1)  = sum1
      ps1(:,jm) = sum2

      sum1 = sum(ps2(:,1 ))/im
      sum2 = sum(ps2(:,jm))/im
      ps2(:,1)  = sum1
      ps2(:,jm) = sum2


      do k=1,nl

! delp = pressure thickness: the psudo-density in a hydrostatic system.

      delp2(:,:,k) = dap(k) + dbk(k)*ps2(:,:)
 
! N-S componenet
 
      do j=j1,j2+1
         d5 = 0.5 * cose(j)
         ymass(:,j,k) = cry(:,j,k)*d5*(delp2(:,j,k)+delp2(:,j-1,k))
      enddo
 
      do j=j1,j2
         dpi(:,j,k) = (ymass(:,j,k)-ymass(:,j+1,k)) *acosp(j)/dlat(j)
      end do
 
! polar caps

      sum1 = -sum(ymass(:,j1  ,k))*rcap
      sum2 =  sum(ymass(:,j2+1,k))*rcap

      dpi(:, 1,k) = sum1
      dpi(:,jm,k) = sum2
 
! E-W component

      do j=j1,j2
      do i=2,im
         pu(i,j,k) = 0.5 * (delp2(i,j,k) + delp2(i-1,j,k))
      enddo
      enddo
 
      do j=j1,j2
         pu(1,j,k) = 0.5 * (delp2(1,j,k) + delp2(im,j,k))
      enddo
 
      do j=j1,j2
      do i=1,im
         xmass(i,j,k) = pu(i,j,k)*crx(i,j,k)
      enddo
      enddo
 
      do j=j1,j2
      do i=1,im-1
         dpi(i,j,k) = dpi(i,j,k) + xmass(i,j,k) - xmass(i+1,j,k)
      enddo
      enddo
 
      do j=j1,j2
         dpi(im,j,k) = dpi(im,j,k) + xmass(im,j,k) - xmass(1,j,k)
      enddo

      enddo ! vertical layer loop 

!****6***0*********0*********0*********0*********0*********0**********72
! Compute Courant number at cell center (upwinding)
!****6***0*********0*********0*********0*********0*********0**********72

! E-W direction

      do k=1,nl 

      do j=j1,j2
      do i=1,im-1
         if(crx(i,j,k)*crx(i+1,j,k) .gt. 0.) then
            if(crx(i,j,k) .gt. 0.) then
            u(i,j,k) = crx(i,j,k)
            else
            u(i,j,k) = crx(i+1,j,k)
            endif
         else
            u(i,j,k) = 0.
         endif
      enddo
      enddo
 
      i=im
      do j=j1,j2
      if(crx(i,j,k)*crx(1,j,k) .gt. 0.) then
         if(crx(i,j,k) .gt. 0.) then
         u(i,j,k) = crx(i,j,k)
         else
         u(i,j,k) = crx(1,j,k)
         endif
      else
         u(i,j,k) = 0.
      endif
      enddo

      enddo ! vertical layer loop 

! N-S direction

      do k=1,nl 

      do j=j1,j2
      do i=1,im
         if(cry(i,j,k)*cry(i,j+1,k) .gt. 0.) then
            if(cry(i,j,k) .gt. 0.) then
            v(i,j,k) = cry(i,j,k)/dlat(j)
            else
            v(i,j,k) = cry(i,j+1,k)/dlat(j)
            endif
         else
            v(i,j,k) = 0.
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

      do j=1,jm

         dg1(1:im) = 0.
         do k=1,nl
            dg1(1:im) = dg1(1:im) + dpi(1:im,j,k)
         enddo
 
! Compute PS2 (PS at n+1) using the hydrostatic assumption.
! Changes (increases) to surface pressure = total column mass convergence
 
         ps2(1:im,j)  = ps1(1:im,j) + dg1(1:im)
 
! compute vertical mass flux from mass conservation principle:
! lower boundary of the first (uppermost) layer 

         w(1:im,j,1) = dpi(1:im,j,1) - dbk(1)*dg1(1:im)
 
         do k=2,nl-1
         w(1:im,j,k) = w(1:im,j,k-1) + dpi(1:im,j,k) - dbk(k)*dg1(1:im)
         enddo

! Earth's surface 

         w(1:im,j,nl) = 0.

      enddo 
 
!MIC$ do all
!MIC$* shared(deform,NL,im,jm,delp,delp1,delp2,DPI,DAP,DBK,PS1,PS2)
!MIC$* private(i,j,k)

      do k=1,nl

         delp1(:,:,k) = dap(k) + dbk(k)*ps1(:,:)
         delp2(:,:,k) = dap(k) + dbk(k)*ps2(:,:)
         delp (:,:,k) = delp1(:,:,k) + dpi(:,:,k)
 
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

!****6***0*********0*********0*********0*********0*********0**********72
! Do transport one tracer at a time.
!****6***0*********0*********0*********0*********0*********0**********72
 
      DO ic=1,nc
 
!MIC$ do all autoscope
!MIC$* shared(q,DQ,delp1,U,V,j1,j2,JS,JN,im,jm,IML,IC,IORD,JORD)
!MIC$* shared(CRX,CRY,PU,xmass,ymass,fx,fy,acosp,rcap,qz)
!MIC$* private(i,j,k,jt,wk,DG2)

      do k=1,nl

! for the polar caps: replace the mixing ratio at the northest (southest)
! latitude with the zonal mean. 

         sum1 = sum(q(:,1 ,k,ic))/im
         sum2 = sum(q(:,jm,k,ic))/im
         q(:,1 ,k,ic) = sum1
         q(:,jm,k,ic) = sum2
      
! Initialize DQ
 
         dq(:,:,k) = q(:,:,k,ic)*delp1(:,:,k)

! E-W advective cross term
      call xadv(im,jm,j1,j2,q(1,1,k,ic),u(1,1,k),js(k),jn(k),iml, &
                wk(1,1,1))

      wk(:,:,1) = q(:,:,k,ic) + 0.5*wk(:,:,1)
 
! N-S advective cross term
      do j=j1,j2
      do i=1,im
      jt = float(j) - v(i,j,k)
      wk(i,j,2) = v(i,j,k) * (q(i,jt,k,ic) - q(i,jt+1,k,ic))
      enddo
      enddo
 
      do j=j1,j2
      do i=1,im
      wk(i,j,2) = q(i,j,k,ic) + 0.5*wk(i,j,2)
      enddo
      enddo

!****6***0*********0*********0*********0*********0*********0**********72
! compute flux in  E-W direction
      call xtp(im,jm,iml,j1,j2,jn(k),js(k),pu(1,1,k),dq(1,1,k), &
               wk(1,1,2),crx(1,1,k),fx(1,1,k),xmass(1,1,k),iord)

! compute flux in  N-S direction
      call ytp(im,jm,j1,j2,acosp,dlat,rcap,dq(1,1,k),wk(1,1,1), &
               cry(1,1,k),dg2,ymass(1,1,k),wk(1,1,3),wk(1,1,4), &
               wk(1,1,5),wk(1,1,6),fy(1,1,k),jord,nud)
!****6***0*********0*********0*********0*********0*********0**********72

      if(ZCROSS) then

! qz is the horizontal advection modified value for input to the
! vertical transport operator FZPPM
! Note: DQ contains only first order upwind contribution.

      do j=1,JM
      do i=1,IM
      qz(i,j,k) = DQ(i,j,k) / delp(i,j,k)
      enddo
      enddo

      else

      do j=1,JM
      do i=1,IM
      qz(i,j,k) = q(i,j,k,IC)
      enddo
      enddo

      endif
 
      enddo ! k-loop
 
!****6***0*********0*********0*********0*********0*********0**********72
! Compute fluxes in the vertical direction
      call FZPPM(qz,fz,IM,JM,NL,DQ,W,delp,KORD)
!****6***0*********0*********0*********0*********0*********0**********72
 
      if( MFCT ) then

      if (debug) write(nud,*) '* mfct on' 
       
! qlow is the low order "monotonic" solution
 
!MIC$ do all
!MIC$* shared(NL,im,jm,j1,jm1,qlow,DQ,delp2)
!MIC$* private(i,j,k)

      DO k=1,NL

      DO j=1,JM
      DO i=1,IM
      qlow(i,j,k) = DQ(i,j,k) / delp2(i,j,k)
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
      call FCT3D(Q(1,1,1,IC),qlow,fx,fy,fz,IM,JM,NL,j1,j2,delp2, &
                 DPI,qz,wk,Qmax,Qmin,DG2,U,V,acosp,dlat,RCAP)
! Note: Q is destroyed!!!
!****6***0*********0*********0*********0*********0*********0**********72
      ENDIF
 
! Final update

!MIC$ do all autoscope
!MIC$* private(i,j,k,sum1,sum2)

      do k=1,NL

      do j=j1,j2
      do i=1,IM
      DQ(i,j,k) = DQ(i,j,k) +  fx(i,j,k) - fx(i+1,j,k)                   &
                            + (fy(i,j,k) - fy(i,j+1,k))*acosp(j)/dlat(j) &
                            +  fz(i,j,k) - fz(i,j,k+1)
      enddo
      enddo
! poles:
      sum1 = fy(IM,j1  ,k)
      sum2 = fy(IM,J2+1,k)

      do i=1,IM-1
      sum1 = sum1 + fy(i,j1  ,k)
      sum2 = sum2 + fy(i,J2+1,k)
      enddo
 
      DQ(1, 1,k) = DQ(1, 1,k) - sum1*RCAP + fz(1, 1,k) - fz(1, 1,k+1)
      DQ(1,JM,k) = DQ(1,JM,k) + sum2*RCAP + fz(1,JM,k) - fz(1,JM,k+1)
 
      do i=2,IM
      DQ(i, 1,k) = DQ(1, 1,k)
      DQ(i,JM,k) = DQ(1,JM,k)
      enddo

      enddo
 
!****6***0*********0*********0*********0*********0*********0**********72
      if(FILL) call qckxyz(DQ,DG2,IM,JM,NL,j1,j2,cosp,acosp,IC)
!****6***0*********0*********0*********0*********0*********0**********72
 
!MIC$ do all
!MIC$* shared(q,IC,NL,j1,im,jm,jm1,DQ,delp2)
!MIC$* private(i,j,k)

!******************************************************************
! finally, convert dq to q
!******************************************************************
! We have two options here:
! - to conserve the total mass of each tracer (cnst=2), use the 
!   surface pressure predicted by the dynamical core
! - to preserve constant mixing ration, use the surface pressure
!   calculated in this subroutine (i.e., the updated ps2)

      select case (cnst)
      case(1)
        q(:,:,:,IC) = DQ(:,:,:) / delp2(:,:,:)
      case(2)
        q(:,:,:,IC) = DQ(:,:,:) / delp2dyn(:,:,:)
      end select

      enddo !tracer loop

      if (debug) then
         write(nud,*) '* Leaving routine tpcore'
         do ic=1,nc
            write(nud,*) '* tracer',ic,'* max. & min. q =', &
                       maxval(q(:,:,:,ic)),minval(q(:,:,:,ic))
         end do
      endif

      RETURN
      END
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine FCT3D(P,plow,fx,fy,fz,im,jm,km,j1,j2,delp,adx,ady, &
                       wk1,Qmax,Qmin,wkx,CRX,CRY,acosp,dlat,RCAP)
!****6***0*********0*********0*********0*********0*********0**********72
 
! MFCT Limiter
! plow: low order solution matrix
! P: current solution matrix
 
      PARAMETER (esl = 1.E-30)
      REAL P(IM,JM,km),CRX(IM,JM,km),CRY(IM,JM,km),plow(IM,JM,km), &
           Qmax(IM,JM,km),Qmin(IM,JM,km),acosp(*),delp(im,jm,km), &
           adx(IM,JM,km),ady(IM,JM,km),fx(IM+1,JM,km), &
           fy(IM,JM,km),fz(im,jm,km+1),wk1(IM,JM,km), &
           wkx(im,jm),wkn(im,jm),dlat(jm)
 
      JM1 = JM-1
 
! Find local min/max of the low-order monotone solution
      call hilo3D(P,im,jm,km,j1,j2,adx,ady,Qmax,Qmin,wkx,wkn)
      call hilo3D(plow,im,jm,km,j1,j2,Qmax,Qmin,wk1,P,wkx,wkn)
! P is destroyed!
 
!     GOTO 123

!MIC$ do all autoscope
!MIC$* shared(im,j1,j2,km,CRX,CRY,adx,ady,Qmax,Qmin)
!MIC$* private(i,j,k,IT,JT,PS1,PS2,PN1,PN2)

      DO k=1,km
      do j=j1,j2
      DO i=1,IM
 
      IT = NINT( float(i) - CRX(i,j,k) )
! Wrap around in E-W
      if(IT .lt. 1) then
            IT = IM + IT
      elseif(IT .GT. IM) then
            IT = IT - IM
      endif
 
      JT = NINT( float(j) - CRY(i,j,k) )
      Qmax(i,j,k) = max(Qmax(i,j,k), adx(IT,JT,k))
      Qmin(i,j,k) = min(Qmin(i,j,k), ady(IT,JT,k))
      enddo
      enddo
 
! Poles:
      PS1 = max(Qmax(1, 1,k), adx(1, 1,k))
      PS2 = min(Qmin(1, 1,k), ady(1, 1,k))
 
      PN1 = max(Qmax(1,JM,k), adx(1,JM,k))
      PN2 = min(Qmin(1,JM,k), ady(1,JM,k))
      DO i=1,IM
      Qmax(i, 1,k) = PS1
      Qmin(i, 1,k) = PS2
 
      Qmax(i,JM,k) = PN1
      Qmin(i,JM,k) = PN2
      enddo
      enddo
 
! 123   continue
! Flux Limiter
 
!MIC$ do all autoscope
!MIC$* shared(adx,ady,fx,fy,fz,plow,Qmax,Qmin,delp)
!MIC$* private(wkx,wkn)
!MIC$* private(i,j,k,ain,aou,bin,bou,cin,cou,btop,bdon)

      DO k=1,km
 
      DO j=j1,j2
      DO i=1,IM
         if(fx(i,j,k) .gt. 0.) then
         Ain = fx(i,j,k)
         Aou = 0.
         else
         Ain = 0.
         Aou = -fx(i,j,k)
         endif
 
         if(fx(i+1,j,k) .gt. 0.) then
         Aou = Aou + fx(i+1,j,k)
         else
         Ain = Ain - fx(i+1,j,k)
         endif
 
         if(fy(i,j,k) .gt. 0.) then
         Bin = fy(i,j,k)
         Bou = 0.
         else
         Bin = 0.
         Bou = -fy(i,j,k)
         endif
 
         if(fy(i,j+1,k) .gt. 0.) then
         Bou = Bou + fy(i,j+1,k)
         else
         Bin = Bin - fy(i,j+1,k)
         endif
 
         if(fz(i,j,k) .gt. 0.) then
         Cin = fz(i,j,k)
         Cou = 0.
         else
         Cin = 0.
         Cou = -fz(i,j,k)
         endif
 
         if(fz(i,j,k+1) .gt. 0.) then
         Cou = Cou + fz(i,j,k+1)
         else
         Cin = Cin - fz(i,j,k+1)
         endif

!HW 2009 - 
!****6***0*********0*********0*********0*********0*********0**********72
         wkx(i,j) = Ain + Bin*acosp(j)/dlat(j) + Cin
         wkn(i,j) = Aou + Bou*acosp(j)/dlat(j) + Cou
!****6***0*********0*********0*********0*********0*********0**********72
!HW 2009 - 
      enddo
      enddo
 
      DO j=j1,j2
      DO i=1,IM
      adx(i,j,k) = delp(i,j,k)*(Qmax(i,j,k)-plow(i,j,k))/(wkx(i,j)+esl)
      ady(i,j,k) = delp(i,j,k)*(plow(i,j,k)-Qmin(i,j,k))/(wkn(i,j)+esl)
      enddo
      enddo
 
! S Pole
      Ain = 0.
      Aou = 0.
      DO i=1,IM
      if(fy(i,j1,k).gt. 0.) then
           Aou = Aou + fy(i,j1,k)
      else
           Ain = Ain + fy(i,j1,k)
      endif
      enddo
      Ain = -Ain * RCAP
      Aou =  Aou * RCAP
 
! add vertical contribution...
 
      i=1
      j=1
      if(fz(i,j,k) .gt. 0.) then
      Cin = fz(i,j,k)
      Cou = 0.
      else
      Cin = 0.
      Cou = -fz(i,j,k)
      endif
 
      if(fz(i,j,k+1) .gt. 0.) then
      Cou = Cou + fz(i,j,k+1)
      else
      Cin = Cin - fz(i,j,k+1)
      endif
 
!****6***0*********0*********0*********0*********0*********0**********72
      btop = delp(1,1,k)*(Qmax(1,1,k)-plow(1,1,k))/(Ain+Cin+esl)
      bdon = delp(1,1,k)*(plow(1,1,k)-Qmin(1,1,k))/(Aou+Cou+esl)
!****6***0*********0*********0*********0*********0*********0**********72
 
      DO i=1,IM
      adx(i,j,k) = btop
      ady(i,j,k) = bdon
      enddo
! N Pole
      J=JM
      Ain = 0.
      Aou = 0.
      DO i=1,IM
      if(fy(i,j2+1,k).gt. 0.) then
           Ain = Ain + fy(i,j2+1,k)
      else
           Aou = Aou + fy(i,j2+1,k)
      endif
      enddo
      Ain =  Ain * RCAP
      Aou = -Aou * RCAP
 
! add vertical contribution...
 
      i=1
      if(fz(i,j,k) .gt. 0.) then
      Cin = fz(i,j,k)
      Cou = 0.
      else
      Cin = 0.
      Cou = -fz(i,j,k)
      endif
 
      if(fz(i,j,k+1) .gt. 0.) then
      Cou = Cou + fz(i,j,k+1)
      else
      Cin = Cin - fz(i,j,k+1)
      endif
 
!****6***0*********0*********0*********0*********0*********0**********72
      btop = delp(1,j,k)*(Qmax(1,j,k)-plow(1,j,k))/(Ain+Cin+esl)
      bdon = delp(1,j,k)*(plow(1,j,k)-Qmin(1,j,k))/(Aou+Cou+esl)
!****6***0*********0*********0*********0*********0*********0**********72
 
      DO i=1,IM
      adx(i,j,k) = btop
      ady(i,j,k) = bdon
      enddo
 
      if(j1 .ne. 2) then
      DO i=1,IM
! SP
      adx(i,2,k) = adx(i,1,k)
      ady(i,2,k) = ady(i,1,k)
! NP
      adx(i,JM1,k) = adx(i,JM,k)
      ady(i,JM1,k) = ady(i,JM,k)
      enddo
      endif
      
      enddo
 
!MIC$ do all autoscope
!MIC$* shared(fz,adx,ady,im,jm,km)
!MIC$* private(i,j,k)

      DO k=1,km
      DO j=j1,j2
      do i=2,IM
      if(fx(i,j,k) .gt. 0.) then
      fx(i,j,k) = min(1.,ady(i-1,j,k),adx(i,j,k))*fx(i,j,k)
      else
      fx(i,j,k) = min(1.,adx(i-1,j,k),ady(i,j,k))*fx(i,j,k)
      endif
      enddo
      enddo
 
! For i=1
      DO j=j1,j2
      if(fx(1,j,k) .gt. 0.) then
      fx(1,j,k) = min(1.,ady(IM,j,k),adx(1,j,k))*fx(1,j,k)
      else
      fx(1,j,k) = min(1.,adx(IM,j,k),ady(1,j,k))*fx(1,j,k)
      endif
      fx(IM+1,j,k) = fx(1,j,k)
      enddo
 
      do j=j1,j2+1
      do i=1,IM
      if(fy(i,j,k) .gt. 0.) then
        fy(i,j,k) = min(1.,ady(i,j-1,k),adx(i,j,k))*fy(i,j,k)
      else
        fy(i,j,k) = min(1.,adx(i,j-1,k),ady(i,j,k))*fy(i,j,k)
      endif
      enddo
      enddo

      if(k .ne. 1) then
      do j=1,jm
      do i=1,im
      if(fz(i,j,k) .gt. 0.) then
        fz(i,j,k) = min(1.,ady(i,j,k-1),adx(i,j,k))*fz(i,j,k)
      else
        fz(i,j,k) = min(1.,adx(i,j,k-1),ady(i,j,k))*fz(i,j,k)
      endif
      enddo
      enddo
      endif

      enddo
 
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine ytp(IMR,JNP,j1,j2,acosp,dlat,RCAP,DQ,P,C,DC2 &
                    ,ymass,fy1,A6,AR,AL,fy2,JORD,nud)
!****6***0*********0*********0*********0*********0*********0**********72
!     do-loops re-written and "dlat" introduced by Hui Wan 
!     for PlanetSimulator (2009-05).

      implicit none

      integer imr,jnp,j1,j2,jord,nud
      REAL P(IMR,JNP),C(IMR,JNP),ymass(IMR,JNP),fy2(IMR,JNP), &
           DC2(IMR,JNP),DQ(IMR,JNP),acosp(JNP),dlat(JNP),rcap

! Work array
      REAL fy1(IMR,JNP),AR(IMR,JNP),AL(IMR,JNP),A6(IMR,JNP)
      real sum1,sum2
      integer i,j,jt,jmr
 
      JMR = JNP - 1
 
      if(JORD.eq.1) then  
!     upwind scheme

         do j=2,jnp
         do i=1,imr
            if ( c(i,j).le.0. ) then
               jt = j
            else
               jt = j-1
            endif
         fy1(i,j) = p(i,jt)
         enddo
         enddo

         do j=j1,j2+1
         do i=1,imr
         fy2(i,j) = 0.
         enddo
         enddo

      else

         call ymist(IMR,JNP,j1,P,DC2)
 
         if(JORD.LE.0 .or. JORD.GE.3) then
!        PPM
!!           call fyppm(C,P,DC2,fy1,fy2,IMR,JNP,j1,j2,A6,AR,AL,JORD)
             write(nud,*) 'PPM algorithm not yet adapted'
             stop

         else
!        van Leer
           do j=2,jnp
           do i=1,imr
              if ( c(i,j).le.0. ) then
                 jt = j
              else
                 jt = j-1
              endif
           fy1(i,j) = p(i,jt)
           fy2(i,j) = (sign(1.,C(i,j))-C(i,j)/dlat(jt))*DC2(i,jt)
           enddo
           enddo

         endif ! van Leer or PPM
      endif    ! upwind or higher order

      do j=2,jnp
      do i=1,imr
      fy1(i,j) = fy1(i,j)*ymass(i,j)
      fy2(i,j) = fy2(i,j)*ymass(i,j)
      enddo
      enddo
 
      DO j=j1,j2
      DO i=1,IMR
      DQ(i,j) = DQ(i,j) + (fy1(i,j) - fy1(i,j+1)) * acosp(j)/dlat(j)
      enddo
      enddo
 
! Poles
      sum1 = fy1(IMR,j1  )
      sum2 = fy1(IMR,J2+1)
      do i=1,IMR-1
      sum1 = sum1 + fy1(i,j1  )
      sum2 = sum2 + fy1(i,J2+1)
      enddo
 
      sum1 = DQ(1,  1) - sum1 * RCAP
      sum2 = DQ(1,JNP) + sum2 * RCAP
      do i=1,IMR
      DQ(i,  1) = sum1
      DQ(i,JNP) = sum2
      enddo
 
!     if(j1.ne.2) then
!     do i=1,IMR
!     DQ(i,  2) = sum1
!     DQ(i,JMR) = sum2
!     enddo
!     endif

      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine  ymist(IMR,JNP,j1,P,DC)
!****6***0*********0*********0*********0*********0*********0**********72
!     re-written by Hui Wan for the PlanetSimulator (2009-05)

      implicit none

!     input and output

      integer imr,jnp,j1
      real p(imr,jnp),dc(imr,jnp)

!     local variables

      integer imh,jmr,i,j
      real r24
      parameter ( r24 = 1./24. )
      real tmp,pmin,pmax
!
! 2nd order version for scalars
!
      imh = imr / 2
      jmr = jnp - 1

      do j=2,jmr 
      do i=1,imr
      tmp = 0.25*(p(i,j+1) - p(i,j-1))
      pmax = max(p(i,j-1),p(i,j),p(i,j+1)) - p(i,j)
      pmin = p(i,j) - min(p(i,j-1),p(i,j),p(i,j+1))
      dc(i,j) = sign(min(abs(tmp),pmin,pmax),tmp)
      enddo
      enddo
!
! Poles:
! Determine slopes in polar caps for scalars!
 
      do i=1,IMH
! South
      tmp = 0.25*(p(i,2) - p(i+imh,2))
      Pmax = max(p(i,2),p(i,1), p(i+imh,2)) - p(i,1)
      Pmin = p(i,1) - min(p(i,2),p(i,1), p(i+imh,2))
      DC(i,1)=sign(min(abs(tmp),Pmax,Pmin),tmp)
! North.
      tmp = 0.25*(p(i+imh,JMR) - p(i,JMR))
      Pmax = max(p(i+imh,JMR),p(i,jnp), p(i,JMR)) - p(i,JNP)
      Pmin = p(i,JNP) - min(p(i+imh,JMR),p(i,jnp), p(i,JMR))
      DC(i,JNP) = sign(min(abs(tmp),Pmax,pmin),tmp)
      enddo
 
! Scalars:
      do i=imh+1,IMR
      DC(i,  1) =  - DC(i-imh,  1)
      DC(i,JNP) =  - DC(i-imh,JNP)
      enddo

      return
      end
 
