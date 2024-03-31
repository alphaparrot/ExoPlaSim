! =======================
! VEGETATION MODULE SIMBA              ! by Axel Kleidon & Pablo Paiewonsky
! =======================

module vegmod
use landmod
implicit none

character(len=80) :: veg_version = 'SIMBA 18-Jun-2014 by Edi'

integer :: jhor

! parameter for soil component

real, parameter :: vws_crit  = 0.25             ! critical soil wetness

! parameter for vegetation component
! all values are in m-k-s units

real, parameter :: q10       = 2.0              ! q10 value for resp

!     * parameters for land surface parameters

real, parameter :: zlaimax   =  9.0               ! max LAI
real, parameter :: zlaimin   =  0.1               ! min LAI for bare soil
real, parameter :: valb_min  =  0.12              ! min albedo for veg
real, parameter :: valb_max  =  0.30              ! max albedo for bare soil
real, parameter :: vsalb_min =  0.20              ! snow albedo for veg
real, parameter :: vwmax_max =  0.5               ! max SWHC for veg
real, parameter :: vwmax_min =  0.05              ! SWHC for bare soil
real, parameter :: vz0_min   =  0.05              ! min roughness for bare soil
real, parameter :: cveg_a    =  2.9               ! forest cover conversion
real, parameter :: cveg_b    =  1.0               ! forest cover conversion
real, parameter :: cveg_c    =  1.570796326794897 ! forest cover conversion
real, parameter :: cveg_d    =  -1.9              ! new forest cover parameter
real, parameter :: cveg_f    =  3.0               ! new soil cover parameter
real, parameter :: cveg_k    =  0.5               ! conversion LAI
real, parameter :: cveg_l    =  0.25              ! LAI-biomass relationship
real, parameter :: ct_crit   =  5.0               ! for temp limit function

!* beta factor defined as GPP  = (1 + beta ln(C/C0))
!  beta  = 0.3-0.6, C0 = 360ppm - Harvey 1989

real, parameter :: co2_ref   = 360.     ! co2 reference concentration in ppm
real, parameter :: co2_comp  = 0.       ! light compensation point in ppm
real, parameter :: co2_sens  = 0.3      ! beta factor

real            :: cveg_e               ! forest cover conversion
real            :: cveg_g               ! soil cover conversion
real            :: zgpp                 ! gross primary production (kg C /m2 /s)
real            :: zgppl                ! light limited gpp (kg C /m2 /s)
real            :: zgppw                ! water limited gpp (kg C /m2 /s)
real            :: znpp                 ! net primary production (kg C /m2 /s)
real            :: zforest              ! forest cover
real            :: zlitter              ! litterfall (kg C /m2 /s)
real            :: znogrow              ! no growth allocation (kg C /m2 /s)
real            :: zvalb(NLIGHTS) = 0.12! albedo
real            :: zvrhs                ! surface conductance
real            :: zres                 ! heterotrophic respiration (kg C /m2 /s)
real            :: zveg                 ! vegetation cover (fract.)
real            :: zvz0                 ! roughness length
real            :: zwmax                ! bucket depth

! namelist variables

integer         :: ncveg     = 1        ! vegetation accelerator
real            :: forgrow   = 1.0      ! 
real            :: rinidagg  = 0.5      ! initial value for dagg(:)
real            :: rinidsc   = 1.0      ! initial value for dsc(:)
real            :: rinidmr   = 2.0      ! initial value for dmr(:)
real            :: rinisoil  = 0.0      ! initial value for dcsoil(:)
real            :: riniveg   = 0.0      ! initial value for dcveg(:)

! dynamically calculated land surface parameters

real, allocatable :: dagg(:)     ! above ground growth
real, allocatable :: dsc(:)      ! stomatal conductance
real, allocatable :: dmr(:)      ! maximum roughness
real, allocatable :: dcsoil(:)   ! organic carbon stored in soil
real, allocatable :: dcveg(:)    ! organic carbon stored in biomass
real, allocatable :: dgrow(:)    ! biomass growth parameter
real, allocatable :: dlai(:)     ! leaf area index (integer value)
real, allocatable :: dvsoil(:)   ! soil cover

! accumulated output variables

real, allocatable :: agpp(:)     ! gross primary production
real, allocatable :: agppl(:)    ! light limited GPP
real, allocatable :: agppw(:)    ! water limited GPP
real, allocatable :: alitter(:)  ! litter production
real, allocatable :: anogrow(:)  ! no growth allocation 
real, allocatable :: anpp(:)     ! net primary production
real, allocatable :: aresh(:)    ! heterotrophic respiration
real, allocatable :: adlai(:)    ! leaf area index
real, allocatable :: adcveg(:)   ! vegetative carbon
real, allocatable :: adcsoil(:)  ! soil carbon

! Non-accumulated output variables

real, allocatable :: dgpp(:)      ! gross primary production
real, allocatable :: dgppl(:)     ! light limited GPP
real, allocatable :: dgppw(:)     ! water limited GPP
real, allocatable :: dlitter(:)   ! litter production
real, allocatable :: dnogrow(:)   ! no growth allocation 
real, allocatable :: dnpp(:)      ! net primary production
real, allocatable :: dresh(:)     ! heterotrophic respiration

real :: zalbsn =   0.           ! snow albedo

!     * output of cumulative fluxes

integer :: ibiomass
integer :: klight
real    :: zlaim                 ! maximum lai
real    :: zvegm                 ! maximum veg cover
real    :: zft
real    :: zsvp                  ! saturation vapor pressure
real    :: zvpd                  ! vapor pressure deficit
real    :: zbeta
real    :: zglacfree             ! glacier free fraction

end module vegmod


!========!
! VEGINI !
!========!

subroutine vegini
use vegmod
implicit none

namelist/vegmod_nl/ncveg,forgrow,rinidagg,rinidsc,rinidmr,rinisoil &
                  ,riniveg

if (mypid == NROOT) then 
   open(12,file=vegmod_namelist)
   read(12,vegmod_nl)
   write(nud,'(/,"***********************************************")')
   write(nud,'("* VEGMOD: ",a35," *")') trim(veg_version)
   write(nud,'("***********************************************")')
   write(nud,'("* Namelist VEGMOD_NL from <",a18,"> *")') &
         vegmod_namelist
   write(nud,'("***********************************************")')
   write(nud,vegmod_nl)
   close(12)
endif

call mpbci(ncveg)
call mpbcr(forgrow)
call mpbcr(rinidagg)
call mpbcr(rinidsc)
call mpbcr(rinidmr)
call mpbcr(rinisoil)
call mpbcr(riniveg)

allocate(dagg(NHOR))    ; dagg(:)    = rinidagg ! above ground growth
allocate(dsc(NHOR))     ; dsc(:)     = rinidsc  ! stomatal conductance
allocate(dmr(NHOR))     ; dmr(:)     = rinidmr  ! maximum roughness
allocate(dcsoil(NHOR))  ; dcsoil(:)  = rinisoil ! organic carbon stored in soil
allocate(dcveg(NHOR))   ; dcveg(:)   = riniveg  ! organic carbon stored in biomass
allocate(dgrow(NHOR))   ; dgrow(:)   = forgrow  ! biomass growth parameter
allocate(dlai(NHOR))    ; dlai(:)    = 0.0      ! leaf area index
allocate(dvsoil(NHOR))  ; dvsoil(:)  = 0.0      ! soil cover
allocate(dgpp(NHOR))    ; dgpp(:)    = 0.0      ! gross primary production
allocate(dgppl(NHOR))   ; dgppl(:)   = 0.0      ! light limited GPP
allocate(dgppw(NHOR))   ; dgppw(:)   = 0.0      ! water limited GPP
allocate(dlitter(NHOR)) ; dlitter(:) = 0.0      ! litter production
allocate(dnogrow(NHOR)) ; dnogrow(:) = 0.0      ! no growth allocation 
allocate(dnpp(NHOR))    ; dnpp(:)    = 0.0      ! net primary production
allocate(dresh(NHOR))   ; dresh(:)   = 0.0      ! heterotrophic respiration
allocate(agpp(NHOR))    ; agpp(:)    = 0.0      ! gross primary production
allocate(agppl(NHOR))   ; agppl(:)   = 0.0      ! light limited GPP
allocate(agppw(NHOR))   ; agppw(:)   = 0.0      ! water limited GPP
allocate(alitter(NHOR)) ; alitter(:) = 0.0      ! litter production
allocate(anogrow(NHOR)) ; anogrow(:) = 0.0      ! no growth allocation 
allocate(anpp(NHOR))    ; anpp(:)    = 0.0      ! net primary production
allocate(aresh(NHOR))   ; aresh(:)   = 0.0      ! heterotrophic respiration
allocate(adlai(NHOR))   ; adlai(:)   = 0.0      ! leaf area index 
allocate(adcveg(NHOR))  ; adcveg(:)  = 0.0      ! vegetative carbon
allocate(adcsoil(NHOR)) ; adcsoil(:) = 0.0      ! soil carbon

call mpsurfgp('dagg'    ,dagg    ,NHOR,1)
call mpsurfgp('dsc'     ,dsc     ,NHOR,1)
call mpsurfgp('dmr'     ,dmr     ,NHOR,1)
call mpsurfgp('dgrow'   ,dgrow   ,NHOR,1)

if (nrestart == 0) then
   where (dls(:) == 0) dforest(:) = 0.0
else ! read variables from restart file
   call mpgetgp('agpp'   ,agpp   ,NHOR,1)
   call mpgetgp('agppl'  ,agppl  ,NHOR,1)
   call mpgetgp('agppw'  ,agppw  ,NHOR,1)
   call mpgetgp('alitter',alitter,NHOR,1)
   call mpgetgp('anogrow',anogrow,NHOR,1)
   call mpgetgp('anpp'   ,anpp   ,NHOR,1)
   call mpgetgp('aresh'  ,aresh  ,NHOR,1)
   call mpgetgp('dcsoil' ,dcsoil ,NHOR,1)
   call mpgetgp('dcveg'  ,dcveg  ,NHOR,1)
   call mpgetgp('adcsoil' ,adcsoil ,NHOR,1)
   call mpgetgp('adcveg'  ,adcveg  ,NHOR,1)
   call mpgetgp('adlai'   ,adlai   ,NHOR,1)
endif ! (nrestart == 0)
return
end subroutine vegini


!========!
! VEGOUT !
!========!

subroutine vegout
use vegmod
implicit none
real :: zfac

! We don't need to do the low-io switch here, because this only gets called when writegp gets called,
! and the default is to write time-averaged data. We will create a separate vegsnapout routine for 
! non-lowio, snapshot, high-cadence writes
zfac = 1.0 / naccuout

! scale accumulated output variables

agpp(:)    = agpp(:)    * zfac
agppl(:)   = agppl(:)   * zfac
agppw(:)   = agppw(:)   * zfac
alitter(:) = alitter(:) * zfac
anogrow(:) = anogrow(:) * zfac
anpp(:)    = anpp(:)    * zfac
aresh(:)   = aresh(:)   * zfac
adlai(:)   = adlai(:)   * zfac
adcveg(:)  = adcveg(:)  * zfac
adcsoil(:) = adcsoil(:) * zfac

! write output variables

call writegp(40,adlai  ,299,0) ! leaf area index

call writegp(40,agpp   ,300,0) ! gross primary production
call writegp(40,anpp   ,301,0) ! net primary production
call writegp(40,agppl  ,302,0) ! light limited GPP 
call writegp(40,agppw  ,303,0) ! water limited GPP
call writegp(40,adcveg ,304,0) ! vegetation carbon
call writegp(40,adcsoil,305,0) ! soil carbon
call writegp(40,anogrow,306,0) ! no growth allocation 
call writegp(40,aresh  ,307,0) ! heterotrophic respiration
call writegp(40,alitter,308,0) ! litter production

! reset accumulated output variables

agpp(:)    = 0.0
agppl(:)   = 0.0
agppw(:)   = 0.0
alitter(:) = 0.0
anogrow(:) = 0.0
anpp(:)    = 0.0
aresh(:)   = 0.0
adlai(:)   = 0.0
adcveg(:)  = 0.0
adcsoil(:) = 0.0

return
end subroutine vegout

!========!
! VEGSNAPOUT !
!========!

subroutine vegsnapout(kcode)
use vegmod
implicit none
integer kcode


! write output variables

call writegp(kcode,dlai  ,299,0) ! leaf area index

call writegp(kcode,dgpp   ,300,0) ! gross primary production
call writegp(kcode,dnpp   ,301,0) ! net primary production
call writegp(kcode,dgppl  ,302,0) ! light limited GPP 
call writegp(kcode,dgppw  ,303,0) ! water limited GPP
call writegp(kcode,dcveg ,304,0) ! vegetation carbon
call writegp(kcode,dcsoil,305,0) ! soil carbon
call writegp(kcode,dnogrow,306,0) ! no growth allocation 
call writegp(kcode,dresh  ,307,0) ! heterotrophic respiration
call writegp(kcode,dlitter,308,0) ! litter production


return
end subroutine vegsnapout

!=========!
! VEGSTOP !
!=========!

subroutine vegstop
use vegmod
implicit none

call mpputgp('agpp'   ,agpp   ,NHOR,1)
call mpputgp('agppl'  ,agppl  ,NHOR,1)
call mpputgp('agppw'  ,agppw  ,NHOR,1)
call mpputgp('alitter',alitter,NHOR,1)
call mpputgp('anogrow',anogrow,NHOR,1)
call mpputgp('anpp'   ,anpp   ,NHOR,1)
call mpputgp('aresh'  ,aresh  ,NHOR,1)
call mpputgp('dcsoil' ,dcsoil ,NHOR,1)
call mpputgp('dcveg'  ,dcveg  ,NHOR,1)
call mpputgp('adcsoil' ,adcsoil ,NHOR,1)
call mpputgp('adcveg'  ,adcveg  ,NHOR,1)
call mpputgp('adlai'   ,adlai   ,NHOR,1)

return
end subroutine vegstop


!=========!
! VEGSTEP !
!=========!

subroutine vegstep
use vegmod
implicit none

! Initialization and copy of puma fields needed

cveg_e = atan(cveg_d) ! forest cover
cveg_g = atan(-cveg_f)   ! soil cover
if (co2 .lt. 1.0e-8) then !No CO2 limit
  zbeta = 0.0
else
  zbeta  = max(0.,1.+co2_sens*log((co2-co2_comp)/(co2_ref-co2_comp)))
endif
  
! Following plasim arrays are used but not modified
! -------------------------------------------------
! dwatc(:) = soil wetness [m]
! dswfl(:) = short wave radiation [W/m2]
! devap(:) = surface evaporation (negative)

do jhor = 1 , NHOR
  if (dls(jhor) > 0.0 .and. dglac(jhor) < 0.90) then ! land cell with < 90 % glacier

    zwmax = dwmax(jhor)   ! Bucket depth
    zglacfree = 1.0 - dglac(jhor) ! glacier free fraction
    ! zsvp = saturation vapor pressure [Pa] (Magnus-Teten)
    ! zvpd = vapor pressure deficit    [Pa] (set to >= 0.01 Pa)

    zsvp = ra1*exp(ra2*(dt(jhor,NLEP)-TMELT)/(dt(jhor,NLEP)-ra4))
    zvpd = max(0.01, zsvp - dp(jhor)*dq(jhor,NLEV)/0.622)

    ! structurally limited leaf area index (eq. 6.6 PlaSim RM)
    ! zlaimax : set to 9.0 - local constant
    ! zlaimin : set to 0.1 - minimal LAI for bare soil

    zlaim = zlaimin + 0.63661977236758 * zlaimax * atan(cveg_l*2.0*dagg(jhor)*dcveg(jhor))
    ! where 0.6366... = 2/pi

    ! structurally limited vegetation cover (eq. 6.5 PlaSim RM)
    ! cveg_k : set to 0.5 (local constant)

    zvegm = 1.0 - exp(-cveg_k * zlaim)

    ! calculate water stress factor (eq. 6.3 PlaSim RM)
    ! water stress = soil wetness / (field capacity * critical soil wetness)

    zveg  = 0.0
    if (zwmax > 0.0) then
       zveg = min(1.0,zvegm,max(0.0,dwatc(jhor) / (zwmax * vws_crit)))
    endif

    ! calculate gross primary productivity

    zft = min(1.0,max(0.0,(dt(jhor,NLEP) - TMELT) / ct_crit))

    ! light limited gpp [kg C / m2 /s]

    zgppl = rlue * zbeta * zft * zveg * dfd(jhor,NLEP)

    ! water limited gpp [kg C / m2 /s]
    ! co2    : initialized with 360.0 [ppmv] - member of namelist $RADPAR
    ! co2conv:                               - member of namelist $LANDPAR

    if (co2 .lt. 1.0e-8) then !In limit of no CO2 there's no growth
      zgppw = 0.
    else
      zgppw = co2conv * dp(jhor) * zveg * (-devap(jhor)) / zvpd * (0.3 * co2)
    endif
    zgppw = min(1.0,max(0.0,zgppw)) ! limit to range (0..1)

    ! combine light & water limitations (eq. 6.1 PlaSim RM)

    zgpp = min(zgppl,zgppw)

    ! net primary production [kg C / m2 /s]
    ! dgrow : initialized from forgrow or read from surface file

    znpp = dgrow(jhor) * 0.5 * zgpp * zglacfree

    ! litterfall [kg C / m2 /s]
    ! tau_veg : initialized to 10 [years] - member of namelist $LANDPAR

    zlitter = dcveg(jhor) / tau_veg

    ! heterotrophic respiration [kg C / m2 /s]
    ! q10 : set to 2.0 (local constant)
    ! tau_soil : initialized to 42 [years] - member of namelist $LANDPAR

    zres = q10**((dt(jhor,NLEP)-TMELT-10.)/10.) * dcsoil(jhor)/tau_soil

    ! no growth allocation [kg C / m2 /s]

    znogrow  = (1.0 - dgrow(jhor)) * 0.5 * zgpp

    ! carbon stored in biomass (ncveg = time accelerator)

    dcveg(jhor) = dcveg(jhor) + (znpp - zlitter) * deltsec * ncveg  !! FL

    ! carbon stored in soil

    dcsoil(jhor) = dcsoil(jhor) + (zlitter - zres) * deltsec * ncveg !! FL

    ! forest cover
    ! dagg : above ground growth factor : initialized from riidagg

    zforest = (atan(dcveg(jhor) - cveg_a) - cveg_e) / (cveg_c - cveg_e)
    zforest = min(1.0,max(0.0,zforest))

    ! soil cover

!   dvsoil(jhor) = atan((2.0*(1.0-dagg(jhor))*dcveg(jhor)-cveg_a)/cveg_b)/cveg_c+cveg_d
!   dvsoil(jhor) = min(1.0,max(0.0,cveg_e * dvsoil(jhor)))
    dvsoil(jhor) = (atan(dcveg(jhor) - cveg_f) - cveg_g) / (cveg_c - cveg_g)
    dvsoil(jhor) = min(1.0,max(0.0,dvsoil(jhor)))

    ! discretization of vegetation state

    if (rnbiocats >= 2.0) then
      ibiomass      = zforest * rnbiocats
      zforest       = min(1.0, real(ibiomass)/(rnbiocats-1.0))
      ibiomass      = dvsoil(jhor) * rnbiocats
      dvsoil(jhor)  = min(1.0, real(ibiomass)/(rnbiocats-1.0))
    endif

    ! derivation of land surface parameters
    
    !We will need to adapt this for input spectra. Something to consider would be using the Qntenna code
    !by Trevor Arp and Nathaniel Gabor (Arp, et al 2020) to pre-compute the absorbance spectrum of 
    !vegetation for a given stellar spectrum. It's a toy model, but it's a pretty insightful one, and
    !the vegetation model here is sort of a toy model too.

    dlai(jhor)  = -log(1.0 - zveg)/cveg_k
    zvz0  = dmr(jhor) * zforest+vz0_min*(1.0-zforest)
    zvalb(:) = valb_min * zveg + valb_max * (1.0 - zveg)
    zwmax = vwmax_max * dvsoil(jhor) + vwmax_min * (1.0 - dvsoil(jhor))
    zvrhs = dsc(jhor) * min(1.0,max(0.0,dwatc(jhor)/(zwmax*vws_crit)))

    ! modify albedo and surface conductance due to snow

    if (dsnow(jhor) > 0.0) then

      ! copy PUMA-2 formulation of snow albedo and mix it with forested albedo
      do klight=1,NLIGHTS
        zalbsn = (albsmax(klight) - albsmin(klight)) * (dt(jhor,NLEP)-263.16) / (TMELT-263.16)
        zalbsn = max(albsmin(klight),min(albsmax(klight),albsmax(klight)-zalbsn))
        zvalb(klight)  = zvalb(klight)+(zalbsn-zvalb(klight))*dsnow(jhor)/(dsnow(jhor)+0.01)
        zvalb(klight)  = vsalb_min * zforest + zvalb(klight) * (1.-zforest)
      enddo
      zvrhs  = 1.0
    endif

    ! interactive coupling

    if (nveg == 2) then
       dz0(jhor)     = sqrt(zvz0*zvz0+dz0climo(jhor)*dz0climo(jhor))
       dwmax(jhor)   = zwmax
       drhs(jhor)    = zvrhs
       dalb(jhor,:)    = zvalb(:)
       dforest(jhor) = zforest
       !!! TODO: CO2 coupling--both local as surface source/sinks and global source/sinks
       ! CO2 flux per dt is area_integral((dresh - dnpp)*dt*mass_conversion)
    endif

  ! accumulate output variables

  agpp(jhor)    = agpp(jhor)    + zgpp
  agppl(jhor)   = agppl(jhor)   + zgppl
  agppw(jhor)   = agppw(jhor)   + zgppw
  alitter(jhor) = alitter(jhor) + zlitter
  anogrow(jhor) = anogrow(jhor) + znogrow
  anpp(jhor)    = anpp(jhor)    + znpp
  aresh(jhor)   = aresh(jhor)   + zres
  adlai(jhor)   = adlai(jhor)   + dlai(jhor)
  adcveg(jhor)  = adcveg(jhor)  + dcveg(jhor)
  adcsoil(jhor)  = adcsoil(jhor)  + dcsoil(jhor)
  
  ! non-accumulated output variables
  
  dgpp(jhor)    = zgpp
  dgppl(jhor)   = zgppl
  dgppw(jhor)   = zgppw
  dlitter(jhor) = zlitter
  dnogrow(jhor) = znogrow
  dnpp(jhor)    = znpp
  dresh(jhor)   = zres

  endif ! (dls(jhor) > 0.0 .and. dglac(jhor) < 0.9)
enddo ! jhor

! Send some arrays to GUI

if (ngui > 0) then
   call guihor("DCVEG"   // char(0),dcveg  ,1,1000.0,0.0)
   call guihor("DVSOIL"  // char(0),dvsoil ,1,1000.0,0.0)
   call guihor("DFOREST" // char(0),dforest,1,1000.0,0.0)
endif

return
end subroutine vegstep
