      module pumamod
      use resmod
!
!     Parameter
!
      parameter(NPRO = NPRO_ATM)        ! Number of processes (resmod)
      parameter(NLEV=NLEV_ATM)
      parameter(NLAT = NLAT_ATM)        ! Number of latitudes (resmod)
      parameter(NLON = NLAT + NLAT)     ! Number of longitudes
      parameter(NLPP = NLAT / NPRO)     ! Latitudes per process
      parameter(NHOR = NLON * NLPP)     ! Horizontal part
      parameter(NUGP = NLON * NLAT)     ! Number of gridpoints
      parameter(NTRU = (NLON-1) / 3)       ! Triangular truncation
      parameter(NTP1 = NTRU + 1)           ! Truncation + 1
      parameter(NRSP =(NTRU+1)*(NTRU+2))   ! No of real global    modes
      parameter(NCSP = NRSP / 2)           ! No of complex global modes
      parameter(NSPP = (NRSP+NPRO-1)/NPRO) ! Modes per process
      parameter(NESP = NSPP * NPRO)        ! Dim of spectral fields
      parameter(NROOT = 0)              ! Master node
      parameter(SOLAR_DAY=86400.)
      parameter(TFREEZE=271.25)
      parameter(CRHOI=920.)
!
      real (kind=8) :: gwd(NHOR)
!
      real :: ccls(NHOR)
      real :: cctaux(NHOR,0:13) = 0.
      real :: cctauy(NHOR,0:13) = 0.
      real :: ccprl(NHOR,0:13)  = 0.
      real :: ccprc(NHOR,0:13)  = 0.
      real :: ccevap(NHOR,0:13) = 0.
      real :: ccroff(NHOR,0:13) = 0.
      real :: cciced(NHOR,0:13) = 0.
      real :: ccsnow(NHOR,0:13) = 0.
      real :: ccsst(NHOR,0:13)  = 0.
      real :: ccswr(NHOR,0:13)  = 0.
      real :: cclwr(NHOR,0:13)  = 0.
      real :: ccshf(NHOR,0:13)  = 0.
      real :: cclhf(NHOR,0:13)  = 0.
!
      real :: cltaux(NHOR) = 0.
      real :: cltauy(NHOR) = 0.
      real :: clprl(NHOR)  = 0.
      real :: clprc(NHOR)  = 0.
      real :: clevap(NHOR) = 0.
      real :: clroff(NHOR) = 0.
      real :: cliced(NHOR) = 0.
      real :: clsnow(NHOR) = 0.
      real :: clsst(NHOR)  = 0.
      real :: clswr(NHOR)  = 0.
      real :: cllwr(NHOR)  = 0.
      real :: clshf(NHOR)  = 0.
      real :: cllhf(NHOR)  = 0.
!
      real :: yls(NHOR)     = 0.
      real :: ytaux(NHOR)   = 0.
      real :: ytauy(NHOR)   = 0.
      real :: yprl(NHOR)    = 0.
      real :: yprc(NHOR)    = 0.
      real :: yevap(NHOR)   = 0.
      real :: ysst(NHOR)    = 0.
      real :: ypme(NHOR)    = 0.
      real :: yroff(NHOR)   = 0.
      real :: yheat(NHOR)   = 0.
      real :: yicesnow(NHOR)= 0.
      real :: yfldo(NHOR)   = 0.
      real :: ymld(NHOR)    = 50.
!
      real :: ols(NHOR)     = 0.
      real :: otaux(NHOR)   = 0.
      real :: otauy(NHOR)   = 0.
      real :: oprl(NHOR)    = 0.
      real :: oprc(NHOR)    = 0.
      real :: oevap(NHOR)   = 0.
      real :: osst(NHOR)    = 0.
      real :: opme(NHOR)    = 0.
      real :: oheat(NHOR)   = 0.
      real :: oroff(NHOR)   = 0.
      real :: oofl(NHOR)    = 0.
      real :: oicesnow(NHOR)= 0.
      real :: oice(NHOR)    = 0.
      real :: osnow(NHOR)   = 0.
!
!     global integer
!
      integer :: ndatim(7) = -1
      integer :: nstep    = 0           ! time step
      integer :: naccuout = 0           ! counter for accumulated output
      integer :: nrestart = 0           ! switch for restart
      integer :: n_run_years = 0 
      integer :: n_run_months = 0
      integer :: n_run_days   = -1
      integer :: n_start_step = 0
      integer :: ntspd = 32
      integer :: n_days_per_month = 30
      integer :: n_days_per_year = 360
      integer :: n_start_year = 1
      integer :: n_start_month = 1
      integer :: nafter = 32
      integer :: naomod = 320
      integer :: nlsg   = 1
      integer :: ngui   = 0
!
!     Parallel Stuff
!
      integer :: mpinfo  = 0
      integer :: mypid   = 0
      integer :: myworld = 0
      integer :: nproc   = NPRO
!
!     needed to compile surfmod
!
      real :: so(NSPP),dls(NHOR)
!
      end module pumamod
!
      program climate
      use pumamod
!
      call clini
!
      mocd  = n_run_months       ! month countdown
      nscd  = n_run_days * ntspd ! step countdown
      if (nrestart == 0) nstep = n_start_step ! timestep since 01-01-0001
!
      call updatim(nstep)  ! set date & time array ndatim
!
      do while (mocd > 0 .or. nscd > 0)  ! main loop
       call clstep
       if(mod(nstep,nafter) == 0) then
        call clout
       endif
       iyea  = ndatim(1)    ! current year
       imon  = ndatim(2)    ! current month
       nstep = nstep + 1
       call updatim(nstep)  ! set date & time array ndatim
       if (imon /= ndatim(2)) then
        mocd = mocd - 1 ! next month
        write(nud,"('Completed month ',I2.2,'-',I4.4)") imon,iyea
       endif
       if (nscd  > 0) nscd = nscd - 1
       if (nscd == 0) exit
      enddo ! month countdown
!
      call clstop
!
      stop
      end
!
!     ================
!     SUBROUTINE CLINI
!     ================
!
      subroutine clini
      use pumamod
!
!     version identifier (date)
!
      character(len=80) :: version = '25.07.2008 by Larry'
!
      integer :: ih(8)
!
      real :: zin(NLON*NLAT)
      real :: zcltaux(NLON*NLAT,0:13) = 0.
      real :: zcltauy(NLON*NLAT,0:13) = 0.
      real :: zclprl(NLON*NLAT,0:13)  = 0.
      real :: zclprc(NLON*NLAT,0:13)  = 0.
      real :: zclevap(NLON*NLAT,0:13) = 0.
      real :: zclroff(NLON*NLAT,0:13) = 0.
      real :: zcliced(NLON*NLAT,0:13) = 0.
      real :: zclsnow(NLON*NLAT,0:13) = 0.
      real :: zclsst(NLON*NLAT,0:13)  = 0.
      real :: zclswr(NLON*NLAT,0:13)  = 0.
      real :: zcllwr(NLON*NLAT,0:13)  = 0.
      real :: zclshf(NLON*NLAT,0:13)  = 0.
      real :: zcllhf(NLON*NLAT,0:13)  = 0.
      real :: zlsm(NLON*NLAT)         = 0.
!
      logical :: lrestart = .false.
!
      character(len=80) :: climatefile='climate.txt'
      character(len=80) :: surface='surface.txt'

      namelist/clpar/ntspd,nafter,n_run_days,n_run_months,n_run_years  &
     &              ,n_start_year,n_start_month,nlsg,climatefile       &
     &              ,surface     
!
!     get process id
!
      call ompi_info(nproc,mypid)
!
!     read and print namelist and distribute it
!
      if (mypid == NROOT) then
         open(12,file='cl_namelist',form='formatted')
         read(12,clpar)
         write(nud,'(/," *******************************************")')
         write(nud,'(" * CLIMATE ",a30," *")') trim(version)
         write(nud,'(" *******************************************")')
         write(nud,'(" * Namelist CLPAR from <cl_namelist> *")')
         write(nud,'(" *******************************************")')
         write(nud,clpar)
         close(12)
         n_start_step = ntspd * (n_start_year     * n_days_per_year     &
                      + (n_start_month-1) * n_days_per_month)
         if(n_run_years > 0) n_run_months=n_run_years*12
      endif 
!
      call mpbci(ntspd)
      call mpbci(nlsg)
      call mpbci(nafter)
      call mpbci(n_days_per_month)
      call mpbci(n_days_per_year)
      call mpbci(n_run_days)
      call mpbci(n_run_months)
      call mpbci(n_run_years)
      call mpbci(n_start_step)
      call mpbci(n_start_year)
      call mpbci(n_start_month)
!
      call calini(n_days_per_month,n_days_per_year,n_start_step,ntspd   &
     &           ,solar_day,mypid)
!
      if(mypid == NROOT) then 
       call restart_ini(lrestart,'climate_restart')
       if (lrestart) then
        nrestart = 1
       endif
      endif
      call mpbci(nrestart)
!
      if (nrestart == 0) then
!
!     new start (read start file)
!
       if (mypid == NROOT) then
        open(20,FILE=trim(climatefile),FORM='formatted')
        do
         read(20,'(8I10)',IOSTAT=iostat) ih(:)
         if (iostat /= 0) exit
         read(20,'(8E12.6)',IOSTAT=iostat) zin(:)
         if (iostat /= 0) exit
         if (ih(1) == 169) then ! Code 169 for SST
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'SST month ',jmon,' found in ',trim(climatefile)
          zclsst(:,jmon) = zin(:)
         elseif (ih(1) == 180) then ! Code 180 for taux
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'TAUX month ',jmon,' found in ',trim(climatefile)
          zcltaux(:,jmon) = zin(:)
         elseif (ih(1) == 181) then ! Code 181 for tauy
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'TAUY month ',jmon,' found in ',trim(climatefile)
          zcltauy(:,jmon) = zin(:)
         elseif (ih(1) == 182) then ! Code 182 for evap
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'EVAP month ',jmon,' found in ',trim(climatefile)
          zclevap(:,jmon) = zin(:)
         elseif (ih(1) == 142) then ! Code 142 for prl
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'PRL month ',jmon,' found in ',trim(climatefile)
          zclprl(:,jmon) = zin(:)
         elseif (ih(1) == 143) then ! Code 143 for prc
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'PRC month ',jmon,' found in ',trim(climatefile)
          zclprc(:,jmon) = zin(:)
         elseif (ih(1) == 160) then ! Code 160 for roff
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'ROFF month ',jmon,' found in ',trim(climatefile)
          zclroff(:,jmon) = zin(:)
         elseif (ih(1) == 211) then ! Code 211 for iced
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'ICED month ',jmon,' found in ',trim(climatefile)
          zcliced(:,jmon) = zin(:)
         elseif (ih(1) == 141) then ! Code 141 for snow
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'SNOW month ',jmon,' found in ',trim(climatefile)
          zclsnow(:,jmon) = zin(:)
         elseif (ih(1) == 176) then ! Code 176 for swr
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'SWR month ',jmon,' found in ',trim(climatefile)
          zclswr(:,jmon) = zin(:)
         elseif (ih(1) == 177) then ! Code 177 for lwr
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'LWR month ',jmon,' found in ',trim(climatefile)
          zcllwr(:,jmon) = zin(:)
         elseif (ih(1) == 146) then ! Code 146 for shf
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'SHF month ',jmon,' found in ',trim(climatefile)
          zclshf(:,jmon) = zin(:)
         elseif (ih(1) == 147) then ! Code 147 for lhf
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'LHF month ',jmon,' found in ',trim(climatefile)
          zcllhf(:,jmon) = zin(:)
         elseif (ih(1) == 172) then ! Code 172 for Land-Sea-Mask
          write(nud,*) 'land-sea-mask found in ',trim(climatefile)
          zlsm(:) = zin(:)
         endif
        enddo
        close(20)
!
!       read sst and ice from surface.txt (if exist) and replace
!
        open(20,FILE=trim(surface),FORM='formatted')
        do
         read(20,'(8I10)',IOSTAT=iostat) ih(:)
         if (iostat /= 0) exit
         read(20,'(8E12.6)',IOSTAT=iostat) zin(:)
         if (iostat /= 0) exit
         if (ih(1) == 169) then ! Code 169 for SST
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'SST month ',jmon,' found in ',trim(surface)
          zclsst(:,jmon) = zin(:)
         elseif (ih(1) == 211) then ! Code 211 for iced
          jmon=ih(3)
          if (jmon > 12) jmon = mod(jmon/100,100)
          write(nud,*) 'ICED month ',jmon,' found in ',trim(surface)
          zcliced(:,jmon) = zin(:)
         elseif (ih(1) == 172) then ! Code 172 for Land-Sea-Mask
          write(nud,*) 'land-sea-mask found in ',trim(surface)
          zlsm(:) = zin(:)
         endif
        enddo
        close(20)
       endif ! mypid
       zclsst(:,0)=zclsst(:,12)
       zclsst(:,13)=zclsst(:,1)
       zcltaux(:,0)=zcltaux(:,12)
       zcltaux(:,13)=zcltaux(:,1)
       zcltauy(:,0)=zcltauy(:,12)
       zcltauy(:,13)=zcltauy(:,1)
       zclprl(:,0)=zclprl(:,12)
       zclprl(:,13)=zclprl(:,1)
       zclprc(:,0)=zclprc(:,12)
       zclprc(:,13)=zclprc(:,1)
       zclevap(:,0)=zclevap(:,12)
       zclevap(:,13)=zclevap(:,1)
       zclroff(:,0)=zclroff(:,12)
       zclroff(:,13)=zclroff(:,1)
       zcliced(:,0)=zcliced(:,12)
       zcliced(:,13)=zcliced(:,1)
       zclsnow(:,0)=zclsnow(:,12)
       zclsnow(:,13)=zclsnow(:,1)
       zclswr(:,0)=zclswr(:,12)
       zclswr(:,13)=zclswr(:,1)
       zcllwr(:,0)=zcllwr(:,12)
       zcllwr(:,13)=zcllwr(:,1)
       zclshf(:,0)=zclshf(:,12)
       zclshf(:,13)=zclshf(:,1)
       zcllhf(:,0)=zcllhf(:,12)
       zcllhf(:,13)=zcllhf(:,1)
!
       call mpscgp(zclsst,ccsst,14)
       call mpscgp(zcltaux,cctaux,14)
       call mpscgp(zcltauy,cctauy,14)
       call mpscgp(zclevap,ccevap,14)
       call mpscgp(zclprl,ccprl,14)
       call mpscgp(zclprc,ccprc,14)
       call mpscgp(zclroff,ccroff,14)
       call mpscgp(zcliced,cciced,14)
       call mpscgp(zclsnow,ccsnow,14)
       call mpscgp(zclswr,ccswr,14)
       call mpscgp(zcllwr,cclwr,14)
       call mpscgp(zclshf,ccshf,14)
       call mpscgp(zcllhf,cclhf,14)
       call mpscgp(zls,ccls,1)
!
!     make clsst >= tfreeze
! 
       do jm=0,13
        ccsst(:,jm)=AMAX1(ccsst(:,jm),TFREEZE)
       enddo
!
!     initialize 
!
       naccuout=0
!
       call clget
!
       yls(:)=real(nint(ccls(:)))
       ytaux(:)=cltaux(:)
       ytauy(:)=cltauy(:)
       yprl(:)=clprl(:)
       yprc(:)=clprc(:)
       yevap(:)=clevap(:)
       ysst(:)=clsst(:)
       ypme(:)=clprl(:)+clprc(:)+clevap(:)
       yroff(:)=clroff(:)
       yicesnow(:)=cliced(:)*CRHOI/1000.+clsnow(:)
       yheat(:)=clswr(:)+cllwr(:)+clshf(:)+cllhf(:)
!
      else

!
!     restart from restart file
!
       if (mypid == NROOT) then
        call get_restart_integer('nstep'   ,nstep)
        call get_restart_integer('naccu'   ,naccuout)
       endif
!
       call mpgetgp('ccls',ccls,NHOR,1)
       call mpgetgp('ccsst',ccsst,NHOR,14)
       call mpgetgp('cctaux',cctaux,NHOR,14)
       call mpgetgp('cctauy',cctauy,NHOR,14)
       call mpgetgp('ccprl',ccprl,NHOR,14)
       call mpgetgp('ccprc',ccprc,NHOR,14)
       call mpgetgp('ccevap',ccevap,NHOR,14)
       call mpgetgp('ccroff',ccroff,NHOR,14)
       call mpgetgp('cciced',cciced,NHOR,14)
       call mpgetgp('ccsnow',ccsnow,NHOR,14)
       call mpgetgp('ccswr',ccswr,NHOR,14)
       call mpgetgp('cclwr',cclwr,NHOR,14)
       call mpgetgp('ccshf',ccshf,NHOR,14)
       call mpgetgp('cclhf',cclhf,NHOR,14)
       call mpgetgp('ols',ols,NHOR,1)
       call mpgetgp('otaux',otaux,NHOR,1)
       call mpgetgp('otauy',otauy,NHOR,1)
       call mpgetgp('oprl',oprl,NHOR,1)
       call mpgetgp('oprc',oprc,NHOR,1)
       call mpgetgp('oevap',oevap,NHOR,1)
       call mpgetgp('osst',osst,NHOR,1)
       call mpgetgp('opme',opme,NHOR,1)
       call mpgetgp('oicesnow',oicesnow,NHOR,1)
       call mpgetgp('oice',oice,NHOR,1)
       call mpgetgp('osnow',osnow,NHOR,1)
       call mpgetgp('oheat',oheat,NHOR,1)
       call mpgetgp('oofl',oofl,NHOR,1)
       call mpgetgp('oroff',oroff,NHOR,1)
       yls(:)=ccls(:)
!
      endif
!
      call mpbci(nstep)
      call mpbci(naccuout)
!
!     open output file
!
      if (mypid == NROOT) then
       open(31,file='cl_output',form='unformatted')
      endif
!
      if (nlsg > 0) then
       call mpgagp(zlsm,yls,1)
       call ntomin(nstep,ndatim(5),ndatim(4),ndatim(3),ndatim(2)        &
     &            ,ndatim(1))
       if(mypid == nroot) then
!
!      check if day_per_year are right (lsg coupling needs 360)
!
        kdpy=n_days_per_year
        if(kdpy .ne. 360.) then
         write(nud,*) '!ERROR! for LSG coupling you need to set '
         write(nud,*) '        n_days_per_year in plasim namelist INP '
         write(nud,*) '        to 360 !'
         write(nud,*) '        at the moment, n_days_per_year= ',kdpy
         write(nud,*) 'Model stoped!'
         stop 'ERROR'
        endif
        call clsgini(ndatim,ntspd,naomod,nlsg,ngui,zlsm,nlon,nlat)
       endif
      endif
!
      if (mypid == NROOT) then
       open(31,file='cl_output',form='unformatted')
      endif
!
      return
      end subroutine clini

!     ===================================================================
!     SUBROUTINE clstep
!     ===================================================================

      subroutine clstep
      use pumamod
!
      real :: zsst(nlon*nlat)
      real :: ztaux(nlon*nlat)
      real :: ztauy(nlon*nlat)
      real :: zpme(nlon*nlat)
      real :: zroff(nlon*nlat)
      real :: zice(nlon*nlat)
      real :: zheat(nlon*nlat)
      real :: zfldo(nlon*nlat)
      real :: zmld(nlon*nlat)
!
!     dbug arrays
!
      real,allocatable :: zprf1(:),zprf2(:),zprf3(:),zprf4(:),zprf5(:)
!
!     get climate
!
      call clget
!
!     set climate (exept sst which is set later according to oceanmod)
!
      yls(:)=ccls(:) 
      ytaux(:)=cltaux(:)
      ytauy(:)=cltauy(:)
      yprl(:)=clprl(:)
      yprc(:)=clprc(:)
      yevap(:)=clevap(:)
      ypme(:)=clprl(:)+clprc(:)+clevap(:) 
      yroff(:)=clroff(:) 
      yicesnow(:)=cliced(:)*CRHOI/1000.+clsnow(:)
      yheat(:)=clswr(:)+cllwr(:)+clshf(:)+cllhf(:)
!
      ols(:)=ccls(:)+ols(:)
      otaux(:)=cltaux(:)+otaux(:)
      otauy(:)=cltauy(:)+otauy(:)
      oprl(:)=clprl(:)+oprl(:)
      oprc(:)=clprc(:)+oprc(:)
      oevap(:)=clevap(:)+oevap(:)
      osst(:)=clsst(:)+osst(:)
      opme(:)=clprl(:)+clprc(:)+clevap(:)+opme(:)
      oroff(:)=clroff(:)+oroff(:)
      oicesnow(:)=cliced(:)*CRHOI/1000.+clsnow(:)+oicesnow(:)
      oice(:)=cliced(:)+oice(:)
      osnow(:)=clsnow(:)+osnow(:)
      oheat(:)=clswr(:)+cllwr(:)+clshf(:)+cllhf(:)+oheat(:)
      naccuout=naccuout+1
!
!     do the lsg coupling
!
      if (nlsg > 0) then
       if(nlsg < 2) call mpgagp(zsst,ysst,1)
       call mpgagp(ztaux,ytaux,1)
       call mpgagp(ztauy,ytauy,1)
       call mpgagp(zpme,ypme,1)
       call mpgagp(zroff,yroff,1)
       call mpgagp(zice,yicesnow,1)
       call mpgagp(zheat,yheat,1)
!!       call mpgagp(zmld,ymld,1)
       call ntomin(nstep,ndatim(5),ndatim(4),ndatim(3),ndatim(2)        &
     &            ,ndatim(1))
       if(mypid == nroot) then
        call clsgstep(ndatim,nstep,zsst,ztaux,ztauy,zpme,zroff,zice     &
     &               ,zheat,zfldo)
       endif
       if(nlsg < 2) then
        if(mod(nstep,naomod) == naomod-1) call mpscgp(zfldo,yfldo,1)
       else
        yfldo(:)=0.
       endif
       oofl(:)=oofl(:)+yfldo(:)
      endif
!
!     if lsg coupling due to ahfl and osst, replace ssts
!
      if(nlsg > 1) then
       if(mod(nstep,naomod) == naomod-1) call mpscgp(zsst,ysst,1)
      endif
!
!     set sst according to oceanmod
!
      ysst(:)=clsst(:)
!
      return
      end subroutine clstep

!     ===================================
!     SUBROUTINE CLSTOP
!     ===================================

      subroutine clstop
      use pumamod
!
!     write restart file
!
      if (mypid == NROOT) then
       call restart_prepare('climate_status')
       call put_restart_integer('nstep',nstep)
       call put_restart_integer('naccu',naccuout)
      endif
!
      call mpputgp('ccls',ccls,NHOR,1)
      call mpputgp('ccsst',ccsst,NHOR,14)
      call mpputgp('cctaux',cctaux,NHOR,14)
      call mpputgp('cctauy',cctauy,NHOR,14)
      call mpputgp('ccprl',ccprl,NHOR,14)
      call mpputgp('ccprc',ccprc,NHOR,14)
      call mpputgp('ccevap',ccevap,NHOR,14)
      call mpputgp('ccroff',ccroff,NHOR,14)
      call mpputgp('cciced',cciced,NHOR,14)
      call mpputgp('ccsnow',ccsnow,NHOR,14)
      call mpputgp('ccswr',ccswr,NHOR,14)
      call mpputgp('cclwr',cclwr,NHOR,14)
      call mpputgp('ccshf',ccshf,NHOR,14)
      call mpputgp('cclhf',cclhf,NHOR,14)
      call mpputgp('ols',ols,NHOR,1)
      call mpputgp('otaux',otaux,NHOR,1)
      call mpputgp('otauy',otauy,NHOR,1)
      call mpputgp('oprl',oprl,NHOR,1)
      call mpputgp('oprc',oprc,NHOR,1)
      call mpputgp('oevap',oevap,NHOR,1)
      call mpputgp('osst',osst,NHOR,1)
      call mpputgp('opme',opme,NHOR,1)
      call mpputgp('oicesnow',oicesnow,NHOR,1)
      call mpputgp('oice',oice,NHOR,1)
      call mpputgp('osnow',osnow,NHOR,1)
      call mpputgp('oheat',oheat,NHOR,1)
      call mpputgp('oofl',oofl,NHOR,1)
      call mpputgp('oroff',oroff,NHOR,1)
!
!     close output file
!
      if (mypid == NROOT) then
       close(31)
      endif
!
!     stop lsg coupling
!
      if(nlsg > 0) then
       if(mypid == NROOT) then
        call clsgstop
       endif
      endif
!
      return
      end subroutine clstop

!     ===================================
!     SUBROUTINE CLOUT
!     ===================================

      subroutine clout
      use pumamod
!
      integer ih(8)
!
      zn=1./real(naccuout)
      ols(:)=ols(:)*zn
      osst(:)=osst(:)*zn
      otaux(:)=otaux(:)*zn
      otauy(:)=otauy(:)*zn
      oprl(:)=oprl(:)*zn
      oprc(:)=oprc(:)*zn
      oevap(:)=oevap(:)*zn
      oroff(:)=oroff(:)*zn
      oice(:)=oice(:)*zn
      oofl(:)=oofl(:)*zn
      oicesnow(:)=oicesnow(:)*zn
      opme(:)=opme(:)*zn
      osnow(:)=osnow(:)*zn
      oheat(:)=oheat(:)*zn
!
      call writegp(31,oprl,142,1)
      call writegp(31,oprc,143,1)
      call writegp(31,oroff,160,1)
      call writegp(31,osst,169,1)
      call writegp(31,ols,172,1)
      call writegp(31,otaux,180,1)
      call writegp(31,otauy,181,1)
      call writegp(31,oevap,182,1)
      call writegp(31,oice,211,1)
      call writegp(31,oofl,222,1)
      call writegp(31,oicesnow,69,1)
      call writegp(31,oheat,263,1)
      call writegp(31,opme,260,1)
      call writegp(31,osnow,141,1)
!
      naccuout=0
      ols(:)=0.
      osst(:)=0.
      otaux(:)=0.
      otauy(:)=0.
      oprl(:)=0.
      oprc(:)=0.
      oevap(:)=0.
      oroff(:)=0.
      oice(:)=0.
      oofl(:)=0.
      oicesnow(:)=0.
      opme(:)=0.
      osnow(:)=0.
      oheat(:)=0.
!
      return
      end subroutine clout

!====================================================================
!      SUBROUTINE CLGET
!====================================================================

      subroutine clget
      use pumamod

!     *********************
!     * get  annual cycle *
!     *********************

!     modified for 14 months version of SST
!     jm2: 0=Dec, 1-12:Months, 13=Jan

      jstep=nstep+1    !advance time step (ts = ts(t)+ dtsdt)
      call ntomin(jstep,nmin,nhour,nday,nmonth,nyear)
!
      jm1 = nmonth
      if (nday > 15) then ! Does not work for other planets
         jm2=jm1+1
      else
         jm2=jm1-1
      endif

!     ****************************************
!*    interpolate
!     ****************************************

       zgw2 = abs((nday-1)*1440+nhour*60+nmin-21600)/43200.
       zgw1 = 1.0 - zgw2
       clsst(:)   = zgw1*ccsst(:,jm1)+zgw2*ccsst(:,jm2)
       cltaux(:)  = zgw1*cctaux(:,jm1)+zgw2*cctaux(:,jm2)
       cltauy(:)  = zgw1*cctauy(:,jm1)+zgw2*cctauy(:,jm2)
       clprl(:)  = zgw1*ccprl(:,jm1)+zgw2*ccprl(:,jm2)
       clprc(:)  = zgw1*ccprc(:,jm1)+zgw2*ccprc(:,jm2)
       clevap(:)  = zgw1*ccevap(:,jm1)+zgw2*ccevap(:,jm2)
       clroff(:)  = zgw1*ccroff(:,jm1)+zgw2*ccroff(:,jm2)
       cliced(:)  = zgw1*cciced(:,jm1)+zgw2*cciced(:,jm2)
       clsnow(:)  = zgw1*ccsnow(:,jm1)+zgw2*ccsnow(:,jm2)
       clswr(:)  = zgw1*ccswr(:,jm1)+zgw2*ccswr(:,jm2)
       cllwr(:)  = zgw1*cclwr(:,jm1)+zgw2*cclwr(:,jm2)
       clshf(:)  = zgw1*ccshf(:,jm1)+zgw2*ccshf(:,jm2)
       cllhf(:)  = zgw1*cclhf(:,jm1)+zgw2*cclhf(:,jm2)

       clsnow(:)=max(clsnow(:),0.)
       cliced(:)=max(cliced(:),0.)
       clprl(:)=max(clprl(:),0.)
       clprc(:)=max(clprc(:),0.)
       clroff(:)=max(clroff(:),0.)
       clevap(:)=min(clevap(:),0.)
       where(cliced(:)==0.) clsnow=0.

       return
       end subroutine clget

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

      subroutine writegp(kunit,pf,kcode,klev)
      use pumamod
      real :: pf(NHOR)
      real :: zf(NUGP)
      integer(kind=4) :: ihead(8)
      real(kind=4)  :: zzf(NUGP)

      istep = nstep
      call ntomin(istep,nmin,nhour,nday,nmonth,nyear)

      ihead(1) = kcode
      ihead(2) = klev
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
      ihead(4) = nmin + 100 * nhour
      ihead(5) = NLON
      ihead(6) = NLAT
      ihead(7) = 1
      ihead(8) = n_days_per_year

      call mpgagp(zf,pf,1)

      if (mypid == NROOT) then
       write(kunit) ihead
       zzf(:)=zf(:)
       write(kunit) zzf
      endif

      return
      end
!
!     dummies
!
      subroutine get_surf_array(yn,pa,k1,k2,kread)
      character (len=*) :: yn ! array name
      real :: pa(k1,k2)       ! array to receive values
      integer :: k1           ! 1. dim of array
      integer :: k2           ! 2. dim of array
      integer :: kread        ! flag to indicate success
      return
      end
      subroutine get_3d_year(yn,pa,kdim,klev,kread)
      character (len=*) :: yn
      real :: pa(kdim,klev,0:13)
      return
      end
      subroutine get_surf_year(yn,pa,kdim,kread)
      character (len=*) :: yn
      integer :: ih(8)
      real :: pa(kdim,0:13)
      return
      end

