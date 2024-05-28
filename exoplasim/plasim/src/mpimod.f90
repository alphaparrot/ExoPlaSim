      module mpimod
      use pumamod
!      include 'mpif.h'
      use mpi

      integer :: mpi_itype = MPI_INTEGER4
      integer :: mpi_rtype = MPI_REAL4
      integer :: mpi_ltype = MPI_LOGICAL

      character(len=80) ynode(NPRO)           ! node names

      end module mpimod
!
!     interface routines to MPI:
!


!     ================
!     SUBROUTINE MPBCI
!     ================

      subroutine mpbci(k) ! broadcast 1 integer
      use mpimod
      integer :: k(*)

      call mpi_bcast(k,1,mpi_itype,NROOT,myworld,mpinfo)

      return
      end subroutine mpbci

!     =================
!     SUBROUTINE MPBCIN
!     =================

      subroutine mpbcin(k,n) ! broadcast n integer
      use mpimod

      integer :: k(n)

      call mpi_bcast(k,n,mpi_itype,NROOT,myworld,mpinfo)

      return
      end subroutine mpbcin

!     ================
!     SUBROUTINE MPBCR
!     ================

      subroutine mpbcr(p) ! broadcast 1 real
      use mpimod
      REAL :: p(*)

      call mpi_bcast(p,1,mpi_rtype,NROOT,myworld,mpinfo)

      return
      end subroutine mpbcr

!     =================
!     SUBROUTINE MPBCRN
!     =================

      subroutine mpbcrn(p,n) ! broadcast n real
      use mpimod
      real :: p(n)

      call mpi_bcast(p,n,mpi_rtype,NROOT,myworld,mpinfo)

      return
      end subroutine mpbcrn

!     ================
!     SUBROUTINE MPBCL
!     ================

      subroutine mpbcl(l) ! broadcast 1 logical
      use mpimod
      logical :: l(*)

      call mpi_bcast(l,1,mpi_ltype,NROOT,myworld,mpinfo)

      return
      end subroutine mpbcl

!     =================
!     SUBROUTINE MPSCIN
!     =================

      subroutine mpscin(k,n) ! scatter n integer
      use mpimod

      integer :: k(n)

      call mpi_scatter(k,n,mpi_itype,k,n,mpi_itype,NROOT,myworld,mpinfo)

      return
      end subroutine mpscin

!     =================
!     SUBROUTINE MPSCRN
!     =================

      subroutine mpscrn(p,n) ! scatter n real
      use mpimod

      real :: p(*)

      call mpi_scatter(p,n,mpi_rtype,p,n,mpi_rtype,NROOT,myworld,mpinfo)

      return
      end subroutine mpscrn

!     =================
!     SUBROUTINE MPSCDN
!     =================

      subroutine mpscdn(p,n) ! scatter n double precision
      use mpimod

      real (kind=8) :: p(*)

      call mpi_scatter(p,n,MPI_REAL8,p,n,MPI_REAL8,NROOT,myworld,mpinfo)

      return
      end subroutine mpscdn

!     =================
!     SUBROUTINE MPSCGP
!     =================

      subroutine mpscgp(pf,pp,klev) ! scatter gridpoint fields
      use mpimod

      real :: pf(NUGP,klev)
      real :: pp(NHOR,klev)

      do jlev = 1 , klev
         call mpi_scatter(pf(:,jlev),NHOR,mpi_rtype,                     &
     &                    pp(:,jlev),NHOR,mpi_rtype,                     &
     &                    NROOT,myworld,mpinfo)
      enddo

      return
      end subroutine mpscgp

!     =================
!     SUBROUTINE MPGAGP
!     =================

      subroutine mpgagp(pf,pp,klev) ! gather gridpoint fields
      use mpimod

      real :: pf(NLON*NLAT,klev)
      real :: pp(NHOR,klev)

      do jlev = 1 , klev
         call mpi_gather(pp(:,jlev),NHOR,mpi_rtype,                      &
     &                   pf(:,jlev),NHOR,mpi_rtype,                      &
     &                   NROOT,myworld,mpinfo)
      enddo

      return
      end subroutine mpgagp

!     ===================
!     SUBROUTINE MPGALLGP
!     ===================

      subroutine mpgallgp(pf,pp,klev) ! gather gritpoint to all
      use mpimod

      real :: pf(NLON*NLAT,klev)
      real :: pp(NHOR,klev)

      do jlev = 1 , klev
         call mpi_allgather(pp(:,jlev),NHOR,mpi_rtype,                   &
     &                      pf(:,jlev),NHOR,mpi_rtype,                   &
     &                      myworld,mpinfo)
      enddo

      return
      end subroutine mpgallgp

!     =================
!     SUBROUTINE MPSCSP
!     =================

      subroutine mpscsp(pf,pp,klev) ! scatter spectral fields
      use mpimod

      real :: pf(NESP,klev)
      real :: pp(NSPP,klev)

      do jlev = 1 , klev
         call mpi_scatter(pf(:,jlev),NSPP,mpi_rtype                      &
     &                   ,pp(:,jlev),NSPP,mpi_rtype                      &
     &                   ,NROOT,myworld,mpinfo)
      enddo

      return
      end subroutine mpscsp

!     =================
!     SUBROUTINE MPGASP
!     =================

      subroutine mpgasp(pf,pp,klev) ! gather spectral fields
      use mpimod

      real :: pf(NESP,klev)
      real :: pp(NSPP,klev)

      do jlev = 1 , klev
         call mpi_gather(pp(:,jlev),NSPP,mpi_rtype                       &
     &                  ,pf(:,jlev),NSPP,mpi_rtype                       &
     &                  ,NROOT,myworld,mpinfo)
      enddo

      return
      end subroutine mpgasp

!     =================
!     SUBROUTINE MPGACS
!     =================

      subroutine mpgacs(pcs) ! gather cross sections
      use mpimod

      real :: pcs(NLAT,NLEV)

      do jlev = 1 , NLEV
         call mpi_gather(pcs(:,jlev),NLPP,mpi_rtype                      &
     &                  ,pcs(:,jlev),NLPP,mpi_rtype                      &
     &                  ,NROOT,myworld,mpinfo)

      enddo
      return
      end subroutine mpgacs

!     ===================
!     SUBROUTINE MPGALLSP
!     ===================

      subroutine mpgallsp(pf,pp,klev) ! gather spectral to all
      use mpimod

      real :: pf(NESP,klev)
      real :: pp(NSPP,klev)

      do jlev = 1 , klev
         call mpi_allgather(pp(:,jlev),NSPP,mpi_rtype                    &
     &                     ,pf(:,jlev),NSPP,mpi_rtype                    &
     &                     ,myworld,mpinfo)
      enddo

      return
      end subroutine mpgallsp

!     ================
!     SUBROUTINE MPSUM
!     ================

      subroutine mpsum(psp,klev) ! sum spectral fields
      use mpimod

      real :: psp(NESP*klev)
      real :: tmp(NESP*klev)

      call mpi_reduce(psp,tmp,NESP*klev,mpi_rtype,MPI_SUM     &
     &               ,NROOT,myworld,mpinfo)
      if (mypid == NROOT) psp = tmp

      return
      end subroutine mpsum

!     ==================
!     SUBROUTINE MPSUMSC
!     ==================

      subroutine mpsumsc(psf,psp,klev) ! sum & scatter spectral
      use mpimod

      real :: psf(NESP,klev)
      real :: psp(NSPP,klev)

      do jlev = 1 , klev
         call mpi_reduce_scatter(psf(:,jlev),psp(:,jlev),nscatsp        &
     &                          ,mpi_rtype,MPI_SUM,myworld,mpinfo)
      enddo

      return
      end subroutine mpsumsc

!     ====================
!     SUBROUTINE MPSUMR
!     ====================

      subroutine mpsumr(pr,kdim) ! sum kdim reals
      use mpimod

      real pr(kdim)
      real tmp(kdim)

      call mpi_reduce(pr,tmp,kdim,mpi_rtype,MPI_SUM,NROOT,myworld,mpinfo)
      if (mypid == NROOT) pr = tmp

      return
      end subroutine mpsumr

!     ====================
!     SUBROUTINE MPSUMBCR
!     ====================

      subroutine mpsumbcr(pr,kdim) ! sum & broadcast kdim reals
      use mpimod

      real pr(kdim)
      real tmp(kdim)

      call mpi_allreduce(pr,tmp,kdim,mpi_rtype,MPI_SUM,myworld,mpinfo)
      pr = tmp

      return
      end subroutine mpsumbcr


!     ==================
!     SUBROUTINE MPABORT
!     ==================

      subroutine mpabort(ym)
      use mpimod

      character (len=* ) :: ym
      character (len=64) :: ystar = ' '
      character (len=64) :: yline = ' '
      character (len=64) :: yabor = 'Program aborted'
      character (len=64) :: ymess = ' '
      character (len=64) :: yhead = ' '

      ilmess = len_trim(ym)
      ilabor = len_trim(yabor)
      ilen   = 60

      do j = 1 , ilen+4
         ystar(j:j) = '*'
         yline(j:j) = '-'
      enddo

      ioff = 2
      if (ilmess < ilen-1) ioff = ioff + (ilen - ilmess) / 2
      ymess(1+ioff:ilmess+ioff) = trim(ym)

      ioff = 2
      if (ilabor < ilen-1) ioff = ioff + (ilen - ilabor) / 2
      yhead(1+ioff:ilabor+ioff) = trim(yabor)

      yline(1:1) = '*'
      ymess(1:1) = '*'
      yhead(1:1) = '*'
      yline(2:2) = ' '
      ymess(2:2) = ' '
      yhead(2:2) = ' '
      j = ilen + 4
      yline(j:j) = '*'
      ymess(j:j) = '*'
      yhead(j:j) = '*'
      j = ilen + 3
      yline(j:j) = ' '
      ymess(j:j) = ' '
      yhead(j:j) = ' '

      if (mypid == NROOT) then
         open (44,file='Abort_Message')
         write(44,'(A)') trim(ystar)
         write(44,'(A)') trim(yhead)
         write(44,'(A)') trim(yline)
         write(44,'(A)') trim(ymess)
         write(44,'(A)') trim(ystar)
         close(44)
         write(nud,'(/,A)') trim(ystar)
         write(nud,'(A)') trim(yhead)
         write(nud,'(A)') trim(yline)
         write(nud,'(A)') trim(ymess)
         write(nud,'(A,/)') trim(ystar)
      call mpi_abort(myworld,mpinfo,mpinfo)
      endif

      stop
      end
      

!     ==================
!     SUBROUTINE MPSTART
!     ==================

      subroutine mpstart(kworld) ! initialization
      use mpimod
      integer :: itest = 0
      real    :: rtest = 0.0
      logical :: ltest = .true.

      if (kind(itest) == 8) mpi_itype = MPI_INTEGER8
      if (kind(rtest) == 8) mpi_rtype = MPI_REAL8

      if (kworld < 0) then
         call mpi_init(mpinfo)
         myworld=MPI_COMM_WORLD
      else
         myworld = kworld
      endif

      call mpi_comm_size(myworld,nproc,mpinfo)
      call mpi_comm_rank(myworld,mypid,mpinfo)

      if (nproc .ne. NPRO .and. mypid == NROOT) then
         write(nud,*)'Compiled for ',NPRO,' nodes'
         write(nud,*)'Running on ',nproc,' nodes'   
         call mpi_abort(myworld,mpinfo,mpinfo)
      endif

      allocate(ympname(npro)) ; ympname(:) = ' '
      call mpi_get_processor_name(ympname(1),ilen,mpinfo)

      call mpi_gather(ympname,80,MPI_CHARACTER,   &   
                      ympname,80,MPI_CHARACTER,   &   
                      NROOT,myworld,mpinfo)

      return
      end subroutine mpstart

!     =================
!     SUBROUTINE MPSTOP
!     =================

      subroutine mpstop
      use mpimod

      call mpi_barrier(myworld,mpinfo)
      call mpi_finalize(mpinfo)

      return
      end subroutine mpstop

!     ===================
!     SUBROUTINE MPREADGP
!     ===================

      subroutine mpreadgp(ktape,p,kdim,klev)
      use mpimod

      real p(kdim,klev)
      real z(NLON*NLAT,klev)

      z = 0.0

      if (mypid == NROOT) read (ktape) z(:,:)

      if (kdim == NHOR) then
       call mpscgp(z,p,klev)
      else
       if (mypid == NROOT) p = z
      endif

      return
      end subroutine mpreadgp

!     ====================
!     SUBROUTINE MPWRITEGP
!     ====================

      subroutine mpwritegp(ktape,p,kdim,klev)
      use mpimod

      real p(kdim,klev)
      real z(NLON*NLAT,klev)

      if (kdim == NHOR) then
       call mpgagp(z,p,klev)
       if (mypid == NROOT) write(ktape) z(1:NLON*NLAT,:)
      else
       if (mypid == NROOT) write(ktape) p(1:NLON*NLAT,:)
      endif

      return
      end subroutine mpwritegp


!     =====================
!     SUBROUTINE MPWRITEGPH
!     =====================

      subroutine mpwritegph(ktape,p,kdim,klev,ihead)
      use mpimod

      real p(kdim,klev)
      real z(NLON*NLAT,klev)
!
      real(kind=4) :: zp(kdim,klev)
      real(kind=4) :: zz(NLON*NLAT,klev)
!

      integer ihead(8)

      if (kdim == NHOR) then
       call mpgagp(z,p,klev)
       if (mypid == NROOT) then
          write(ktape) ihead
          zz(:,:)=z(:,:)
          write(ktape) zz(1:NLON*NLAT,:)
       endif
      else
       if (mypid == NROOT) then
          write(ktape) ihead
          zp(:,:)=p(:,:)
          write(ktape) zp(1:NLON*NLAT,:)
       endif
      endif

      return
      end subroutine mpwritegph

!     ===================
!     SUBROUTINE MPREADSP
!     ===================

      subroutine mpreadsp(ktape,p,kdim,klev)
      use mpimod

      real p(kdim,klev)
      real z(NESP,klev)

      z = 0.0
      if (mypid == NROOT) read(ktape) ((z(i,j),i=1,NRSP),j=1,klev)
      if (kdim == NSPP) then
         call mpscsp(z,p,klev)
      else
         if (mypid == NROOT) p = z
      endif

      return
      end subroutine mpreadsp

!     ====================
!     SUBROUTINE MPWRITESP
!     ====================

      subroutine mpwritesp(ktape,p,kdim,klev)
      use mpimod

      real p(kdim,klev)
      real z(NESP,klev)

      if (kdim == NSPP) then
         call mpgasp(z,p,klev)
         if (mypid == NROOT) write(ktape) ((z(i,j),i=1,NRSP),j=1,klev)
      else
         if (mypid == NROOT) write(ktape) ((z(i,j),i=1,NRSP),j=1,klev)
      endif

      return
      end subroutine mpwritesp


!     ===================
!     SUBROUTINE OMPI_INFO
!     ===================

      subroutine ompi_info(nprocess,npid)    ! get nproc and pid
      use mpimod

      myworld=MPI_COMM_WORLD

      call mpi_comm_size(myworld,nprocess,mpinfo)
      call mpi_comm_rank(myworld,npid,mpinfo)

      return
      end subroutine ompi_info


!     ==================
!     SUBROUTINE MPGETSP
!     ==================

      subroutine mpgetsp(yn,p,kdim,klev)
      use mpimod

      character (len=*) :: yn
      real :: p(kdim,klev)
      real :: z(NESP,klev)

      z(:,:) = 0.0
      if (mypid == NROOT) call get_restart_array(yn,z,NRSP,NESP,klev)
      call mpscsp(z,p,klev)

      return
      end subroutine mpgetsp


!     ==================
!     SUBROUTINE MPGETGP
!     ==================

      subroutine mpgetgp(yn,p,kdim,klev)
      use mpimod

      character (len=*) :: yn
      real :: p(kdim,klev)
      real :: z(NUGP,klev)

      if (mypid == NROOT) call get_restart_array(yn,z,NUGP,NUGP,klev)
      call mpscgp(z,p,klev)

      return
      end subroutine mpgetgp


!     ==================
!     SUBROUTINE MPPUTSP
!     ==================

      subroutine mpputsp(yn,p,kdim,klev)
      use mpimod

      character (len=*) :: yn
      real :: p(kdim,klev)
      real :: z(NESP,klev)

      call mpgasp(z,p,klev)
      if (mypid == NROOT) call put_restart_array(yn,z,NRSP,NESP,klev)

      return
      end subroutine mpputsp


!     ==================
!     SUBROUTINE MPPUTGP
!     ==================

      subroutine mpputgp(yn,p,kdim,klev)
      use mpimod

      character (len=*) :: yn
      real :: p(kdim,klev)
      real :: z(NUGP,klev)

      call mpgagp(z,p,klev)
      if (mypid == NROOT) call put_restart_array(yn,z,NUGP,NUGP,klev)

      return
      end subroutine mpputgp


!     ===================
!     SUBROUTINE MPSURFGP
!     ===================

      subroutine mpsurfgp(yn,p,kdim,klev)
      use mpimod

      character (len=*) :: yn
      real :: p(kdim,klev)
      real :: z(NUGP,klev)
      integer :: iread(1)

      if (mypid == NROOT) call get_surf_array(yn,z,NUGP,klev,iread)
      call mpbci(iread)
      if (iread(1) == 1) call mpscgp(z,p,klev)

      return
      end subroutine mpsurfgp


!     ===================
!     SUBROUTINE MPMAXVAL
!     ===================

      subroutine mpmaxval(p,kdim,klev,pmax)
      use mpimod

      real :: p(kdim,klev)
      real :: pmax(1)
      real zmax(1)

      zmax(1) = maxval(p(:,:))
      call mpi_allreduce(zmax,pmax,1,mpi_rtype,MPI_MAX,myworld,mpinfo)

      return
      end subroutine mpmaxval


!     ===================
!     SUBROUTINE MPSUMVAL
!     ===================

      subroutine mpsumval(p,kdim,klev,psum)
      use mpimod

      real :: p(kdim,klev)
      real :: psum(1)
      real :: zsum(1)

      zsum = sum(p(:,:))
      call mpi_allreduce(zsum,psum,1,mpi_rtype,MPI_SUM,myworld,mpinfo)

      return
      end subroutine mpsumval

      subroutine mrsum(k) ! sum up 1 integer
      return
      end 

      subroutine mrbci(k) ! broadcast 1 integer
      return
      end 

      subroutine mrdiff(p,d,n)
      real :: p(n)
      real :: d(n)
      return
      end 

      subroutine mrdimensions ! used in mpimod_multi.f90
      return
      end 

