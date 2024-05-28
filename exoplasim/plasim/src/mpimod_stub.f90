!     =======================================
!     mpimod_dummy.f90
!     ----------------
!     This module replaces <mpimod.f90> for 
!     single CPU runs
!
!     The module is shared by PUMA and PlaSim
!     =======================================

      subroutine mpbci(k) ! broadcast 1 integer
      return
      end

      subroutine mpbcin(k,n) ! broadcast n integer
      integer :: k(n)
      return
      end

      subroutine mpbcr(p) ! broadcast 1 real
      return
      end

      subroutine mpbcrn(p,n) ! broadcast n real
      return
      end

      subroutine mpbcl(k) ! broadcast 1 logical
      return
      end

      subroutine mpscin(k,n) ! scatter n integer
      return
      end

      subroutine mpscrn(p,n) ! scatter n real
      return
      end

      subroutine mpscdn(p,n) ! scatter n double precision
      return
      end

      subroutine mpscsp(pf,pp,klev) ! scatter spectral fields
      use pumamod
      real pf(NESP,klev)
      real pp(NSPP,klev)
      pp(1:NSPP,1:klev) = pf(1:NSPP,1:klev)
      return
      end

      subroutine mpscgp(pf,pp,klev) ! scatter gridpoint fields
      use pumamod
      real pf(NLON*NLAT,klev)
      real pp(NHOR,klev)
      pp(1:NHOR,1:klev) = pf(1:NHOR,1:klev)
      return
      end

      subroutine mpgasp(pf,pp,klev) ! gather spectral fields
      use pumamod
      real pf(NESP,klev)
      real pp(NSPP,klev)
      pf(1:NSPP,1:klev) = pp(1:NSPP,1:klev)
      return
      end

      subroutine mpgagp(pf,pp,klev) ! gather gridpoint fields
      use pumamod
      real pf(NHOR,klev)
      real pp(NHOR,klev)
      pf = pp
      return
      end

      subroutine mpgacs(pcs) ! gather cross sections
      return
      end

      subroutine mpgallsp(pf,pp,klev) ! gather spectral to all
      use pumamod
      real pf(NESP,klev)
      real pp(NSPP,klev)
      pf(1:NSPP,1:klev) = pp(1:NSPP,1:klev)
      return
      end

      subroutine mpsum(psp,klev) ! sum spectral fields
      return
      end

      subroutine mpsumsc(psf,psp,klev) ! sum & scatter spectral
      use pumamod
      real psf(NESP,klev)
      real psp(NSPP,klev)
      psp(1:NSPP,1:klev) = psf(1:NSPP,1:klev)
      return
      end

      subroutine mpsumr(pr,kdim) ! sum kdim reals
      return
      end subroutine mpsumr

      subroutine mpsumbcr(pr,kdim) ! sum & broadcast kdim reals
      return
      end

!     ==================
!     SUBROUTINE MPABORT
!     ==================

      subroutine mpabort(ym)

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

      stop
      end


      subroutine mpstart(kworld) ! initialization
      use pumamod
      if (NPRO > 1) then
        write(nud,*)'error : scalar version compiled with NPRO > 1!'
        stop
      endif
      return
      end

      subroutine mpstop
      return
      end

      subroutine mpreadsp(ktape,p,kdim,klev)
      real p(kdim,klev)
      read (ktape) p
      return
      end

      subroutine mpreadgp(ktape,p,kdim,klev)
      real p(kdim,klev)
      read (ktape) p
      return
      end

      subroutine mpwritesp(ktape,p,kdim,klev)
      real p(kdim,klev)
      write (ktape) p
      return
      end

      subroutine mpwritegp(ktape,p,kdim,klev)
      real p(kdim,klev)
      write (ktape) p
      return
      end

      subroutine mpwritegph(ktape,p,kdim,klev,ihead)
      real p(kdim,klev)
!
      real(kind=4) zp(kdim,klev)
!
      integer ihead(8)
      write(ktape) ihead
      zp(:,:)=p(:,:)
      write(ktape) zp

      return
      end


      subroutine ompi_info(nprocess,pid)    ! get nproc and pid
      integer nprocess, pid
      nprocess = 1
      pid = 0
      return
      end subroutine ompi_info


      subroutine mpgetsp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call get_restart_array(yn,p,kdim,kdim,klev)
      return
      end subroutine mpgetsp


      subroutine mpgetgp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call get_restart_array(yn,p,kdim,kdim,klev)
      return
      end subroutine mpgetgp


      subroutine mpputsp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call put_restart_array(yn,p,kdim,kdim,klev)
      return
      end subroutine mpputsp


      subroutine mpputgp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call put_restart_array(yn,p,kdim,kdim,klev)
      return
      end subroutine mpputgp


      subroutine mpsurfgp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call get_surf_array(yn,p,kdim,klev,iread)
      return
      end subroutine mpsurfgp


      subroutine mpmaxval(p,kdim,klev,pmax)
      real :: p(kdim,klev)
      pmax = maxval(p(:,:))
      return
      end subroutine mpmaxval


      subroutine mpsumval(p,kdim,klev,psum)
      real :: p(kdim,klev)
      psum = sum(p(:,:))
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

      subroutine mrdimensions
      return
      end 



