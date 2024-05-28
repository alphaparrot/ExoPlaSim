!****6***0*********0*********0*********0*********0*********0**********72
      subroutine hilo3D(P,im,jm,km,j1,j2,Pmax,Pmin,Qmax,Qmin,bt,bd)
!****6***0*********0*********0*********0*********0*********0**********72
 
      REAL P(IM,JM,km),Pmax(IM,JM,km),Pmin(IM,JM,km),        &
           Qmax(IM,JM,km),Qmin(IM,JM,km),bt(im,jm),bd(im,jm)
 
!MIC$ do all autoscope
!MIC$* shared(P,Qmax,Qmin,im,jm,j1,j2)
!MIC$* private(k,bt,bd)

      DO k=1,km
      call hilo(P(1,1,k),im,jm,j1,j2,Qmax(1,1,k),Qmin(1,1,k),bt,bd)
      enddo
 
      km1 = km-1
      km2 = km-2
 
!MIC$ do all autoscope
!MIC$* shared(Pmax,Pmin,Qmax,Qmin,im,jm,km,km1,km2)
!MIC$* private(i,j)

      DO j=1,jm
      DO i=1,im
! k=1 and k=km
      Pmax(i,j, 1) = max(Qmax(i,j,  2),Qmax(i,j, 1))
      Pmin(i,j, 1) = min(Qmin(i,j,  2),Qmin(i,j, 1))
      Pmax(i,j,km) = max(Qmax(i,j,km1),Qmax(i,j,km))
      Pmin(i,j,km) = min(Qmin(i,j,km1),Qmin(i,j,km))
! k=2 and k=km1
      Pmax(i,j,  2) = max(Qmax(i,j,  3),Pmax(i,j, 1))
      Pmin(i,j,  2) = min(Qmin(i,j,  3),Pmin(i,j, 1))
      Pmax(i,j,km1) = max(Qmax(i,j,km2),Pmax(i,j,km))
      Pmin(i,j,km1) = min(Qmin(i,j,km2),Pmin(i,j,km))
      enddo
      enddo
 
!MIC$ do all autoscope
!MIC$* shared(Pmax,Pmin,Qmax,Qmin,im,jm,km,km1,km2)
!MIC$* private(i,j,k)

      DO k=3,km2
      DO j=1,jm
      DO i=1,im
      Pmax(i,j,k) = max(Qmax(i,j,k-1),Qmax(i,j,k),Qmax(i,j,k+1))
      Pmin(i,j,k) = min(Qmin(i,j,k-1),Qmin(i,j,k),Qmin(i,j,k+1))
      enddo
      enddo
      enddo
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine hilo(q,im,jm,j1,j2,qmax,qmin,bt,bd)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL q(IM,JM),Qmax(IM,JM),Qmin(IM,JM),bt(IM,*),bd(IM,*)
 
! y-sweep
      DO j=j1,j2
      DO i=1,IM
      bt(i,j) = max(q(i,j-1),q(i,j),q(i,j+1))
      bd(i,j) = min(q(i,j-1),q(i,j),q(i,j+1))
      enddo
      enddo
!
! x-sweep
      IM1 = IM-1
      DO j=j1,j2
      DO i=2,IM1
      Qmax(i,j) = max(bt(i-1,j),bt(i,j),bt(i+1,j))
      Qmin(i,j) = min(bd(i-1,j),bd(i,j),bd(i+1,j))
      enddo
      enddo
!
      DO j=j1,j2
!     i = 1
      Qmax(1,j) = max(bt(IM,j),bt(1,j),bt(2,j))
      Qmin(1,j) = min(bd(IM,j),bd(1,j),bd(2,j))
!     i = IM
      Qmax(IM,j) = max(bt(IM1,j),bt(IM,j),bt(1,j))
      Qmin(IM,j) = min(bd(IM1,j),bd(IM,j),bd(1,j))
      enddo
!
! N. Pole:
      Pmax = q(1,JM)
      Pmin = q(1,JM)
      do i=1,IM
      if(q(i,j2) .gt. Pmax) then
            Pmax = q(i,j2)
      elseif(q(i,j2) .lt. Pmin) then
            Pmin = q(i,j2)
      endif
      enddo
!
      do i=1,IM
      Qmax(i,JM) = Pmax
      Qmin(i,JM) = Pmin
      enddo
!
! S. Pole:
      Pmax = q(1,1)
      Pmin = q(1,1)
      do i=1,IM
      if(q(i,j1) .gt. Pmax) then
            Pmax = q(i,j1)
      elseif(q(i,j1) .lt. Pmin) then
            Pmin = q(i,j1)
      endif
      enddo
!
      do i=1,IM
      Qmax(i,1) = Pmax
      Qmin(i,1) = Pmin
      enddo
!
      if(j1 .ne. 2) then
      JM1 = JM-1
      do i=1,IM
      Qmax(i,2) = Qmax(i,1)
      Qmin(i,2) = Qmin(i,1)
!
      Qmax(i,JM1) = Qmax(i,JM)
      Qmin(i,JM1) = Qmin(i,JM)
      enddo
      endif
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine FZPPM(P,fz,IMR,JNP,NL,DQ,WZ,delp,KORD)
!****6***0*********0*********0*********0*********0*********0**********72
 
      PARAMETER ( R23 = 2./3., R3 = 1./3.)
      REAL WZ(IMR,JNP,NL),P(IMR,JNP,NL),DQ(IMR,JNP,NL), &
           fz(IMR,JNP,NL+1),delp(IMR,JNP,NL)
! local 2d arrays
      REAL AR(IMR,NL),AL(IMR,NL),A6(IMR,NL),delq(IMR,NL),DC(IMR,NL)
 
      km = NL
      km1 = NL-1
      LMT = max(KORD - 3, 0)
 
! find global min/max
 
      call vmax1d(P(1,1, 1),Tmax,Tmin,IMR*JNP)
      call vmax1d(P(1,1,NL),Bmax,Bmin,IMR*JNP)

!MIC$ do all autoscope
!MIC$* shared(LMT,Tmax,Tmin,Bmax,Bmin,JNP,IMR)
!MIC$* shared(fz,DQ,WZ)
!MIC$* private(i,j,k,c1,c2,tmp,qmax,qmin,A1,A2,d1,d2,qm,dp,c3)
!MIC$* private(cmax,cmin,DC,delq,AR,AL,A6,CM,CP)

      do j=1,JNP

      do k=2,km
      do i=1,IMR
      A6(i,k) = delp(i,j,k-1) + delp(i,j,k)
      enddo
      enddo

      do k=1,km1
      do i=1,IMR
      delq(i,k) = P(i,j,k+1) - P(i,j,k)
      enddo
      enddo
 
      DO k=2,km1
      DO I=1,IMR
      c1 = (delp(i,j,k-1)+0.5*delp(i,j,k))/A6(i,k+1)
      c2 = (delp(i,j,k+1)+0.5*delp(i,j,k))/A6(i,k)
      tmp = delp(i,j,k)*(c1*delq(i,k) + c2*delq(i,k-1)) &
           / (A6(i,k)+delp(i,j,k+1))
      Qmax = max(P(i,j,k-1),P(i,j,k),P(i,j,k+1)) - P(i,j,k)
      Qmin = P(i,j,k) - min(P(i,j,k-1),P(i,j,k),P(i,j,k+1))
      DC(i,k) = sign(min(abs(tmp),Qmax,Qmin), tmp)
      enddo
      enddo
 
!****6***0*********0*********0*********0*********0*********0**********72
! Compute the first guess at cell interface
! First guesses are required to be continuous.
!****6***0*********0*********0*********0*********0*********0**********72
 
! Interior.
 
      DO k=3,km1
      DO i=1,IMR
      c1 = delq(i,k-1)*delp(i,j,k-1) / A6(i,k)
      A1 = A6(i,k-1) / (A6(i,k) + delp(i,j,k-1))
      A2 = A6(i,k+1) / (A6(i,k) + delp(i,j,k))
      AL(i,k) = P(i,j,k-1) + c1 + 2./(A6(i,k-1)+A6(i,k+1)) *  &
                ( delp(i,j,k  )*(c1*(A1 - A2)+A2*DC(i,k-1)) - &
                                delp(i,j,k-1)*A1*DC(i,k  ) )
      enddo
      enddo
 
! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      DO i=1,IMR
      d1 = delp(i,j,1)
      d2 = delp(i,j,2)
      qm = (d2*P(i,j,1)+d1*P(i,j,2)) / (d1+d2)
      dp = 2.*(P(i,j,2)-P(i,j,1)) / (d1+d2)
      c1 = 4.*(AL(i,3)-qm-d2*dp) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
      c3 = dp - 0.5*c1*(d2*(5.*d1+d2)-3.*d1**2)
      AL(i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
      AL(i,1) = d1*(2.*c1*d1**2-c3) + AL(i,2)
      DC(i,1) =  P(i,j,1) - AL(i,1)
! No over- and undershoot condition
      AL(i,1) = max(Tmin,AL(i,1))
      AL(i,1) = min(Tmax,AL(i,1))
      Cmax = max(P(i,j,1), P(i,j,2))
      Cmin = min(P(i,j,1), P(i,j,2))
      AL(i,2) = max(Cmin,AL(i,2))
      AL(i,2) = min(Cmax,AL(i,2))
      enddo
      
 
! Bottom
      DO i=1,IMR
      d1 = delp(i,j,km )
      d2 = delp(i,j,km1)
      qm = (d2*P(i,j,km)+d1*P(i,j,km1)) / (d1+d2)
      dp = 2.*(P(i,j,km1)-P(i,j,km)) / (d1+d2)
	c1 = 4.*(AL(i,km1)-qm-d2*dp) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
	c3 = dp - 0.5*c1*(d2*(5.*d1+d2)-3.*d1**2)
      AL(i,km) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
      AR(i,km) = d1*(2.*c1*d1**2-c3) + AL(i,km)
      DC(i,km) = AR(i,km) -  P(i,j,km)
! No over- and undershoot condition
      Cmax = max(P(i,j,km), P(i,j,km1))
      Cmin = min(P(i,j,km), P(i,j,km1))
      AL(i,km) = max(Cmin,AL(i,km))
      AL(i,km) = min(Cmax,AL(i,km))
      AR(i,km) = max(Bmin,AR(i,km))
      AR(i,km) = min(Bmax,AR(i,km))
      enddo

      do k=1,km1
      do i=1,IMR
      AR(i,k) = AL(i,k+1)
      enddo
      enddo
 
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
 
      do k=1,NL
      do i=1,IMR
      A6(i,k) = 3.*(P(i,j,k)+P(i,j,k) - (AL(i,k)+AR(i,k)))
      enddo
      enddo
 
! Top
      call lmtppm(DC(1,1),A6(1,1),AR(1,1),AL(1,1),P(1,j,1), &
                  IMR,0)
! Bottom
      call lmtppm(DC(1,NL),A6(1,NL),AR(1,NL),AL(1,NL), &
                  P(1,j,NL),IMR,0)
! Interior.
      if(LMT.LE.2) then
      do k=2,NL-1
      call lmtppm(DC(1,k),A6(1,k),AR(1,k),AL(1,k),P(1,j,k), &
                  IMR,LMT)
      enddo
      endif
 
      DO k=2,NL
      DO i=1,IMR
      IF(WZ(i,j,k-1).GT.0.) then
             CM = WZ(i,j,k-1) / delp(i,j,k-1)
      DC(i,k) = P(i,j,k-1)
      fz(i,j,k) = AR(i,k-1)+0.5*CM*(AL(i,k-1)-AR(i,k-1)+ &
                            A6(i,k-1)*(1.-R23*CM))
      else
             CP = WZ(i,j,k-1) / delp(i,j,k)
      DC(i,k) = P(i,j,k)
      fz(i,j,k) = AL(i,k)+0.5*CP*(AL(i,k)-AR(i,k)- &
                          A6(i,k)*(1.+R23*CP))
      endif
      enddo
      enddo
 
      DO k=2,NL
      DO i=1,IMR
      fz(i,j,k) = WZ(i,j,k-1) * (fz(i,j,k) - DC(i,k))
      DC(i,k) = WZ(i,j,k-1) * DC(i,k)
      enddo
      enddo

      do i=1,IMR
      fz(i,j,   1) = 0.
      fz(i,j,NL+1) = 0.
      DQ(i,j, 1) = DQ(i,j, 1) - DC(i, 2)
      DQ(i,j,NL) = DQ(i,j,NL) + DC(i,NL)
      enddo
 
      do k=2,km1
      do i=1,IMR
      DQ(i,j,k) = DQ(i,j,k) + DC(i,k) - DC(i,k+1)
      enddo
      enddo
 
      enddo
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine xtp(IMR,JNP,IML,j1,j2,JN,JS,PU,DQ,Q,C,fx2,xmass,IORD)
!****6***0*********0*********0*********0*********0*********0**********72

      REAL C(IMR,*),DC(-IML:IMR+IML+1),xmass(IMR,JNP), &
           fx1(IMR+1),DQ(IMR,JNP),qtmp(-IML:IMR+1+IML)
      REAL PU(IMR,JNP),Q(IMR,JNP),fx2(IMR+1,JNP)
      INTEGER ISAVE(IMR)
!
      IMP = IMR + 1
!
! van Leer at high latitudes
      jvan = max(1,JNP/20)
      j1vl = j1+jvan
      j2vl = j2-jvan
!
      do j=j1,j2
!
      do i=1,IMR
      qtmp(i) = q(i,j)
      enddo
!
      if(j.ge.JN .or. j.le.JS) then
  !
!****6***0*********0*********0*********0*********0*********0**********72
! *** Conservative (flux-form) Semi-Lagrangian transport ***
!****6***0*********0*********0*********0*********0*********0**********72
!    
        do i=-IML,0
        qtmp(i)     = q(IMR+i,j)
        qtmp(IMP-i) = q(1-i,j)
        enddo
    !
        IF(IORD.eq.1 .or. j.eq.j1 .or. j.eq.j2) THEN
            DO i=1,IMR
            itmp = INT(c(i,j))
            ISAVE(i) = i - itmp
            iu = i - c(i,j)
            fx1(i) = (c(i,j) - itmp)*qtmp(iu)
            ENDDO
        !
        ! Zero high order contribution
            DO i=1,IMR
            fx2(i,j) = 0.
            enddo
        ELSE
            call xmist(IMR,IML,Qtmp,DC)
        !
            do i=-IML,0
            DC(i)     = DC(IMR+i)
            DC(IMP-i) = DC(1-i)
            enddo
        !
            DO i=1,IMR
            itmp = INT(c(i,j))
            rut  = c(i,j) - itmp
            ISAVE(i) = i - itmp
            iu = i - c(i,j)
            fx1(i  ) = rut*qtmp(iu)
            fx2(i,j) = rut*DC(iu)*(sign(1.,rut) - rut)
            ENDDO
        ENDIF
    !
        do i=1,IMR
        IF(c(i,j).GT.1.) then
    !DIR$ NOVECTOR
            do ist = ISAVE(i),i-1
            fx1(i) = fx1(i) + qtmp(ist)
            enddo
        elseIF(c(i,j).LT.-1.) then
    !DIR$ NOVECTOR
            do ist = i,ISAVE(i)-1
            fx1(i) = fx1(i) - qtmp(ist)
            enddo
        endif
        enddo
    !DIR$ VECTOR
        do i=1,IMR
        fx1(i)   = PU(i,j)*fx1(i)
        fx2(i,j) = PU(i,j)*fx2(i,j)
        enddo
      
      else
!****6***0*********0*********0*********0*********0*********0**********72
! *** Eulerian ***
!****6***0*********0*********0*********0*********0*********0**********72
!
        qtmp(0)     = q(IMR,J)
        qtmp(-1)    = q(IMR-1,J)
        qtmp(IMP)   = q(1,J)
        qtmp(IMP+1) = q(2,J)
    !
        IF(IORD.eq.1 .or. j.eq.j1 .or. j.eq.j2) THEN
            DO i=1,IMR
            iu = float(i) - c(i,j)
            fx1(i) = qtmp(iu)
            ENDDO
        !
        ! Zero high order contribution
            DO i=1,IMR
            fx2(i,j) = 0.
            enddo
        ELSE
            call xmist(IMR,IML,Qtmp,DC)
            DC(0) = DC(IMR)
        !
            if(IORD.eq.2 .or. j.le.j1vl .or. j.ge.j2vl) then
                DO i=1,IMR
                iu = float(i) - c(i,j)
                fx1(i  ) = qtmp(iu)
                fx2(i,j) = DC(iu)*(sign(1.,c(i,j))-c(i,j))
                ENDDO
            else
                call fxppm(IMR,IML,C(1,j),Qtmp,DC,fx1,fx2(1,j),IORD)
            endif
    !
        ENDIF
    !
        DO i=1,IMR
        fx1(i  ) = fx1(i  )*xmass(i,j)
        fx2(i,j) = fx2(i,j)*xmass(i,j)
        ENDDO
!
      endif
!
!
!****6***0*********0*********0*********0*********0*********0**********72
!
      fx1(IMP  ) = fx1(1  )
      fx2(IMP,j) = fx2(1,j)
! Update using low order fluxes.
      DO i=1,IMR
      DQ(i,j) =  DQ(i,j) + fx1(i)-fx1(i+1)
      ENDDO
!
!****6***0*********0*********0*********0*********0*********0**********72
!
      enddo
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine fxppm(IMR,IML,UT,P,DC,fx1,fx2,IORD)
!****6***0*********0*********0*********0*********0*********0**********72
      PARAMETER ( R3 = 1./3., R23 = 2./3. )
      REAL UT(*),fx1(*),P(-IML:IMR+IML+1),DC(-IML:IMR+IML+1)
      REAL AR(0:IMR),AL(0:IMR),A6(0:IMR),fx2(*)
 
      LMT = IORD - 3
 
      DO i=1,IMR
      AL(i) = 0.5*(p(i-1)+p(i)) + (DC(i-1) - DC(i))*R3
      ENDDO
 
      do i=1,IMR-1
      AR(i) = AL(i+1)
      enddo
      AR(IMR) = AL(1)
 
      do i=1,IMR
      A6(i) = 3.*(p(i)+p(i)  - (AL(i)+AR(i)))
      enddo
 
      if(LMT.LE.2) call lmtppm(DC(1),A6(1),AR(1),AL(1),P(1),IMR,LMT)
 
      AL(0) = AL(IMR)
      AR(0) = AR(IMR)
      A6(0) = A6(IMR)
 
! Abs(UT(i)) < 1
      DO i=1,IMR
      IF(UT(i).GT.0.) then
      fx1(i) = P(i-1)
      fx2(i) = AR(i-1) + 0.5*UT(i)*(AL(i-1) - AR(i-1) + &
                             A6(i-1)*(1.-R23*UT(i)) )
      else
      fx1(i) = P(i)
      fx2(i) = AL(i) - 0.5*UT(i)*(AR(i) - AL(i) + &
                           A6(i)*(1.+R23*UT(i)))
      endif
      enddo
 
      DO i=1,IMR
      fx2(i) = fx2(i) - fx1(i)
      enddo
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine fyppm(C,P,DC,fy1,fy2,IMR,JNP,j1,j2,A6,AR,AL,JORD)
!****6***0*********0*********0*********0*********0*********0**********72
      PARAMETER ( R3 = 1./3., R23 = 2./3. )
      REAL C(IMR,*),fy1(IMR,*),P(IMR,*),DC(IMR,*),fy2(IMR,JNP)
      REAL AR(IMR,JNP),AL(IMR,JNP),A6(IMR,JNP)
 
      IMH = IMR / 2
      JMR = JNP - 1
      j11 = j1-1
      IMJM1 = IMR*(J2-J1+2)
      len   = IMR*(J2-J1+3)
      LMT = JORD - 3
 
      DO i=1,IMR*JMR
      AL(i,2) = 0.5*(p(i,1)+p(i,2)) + (DC(i,1) - DC(i,2))*R3
      AR(i,1) = AL(i,2)
      ENDDO
!
! Poles:
!
      DO i=1,IMH
      AL(i,1) = AL(i+IMH,2)
      AL(i+IMH,1) = AL(i,2)
!
      AR(i,JNP) = AR(i+IMH,JMR)
      AR(i+IMH,JNP) = AR(i,JMR)
      enddo
!
      do i=1,len
      A6(i,j11) = 3.*(p(i,j11)+p(i,j11)  - (AL(i,j11)+AR(i,j11)))
      enddo
 
      if(LMT.le.2) call lmtppm(DC(1,j11),A6(1,j11),AR(1,j11), &
                               AL(1,j11),P(1,j11),len,LMT)
 
      DO i=1,IMJM1
      IF(C(i,j1).GT.0.) then
      fy1(i,j1) = P(i,j11)
      fy2(i,j1) = AR(i,j11) + 0.5*C(i,j1)*(AL(i,j11) - AR(i,j11) + &
                               A6(i,j11)*(1.-R23*C(i,j1)) )
      else
      fy1(i,j1) = P(i,j1)
      fy2(i,j1) = AL(i,j1) - 0.5*C(i,j1)*(AR(i,j1) - AL(i,j1) + &
                              A6(i,j1)*(1.+R23*C(i,j1)))
      endif
      ENDDO
 
      DO i=1,IMJM1
      fy2(i,j1) = fy2(i,j1) - fy1(i,j1)
      ENDDO
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine xadv(IMR,JNP,j1,j2,p,UA,JS,JN,IML,adx)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL p(IMR,JNP),adx(IMR,JNP),qtmp(-IMR:IMR+IMR),UA(IMR,JNP)
 
      JMR = JNP-1
      outer: do j=j1,j2
      if(J.GT.JS  .and. J.LT.JN) EXIT outer
 
      do i=1,IMR
      qtmp(i) = p(i,j)
      enddo
 
      do i=-IML,0
      qtmp(i)       = p(IMR+i,j)
      qtmp(IMR+1-i) = p(1-i,j)
      enddo
 
      DO i=1,IMR
      iu = UA(i,j)
      ru = UA(i,j) - iu
      iiu = i-iu
      if(UA(i,j).GE.0.) then
      adx(i,j) = qtmp(iiu)+ru*(qtmp(iiu-1)-qtmp(iiu))
      else
      adx(i,j) = qtmp(iiu)+ru*(qtmp(iiu)-qtmp(iiu+1))
      endif
      enddo
 
      do i=1,IMR
      adx(i,j) = adx(i,j) - p(i,j)
      enddo
      enddo outer
 
! Eulerian upwind
 
      do j=JS+1,JN-1
 
      do i=1,IMR
      qtmp(i) = p(i,j)
      enddo
 
      qtmp(0)     = p(IMR,J)
      qtmp(IMR+1) = p(1,J)
 
      DO i=1,IMR
      IP = i - UA(i,j)
      adx(i,j) = UA(i,j)*(qtmp(ip)-qtmp(ip+1))
      enddo
      enddo
 
      if(j1.ne.2) then
      do i=1,IMR
      adx(i,  2) = 0.
      adx(i,JMR) = 0.
      enddo
      endif

! set cross term due to x-adv at the poles to zero.
      do i=1,IMR
      adx(i,  1) = 0.
      adx(i,JNP) = 0.
      enddo
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine xmist(IMR,IML,P,DC)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL P(-IML:IMR+1+IML),DC(-IML:IMR+1+IML)
!
! 2nd order version.
!
      do i=1,IMR
      tmp = 0.25*(p(i+1) - p(i-1))
      Pmax = max(P(i-1), p(i), p(i+1)) - p(i)
      Pmin = p(i) - min(P(i-1), p(i), p(i+1))
      DC(i) = sign(min(abs(tmp),Pmax,Pmin), tmp)
      enddo
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine lmtppm(DC,A6,AR,AL,P,IM,LMT)
!****6***0*********0*********0*********0*********0*********0**********72
!
! A6 =  CURVATURE OF THE TEST PARABOLA
! AR =  RIGHT EDGE VALUE OF THE TEST PARABOLA
! AL =  LEFT  EDGE VALUE OF THE TEST PARABOLA
! DC =  0.5 * MISMATCH
! P  =  CELL-AVERAGED VALUE
! IM =  VECTOR LENGTH
!
! OPTIONS:
!
! LMT = 0: FULL MONOTONICITY
! LMT = 1: SEMI-MONOTONIC CONSTRAINT (NO UNDERSHOOTS)
! LMT = 2: POSITIVE-DEFINITE CONSTRAINT
!
      PARAMETER ( R12 = 1./12. )
      REAL A6(IM),AR(IM),AL(IM),P(IM),DC(IM)
 
      if(LMT.eq.0) then
! Full constraint
      do i=1,IM
      if(DC(i).eq.0.) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      else
      da1  = AR(i) - AL(i)
      da2  = da1**2
      A6DA = A6(i)*da1
      if(A6DA .lt. -da2) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      elseif(A6DA .gt. da2) then
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
      endif
      enddo

      else if(LMT.eq.1) then
! Semi-monotonic constraint
      outer1: do i=1,IM
      if(abs(AR(i)-AL(i)) .GE. -A6(i)) exit outer1
      if(p(i).lt.AR(i) .and. p(i).lt.AL(i)) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      elseif(AR(i) .gt. AL(i)) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      else
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
      enddo outer1
      
      else if(LMT.eq.2) then
      outer2: do i=1,IM
      if(abs(AR(i)-AL(i)) .GE. -A6(i)) exit outer2
      fmin = p(i) + 0.25*(AR(i)-AL(i))**2/A6(i) + A6(i)*R12
      if(fmin.ge.0.) exit outer2
      if(p(i).lt.AR(i) .and. p(i).lt.AL(i)) then
            AR(i) = p(i)
            AL(i) = p(i)
            A6(i) = 0.
      elseif(AR(i) .gt. AL(i)) then
            A6(i) = 3.*(AL(i)-p(i))
            AR(i) = AL(i) - A6(i)
      else
            A6(i) = 3.*(AR(i)-p(i))
            AL(i) = AR(i) - A6(i)
      endif
      enddo outer2
      endif
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine cosa(cosp,cose,JM,PI,DP)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL cosp(*),cose(*),sine(JM)
 
      do j=2,JM
         ph5  = -0.5*PI + (float(j-1)-0.5)*DP
         sine(j) = SIN(ph5)
      enddo
 
      do J=2,JM-1
         cosp(J) = (sine(j+1)-sine(j))/DP
      enddo
 
      cosp( 1) = 0.
      cosp(JM) = 0.
 
! Define cosine at edges..

      do j=2,JM
         cose(j) = 0.5 * (cosp(j-1)+cosp(j))
      enddo
      cose(1) = cose(2)
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine cosc(cosp,cose,JNP,PI,DP)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL cosp(*),cose(*)
      phi = -0.5*PI
      do j=2,JNP-1
        phi  =  phi + DP
        cosp(j) = cos(phi)
      enddo
      cosp(  1) = 0.
      cosp(JNP) = 0.
 
      do j=2,JNP
        cose(j) = 0.5*(cosp(j)+cosp(j-1))
      enddo
 
      do j=2,JNP-1
        cosp(j) = 0.5*(cose(j)+cose(j+1))
      enddo
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine vmax1d(a,Pmax,Pmin,n)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL a(*)
      Pmax = a(n)
      Pmin = a(n)
      do i=1,n-1
        Pmax = amax1(Pmax,a(i))
        Pmin = amin1(Pmin,a(i))
      enddo
      return
      end

!****6***0*********0*********0*********0*********0*********0**********72
      SUBROUTINE qckxyz(Q,qtmp,IMR,JNP,NLAY,j1,j2,cosp,acosp,IC)
!****6***0*********0*********0*********0*********0*********0**********72
 
      PARAMETER ( tiny = 1.E-30 )
      PARAMETER ( kmax = 200 )
      REAL Q(IMR,JNP,NLAY),qtmp(IMR,JNP),cosp(*),acosp(*)
      integer IP(kmax)
 
      NLM1 = NLAY-1

! Do horizontal filling.

!MIC$ do all autoscope
!MIC$* private(i,j,L,qtmp)

      do L=1,NLAY
         call filns(q(1,1,L),IMR,JNP,j1,j2,cosp,acosp,ip(L),tiny)
         if(ip(L).ne.0) &
         call filew(q(1,1,L),qtmp,IMR,JNP,j1,j2,ip(L),tiny)
     enddo

      ipz = 0
      do L=1,NLAY
      if(ip(L) .ne. 0) then
         ipz = L
         go to 111
      endif
      enddo
      return

111   continue

      if(ipz .eq. 0) return

      if(ipz .eq. 1) then
         lpz = 2
      else
         lpz = ipz
      endif

! Do vertical filling.

!MIC$ do all autoscope
!MIC$* private(i,j,L,qup,qly,dup)

      do j=j1,j2

      if(ipz .eq. 1) then
! Top layer
      do i=1,IMR
      if(Q(i,j,1).LT.0.) then
          Q(i,j,2) = Q(i,j,2) + Q(i,j,1)
          Q(i,j,1) = 0.
      endif
      enddo
      endif
 
      DO L = lpz,NLM1
      do i=1,IMR
      IF( Q(i,j,L).LT.0.) THEN
! From above
          qup =  Q(i,j,L-1)
          qly = -Q(i,j,L)
          dup  = min(qly,qup)
          Q(i,j,L-1) = qup - dup
          Q(i,j,L  ) = dup-qly
! Below
          Q(i,j,L+1) = Q(i,j,L+1) + Q(i,j,L)
          Q(i,j,L)   = 0.
      ENDIF
      enddo
      enddo
 
! BOTTOM LAYER
      L = NLAY
      do i=1,IMR
      IF( Q(i,j,L).LT.0.) THEN
 
! From above
 
          qup = Q(i,j,NLM1)
          qly = -Q(i,j,L)
          dup = min(qly,qup)
          Q(i,j,NLM1) = qup - dup

! From "below" the surface.
          Q(i,j,L) = 0.
      ENDIF
      enddo
      enddo

      RETURN
      END
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine filns(q,IMR,JNP,j1,j2,cosp,acosp,ipy,tiny)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL q(IMR,*),cosp(*),acosp(*)
      LOGICAL first
      DATA first /.true./
      SAVE cap1
 
      if(first) then
      DP = 4.*ATAN(1.)/float(JNP-1)
      cap1 = IMR*(1.-COS((j1-1.5)*DP))/DP
      first = .false.
      endif
 
      ipy = 0
      do j=j1+1,j2-1
      DO i=1,IMR
      IF(q(i,j).LT.0.) THEN
      ipy =  1
      dq  = - q(i,j)*cosp(j)
! North
      dn = q(i,j+1)*cosp(j+1)
      d0 = max(0.,dn)
      d1 = min(dq,d0)
      q(i,j+1) = (dn - d1)*acosp(j+1)
      dq = dq - d1
! South
      ds = q(i,j-1)*cosp(j-1)
      d0 = max(0.,ds)
      d2 = min(dq,d0)
      q(i,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      endif
      enddo
      enddo
 
      do i=1,imr
      IF(q(i,j1).LT.0.) THEN
      ipy =  1
      dq  = - q(i,j1)*cosp(j1)
! North
      dn = q(i,j1+1)*cosp(j1+1)
      d0 = max(0.,dn)
      d1 = min(dq,d0)
      q(i,j1+1) = (dn - d1)*acosp(j1+1)
      q(i,j1) = (d1 - dq)*acosp(j1) + tiny
      endif
      enddo
 
      j = j2
      do i=1,imr
      IF(q(i,j).LT.0.) THEN
      ipy =  1
      dq  = - q(i,j)*cosp(j)
! South
      ds = q(i,j-1)*cosp(j-1)
      d0 = max(0.,ds)
      d2 = min(dq,d0)
      q(i,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      endif
      enddo
 
! Check Poles.
      if(q(1,1).lt.0.) then
      dq = q(1,1)*cap1/float(IMR)*acosp(j1)
      do i=1,imr
      q(i,1) = 0.
      q(i,j1) = q(i,j1) + dq
      if(q(i,j1).lt.0.) ipy = 1
      enddo
      endif
 
      if(q(1,JNP).lt.0.) then
      dq = q(1,JNP)*cap1/float(IMR)*acosp(j2)
      do i=1,imr
      q(i,JNP) = 0.
      q(i,j2) = q(i,j2) + dq
      if(q(i,j2).lt.0.) ipy = 1
      enddo
      endif
 
      return
      end
 
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine filew(q,qtmp,IMR,JNP,j1,j2,ipx,tiny)
!****6***0*********0*********0*********0*********0*********0**********72
      REAL q(IMR,*),qtmp(JNP,IMR)
 
      ipx = 0
! Copy & swap direction for vectorization.
      do i=1,imr
      do j=j1,j2
      qtmp(j,i) = q(i,j)
      enddo
      enddo
 
      do i=2,imr-1
      do j=j1,j2
      if(qtmp(j,i).lt.0.) then
      ipx =  1
! west
      d0 = max(0.,qtmp(j,i-1))
      d1 = min(-qtmp(j,i),d0)
      qtmp(j,i-1) = qtmp(j,i-1) - d1
      qtmp(j,i) = qtmp(j,i) + d1
! east
      d0 = max(0.,qtmp(j,i+1))
      d2 = min(-qtmp(j,i),d0)
      qtmp(j,i+1) = qtmp(j,i+1) - d2
      qtmp(j,i) = qtmp(j,i) + d2 + tiny
      endif
      enddo
      enddo
 
      i=1
      do j=j1,j2
      if(qtmp(j,i).lt.0.) then
      ipx =  1
! west
      d0 = max(0.,qtmp(j,imr))
      d1 = min(-qtmp(j,i),d0)
      qtmp(j,imr) = qtmp(j,imr) - d1
      qtmp(j,i) = qtmp(j,i) + d1
! east
      d0 = max(0.,qtmp(j,i+1))
      d2 = min(-qtmp(j,i),d0)
      qtmp(j,i+1) = qtmp(j,i+1) - d2
 
      qtmp(j,i) = qtmp(j,i) + d2 + tiny
      endif
      enddo
      
      i=IMR
      do j=j1,j2
      if(qtmp(j,i).lt.0.) then
      ipx =  1
! west
      d0 = max(0.,qtmp(j,i-1))
      d1 = min(-qtmp(j,i),d0)
      qtmp(j,i-1) = qtmp(j,i-1) - d1
      qtmp(j,i) = qtmp(j,i) + d1
! east
      d0 = max(0.,qtmp(j,1))
      d2 = min(-qtmp(j,i),d0)
      qtmp(j,1) = qtmp(j,1) - d2
 
      qtmp(j,i) = qtmp(j,i) + d2 + tiny
      endif
      enddo
 
      if(ipx.ne.0) then
      do j=j1,j2
      do i=1,imr
      q(i,j) = qtmp(j,i)
      enddo
      enddo
      else

! Pole
      if(q(1,1).lt.0. .or. q(1,JNP).lt.0.) ipx = 1
      endif
      return
      end    
