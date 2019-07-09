* Source file MATERIAL.FOR |||||||||||||||||||||||||||||||||||||||||||||

*     iModel = 0: van Genuchten
*              1: modified van Genuchten (Vogel and Cislerova)
*              2: Brooks and Corey
*              3: van Genuchte with air entry value of 2 cm
*              4: log-normal (Kosugi)
*              5: dual-porosity (Durner)
*             10: fractal model (Shlomo Orr)
*     lAltern = VG model with alfa and n different for retention curve 
*               and hydraulic conductivity function. (iModel=1)
*               One needs to comment out line 250 in input2.for
*     lUnBound = Unbound n and m in VG-M function.
*                one needs to uncheck m=Par(6) in other routines as well.

      real function FK(iModel,h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,Ks,Kr,Kk,Lambda,n2,m2,mn,mt
      real h,Par(10)
      integer PPar
      logical lAltern,lUnBound

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      Ks=amax1(Par(5),1.e-37)
      BPar=Par(6)
      if(iModel.eq.0.or.iModel.eq.1.or.iModel.eq.3) then    ! VG and modified VG
*        BPar=.5d0
        PPar=2
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
          Qk=Qs
          Kk=Ks
        else if(iModel.eq.1) then
          lAltern=.false.  
          if(.not.lAltern) then
            Qm=Par(7)
            Qa=Par(8)
            Qk=Par(9)
            Kk=Par(10)
          else
            Qm=Qs
            Qa=Qr
            Qk=Qs
            Kk=Ks
            Alfa=Par(9)
            n=Par(10)
          end if
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
        lUnBound=.false.
        if(lUnBound) then
          m=Par(6)
          BPar=0.5d0
        end if
        HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
        HH=max(dble(h),HMin)
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Qeek=dmin1((Qk-Qa)/(Qm-Qa),Qees)
        Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
        Hk=-1.d0/Alfa*(Qeek**(-1.d0/m)-1.d0)**(1.d0/n)
        if(dble(h).lt.Hk) then
          if(.not.lUnBound) then            ! m=1-1/n
            Qee=(1.d0+(-Alfa*HH)**n)**(-m)
            Qe =(Qm-Qa)/(Qs-Qa)*Qee
            Qek=(Qm-Qa)/(Qs-Qa)*Qeek
            FFQ =1.d0-(1.d0-Qee **(1.d0/m))**m
            FFQk=1.d0-(1.d0-Qeek**(1.d0/m))**m
            if(FFQ.le.0.d0) FFQ=m*Qee**(1.d0/m)
            Kr=(Qe/Qek)**Bpar*(FFQ/FFQk)**PPar*Kk/Ks
            if(iModel.eq.0) Kr=Qe**Bpar*(FFQ)**PPar
*           Gardner's model {K=Ks*[exp(-a*h)]}
c            Kr=dexp(-BPar*dble(h))
            FK=sngl(max(Ks*Kr,1.d-37))
          else                              ! unbounded n and m
            mn=m*n
            mt=1.d0
*           Calculate complete Beta function
            AA=m+mt/n
            BB=1.d0-mt/n
            if(BB.le.0.004d0) then
              write(*,*) 1010
1010          format(/5x,'Parameter N is too small, not executed!')
              stop
            end if
            Beta=Gamma(AA)*Gamma(BB)/Gamma(m+1.)
            WCL=dmax1(2.d0/(2.d0+m),0.2d0)
            dlgKs=dlog10(Ks)
*           Water content
            AX=Alfa*(-h)
            if(AX.lt.1.d-20) then
              Qe=1.0d0
            else
              EX=n*dlog10(AX)
              if(EX.lt.-10.d0) then
                Qe=1.0d0
              else if(EX.lt.10.d0) then
                Qe=(1.+AX**n)**(-m)
              else 
                EX=m*EX
                if(EX.lt.30.d0) then
                  Qe=AX**(-m*n)
                else
                  Qe=0.0d0
                end if
              end if
            end if
*           Conductivity
            if(Qe.le.1.d-10) then
              FK=1.d-37
            else if(Qe.gt.0.999999d0) then
              FK=Ks
            else
              dlgW=dlog10(Qe)
              dlg2=3.0d0-mt+BPar+2.0d0/mn
              dlgC=dlg2*dlgW+dlgKs
              if(dlgC.gt.-37.d0.and.dlgW.gt.(-15.d0*m)) then
                dw=Qe**(1.d0/m)
                if(dw.lt.1.d-06) then
                  dlg1=(3.0d0-mt)*dlog10(n/(Beta*(mn+mt)))
                  dlgC=dlgC+dlg1
                  FK=10.**dlgC
                  return
                end if
                if(Qe-WCL.le.0.d0) then
                  Term=BInc(dw,AA,BB,Beta)
                else
                  Term=1.d0-BInc(1.d0-dw,BB,AA,Beta)
                end if
                Kr=Qe**BPar*Term
                if(mt.lt.1.5d0) Kr=Kr*Term
                dlgC=dlog10(Kr)+dlgKs
              end if
              dlgC=dmax1(-37.d0,dlgC)
              FK=10.**dlgC
            end if
          end if
        end if
        if(dble(h).ge.Hk.and.dble(h).lt.Hs) then
          Kr=(1.d0-Kk/Ks)/(Hs-Hk)*(dble(h)-Hs)+1.d0
          FK=sngl(Ks*Kr)
        end if
        if(dble(h).ge.Hs) FK=sngl(Ks)
      else if(iModel.eq.2) then                  ! Brooks and Cores
*        BPar=1.d0
        Lambda=2.d0  !  !=2 for Mualem Model, =1.5 for Burdine model
        Hs=-1.d0/Alfa
        if(h.lt.Hs) then
          Kr=1.d0/(-Alfa*h)**(n*(BPar+Lambda)+2.d0)
          FK=sngl(max(Ks*Kr,1.d-37))
        else
          FK=sngl(Ks)
        end if
      else if(iModel.eq.4) then                  ! Log-normal model
        Hs=0.d0
        if(h.lt.Hs) then
          Qee=qnorm(dlog(-h/Alfa)/n)
          t=qnorm(dlog(-h/Alfa)/n+n)
          Kr=Qee**Bpar*t*t
          FK=sngl(max(Ks*Kr,1.d-37))
        else
          FK=sngl(Ks)
        end if
      else if(iModel.eq.5) then                  ! Dual-porosity model
        w2=Par(7)
        Alfa2=Par(8)
        n2=Par(9)
        m =1.d0-1.d0/n
        m2=1.d0-1.d0/n2
        w1=1.d0-w2
        Sw1=w1*(1.d0+(-Alfa *h)**n )**(-m )
        Sw2=w2*(1.d0+(-Alfa2*h)**n2)**(-m2)
        Qe=Sw1+Sw2
        Sv1=(-Alfa *h)**(n -1)
        Sv2=(-Alfa2*h)**(n2-1)
        Sk1=w1*Alfa *(1.d0-Sv1*(1.d0+(-Alfa *h)**n )**(-m ))
        Sk2=w2*Alfa2*(1.d0-Sv2*(1.d0+(-Alfa2*h)**n2)**(-m2))
        rNumer=Sk1+Sk2
        rDenom=w1*Alfa+w2*Alfa2
        if(rDenom.ne.0.) Kr=Qe**BPar*(rNumer/rDenom)**2
        FK=sngl(max(Ks*Kr,1.d-37))
      else if(iModel.eq.-1) then                  ! Shlomo Orr
        ha=Alfa
        D=n
        Kr=1.
        if(-h.gt.ha)
     !    Kr=(1.-(1.-(-ha/h)**(3.-D))/(1.-Qr))**(D/(3.-D))
        FK=sngl(max(Ks*Kr,1.d-37))
      end if

      return
      end

************************************************************************

      real function FC(iModel,h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,n2,m2
      real h,Par(10)

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      if(iModel.eq.0.or.iModel.eq.1.or.iModel.eq.3) then
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
        else if(iModel.eq.1) then
          Qm=Par(7)
          Qa=Par(8)
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
c        m=Par(6)
        HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
        HH=max(dble(h),HMin)
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
        if(dble(h).lt.Hs) then
          C1=(1.d0+(-Alfa*HH)**n)**(-m-1.d0)
          C2=(Qm-Qa)*m*n*(Alfa**n)*(-HH)**(n-1.d0)*C1
          FC=sngl(max(C2,1.d-37))
          return
        else
          FC=0.0
        end if
      else if(iModel.eq.2) then
        Hs=-1.d0/Alfa
        if(h.lt.Hs) then
          C2=(Qs-Qr)*n*Alfa**(-n)*(-h)**(-n-1.d0)
          FC=sngl(max(C2,1.d-37))
        else
          FC=0.0
        end if
      else if(iModel.eq.4) then
        Hs=0.d0
        if(h.lt.Hs) then
          t=exp(-1.d0*(dlog(-h/Alfa))**2.d0/(2.d0*n**2.d0))
          C2=(Qs-Qr)/(2.d0*3.141592654)**0.5d0/n/(-h)*t
          FC=sngl(max(C2,1.d-37))
        else
          FC=0.0
        end if
      else if(iModel.eq.5) then
        w2=Par(7)
        Alfa2=Par(8)
        n2=Par(9)
        m =1.d0-1.d0/n
        m2=1.d0-1.d0/n2
        w1=1.d0-w2
        C1a=(1.d0+(-Alfa *h)**n )**(-m -1.d0)
        C1b=(1.d0+(-Alfa2*h)**n2)**(-m2-1.d0)
        C2a=(Qs-Qr)*m *n *(Alfa **n )*(-h)**(n -1.d0)*C1a*w1
        C2b=(Qs-Qr)*m2*n2*(Alfa2**n2)*(-h)**(n2-1.d0)*C1b*w2
        FC=C2a+C2b
      else if(iModel.eq.-1) then                  ! Shlomo Orr
        ha=Alfa
        D=n
        if(-h.lt.ha) then
          C1=0.
        else
          C1=-1.*ha**(3.-D)*(D-3.)*(-h)**(D-4.)
        end if
        FC=sngl(max(C1,1.d-37))
      end if

      return
      end

************************************************************************

      real function FQ(iModel,h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,n2,m2
      real h,Par(10)

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      if(iModel.eq.0.or.iModel.eq.1.or.iModel.eq.3) then
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
        else if(iModel.eq.1) then
          Qm=Par(7)
          Qa=Par(8)
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
c        m=Par(6)
        HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
        HH=max(dble(h),HMin)
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
        if(dble(h).lt.Hs) then
          Qee=(1.d0+(-Alfa*HH)**n)**(-m)
          FQ=sngl(max(Qa+(Qm-Qa)*Qee,1.d-37))
          return
        else
          FQ=sngl(Qs)
        end if
      else if(iModel.eq.2) then
        Hs=-1.d0/Alfa
        if(h.lt.Hs) then
          Qee=(-Alfa*h)**(-n)
          FQ=sngl(max(Qr+(Qs-Qr)*Qee,1.d-37))
        else
          FQ=sngl(Qs)
        end if
      else if(iModel.eq.4) then
        Hs=0.d0
        if(h.lt.Hs) then
          Qee=qnorm(dlog(-h/Alfa)/n)
          FQ=sngl(max(Qr+(Qs-Qr)*Qee,1.d-37))
        else
          FQ=sngl(Qs)
        end if
      else if(iModel.eq.5) then
        w2=Par(7)
        Alfa2=Par(8)
        n2=Par(9)
        m =1.d0-1.d0/n
        m2=1.d0-1.d0/n2
        w1=1.d0-w2
        Sw1=w1*(1.d0+(-Alfa *h)**n )**(-m )
        Sw2=w2*(1.d0+(-Alfa2*h)**n2)**(-m2)
        Qe=Sw1+Sw2
        FQ=sngl(max(Qr+(Qs-Qr)*Qe,1.d-37))
      else if(iModel.eq.-1) then                  ! Shlomo Orr
        ha=Alfa
        D=n
        FQ=Qs
        if(ha.gt.0.) FQ=max(min(Qs,Qs+(-h/ha)**(D-3.)-1.),Qr)
      end if

      return
      end

************************************************************************

      real function FH(iModel,Qe,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,n2,m2
      real Qe,Par(10)

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      if(iModel.eq.0.or.iModel.eq.1.or.iModel.eq.3) then
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
        else if(iModel.eq.1) then
          Qm=Par(7)
          Qa=Par(8)
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
c        m=Par(6)
        HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
        QeeM=(1.d0+(-Alfa*HMin)**n)**(-m)
        Qee=dmin1(dmax1(Qe*(Qs-Qa)/(Qm-Qa),QeeM),.999999999999999d0)
        FH=sngl(max(-1.d0/Alfa*(Qee**(-1.d0/m)-1.d0)**(1.d0/n),-1.d37))
      else if(iModel.eq.2) then
        FH=sngl(max(-1.d0/Alfa*max(Qe,1.e-10)**(-1.d0/n),-1.d37))
      else if(iModel.eq.4) then
        if(Qe.gt.0.9999) then
          FH=0.0
        else if(Qe.lt.0.00001) then
          FH=-1.e+8
        else
          y=Qe*2.d0
          if(y.lt.1.) p=sqrt(-dlog(y/2.d0))
          if(y.ge.1.) p=sqrt(-dlog(1-y/2.d0))
          x=p-(1.881796+0.9425908*p+0.0546028*p**3)/
     !       (1.+2.356868*p+0.3087091*p**2+0.0937563*p**3+0.021914*p**4)
          if(y.ge.1.) x=-x
          FH=sngl(-Alfa*exp(sqrt(2.)*n*x))
        end if
      else if(iModel.eq.5) then
        w2=Par(7)
        Alfa2=Par(8)
        n2=Par(9)
        m =1.d0-1.d0/n
        m2=1.d0-1.d0/n2
        w1=1.d0-w2
        Qee=Qe
        if(Qee.gt.0.9999d0) then
          FH=0.0
        else if(Qee.lt.0.00001d0) then
          FH=-1.e+8
        else
          h=xMualem(Qee,Par,10)
          FH=sngl(max(h,-1.d37))
        end if
      else if(iModel.eq.-1) then                  ! Shlomo Orr
        ha=Alfa
        D=n
        th=Qr+(Qs-Qr)*Qe
        h=0.
        if(D.ne.3.) h=-ha*(1.-Qs+th)**(1./(D-3.))
        FH=sngl(max(h,-1.d37))
      end if

      return
      end

************************************************************************

      real function FS(iModel,h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,n2,m2
      real h,Par(10)

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      if(iModel.eq.0.or.iModel.eq.1.or.iModel.eq.3) then
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
        else if(iModel.eq.1) then
          Qm=Par(7)
          Qa=Par(8)
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
c        m=Par(6)
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
        if(h.lt.Hs) then
          HMin=-1.d300**(1./n)/max(Alfa,1.d0)
          HH=max(dble(h),HMin)
          Qee=(1.d0+(-Alfa*HH)**n)**(-m)
          Qe=Qee*(Qm-Qa)/(Qs-Qa)
          FS=sngl(max(Qe,1.d-37))
        else
          FS=1.0
        end if
      else if(iModel.eq.2) then
        Hs=-1.d0/Alfa
        if(h.lt.Hs) then
          Qe=(-Alfa*h)**(-n)
          FS=sngl(max(Qe,1.d-37))
        else
          FS=1.0
        end if
      else if(iModel.eq.4) then
        Hs=0.d0
        if(h.lt.Hs) then
          Qee=qnorm(dlog(-h/Alfa)/n)
          FS=sngl(max(Qee,1.d-37))
        else
          FS=1.0
        end if
      else if(iModel.eq.5) then
        w2=Par(7)
        Alfa2=Par(8)
        n2=Par(9)
        m =1.d0-1.d0/n
        m2=1.d0-1.d0/n2
        w1=1.d0-w2
        Sw1=w1*(1.d0+(-Alfa *h)**n )**(-m )
        Sw2=w2*(1.d0+(-Alfa2*h)**n2)**(-m2)
        Qe=Sw1+Sw2
        FS=sngl(max(Qe,1.d-37))
      else if(iModel.eq.-1) then                  ! Shlomo Orr
        ha=Alfa
        D=n
        FS=1.
        if(ha.gt.0.) 
     !    FS=max(min(1.,(Qs+(-h/ha)**(D-3.)-1.-Qr)/(Qs-Qr)),0.)
      end if

      return
      end

************************************************************************

      real function FKQ(iModel,th,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,Ks,Kr,Kk
      real th,Par(10)
      integer PPar

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      Ks=Par(5)
      BPar=Par(6)
      if(iModel.eq.0.or.iModel.eq.1.or.iModel.eq.3) then     ! VG and modified VG
        PPar=2
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
          Qk=Qs
          Kk=Ks
        else if(iModel.eq.1) then
          Qm=Par(7)
          Qa=Par(8)
          Qk=Par(9)
          Kk=Par(10)
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Qeek=dmin1((Qk-Qa)/(Qm-Qa),Qees)
        if(dble(th).lt.Qk) then
          Qee=(dble(th)-Qa)/(Qm-Qa)
          Qe =(Qm-Qa)/(Qs-Qa)*Qee
          Qek=(Qm-Qa)/(Qs-Qa)*Qeek
          FFQ =1.d0-(1.d0-Qee **(1.d0/m))**m
          FFQk=1.d0-(1.d0-Qeek**(1.d0/m))**m
          if(FFQ.le.0.d0) FFQ=m*Qee**(1.d0/m)
          Kr=(Qe/Qek)**Bpar*(FFQ/FFQk)**PPar*Kk/Ks
          FKQ=sngl(max(Ks*Kr,1.d-37))
        end if
        if(dble(th).ge.Qs) FKQ=sngl(Ks)
      else if(iModel.eq.-1) then                  ! Shlomo Orr
        D=n
        Kr=1.
        Qx=0.
        if(D.ne.3.) Qx=Qr+2.*(1.-Qs)/(D/(3.-D)-2.)
        if(th.gt.Qx) then
          S=th/Qs
          if(D.ne.3.) Kr=(1.-Qs*(1-S)/(1.-Qr))**(D/(3.-D))
        else
          Sx=Qx/Qs
          if(D.ne.3.) then
            Kr=(1.-Qs*(1-Sx)/(1.-Qr))**(D/(3.-D))
            if(Qx.gt.Qr) Kr=Kr*(th-Qr)**2/(Qx-Qr)**2
          end if
        end if
        FKQ=sngl(max(Ks*Kr,1.d-37))
      end if

      return
      end

************************************************************************

      real function FKS(iModel,S,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,Ks,Kr,Kk
      real S,th,Par(10)
      integer PPar

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      Ks=Par(5)
      BPar=Par(6)
      if(iModel.eq.0.or.iModel.eq.1.or.iModel.eq.3) then   ! VG and modified VG
        PPar=2
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
          Qk=Qs
          Kk=Ks
        else if(iModel.eq.1) then
          Qm=Par(7)
          Qa=Par(8)
          Qk=Par(9)
          Kk=Par(10)
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Qeek=dmin1((Qk-Qa)/(Qm-Qa),Qees)
        th=S*(Qs-Qr)+Qr
        if(dble(th).lt.Qk) then
          Qee=(dble(th)-Qa)/(Qm-Qa)
          Qe =(Qm-Qa)/(Qs-Qa)*Qee
          Qek=(Qm-Qa)/(Qs-Qa)*Qeek
          FFQ =1.d0-(1.d0-Qee **(1.d0/m))**m
          FFQk=1.d0-(1.d0-Qeek**(1.d0/m))**m
          if(FFQ.le.0.d0) FFQ=m*Qee**(1.d0/m)
          Kr=(Qe/Qek)**Bpar*(FFQ/FFQk)**PPar*Kk/Ks
          FKS=sngl(max(Ks*Kr,1.d-37))
        end if
        if(dble(th).ge.Qs) FKS=sngl(Ks)
      else if(iModel.eq.-1) then                  ! Shlomo Orr
        D=n
        Kr=1.
        Qx=0.
        th=S*(Qs-Qr)+Qr
        if(D.ne.3.) Qx=Qr+2.*(1.-Qs)/(D/(3.-D)-2.)
        if(th.gt.Qx) then
          Se=th/Qs
          if(D.ne.3.) Kr=(1.-Qs*(1-Se)/(1.-Qr))**(D/(3.-D))
        else
          Sx=Qx/Qs
          if(D.ne.3.) then
            Kr=(1.-Qs*(1-Sx)/(1.-Qr))**(D/(3.-D))
            if(Qx.gt.Qr) Kr=Kr*(th-Qr)**2/(Qx-Qr)**2
          end if
        end if
        FKQ=sngl(max(Ks*Kr,1.d-37))
      end if

      return
      end

************************************************************************

      double precision function qnorm(x)

      implicit double precision(A-H,O-Z)

      z=abs(x/2.**0.5)
      t=1./(1.+0.5*z)
      erfc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
     !     t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
     !     t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))
      if(x.lt.0.) erfc=2.-erfc
      qnorm=erfc/2.

      return
      end

************************************************************************

*     Evaluate h for given theta_e for dual-porosity function

      real*8 function xMualem(SE,Par,NPar)

      implicit real*8(A-H,O-Z)
      real Par
      dimension Par(NPar)

      x1=-1.e-6
      x2=-1.e+6
      call ZBRAK(X1,X2,XB1,XB2,SE,Par,NPar)
      hhh=ZBRENT(XB1,XB2,SE,Par,NPar)
      xMualem=hhh    ! for calculation hh  
      if(hhh.ne.0.) then
c        xMualem=1./hhh ! for integration
      else
        PAUSE 'xMualem: h is equal to zero!'
      end if
      return
      end

************************************************************************

      function DoublePor(hh,SE,Par,NPar)

*     Double porosity function - for evaluation of h for given theta_e

      implicit real*8(A-H,O-Z)
      real Par
      dimension Par(NPar)

      wcr=Par(1)
      wcs=Par(2)
      Alpha=Par(3)
      rn=Par(4)
      rm=1.-1./rn
      w2=Par(7)
      w1=1.-W2
      Alpha2=Par(8)
      rn2=Par(9)
      rm2=1.-1./rn2

      Sw1=w1*(1.+(-Alpha *hh)**rn )**(-rm )
      Sw2=w2*(1.+(-Alpha2*hh)**rn2)**(-rm2)
      rwc=Sw1+Sw2
      DoublePor=SE-rwc

      return
      end

************************************************************************

*     Bracketing of the root, Numerical recepies (345)

      subroutine ZBRAK(X1,X2,XB1,XB2,SE,Par,NPar)

      implicit real*8(A-H,O-Z)
      real Par
      dimension Par(NPar)

      NB=1
      NBB=NB
      NB=0
      n=1000

      dlh=(dlog10(-X2)-dlog10(-X1))/(N-1)
      FP=DoublePor(X1,SE,Par,NPar)
      do 11 i=1,n
        dx2=dlog10(-X1)+(i  )*dlh
        X2=-10**dx2
        FC=DoublePor(X2,SE,Par,NPar)
        if(FC*FP.lt.0.) then
          XB1=X1
          XB2=X2
          return
        end if
        FP=FC
        X1=X2
        if(NBB.eq.NB) return
11    continue

      return
      end

************************************************************************

*     Brent method of finding root that lies between x1 and x2, 
*     Numerical recepies (354)

      real*8 function ZBRENT(X1,X2,SE,Par,NPar)

      implicit real*8(A-H,O-Z)
      parameter (ITMAX=100,EPS=3.E-8,TOL=1.e-6)
      real Par
      dimension Par(NPar)

      A=X1
      B=X2
      FA=DoublePor(A,SE,Par,NPar)
      FB=DoublePor(B,SE,Par,NPar)
      IF(FB*FA.GT.0.) PAUSE 'Root must be bracketed for ZBRENT.'
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRENT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=DoublePor(B,SE,Par,NPar)
11    CONTINUE
      PAUSE 'ZBRENT exceeding maximum iterations.'
      ZBRENT=B

      RETURN
      END

************************************************************************

      function Gamma(Z)

*     Purpose:  To calculate the Gamma function for positive Z

      implicit real*8 (A-H,O-Z)
      if(Z.lt.33.) goto 11
      Gamma=1.d36
      return
11    X=Z
      Gamma=1.0
      if(X-2.0) 14,14,13
12    if(X-2.0) 16,16,13
13    X=X-1.0
      Gamma=Gamma*X
      goto 12
14    if(X-1.0) 15,17,16
15    Gamma=Gamma/X
      X=X+1.0
16    Y=X-1.0
      FY=1.0-Y*(.5771017-Y*(.985854-Y*(.8764218-Y*(.8328212-Y*(.5684729-
     !Y*(.2548205-.0514993*Y))))))
      Gamma=Gamma*FY
17    return
      end

************************************************************************

      function BInc(X,A,B,Beta)

*     Purpose: To calculate the incomplete Beta-function

      implicit real*8 (A-H,O-Z)
      dimension T(200)
      data NT/10/
      NT1=NT+1
      T(1)=-(A+B)*X/(A+1.0)
      do 11 i=2,NT,2
        Y=float(i/2)
        Y2=float(i)
        T(i)=Y*(B-Y)*X/((A+Y2-1.0)*(A+Y2))
        T(i+1)=-(A+Y)*(A+B+Y)*X/((A+Y2)*(A+Y2+1.0))
11     continue
      BInc=1.0
      do 12 i=1,NT
        k=NT1-i
        BInc=1.+T(k)/BInc
12     continue
      BInc=X**A*(1.-X)**B/(BInc*A*Beta)
      return
      end

************************************************************************

      subroutine qromb(a,b,ss,iModel,Par)
      integer JMAX,JMAXP,K,KM,iModel
      real a,b,ss,EPS,Par(11)
      parameter (eps=1.e-6, jMax=20, jMaxP=jMax+1, K=5, KM=K-1)
      integer j
      real dss,h(JMAXP),s(JMAXP)

      h(1)=1.
      do 11 j=1,jMax
        call trapzd(a,b,s(j),j,iModel,Par)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.eps*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      pause 'too many steps in qromb'
      end

***********************************************************************

      subroutine trapzd(a,b,s,n,iModel,Par)
      integer n
      real a,b,s
      integer it,j,iModel
      real del,sum,tnm,x,Par(11)

      if (n.eq.1) then
        s=0.5*(b-a)*(FH(iModel,a,Par)+FH(iModel,b,Par))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+FH(iModel,x,Par)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      end

************************************************************************

      subroutine polint(xa,ya,n,x,y,dy)
      integer n,NMAX
      real dy,x,y,xa(n),ya(n)
      parameter (NMAX=10)
      integer i,m,ns
      real den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
