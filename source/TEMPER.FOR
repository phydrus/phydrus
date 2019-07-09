* Source file TEMPER.FOR |||||||||||||||||||||||||||||||||||||||||||||||

*     Calculation of heat transport

      subroutine Temper(N,NMat,x,dt,t,MatNum,TempO,TempN,TPar,Ampl,B,D,
     !                  E,F,vOld,vNew,ThOld,ThNew,Cap,Cond,Sink,tPeriod,
     !                  kTopT,tTop,kBotT,tBot,lVapor,ThVOld,ThVNew,
     !                  vVOld,vVNew,g0,lEnBal,HeatFl,xConv,tConv,dtMaxT,
     !                  iCampbell,iTemp)

      logical lVapor,lEnBal
      double precision B,D,E,F,t
      dimension x(N),MatNum(N),TempO(N),TempN(N),TPar(10,NMat),B(N),
     !          D(N),E(N),F(N),vOld(N),vNew(N),ThOld(N),ThNew(N),Cap(N),
     !          Cond(N),Sink(N),ThVOld(N),ThVNew(N),vVOld(N),vVNew(N),
     !          g0(N)

      Cw=TPar(9,1)
      PI=3.141592654
      tTopA=tTop
      if(tPeriod.gt.0.) tTopA=tTop+Ampl*sin(2.*PI*sngl(t)/tPeriod-7.*PI/
     !                                      12.)

*     Upper Flux BC
      if(kTopT.lt.0.and..not.lEnBal) then
        HeatFl=Cw*tTopA*(vNew(N)+vOld(N))/2.
        if(lVapor) then
          Cv=1.8e+6/xConv/tConv/tConv
          xLat=xLatent(TempN(N))/xConv/tConv/tConv
          HeatFl=HeatFl+Cv*tTopA*(vVNew(N)+vVOld(N))/2.+
     !                  xLat*(vVNew(N)+vVOld(N))/2.
        end if
      end if

      do 11 Level=1,2
        if(Level.eq.1) then
          call TempCap(N,NMat,MatNum,ThOld,vOld,Cap,Cond,TPar,iCampbell,
     !                 xConv,tConv)
        else
          call TempCap(N,NMat,MatNum,ThNew,vNew,Cap,Cond,TPar,iCampbell,
     !                 xConv,tConv)
        end if
        call SetUpHeat(N,Level,lEnBal,kBotT,kTopT,tBot,tTopA,HeatFl,dt,
     !                 x,Cw,B,D,E,F,TempO,Cap,Cond,vNew,vOld,Sink)
11    continue

*     Adjust matrix for vapor flow effects
      if(lVapor) call TempAdj(N,x,dt,TempN,TempO,B,D,E,F,vVOld,vVNew,g0,
     !                        kTopT,kBotT,tTopA,tBot,ThVNew,ThVOld,
     !                        xConv,tConv)

*     Solve matrix equation
      call BanSol(N,B,D,E,F)
      do 12 i=1,N
        TempN(i)=sngl(F(i))
12    continue

*     Max time step
      dtMaxT=1.e+30
      Cv=1.8e+6/xConv/tConv/tConv
      if(lEnBal) then
        if(iTemp.eq.0) then
          dtMaxT=0.
          dtMaxT=abs(Cap(N)*(x(N)-x(N-1))/
     !            (max(0.1,TempN(N))*(Cw*vNew(N)+Cv*vVNew(N))))
        end if
        dTempMax=2.
        dTemp=abs(TempN(N)-TempO(N))
        dtMaxT=min(dtMaxT,dTempMax/max(0.1,dTemp)*dt)
      end if

      return
      end

************************************************************************

      subroutine SetUpHeat(N,Level,lEnBal,kBotT,kTopT,tBot,tTopA,HeatFl,
     !                     dt,x,Cw,B,D,E,F,TempO,Cap,Cond,vNew,vOld,
     !                     Sink)

      logical lEnBal
      double precision B,D,E,F
      dimension x(N),TempO(N),B(N),D(N),E(N),F(N),vOld(N),vNew(N),
     !          Cap(N),Cond(N),Sink(N)

      dx=x(2)-x(1)
      if(kBotT.gt.0) then
        D(1)=1.
        E(1)=0.
        F(1)=tBot
      else if(kBotT.lt.0) then
        if(Level.eq.2) then
          D(1)=dx/2./dt*Cap(1)+(Cond(1)+Cond(2))/dx/4.+
     !         Cw*(2.*vNew(1)+vNew(2))/12.+
     !         dx/24.*Cw*(3.*Sink(1)+Sink(2))
          E(1)=-(Cond(1)+Cond(2))/4./dx+
     !         Cw*(2.*vNew(2)+vNew(1))/12.+dx/24.*Cw*(Sink(1)+Sink(2))
        else
          F(1)=TempO(1)*(dx/2./dt*Cap(1)-
     !                   (Cond(1)+Cond(2))/dx/4.-
     !                   Cw*(2.*vOld(1)+vOld(2))/12.-
     !                   dx/24.*Cw*(3.*Sink(1)+Sink(2)))+
     !         TempO(2)*((Cond(1)+Cond(2))/4./dx-
     !                   Cw*(2.*vOld(2)+vOld(1))/12.-
     !                   dx/24.*Cw*(Sink(1)+Sink(2)))+
     !                   tBot*Cw*(vNew(1)+vOld(1))/2.
        end if
      else
        D(1)=-1.
        E(1)=1.
        F(1)=0.
      end if
      do 12 i=2,N-1
        dxA=x(i)-x(i-1)
        dxB=x(i+1)-x(i)
        dx=(x(i+1)-x(i-1))/2.
        if(Level.eq.2) then
          B(i)=-(Cond(i)+Cond(i-1))/4./dxA-Cw*(vNew(i)+
     !         2.*vNew(i-1))/12.+
     !         dxA/24.*Cw*(Sink(i-1)+Sink(i))
          D(i)=(Cond(i-1)+Cond(i))/4./dxA+(Cond(i)+Cond(i+1))/4./dxB+
     !         dx/dt*Cap(i)+
     !         Cw*(vNew(i+1)-vNew(i-1))/12.+
     !         dxA/24.*Cw*(Sink(i-1)+3.*Sink(i))+
     !         dxB/24.*Cw*(3.*Sink(i)+Sink(i+1))
          E(i)=-(Cond(i)+Cond(i+1))/4./dxB+
     !         Cw*(2.*vNew(i+1)+vNew(i))/12.+
     !         dxB/24*Cw*(Sink(i+1)+Sink(i))
        else
          F(i)=TempO(i-1)*((Cond(i)+Cond(i-1))/4./dxA+
     !                     Cw*(vOld(i)+2.*vOld(i-1))/12.-
     !                     dxA/24.*Cw*(Sink(i-1)+Sink(i)))+
     !         TempO(i)*(-Cw*(vOld(i+1)-vOld(i-1))/12.+
     !                   dx/dt*Cap(i)-
     !                   (Cond(i+1)+Cond(i))/4./dxB-
     !                   (Cond(i)+Cond(i-1))/4./dxA-
     !                   dxA/24.*Cw*(Sink(i-1)+3.*Sink(i))-
     !                   dxB/24.*Cw*(3.*Sink(i)+Sink(i+1)))+
     !         TempO(i+1)*((Cond(i+1)+Cond(i))/4./dxB-
     !                     Cw*(2.*vOld(i+1)+vOld(i))/12.-
     !                     dxB/24.*Cw*(Sink(i+1)+Sink(i)))
        end if
12    continue
      if(kTopT.gt.0) then
        B(N)=0.
        D(N)=1.
        F(N)=tTopA
      else if(kTopT.lt.0) then
        dx=x(N)-x(N-1)
        if(Level.eq.2) then
          B(N)=-(Cond(N)+Cond(N-1))/4./dx-
     !         Cw*(vNew(N)+2.*vNew(N-1))/12.+
     !         dx/24.*Cw*(Sink(N-1)+Sink(N))
          D(N)=dx/2./dt*Cap(N)+
     !         (Cond(N-1)+Cond(N))/4./dx-
     !         Cw*(2.*vNew(N)+vNew(N-1))/12.+
     !         dx/24.*Cw*(Sink(N-1)+3.*Sink(N))
        else
          F(N)=TempO(N-1)*((Cond(N)+Cond(N-1))/4./dx+
     !                     Cw*(vOld(N)+2.*vOld(N-1))/12.-
     !                     dx/24.*Cw*(Sink(N-1)+Sink(N)))+
     !         TempO(N)*(dx/2./dt*Cap(N)-
     !                   (Cond(N-1)+Cond(N))/4./dx+
     !                   Cw*(2.*vOld(N)+vOld(N-1))/12.-
     !                   dx/24.*Cw*(Sink(N-1)+3.*Sink(N)))
          if(.not.lEnBal) then
            F(N)=F(N)-HeatFl
            F(N)=F(N)-min(0.,HeatFl)
          else
            F(N)=F(N)-HeatFl
          end if
        end if
      end if

      return
      end

************************************************************************

      subroutine TempAdj(N,x,dt,TempN,TempO,B,D,E,F,vVOld,vVNew,g0,
     !                   kTopT,kBotT,tTopA,tBot,ThVNew,ThVOld,xConv,
     !                   tConv)

      double precision B,D,E,F
      real Lat
      dimension x(N),TempN(N),TempO(N),B(N),D(N),E(N),F(N),vVOld(N),
     !          vVNew(N),ThVNew(N),ThVOld(N),g0(N)

*     Cv  - volumetric specific heat of vapor [J/m3/K,kg/m/s2/K]
*     Lat - volumetric latent heat of vaporization of water [J/m3,kg/m/s2]

      Cv=1.8e+6/xConv/tConv/tConv

      do 13 iLevel=1,2
        do 11 i=1,N
          if(iLevel.eq.1) then
            if(i.eq.1) then 
              vVGrad=(vVOld(i+1)-vVOld(i))/(x(i+1)-x(i))
            else if(i.eq.N) then
              vVGrad=(vVOld(i)-vVOld(i-1))/(x(i)-x(i-1))
            else
              vVGrad=(vVOld(i+1)-vVOld(i-1))/(x(i+1)-x(i-1))*2.
            end if
            Lat=xLatent(TempO(i))/xConv/tConv/tConv
          else
            if(i.eq.1) then 
              vVGrad=(vVNew(i+1)-vVNew(i))/(x(i+1)-x(i))
            else if(i.eq.N) then
              vVGrad=(vVNew(i)-vVNew(i-1))/(x(i)-x(i-1))
            else
              vVGrad=(vVNew(i+1)-vVNew(i-1))/(x(i+1)-x(i-1))*2.
            end if
            Lat=xLatent(TempN(i))/xConv/tConv/tConv
          end if
          ThVGrad=(ThVNew(i)-ThVOld(i))/dt
          g0(i)=-Lat*(vVGrad+ThVGrad)
11      continue

        dx=x(2)-x(1)
        if(kBotT.lt.0) then
          if(iLevel.eq.1) then
            F(1)=F(1)+TempO(1)*(-Cv*(2.*vVOld(1)+vVOld(2))/12.)+
     !                TempO(2)*(-Cv*(2.*vVOld(2)+vVOld(1))/12.)+
     !                dx/12.*(2.*g0(1)+g0(2))+
     !                tBot*Cv*(vVNew(1)+vVOld(1))/2.
          else
            D(1)=D(1)+Cv*(2.*vVNew(1)+vVNew(2))/12.
            E(1)=E(1)+Cv*(2.*vVNew(2)+vVNew(1))/12.
            F(1)=F(1)+dx/12.*(2.*g0(1)+g0(2))
          end if
        end if
        do 12 i=2,N-1
          dxA=x(i)-x(i-1)
          dxB=x(i+1)-x(i)
          dx=(x(i+1)-x(i-1))/2.
          if(iLevel.eq.1) then
            F(i)=F(i)+TempO(i-1)*( Cv*(vVOld(i) +2.*vVOld(i-1))/12.)+
     !                TempO(i)  *(-Cv*(vVOld(i+1)  -vVOld(i-1))/12.)+
     !                TempO(i+1)*(-Cv*(2.*vVOld(i+1) +vVOld(i))/12.)+
     !             dxA*(g0(i-1)+2.*g0(i))/12.+dxB*(2.*g0(i)+g0(i+1))/12.
          else
            B(i)=B(i)-Cv*(vVNew(i)+2.*vVNew(i-1))/12.
            D(i)=D(i)+Cv*(vVNew(i+1)-vVNew(i-1))/12.
            E(i)=E(i)+Cv*(2.*vVNew(i+1)+vVNew(i))/12.
            F(i)=F(i)+dxA*(g0(i-1)+2.*g0(i))/12.+
     !                dxB*(2.*g0(i)+g0(i+1))/12.
          end if
12      continue
        if(kTopT.lt.0) then
          dx=x(N)-x(N-1)
          if(iLevel.eq.1) then
            F(N)=F(N)+TempO(N-1)*(Cv*(vVOld(N)+2.*vVOld(N-1))/12.)+
     !                TempO(N)  *(Cv*(2.*vVOld(N)+vVOld(N-1))/12.)+
     !                dx/12.*(g0(N-1)+2.*g0(N))
            if(kTopT.eq.-1) F(N)=F(N)-tTopA*Cv*(vVNew(N)+vVOld(N))/2.
          else
            B(N)=B(N)-Cv*(vVNew(N)+2.*vVNew(N-1))/12.
            D(N)=D(N)-Cv*(2.*vVNew(N)+vVNew(N-1))/12.
            F(N)=F(N)+dx/12.*(g0(N-1)+2.*g0(N))
          end if
        end if
13    continue

      return
      end

************************************************************************

      subroutine TempCap(N,NMat,MatNum,Theta,Veloc,Cap,Cond,TPar,
     !                   iCampbell,xConv,tConv)

      dimension MatNum(N),Theta(N),Veloc(N),Cap(N),Cond(N),
     !          TPar(10,NMat)

      do 11 i=1,N
        M=MatNum(i)
        th=Theta(i)
        v=Veloc(i)
        Cap(i)=TPar(7,M)*TPar(1,M)+TPar(8,M)*TPar(2,M)+TPar(9,M)*th

        if(Cap(i).eq.0.) then
          write(*,*) 'Heat capacity is equal to zero'
          read(*,*)
          stop
        end if
        Cond(i)=amax1(0.,TPar(4,M)+TPar(5,M)*th+TPar(6,M)*sqrt(th))
        if(iCampbell.eq.1) then
*         TPar(1,M) - the volume fraction of solids
*         TPar(4,M) - the volume fraction of quartz
*         TPar(5,M) - the volume fraction of other minerals
*         TPar(6,M) - the volume fraction of clay
          AA=(0.57+1.73*TPar(4,M)+0.93*TPar(5,M))/(1.-0.74*TPar(4,M)-
     !       0.49*TPar(5,M))-2.8*TPar(1,M)*(1.-TPar(1,M))
          BB=2.8*TPar(1,M)
          xc=amax1(0.005,TPar(6,M))
          CC=1.+2.6/sqrt(xc)
          DD=0.03+0.7*TPar(1,M)**2
          EE=4.
          XLamb=AA+BB*th-(AA-DD)*exp(-(CC*th)**EE)
          Cond(i)=amax1(0.,XLamb*xConv/tConv/tConv/tConv)
        end if

        Cond(i)=Cond(i)+TPar(9,M)*TPar(3,M)*abs(v)
11    continue
      return
      end

************************************************************************

      real function xLatent(Temp)

*     Function calculating the volumetric latent heat of vaporization of
*     water [J/m3,kg/m/s2]
*     row   - density of soil water [kg/m3]
*     xLw   - latent heat of vaporization of water [J/kg,ML2/T2]

      row=(1.-7.37e-6*(Temp-4.)**2+3.79e-8*(Temp-4.)**3)*1000.
      xLw=2.501e+06-2369.2*Temp
      xLatent=row*xLw

      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||