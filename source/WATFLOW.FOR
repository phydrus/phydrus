* Source file WATFLOW.FOR ||||||||||||||||||||||||||||||||||||||||||||||

      subroutine WatFlow(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,
     !                   hOld,MatNum,ParD,ParW,Con,Cap,ConSat,Ah,AK,ATh,
     !                   hSat,hTemp,KodTop,KodBot,rTop,rBot,CosAlf,t,dt,
     !                   x,Sink,P,R,S,FreeD,SeepF,qGWLF,Aqh,Bqh,GWL0L,
     !                   hTop,hBot,hCritA,hCritS,WLayer,Iter,ItCum,
     !                   TopInf,KTOld,KBOld,TolTh,TolH,MaxIt,dtMin,tOld,
     !                   dtOpt,ConvgF,TheTab,ThNew,ThOld,thr,ths,lWTDep,
     !                   TempN,Kappa,KappaO,AThS,ThRR,ConO,ConR,AKS,AhW,
     !                   AThW,AKW,iHyst,iModel,qDrain,zBotDr,BaseGW,
     !                   rSpacing,iPosDr,rKhTop,rKhBot,rKvTop,rKvBot,
     !                   Entres,WetPer,zInTF,GeoFac,lTable,lVapor,xConv,
     !                   tConv,ConLT,ConVT,ConVh,TauW,ThEq,ThVNew,
     !                   ThVOld,nTabMod,iDualPor,ThNewIm,ThOldIm,SinkIm,
     !                   vTop,TempO,iTemp,WTransf,lDensity,Conc,NSD,
     !                   iEnhanc,lCentrif,Radius,hSeep)

      logical ConvgF,ItCrit,FreeD,qGWLF,TopInf,WLayer,SeepF,lWTDep,
     !        qDrain,lTable,lVapor,lDensity,lCentrif
      double precision P,R,S,PB,RB,SB,PT,RT,ST,rMin,t,tOld
      dimension x(NumNP),hNew(NumNP),hOld(NumNP),hTemp(NumNP),thr(NMat),
     !          MatNum(NumNP),ParD(11,NMat),Sink(NumNP),Con(NumNP),
     !          ConTab(NTabD,NMat),CapTab(NTabD,NMat),hTab(NTabD,NMat),
     !          Cap(NumNP),ConSat(NMat),P(NumNP),R(NumNP),S(NumNP),
     !          hSat(NMat),Ah(NumNP),AK(NumNP),ATh(NumNP),ths(NMat),
     !          TheTab(NTabD,NMat),ThNew(NumNP),ThOld(NumNP),
     !          TempN(NumNP),ParW(11,NMat),ConO(NumNP),Kappa(NumNP),
     !          AThS(NumNP),ThRR(NumNP),ConR(NumNP),AKS(NumNP),
     !          AhW(NMat),AThW(NMat),KappaO(NumNP),AKW(NMat),NTab(NMat),
     !          ConLT(NumNP),ConVT(NumNP),ConVh(NumNP),ThEq(NumNP),
     !          ThVNew(NumNP),ThVOld(NumNP),SinkIm(NumNP),TempO(NumNP),
     !          ThNewIm(NumNP),ThOldIm(NumNP),Conc(NSD,NumNP)

      rMax=1.e+10
      rMin=1.d-100

*     Nonequilibrium transport [Ross and Smettem, 2001]
      Rate=1.
      if(TauW.gt.0.) Rate=amin1(1.,amax1(0.000001,1.-exp(-dt/TauW)))

*     Dual porosity mass transfer
      if(iDualPor.gt.0)
     !  call DualPor(NumNP,NMat,MatNum,iDualPor,ThOld,ths,thr,ThNewIm,
     !               ThOldIm,ParD,SinkIm,dt,iModel,hNew,hCritA,x,
     !               WTransf)

11    continue

      Iter=0
      ConvgF=.true.
*     End of ponding (works for both BC and VG)
      if(WLayer.and.hNew(NumNP).gt.0..and.hNew(NumNP).lt.0.00005*xConv.
     !   and.rTop.ge.0.) then
        hNew(NumNP)=FH(iModel,0.9999,ParD(1,MatNum(NumNP)))
        hOld(NumNP)=hNew(NumNP)
        hTemp(NumNP)=hNew(NumNP)
      end if

12    continue

*     Generate terms of matrix equation and solve by Gauss elimination
      if(iHyst.ne.3) then
      call SetMat(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,MatNum,
     !            ParD,Con,Cap,ConSat,Ah,AK,ATh,hSat,hTemp,TheTab,ThEq,
     !            thr,ths,lWTDep,TempN,Iter,ConO,Kappa,AThS,ThRR,ConR,
     !            AKS,AhW,AThW,AKW,iModel,lTable,lVapor,ThVNew,ConLT,
     !            ConVT,ConVh,xConv,tConv,nTabMod,hCritA,lDensity,Conc,
     !            NSD,iEnhanc)
      if(Iter.eq.2.and.iHyst.gt.0)
     !call Hyster(NumNP,NMat,hOld,MatNum,ParD,ParW,ThNew,ThOld,Kappa,
     !            AThS,ThRR,ConO,ConR,AKS,KappaO,Ah,AK,iHyst,iModel,
     !            TolTh)
      else ! Bob Lenhard hysteresis
        call Hyst(NumNP,NMat,ParD,ParW,MatNum,Kappa,hNew,hOld,ThEq,
     !            Con,Cap,0,2)
      end if
      call Reset (NumNP,rTop,rBot,CosAlf,dt,x,hOld,Con,Cap,WLayer,hNew,
     !            Sink,P,R,S,PB,RB,SB,PT,RT,ST,FreeD,qGWLF,Aqh,Bqh,
     !            GWL0L,ThNew,ThOld,vTop,qDrain,zBotDr,BaseGW,rSpacing,
     !            iPosDr,rKhTop,rKhBot,rKvTop,rKvBot,Entres,WetPer,
     !            zInTF,GeoFac,ThEq,Rate,lVapor,lWTDep,TempN,ConLT,
     !            ConVT,ConVh,ThVNew,ThVOld,iDualPor,SinkIm,lDensity,
     !            Conc,NSD,lCentrif,Radius)

      call Shift (NumNP,KodTop,rTop,rBot,hTop,hBot,hCritA,CosAlf,WLayer,
     !            Con,hNew,x,TopInf,KodBot,SeepF,ThNew,ThOld,Sink,dt,
     !            lVapor,lWTDep,ConLT,TempN,ConVh,ConVT,iDualPor,SinkIm,
     !            lDensity,Conc,NSD,lCentrif,Radius,hSeep)
      do 13 i=1,NumNP
        hTemp(i)=hNew(i)
13    continue
      call Gauss (NumNP,KodTop,KodBot,hTop,hBot,hNew,P,R,S,PB,RB,SB,PT,
     !            RT,ST,rMin)
      do 17 i=1,NumNP
        if(abs(hNew(i)).gt.rMax) hNew(i)=sign(rMax,hNew(i))
        if(abs(KodTop).eq.4.and.hNew(i).lt.hCritA.and.i.eq.NumNP)
     !                        hNew(i)=hCritA
        if(abs(KodTop).eq.4.and.hNew(i).lt.hCritA.and.i.gt.NumNP*9/10.
     !     and.Sink(i).le.0.) hNew(i)=hCritA
17    continue
      Iter =Iter+1
      ItCum=ItCum+1

*     Test for convergence
      ItCrit=.true.
      do 14 i=1,NumNP
        m=MatNum(i)
        EpsTh=0.
        EpsH=0.
        if(hTemp(i).lt.hSat(m).and.hNew(i).lt.hSat(m)) !.and.TauW.eq.0.
     !                                                              then
          Th=ThNew(i)+Cap(i)*(hNew(i)-hTemp(i))/(ths(m)-thr(m))/ATh(i)*
     !       Rate
          EpsTh=abs(ThNew(i)-Th)
c          if(TauW.gt.0) EpsH=abs(hNew(i)-hTemp(i))-abs(0.05*hNew(i))
        else
          EpsH=abs(hNew(i)-hTemp(i))
        end if
        if(EpsTh.gt.TolTh.or.EpsH.gt.TolH.or.abs(hNew(i)).gt.rMax*0.999)
     !                                                              then
          ItCrit=.false.
          if(abs(hNew(i)).gt.rMax*0.999) Iter=MaxIt
          goto 15
        end if
14    continue
15    continue
      if(.not.ItCrit.or.(Iter.le.1.or.(Iter.le.2.and.iHyst.gt.0))) then
        if(Iter.lt.MaxIt) then
          goto 12
        else if(dt.le.dtMin) then
          ConvgF=.false.
          write(*,*) ' The numerical solution has not converged ! '
          return
        else
          do 16 i=1,NumNP
            if(iHyst.gt.0) Kappa(i)=KappaO(i)
            hNew(i) =hOld(i)
            hTemp(i)=hOld(i)
            TempN(i)=TempO(i)
16        continue
          KodTop=KTOld
          KodBot=KBOld
          dt=amax1(dt/3,dtMin)
          dtOpt=dt
          t=tOld+dt
          if(TauW.gt.0.) Rate=amin1(1.,amax1(0.000001,1.-exp(-dt/TauW)))
          iTemp=0
          goto 11
        end if
      end if
      if(ItCrit) then
        do 18 i=1,NumNP
          ThNew(i)=ThNew(i)+Cap(i)*(hNew(i)-hTemp(i))*Rate
18      continue
      end if
      if(Wlayer) then
        if(hNew(NumNP).gt.hCritS) then
          KodTop=4
          hTop=hCritS
        end if
      end if

      if(iHyst.eq.3)
     !  call Hyst(NumNP,NMat,ParD,ParW,MatNum,Kappa,hNew,hOld,ThNew,
     !            Con,Cap,0,3)

      return
      end

************************************************************************

      subroutine Reset(N,rTop,rBot,CosAlf,dt,x,hOld,Con,Cap,WLayer,hNew,
     !                 Sink,P,R,S,PB,RB,SB,PT,RT,ST,FreeD,qGWLF,Aqh,Bqh,
     !                 GWL0L,ThNew,ThOld,vTop,qDrain,zBotDr,BaseGW,
     !                 rSpacing,iPosDr,rKhTop,rKhBot,rKvTop,rKvBot,
     !                 Entres,WetPer,zInTF,GeoFac,ThEq,Rate,lVapor,
     !                 lWTDep,Temp,ConLT,ConVT,ConVh,ThVNew,ThVOld,
     !                 iDualPor,SinkIm,lDensity,Conc,NSD,lCentrif,
     !                 Radius)

      logical WLayer,FreeD,qGWLF,qDrain,lVapor,lWTDep,lDensity,lCentrif,
     !        lGeom
      double precision P,R,S,PB,RB,SB,PT,RT,ST,A2,A3,B,F2
      dimension x(N),hOld(N),hNew(N),P(N),R(N),S(N),Con(N),Cap(N),
     !          Sink(N),ThNew(N),ThOld(N),ThEq(N),Temp(N),ConLT(N),
     !          ConVT(N),ConVh(N),ThVNew(N),ThVOld(N),SinkIm(N),
     !          Conc(NSD,N)

      lGeom=.false.   ! Arithmetic average (false), geometric average (true)
      fRE=1.
      Grav=CosAlf
      do 10 i=1,N
        ThNew(i)=ThOld(i)+(ThEq(i)-ThOld(i))*Rate
10    continue

*     Finite differences

*     Bottom BC
      dxB=x(2)-x(1)
      dx=dxB/2.
      ConB=(Con(1)+Con(2))/2.              ! Arithmetic average
      if(lGeom) ConB=(Con(1)*Con(2))**0.5  ! Geometric average
      if(lCentrif) Grav=CosAlf*(Radius+abs((x(1)+x(2))/2.))
      B=ConB*Grav
      if(lDensity) then
        fRE=fRo(1,(Conc(1,1)+Conc(1,2))/2.)
        B=ConB*Grav*fRE
        fRE=fRo(1,Conc(1,1))
      end if
      if(lVapor) ConB=ConB+(ConVh(1)+ConVh(2))/2.
      S(1)=-ConB/dxB
      if(FreeD) rBot=-ConB*Grav*fRE
      F2=Cap(1)*dx/dt*fRE*Rate
      RB=ConB/dxB+F2
      SB=-ConB/dxB
      if(qGWLF) rBot=Fqh(hNew(1)-GWL0L,Aqh,Bqh)
      if(qDrain) rBot=FqDrain(x(1)+hNew(1),zBotDr,BaseGW,rSpacing,
     !                        iPosDr,rKhTop,rKhBot,rKvTop,rKvBot,
     !                        Entres,WetPer,zInTF,GeoFac)
      PB=B-Sink(1)*dx+F2*hNew(1)-(ThNew(1)-ThOld(1))*dx/dt*fRE+
     !   rBot
      if(iDualPor.gt.0) PB=PB-SinkIm(1)*dx
      if(lVapor.or.lWTDep) then
        ConTB=0.
        if(lVapor) ConTB=ConTB+(ConVT(1)+ConVT(2))/2.
        if(lWTDep) ConTB=ConTB+(ConLT(1)+ConLT(2))/2.
        PB=PB+ConTB*(Temp(2)-Temp(1))/dxB-(ThVNew(1)-ThVOld(1))*dx/dt
      end if
      do 11 i=2,N-1
        dxA=x(i)-x(i-1)
        dxB=x(i+1)-x(i)
        dx=(dxA+dxB)/2.
        ConA=(Con(i)+Con(i-1))/2.
        ConB=(Con(i)+Con(i+1))/2.
        if(lGeom) ConA=(Con(i)*Con(i-1))**0.5
        if(lGeom) ConB=(Con(i)*Con(i+1))**0.5
        if(lCentrif) Grav=CosAlf*(Radius+abs(x(i)))
        B=(ConA-ConB)*Grav
        if(lDensity) then
          B=(ConA*fRo(1,(Conc(1,i)+Conc(1,i-1))/2.)-
     !       ConB*fRo(1,(Conc(1,i)+Conc(1,i+1))/2.))*Grav
          fRE=fRo(1,Conc(1,i))
        end if
        if(lCentrif) B=B+CosAlf*Con(i)*dx
        if(lVapor) then
          ConA=ConA+(ConVh(i)+ConVh(i-1))/2.
          ConB=ConB+(ConVh(i)+ConVh(i+1))/2.
        end if
        A2=ConA/dxA+ConB/dxB
        A3=-ConB/dxB
        F2=Cap(i)*dx/dt*fRE*Rate
        R(i)=A2+F2
        P(i)=F2*hNew(i)-(ThNew(i)-ThOld(i))*dx/dt*fRE-B-Sink(i)*dx
        if(iDualPor.gt.0) P(i)=P(i)-SinkIm(i)*dx
        S(i)=A3
        if(lVapor.or.lWTDep) then
          ConTA=0.
          if(lVapor) ConTA=ConTA+(ConVT(i)+ConVT(i-1))/2.
          if(lWTDep) ConTA=ConTA+(ConLT(i)+ConLT(i-1))/2.
          ConTB=0.
          if(lVapor) ConTB=ConTB+(ConVT(i)+ConVT(i+1))/2.
          if(lWTDep) ConTB=ConTB+(ConLT(i)+ConLT(i+1))/2.
          P(i)=P(i)+ConTB*(Temp(i+1)-Temp(i))/dxB-
     !              ConTA*(Temp(i)-Temp(i-1))/dxA-
     !              (ThVNew(i)-ThVOld(i))*dx/dt
        end if
11    continue

*     Top BC
      dxA=x(N)-x(N-1)
      dx=dxA/2.
      ConA=(Con(N)+Con(N-1))/2.
      if(lGeom) ConA=(Con(N)*Con(N-1))**0.5
      if(lCentrif) Grav=CosAlf*(Radius+abs((x(N)+x(N-1))/2.))
      B=ConA*Grav
      if(lVapor) ConA=ConA+(ConVh(N)+ConVh(N-1))/2.
      if(lDensity) then
        B=ConA*Grav*(fRo(1,(Conc(1,N)+Conc(1,N-1))/2.))
        fRE=fRo(1,Conc(1,N))
      end if
      F2=Cap(N)*dx/dt*fRE*Rate
      RT=ConA/dxA+F2
      ST=-ConA/dxA
      PT=F2*hNew(N)-(ThNew(N)-ThOld(N))*dx/dt*fRE-Sink(N)*dx-B
      if(iDualPor.gt.0) PT=PT-SinkIm(N)*dx
      if(lVapor.or.lWTDep) then
        ConTA=0.
        if(lVapor) ConTA=ConTA+(ConVT(N)+ConVT(N-1))/2.
        if(lWTDep) ConTA=ConTA+(ConLT(N)+ConLT(N-1))/2.
        PT=PT-ConTA*(Temp(N)-Temp(N-1))/dxA-(ThVNew(N)-ThVOld(N))*dx/dt
      end if
      vTop=-sngl(ST)*hNew(N-1)-sngl(RT)*hNew(N)+sngl(PT)
      PT=PT-rTop
      if(WLayer) then
        if(hNew(N).gt.0.) then
          RT=RT+1./dt
          PT=PT+amax1(hOld(N),0.)/dt
        else
          PT=PT+amax1(hOld(N),0.)/dt
        end if
      end if
      return
      end

************************************************************************

      subroutine Gauss(N,KodTop,KodBot,hTop,hBot,hNew,P,R,S,PB,RB,SB,PT,
     !                 RT,ST,rMin)

      double precision P,R,S,PB,RB,SB,PT,RT,ST,rMin
      dimension hNew(N),P(N),R(N),S(N)

*     Forward
      if(KodBot.ge.0) then
        P(2)=P(2)-S(1)*hBot
      else
        if(dabs(RB).lt.rMin) RB=rMin
        P(2)=P(2)-PB*S(1)/RB
        R(2)=R(2)-SB*S(1)/RB
      end if
      do 11 i=3,N-1
        if(dabs(R(i-1)).lt.rMin) R(i-1)=rMin
        P(i)=P(i)-P(i-1)*S(i-1)/R(i-1)
        R(i)=R(i)-S(i-1)*S(i-1)/R(i-1)
11    continue
      if(KodTop.gt.0) then
        P(N-1)=P(N-1)-S(N-1)*hTop
      else
        if(dabs(R(N-1)).lt.rMin) R(N-1)=rMin
        P(N)=PT-P(N-1)*ST/R(N-1)
        R(N)=RT-S(N-1)*ST/R(N-1)
      end if

*     Back
      if(dabs(R(N-1)).lt.rMin) R(N-1)=rMin
      if(KodTop.gt.0) then
        hNew(N)=hTop
        hNew(N-1)=sngl(P(N-1)/R(N-1))
      else
        hNew(N)=sngl(P(N)/R(N))
        hNew(N-1)=sngl((P(N-1)-S(N-1)*hNew(N))/R(N-1))
      end if
      do 12 i=N-2,2,-1
        if(dabs(R(i)).lt.rMin) R(i)=rMin
        hNew(i)=sngl((P(i)-S(i)*hNew(i+1))/R(i))
12    continue
      if(KodBot.ge.0) then
        hNew(1)=hBot
      else
        if(dabs(RB).lt.rMin) RB=rMin
        hNew(1)=sngl((PB-SB*hNew(2))/RB)
      end if
      do 13 i=1,N
13    continue        
      return
      end

************************************************************************

      subroutine Shift(N,KodTop,rTop,rBot,hTop,hBot,hCritA,CosAlf,
     !                 WLayer,Con,hNew,x,TopInf,KodBot,SeepF,ThNew,
     !                 ThOld,Sink,dt,lVapor,lWTDep,ConLT,Temp,ConVh,
     !                 ConVT,iDualPor,SinkIm,lDensity,Conc,NSD,lCentrif,
     !                 Radius,hSeep)

      dimension Con(N),hNew(N),x(N),ThNew(N),ThOld(N),Sink(N),ConLT(N),
     !          Temp(N),ConVh(N),ConVT(N),SinkIm(N),ConC(NSD,N)
      logical WLayer,TopInf,SeepF,lVapor,lWTDep,lDensity,lCentrif

      fRE=1.
      Grav=CosAlf

*     Seepage face at the bottom
      if(SeepF) then
        dx=x(2)-x(1)
        if(lDensity) fRE=fRo(1,Conc(1,1))
        if(lCentrif) Grav=CosAlf*(Radius+abs((x(2)+x(1))/2.))
        vBot=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/dx+Grav*fRE)-
     !         dx/2.*fRE*((ThNew(1)-ThOld(1))/dt+Sink(1))
        if(KodBot.ge.0) then
          if(vBot.gt.0.) then
            KodBot=-2
            rBot=0.
          end if
        else
          if(hNew(1).ge.hSeep) then
            KodBot=2
            hBot=hSeep
          end if
        end if
      end if

*     Atmospheric boundary condition
      if(TopInf.and.(abs(KodTop).eq.4.or.
     !              (abs(KodTop).eq.1.and.rTop.gt.0.))) then
        if(KodTop.gt.0) then
          M=N-1
          dx=(x(N)-x(M))
          if(lDensity) fRE=fRo(1,Conc(1,N))
          if(lCentrif) Grav=CosAlf*(Radius+abs((x(N)+x(M))/2.))
          vTop=-(Con(N)+Con(M))/2.*((hNew(N)-hNew(M))/dx+Grav*fRE)-
     !           (ThNew(N)-ThOld(N))*fRE*dx/2./dt-Sink(N)*dx/2.
          if(iDualPor.gt.0) vTop=vTop-SinkIm(N)*dx/2.
          if(lWTDep) vTop=vTop-
     !              (ConLT(N)+ConLT(M))/2.*(Temp(N)-Temp(M))/dx
          if(lVapor) vTop=vTop-
     !             (ConVh(N)+ConVh(M))/2.*(hNew(N)-hNew(M))/dx-
     !             (ConVT(N)+ConVT(M))/2.*(Temp(N)-Temp(M))/dx
          if(abs(vTop).gt.abs(rTop).or.vTop*rTop.le.0.) then 
            if(abs(KodTop).eq.4) KodTop=-4
          end if
          if(KodTop.eq.4.and.hNew(N).le.0.99*hCritA.and.rTop.lt.0.)
     !      KodTop=-4
        else
          if(.not.WLayer) then
            if(hNew(N).gt.0.) then
              if(abs(KodTop).eq.4) KodTop=4
              if(abs(KodTop).eq.1) KodTop=1
              hTop=0.
            end if
          end if
          if(hNew(N).le.hCritA) then
            if(abs(KodTop).eq.4) KodTop=4
            if(abs(KodTop).eq.1) KodTop=1
            hTop=hCritA
          end if
        end if
      end if
      return
      end

************************************************************************

      subroutine SetMat(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,
     !                  MatNum,ParD,Con,Cap,ConSat,Ah,AK,ATh,hSat,hTemp,
     !                  TheTab,theta,thr,ths,lWTDep,Temp,Iter,ConO,
     !                  Kappa,AThS,ThRR,ConR,AKS,AhW,AThW,AKW,iModel,
     !                  lTable,lVapor,ThetaV,ConLT,ConVT,ConVh,xConv,
     !                  tConv,nTabMod,hCritA,lDensity,Conc,NSD,iEnhanc)

      logical lWTDep,lTable,lVapor,lDensity
      dimension hTab(NTabD,NMat),ConTab(NTabD,NMat),CapTab(NTabD,NMat),
     !          hNew(NumNP),MatNum(NumNP),ParD(11,NMat),Con(NumNP),
     !          Cap(NumNP),ConSat(NMat),Ah(NumNP),AK(NumNP),
     !          ATh(NumNP),hSat(NMat),hTemp(NumNP),TheTab(NTabD,NMat),
     !          theta(NumNP),thr(NMat),ths(NMat),Temp(NumNP),
     !          ConO(NumNP),Kappa(NumNP),AThS(NumNP),ThRR(NUmNP),
     !          ConR(NumNP),AKS(NumNP),AhW(NMat),AThW(NMat),AKW(NMat),
     !          NTab(NMat),ConLT(NumNP),ConVT(NumNP),ConVh(NumNP),
     !          ThetaV(NumNP),Conc(NSD,NumNP)

      if(iModel.lt.nTabMod) then
        alh1=alog10(-hTab(1,1))
        dlh =(alog10(-hTab(NTab(1),1))-alh1)/(NTab(1)-1)
      end if
      do 11 i=1,NumNP
        AT=1.
        BT=1.
        if(lWTDep) then       ! Temperature dependence
          TempR=20.
          AT=(75.6-0.1425*Temp(i)-2.38e-4*Temp(i)**2)/              ! Surface tension
     !       (75.6-0.1425*TempR  -2.38e-4*TempR  **2)
          BT= (1.787-0.007*TempR  )/(1.+0.03225*TempR  )/           ! Dynamic viscosity
     !       ((1.787-0.007*Temp(i))/(1.+0.03225*Temp(i)))*
     !        (1.-7.37e-6*(Temp(i)-4.)**2+3.79e-8*(Temp(i)-4.)**3)/ ! Density of soil water
     !        (1.-7.37e-6*(TempR  -4.)**2+3.79e-8*(TempR  -4.)**3)
        end if
        if(lDensity) BT=BT*fRo(1,Conc(1,i))/fRo(2,Conc(1,i))        ! Bulk density/dynamic viscosity
        M=MatNum(i)
        if(Kappa(i).eq.-1) then
          hi1=amin1(hSat(M),hTemp(i)/Ah(i)/AT)
          hi2=amin1(hSat(M), hNew(i)/Ah(i)/AT)
        else if(Kappa(i).eq.+1) then
          hi1=amin1(hSat(M),hTemp(i)/Ah(i)/AhW(M)/AT)
          hi2=amin1(hSat(M), hNew(i)/Ah(i)/AhW(M)/AT)
        end if
        hiM=0.1*hi1+0.9*hi2
        if(iModel.lt.nTabMod) then          ! Conductivity
          if(hi1.ge.hSat(M).and.hi2.ge.hSat(M)) then
            Coni=ConSat(M)
          else if(hiM.gt.hTab(NTab(1),1).and.hiM.le.hTab(1,1).and. 
     !                                                      lTable) then
            iT=int((alog10(-hiM)-alh1)/dlh)+1
            dh=(hiM-hTab(iT,1))/(hTab(iT+1,1)-hTab(iT,1))
            Coni=ConTab(iT,M)+(ConTab(iT+1,M)-ConTab(iT,M))*dh
          else
            Coni=FK(iModel,hiM,ParD(1,M))
          end if
        else if(iModel.eq.nTabMod) then     ! Tables
          if(hi1.ge.hSat(M).and.hi2.ge.hSat(M)) then
            Coni=ConSat(M)
          else if(hiM.ge.hTab(NTab(M),M).and.hiM.le.hTab(1,M)) then
            iT=1
            do 12 j=1,NTab(M)-1
              if(hiM.ge.hTab(j+1,M).and.hiM.lt.hTab(j,M)) iT=j
12          continue
            dh=(hiM-hTab(iT,M))/(hTab(iT+1,M)-hTab(iT,M))
            Coni=ConTab(iT,M)+(ConTab(iT+1,M)-ConTab(iT,M))*dh
          else
            if(hiM.gt.hTab(1,M)) then
              Coni=ConTab(1,M)+(ConSat(M)-ConTab(1,M))*
     !                                         (hTab(1,M)-hi2)/hTab(1,M)
            else if(hiM.lt.hTab(NTab(M),M)) then
              Coni=ConTab(NTab(M),M)-ConTab(NTab(M),M)*
     !                       (alog10(-hiM)-alog10(-hTab(NTab(M),M)))/
     !                                (10.-alog10(-hTab(NTab(M),M)))
            end if
          end if
        end if
        if(iModel.lt.nTabMod) then         ! Capacity and water content
          if(hiM.ge.hSat(M)) then
            Capi=0.
            Thei=ths(M)
          else if(hiM.ge.hTab(NTab(1),1).and.hiM.le.hTab(1,1).and.
     !                                                      lTable) then
            iT=int((alog10(-hiM)-alh1)/dlh)+1
            dh=(hiM-hTab(iT,1))/(hTab(iT+1,1)-hTab(iT,1))
            Capi=CapTab(iT,M)+(CapTab(iT+1,M)-CapTab(iT,M))*dh
            Thei=TheTab(iT,M)+(TheTab(iT+1,M)-TheTab(iT,M))*dh
          else
            Capi=FC(iModel,hiM,ParD(1,M))
            Thei=FQ(iModel,hiM,ParD(1,M))
          end if
        else if(iModel.eq.nTabMod) then     ! Tables
          if(hi2.ge.hSat(M)) then
            Capi=0
            Thei=ths(M)
          else if(hi2.ge.hTab(NTab(M),M).and.hi2.le.hTab(1,M)) then
            iT=1
            do 13 j=1,NTab(M)-1
              if(hi2.ge.hTab(j+1,M).and.hi2.le.hTab(j,M)) iT=j
13          continue
            dh=(hi2-hTab(iT,M))/(hTab(iT+1,M)-hTab(iT,M))
            Capi=CapTab(iT,M)+(CapTab(iT+1,M)-CapTab(iT,M))*dh
            Thei=TheTab(iT,M)+(TheTab(iT+1,M)-TheTab(iT,M))*dh
          else
            if(hi2.gt.hTab(1,M)) then
              Capi=CapTab(1,M)*hi2/hTab(1,M)
              Thei=TheTab(1,M)+(ths(M)-TheTab(1,M))*
     !                                         (hTab(1,M)-hi2)/hTab(1,M)
            else if(hi2.lt.hTab(NTab(M),M)) then
              Capi=CapTab(NTab(M),M)-CapTab(NTab(M),M)*
     !                       (alog10(-hi2)-alog10(-hTab(NTab(M),M)))/
     !                                 (6.-alog10(-hTab(NTab(M),M)))
              Thei=thr(M)+(TheTab(NTab(M),M)-thr(M))*(hi2+1.e+6)/
     !                                   (hTab(NTab(M),M)+1.e+6)
            end if
          end if
        end if
        if(Kappa(i).eq.-1) then  ! Drying
          Con(i)=Coni*AK(i)*BT*AKS(i)
          Cap(i)=Capi*ATh(i)*AThS(i)/Ah(i)/AT
          theta(i)=thr(M)+(Thei-thr(M))*ATh(i)*AThS(i)
        else                     ! Wetting
          Con(i)=ConR(i)+Coni*AK(i)*BT*AKS(i)*AKW(M)
          Cap(i)=Capi*ATh(i)*AThS(i)*AThW(M)/Ah(i)/AhW(M)/AT
          theta(i)=ThRR(i)+AThW(M)*ATh(i)*AThS(i)*(Thei-thr(M))
        end if
        if(Iter.eq.0) ConO(i)=Con(i)
11    continue
      if(lVapor.or.lWTDep) then  ! Vapor flow
        call ConVapor(NumNP,NMat,MatNum,hNew,Temp,Con,Theta,ths,ConLT,
     !                ConVT,ConVh,xConv,tConv,lVapor,hCritA,iEnhanc)
      if(lVapor)
     !    call VaporContent(NumNP,NMat,Theta,ThetaV,Temp,hNew,MatNum,
     !                      ths,xConv)
      end if
      return
      end

************************************************************************

      real function Fqh(GWL,Aqh,Bqh)
      Fqh=Aqh*exp(Bqh*abs(GWL))
      return
      end

************************************************************************

      function FqDrain(GWL,zBotDr,BaseGW,rSpacing,iPosDr,KhTop,KhBot,
     !                 KvTop,KvBot,Entres,WetPer,zInTF,GeoFac)

*    -------------------------------------------------------------------
*     Purpose: determines the drainage flux
*     Based on the SWAP model by van Dam et al. 1997

*     iPosDr   kod for the position of the drain.......................I
*              =1: Homogeneous profile, drain on top of impervious layer
*              =2: Homogeneous profile, drain above impervious layer
*              =3: Heterogeneous profile, drain at interface between
*                                         both soil layers
*              =4: Heterogeneous profile, drain in bottom layer
*              =5: Heterogeneous profile, drain in top layer
*     Input    Calculation: GWL,Pond
*              Common: zBotDr,rSpacing,Entres
*              iPosDr=1: KhTop
*              iPosDr=2: BaseGW,KhTop,WetPer
*              iPosDr=3: BaseGW,KhTop,KhBot,WetPer
*              iPosDr=4: BaseGW,KvTop,KvBot,KhBot,WetPer,zInTF
*              iPosDr=5: BaseGW,KhTop,KvTop,KhBot,WetPer,zInTF,GeoFac
*     GWL      ground water level...................................I(s)
*     zBotDr   coordinate of the bottom of the drainage................I
*     Pond     Ponding (cm).........................................I(s)
*     BaseGW   coordinate of the impervious layer cm (2,3,4,5).........I
*     rSpacing drain spacing (1,2,3,4,5)...............................I
*     KhTop    horizontal saturated hydraulic conductivity above drain
*              (cm/d) (1,2,3,5)........................................I
*     KhBot    horizontal saturated hydraulic conductivity below drain
*              (cm/d) (3,4,5)..........................................I
*     KvTop    vertical saturated hydraulic conductivity above drain
*              (cm/d) (4,5)............................................I
*     KvBot    vertical saturated hydraulic conductivity below drain
*              (cm/d) (4)..............................................I
*     Entres   entrance resistance into the drain and/or ditches (d)
*              (1,2,3,4,5).............................................I
*     WetPer   wet perimeter (cm) (2,3,4,5)............................I
*     zInTF    level of the transition between the upper and lower
*              soil layer (cm) (4,5)...................................I
*     GeoFac   geometry factor (5).....................................I
*     FqDrain  drainage flux (cm/d)....................................O

*     dh       hydraulic difference between drain and the middle of
*              the spacing
*     zImp     adjusted coordinate of the impervious layer cm (2,3,4,5)
*     dBot     depth to the impervious layer below the drain
*     EqD      equilvalent depth (cm)   
*     TotRes   total drainage resistance
*     RVer     vertical drainage resistance
*     RHor     horizontal drainage resistance
*     RRad     radial drainage resistance
*     x        typical length variable 
*     ------------------------------------------------------------------  
*     global variables
      integer iPosDr
      real GWL,zBotDr,BaseGW,rSpacing,FqDrain,KhTop,KhBot,KvTop,KvBot
      real Entres,WetPer,zInTF,GeoFac

*     local variables
      integer i
      real dh,zImp,dBot,pi,TotRes,x,fx,EqD,RVer,RHor,RRad
      parameter(pi=3.14159)

*     drainage flux calculated according to Hooghoudt or Ernst
      dh=GWL-zBotDr

*     contributing layer below drains limited to 1/4 L
      if(iPosDr.gt.1) then
        zImp=max(BaseGW,zBotDr-0.25*rSpacing)
        dBot=(zBotDr-zImp)  
        if(dBot.lt.0.0) STOP 'Error - Bocodrb: dBot negative'
      end if

*     no infiltration allowed
      if(dh.lt.1.0e-10) then
        FqDrain=0.0
        return
      end if

*     case 1: homogeneous, on top of impervious layer
      if(iPosDr.eq.1) then
*     calculation of drainage resistance and drainage flux
        TotRes=rSpacing*rSpacing/(4.*KhTop*abs(dh))+Entres                ! Eq. 8.7

*     case 2,3: in homogeneous profile or at interface of 2 layers
      else if(iPosDr.eq.2.or.iPosDr.eq.3) then
*       calculation of equivalent depth
        x=2.*pi*dBot/rSpacing                                             ! Eq. 8.9
        if(x.gt.0.5) then
          fx=0.0
          do 10 i=1,5,2
            fx=fx+(4.*exp(-2.*i*x))/(i*(1.0-exp(-2.*i*x)))                ! Eq. 8.13
10        continue
          EqD=pi*rSpacing/8./(log(rSpacing/WetPer)+fx)                    ! Eq. 8.12
        else
          if(x.lt.1.0E-6) then
            EqD=dBot
          else
            fx=pi**2/(4.*x)+log(x/(2*pi))                                 ! Eq. 8.11
            EqD=pi*rSpacing/8./(log(rSpacing/WetPer)+fx)                  ! Eq. 8.12
          end if
        end if
        if(EqD.gt.dBot) EqD=dBot

*       calculation of drainage resistance & drainage flux
        if(iPosDr.eq.2) then
          TotRes=rSpacing*rSpacing/(8.*KhTop*EqD+4*KhTop*abs(dh))+Entres  ! Eq. 8.8
        else if(iPosDr.eq.3) then
          TotRes=rSpacing*rSpacing/(8.*KhBot*EqD+4*KhTop*abs(dh))+Entres  ! Eq. 8.14
        end if

*     case 4: drain in bottom layer
      else if(iPosDr.eq.4) then
        if(zBotDr.gt.zInTF) stop 'Error - check zInTF and zBotDr' 
        RVer=max(GWL-zInTF,0.)/KvTop +(min(zInTF,GWL)-zBotDr)/KvBot       ! Eq. 8.16
        RHor=rSpacing*rSpacing/(8.*KhBot*dBot)                            ! Eq. 8.17
        RRad=rSpacing/(pi*sqrt(KhBot*KvBot))*log(dBot/WetPer)             ! Eq. 8.18
        TotRes=RVer+RHor+RRad+Entres                                      ! Eq. 8.15

*     case 5: drain in top layer
      else if(iPosDr.eq.5) then
        if(zBotDr.lt.zInTF) stop 'Error - check zInTF and zBotDr'
        RVer=(GWL-zBotDr)/KvTop                                           ! Eq. 8.19
        RHor=rSpacing*rSpacing/(8.*KhTop*(zBotDr-zInTF)+                  ! Eq. 8.20
     !                          8.*KhBot*(zInTF-zImp))
        RRad=rSpacing/(pi*sqrt(KhTop*KvTop))*log((GeoFac*                 ! Eq. 8.21
     !       (zBotDr-zInTF))/WetPer)
        TotRes=RVer+RHor+RRad+Entres                                      ! Eq. 8.15
      end if
      FqDrain=-dh/TotRes                                                  ! Eq. 8.6

      return
      end

************************************************************************

*     To calculate the velocities

      subroutine Veloc(N,hNew,Con,x,CosAlf,v,ThNew,ThOld,Sink,dt,lVapor,
     !                 lWTDep,ConLT,ConVT,ConVh,Temp,vV,ThVNew,ThVOld,
     !                 lDensity,Conc,NSD,lCentrif,Radius)

      logical lVapor,lWTDep,lDensity,lCentrif
      dimension hNew(N),x(N),Con(N),v(N),ConLT(N),ConVT(N),ConVh(N),
     !          Temp(N),vV(N),Conc(NSD,N),ThNew(N),ThOld(N),Sink(N),
     !          ThVNew(N),ThVOld(N)

      fRE=1.
      Grav=CosAlf
      M=N-1
      dxN=x(N)-x(M)
      if(lDensity) fRE=(fRo(1,Conc(1,N))+fRo(1,Conc(1,M)))/2.
      if(lCentrif) Grav=CosAlf*(Radius+abs((x(N)+x(M))/2.))
      v(N)=-(Con(N)+Con(M))/2.*((hNew(N)-hNew(M))/dxN+fRE*Grav)-
     !     dxN/2.*(fRE*(ThNew(N)-ThOld(N))/dt+Sink(N))
      if(lWTDep) v(N)=v(N)-(ConLT(N)+ConLT(M))/2.*(Temp(N)-Temp(M))/dxN
      vV(N)=0.
      if(lVapor) vV(N)=-(ConVh(N)+ConVh(M))/2.*(hNew(N)-hNew(M))/dxN-
     !                  (ConVT(N)+ConVT(M))/2.*(Temp(N)-Temp(M))/dxN-
     !                  dxN/2.*(ThVNew(N)-ThVOld(N))/dt
      do 11 i=2,N-1
        dxA=x(i+1)-x(i)
        dxB=x(i)-x(i-1)
        if(lDensity) fRE=(fRo(1,Conc(1,i))+fRo(1,Conc(1,i+1)))/2.
        if(lCentrif) Grav=CosAlf*(Radius+abs((x(i+1)+x(i))/2.))
        vA=-(Con(i)+Con(i+1))/2.*((hNew(i+1)-hNew(i))/dxA+fRE*Grav)
        if(lDensity) fRE=(fRo(1,Conc(1,i))+fRo(1,Conc(1,i-1)))/2.
        if(lCentrif) Grav=CosAlf*(Radius+abs((x(i)+x(i-1))/2.))
        vB=-(Con(i)+Con(i-1))/2.*((hNew(i)-hNew(i-1))/dxB+fRE*Grav)
        v(i)=(vA*dxB+vB*dxA)/(dxA+dxB)
        if(lWTDep) then
          vTA=-(ConLT(i)+ConLT(i+1))/2.*(Temp(i+1)-Temp(i))/dxA
          vTB=-(ConLT(i)+ConLT(i-1))/2.*(Temp(i)-Temp(i-1))/dxB
          v(i)=v(i)+(vTA*dxB+vTB*dxA)/(dxA+dxB)
        end if
        vV(i)=0.
        if(lVapor) then
          vVA=-(ConVh(i)+ConVh(i+1))/2.*(hNew(i+1)-hNew(i))/dxA
          vVB=-(ConVh(i)+ConVh(i-1))/2.*(hNew(i)-hNew(i-1))/dxB
          vVA=vVA-(ConVT(i)+ConVT(i+1))/2.*(Temp(i+1)-Temp(i))/dxA
          vVB=vVB-(ConVT(i)+ConVT(i-1))/2.*(Temp(i)-Temp(i-1))/dxB
          vV(i)=(vVA*dxB+vVB*dxA)/(dxA+dxB)
        end if
11    continue
      dx1=x(2)-x(1)
      if(lDensity) fRE=(fRo(1,Conc(1,2))+fRo(1,Conc(1,1)))/2.
      if(lCentrif) Grav=CosAlf*(Radius+abs((x(2)+x(1))/2.))
      v(1)=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/dx1+fRE*Grav)+
     !      dx1/2.*(fRE*(ThNew(1)-ThOld(1))/dt+Sink(1))
      if(lWTDep) v(1)=v(1)-(ConLT(1)+ConLT(2))/2.*(Temp(2)-Temp(1))/dx1
      vV(1)=0.
      if(lVapor) vV(1)=   -(ConVh(1)+ConVh(2))/2.*(hNew(2)-hNew(1))/dx1-
     !                     (ConVT(1)+ConVT(2))/2.*(Temp(2)-Temp(1))/dx1+
     !                      dx1/2.*(ThVNew(1)-ThVOld(1))/dt
      return
      end

************************************************************************

      subroutine Hyster(NumNP,NMat,hOld,MatNum,ParD,ParW,ThNew,ThOld,
     !                  Kappa,AThS,ThRR,ConO,ConR,AKS,KappaO,Ah,AK,
     !                  iHyst,iModel,TolTh)

      real KsD,KsW,Ks,KW

      dimension MatNum(NumNP),ThOld(NumNP),hOld(NumNP),ParD(11,NMat),
     !          ThNew(NumNP),Kappa(NumNP),AThS(NumNP),ThRR(NUmNP),
     !          ConO(NumNP),ConR(NumNP),AKS(NUmNP),KappaO(NumNP),
     !          Ah(NumNP),AK(NumNP),ParW(11,NMat)

      do 11 i=1,NumNP

*       Check for reversal
        KappaO(i)=Kappa(i)
        if((ThNew(i)-ThOld(i))*Kappa(i).ge.-TolTh/1.) goto 11
        Kappa(i)=-Kappa(i)
        m=MatNum(i)
        Thr =ParD(1,m)
        ThsD=ParD(2,m)
        ThsW=ParW(2,m)
        KsD =ParD(5,m)
        KsW =ParW(5,m)

*       Update Ths and Ks for wetting scanning curve
        if(Kappa(i).eq.1) then
          if(ThsW.ge.0.999*ThsD) then
            Ths=ThsD
          else
            RR=1./(ThsD-ThsW)-1./(ThsD-Thr)                    ! Eq. 8
            Ths=ThsD-(ThsD-ThOld(i))/(1.+RR*(ThsD-ThOld(i)))   ! Eq. 8
          end if
          if(KsW.ge.0.999*KsD) then
            Ks=KsD
          else
            RR=1./(KsD-KsW)-1./KsD                             ! Eq. 13
            Ks=KsD-(KsD-ConO(i))/(1.+RR*(KsD-ConO(i)))         ! Eq. 13
          end if
        end if

*       Update parameters for scanning curve
        if(Kappa(i).eq.1) then ! Wetting
          AThS(i)=1.
          SeW=FS(iModel,hOld(i)/Ah(i),ParW(1,m))
          if(SeW.lt.0.999) AThS(i)=(ThOld(i)-Ths)/(1.-SeW)/(Thr-ThsW) ! Eq. 7
          ThRR(I)=Ths-AThS(i)*(ThsW-Thr)                       ! Eq. 7a
          AKS(i)=1.
          ConR(i)=0.
          if(iHyst.eq.2) then
            KW=AK(i)*FK(iModel,hOld(i)/Ah(i),ParW(1,m))
            if(KW.lt.0.999*KsW) AKS(i)=(ConO(i)-Ks)/(KW-KsW)   ! Eq. 12
            ConR(i)=Ks-AKS(i)*KsW                              ! Eq. 12a
          end if
        else ! Drying
          AThS(i)=(ThOld(i)-Thr)/FS(iModel,hOld(i)/Ah(i),ParD(1,m))/
     !            (ThsD-Thr)                                   ! Eq. 5
          ThRR(i)=Thr
          AKS(i)=1.
          ConR(i)=0.
          if(iHyst.eq.2) AKS(i)=ConO(i)/FK(iModel,hOld(i)/Ah(i),
     !                          ParD(1,m))/AK(i)               ! Eq. 10
        end if
11    continue
      return
      end

************************************************************************

*     To calculate isothermal vapor hydraulic conductivity, and 
*     thermal vapor and liquid hydraulic conductivities

      subroutine ConVapor(N,NMat,MatNum,hNew,Temp,Con,Theta,ths,ConLT,
     !                    ConVT,ConVh,xConv,tConv,lVapor,hCritA,iEnhanc)

*     ConVh - Conductivity for vapor phase due to gradient of h [m/s]
*     ConVT - Conductivity for vapor phase due to gradient of T [m2/s/K]
*     ConLT - Conductivity for liquid phase due to gradient of T [m2/s/K]
*     Gwt   - Gain factor [-]
*     Gamma - surface tension [N/m,J/m2], [g/s2]
*     dGamma - derivative of surface tension versus temperature [g/s2/K]
*     Diff0 - diffusivity of water vapor in air [m2/s] (2.12e-5)
*     Tau   - tortuosity [-]
*     row   - density of soil water [kg/m3]
*     rovs  - saturated vapor density [kg/m3]
*     drovs - derivative of saturated vapor density versus temp [kg/m3/K]
*     g     - gravitational acceleration [m/s2] (9.81)
*     xMol  - molecular weight of water [kg/mol] (0.018015)
*     R     - universal gas constant [J/mol/K] (8.314)
*     Hr    - relative humidity [-]
*     eta   - enhancement factor [-]
*     fc    - mass fraction of clay in soil (0.02)
*     T     - temperature [C]

      logical lVapor,lLimit
      dimension hNew(N),Temp(N),Con(N),Theta(N),MatNum(N),ths(NMat),
     !          ConLT(N),ConVT(N),ConVh(N)

      data Gwt, Diff0 ,   g ,   xMol  ,   R  , fc  ,Gamma0
     !    / 7.,2.12e-5, 9.81, 0.018015, 8.314, 0.02, 71.89/

      do 11 i=1,N
        h=hNew(i)/xConv                         ! Conversion to m
        ConLh=Con(i)/xConv*tConv                ! Conversion to m/s
        lLimit=.false.
        if(hNew(N).lt.0.99*hCritA)
     !    lLimit=.true.
        T=Temp(i)
        M=MatNum(i)
        ThetaS=ths(M)

        Gamma=75.6-0.1425*T-0.000238*T*T
        dGamma=-0.1425-0.000479*T
        ConLT(i)=ConLh*h*Gwt*dGamma/Gamma0

        if(lVapor) then
          TKelv=T+273.15
          DiffT=Diff0*(TKelv/273.15)**2
          ThetaA=ThetaS-Theta(i)
          Tau=ThetaA**(7./3.)/ThetaS**2  ! Millington & Quirk
          Diff=Tau*ThetaA*DiffT
          row=(1.-7.37e-6*(T-4.)**2+3.79e-8*(T-4.)**3)*1000.
          rovs=0.001*exp(31.3716-6014.79/TKelv-0.00792495*Tkelv)/TKelv
          Hr=exp(h*xMol*g/R/TKelv)
          if(lLimit) Hr=0.000001           ! ##runs fast
          ConVh(i)=Diff/row*rovs*xMol*g/R/TKelv*Hr

          TKelv1=TKelv+1.
          rovs1=0.001*exp(31.3716-6014.79/TKelv1-0.00792495*TKelv1)/
     !                                                           TKelv1
          drovs=rovs1-rovs
          eta=1.
          if(iEnhanc.eq.1) then
            eta=9.5+3.*Theta(i)/ThetaS-
     !                  8.5*exp(-((1.+2.6/sqrt(fc))*Theta(i)/ThetaS)**4)
          end if
          ConVT(i)=Diff/row*eta*Hr*drovs

*         Conversions to HYDRUS units
          ConVh(i)=ConVh(i)*xConv      /tConv
          ConVT(i)=ConVT(i)*xConv*xConv/tConv
        end if
        ConLT(i)=ConLT(i)*xConv*xConv/tConv
11    continue

      return
      end

*************************************************************************

      subroutine VaporContent(NumNP,NMat,Theta,ThetaV,Temp,hNew,MatNum,
     !                        ths,xConv)

      dimension Theta(NumNP),ThetaV(NumNP),Temp(NumNP),hNew(NumNP),
     !          MatNum(NumNP),ths(NMat)

*     g     - gravitational acceleration [m/s2] (9.81)
*     xMol  - molecular weight of water [kg/mol] (0.018015)
*     R     - universal gas constant [J/mol/K] (8.314)
*     ThetaV- volumetric vapor content expressed as an equivalent water content
*     row   - density of soil water [kg/m3]
*     rovs  - saturated vapor density [kg/m3]
*     rov   - vapor density [kg/m3]
*     Hr    - relative humidity [-]

      g=9.81
      xMol=0.018015
      R=8.314

      do 11 i=1,NumNP
        h=hNew(i)/xConv     ! Conversion to m
        M=MatNum(i)
        TKelv=Temp(i)+273.15
        rovs=0.001*exp(31.3716-6014.79/TKelv-0.00792495*Tkelv)/TKelv
        Hr=exp(h*xMol*g/R/TKelv)
        rov=rovs*Hr
        row=(1.-7.37e-6*(Temp(i)-4.)**2+3.79e-8*(Temp(i)-4.)**3)*1000.
        ThetaV(i)=rov*(ths(M)-Theta(i))/row
11    continue

      return
      end

************************************************************************

*     To calculate isothermal vapor hydraulic conductivity

      real function ConVh(hNew,Theta,ths,xConv,tConv)

*     ConVh - Conductivity for vapor phase due to gradient of h [m/s]
*     Diff0 - diffusivity of water vapor in air [m2/s] (2.12e-5)
*     Tau   - tortuosity [-]
*     row   - density of soil water [kg/m3]
*     rovs  - saturated vapor density [kg/m3]
*     g     - gravitational acceleration [m/s2] (9.81)
*     xMol  - molecular weight of water [kg/mol] (0.018015)
*     R     - universal gas constant [J/mol/K] (8.314)
*     Hr    - relative humidity [-]
*     Temp  - temperature [C]

      data Diff0 ,   g ,   xMol  ,   R
     !    /2.12e-5, 9.81, 0.018015, 8.314/

      Temp=20.   
      h=hNew/xConv                       ! Conversion to m
      TKelv=Temp+273.15
      DiffT=Diff0*(TKelv/273.15)**2
      ThetaA=ths-Theta
      if(ThetaA.gt.0.) Tau=ThetaA**(7./3.)/ths**2         ! Millington & Quirk
      Diff=Tau*ThetaA*DiffT
      row=(1.-7.37e-6*(Temp-4.)**2+3.79e-8*(Temp-4.)**3)*1000.
      rovs=0.001*exp(31.3716-6014.79/TKelv-0.00792495*Tkelv)/TKelv
      Hr=exp(h*xMol*g/R/TKelv)
      ConVh=Diff/row*rovs*xMol*g/R/TKelv*Hr
      ConVh=ConVh*xConv/tConv            ! Conversions to HYDRUS units
      return
      end

*************************************************************************

      real function fRo(iKod,Conc)


*     Ratio of bulk densities (dynamic viscosities) at given and zero 
*     concentrations

      fRo=1.

      if(iKod.eq.1) then            ! bulk density
        a1=1.
        a2=0.75
        a3=0.
        a4=0.
      else if(iKod.eq.2) then       ! dynamic viscosity
        a1=1.
        a2=0.
        a3=0.
        a4=0.
      end if
      fRo=a1+a2*Conc+a3*Conc**2+a4*Conc**3

      return
      end

*************************************************************************

      subroutine DualPor(N,NMat,MatNum,iDualPor,ThOld,ths,thr,ThNewIm,
     !                   ThOldIm,ParD,SinkIm,dt,iModel,hNew,hCritA,x,
     !                   WTransf)

      dimension MatNum(N),ThOld(N),thr(NMat),ths(NMat),ThNewIm(N),
     !          ParD(11,NMat),SinkIm(N),hNew(N),ThOldIm(N),x(N),Par(10)

      WTransf=0.
      do 11 i=1,N
        M=MatNum(i)
        SeIm=min(1.,(ThOldIm(i)-ParD(7,M))/(ParD(8,M)-ParD(7,M)))
        if(iDualPor.eq.1) then                   ! Water Content driven exchange
          Se=min(1.,(ThOld(i)-thr(M))/(ths(M)-thr(M)))
          SinkIm(i)=ParD(9,M)*(Se-SeIm)
          DeltaTh=(Se-SeIm)/(ths(M)-thr(M)+ParD(8,M)-ParD(7,M))*
     !            (ths(M)-thr(M))*(ParD(8,M)-ParD(7,M))
        else if(iDualPor.eq.2) then              ! Pressure head driven exchange
          h=FH(iModel,SeIm,ParD(7,M))
          Par(1)=ParD(7,M)                     ! Variable coeff. (function of h)
          Par(2)=ParD(8,M)
          Par(3)=ParD(9,M)
          Par(4)=ParD(10,M)
          Par(5)=ParD(11,M)
          Par(6)=ParD(6,M)
          CondM=FK(iModel,      h,Par)
          CondF=FK(iModel,hNew(i),Par)
          SinkIm(i)=0.5*(CondM+CondF)*(hNew(i)-h)
          if(i.eq.N.and.abs(hCritA-hNew(N)).lt.-0.001*hCritA.and.
     !                  abs(hCritA-h).lt.-0.01*hCritA) SinkIm(i)=0.
        end if
        if(SinkIm(i).gt.0.) then
          TrMaxIm=(ParD(8,M)-ThOldIm(i))/dt
          if(iDualPor.eq.1) TrMaxIm=min(DeltaTh/dt,TrMaxIm)
          if(SinkIm(i).gt.TrMaxIm) SinkIm(i)=TrMaxIm
        end if
        if(SinkIm(i).lt.0.) then
          TrMaxIm=-(ths(M)-ThOld(i))/dt
          if(iDualPor.eq.1) TrMaxIm=max(DeltaTh/dt,TrMaxIm)
          if(SinkIm(i).lt.TrMaxIm) SinkIm(i)=TrMaxIm
        end if

        ThNewIm(i)=max(min(ThOldIm(i)+SinkIm(i)*dt,ParD(8,M)),ParD(7,M))

        if(i.ge.2) 
     !    WTransf=WTransf+(SinkIm(i-1)+SinkIm(i))/2.*(x(i)-x(i-1))
11    continue

      return
      end

************************************************************************

      subroutine Update(NumNP,lWat,lChem,lTemp,lVapor,iDualPor,lExtrap,
     !                  dt,dtOld,hTemp,hNew,hOld,ThOld,ThNew,vOld,vNew,
     !                  ThVOld,ThVNew,vVOld,vVNew,ThOldIm,ThNewIm,TempO,
     !                  TempN,rTop,xConv,ConSMax,KodTop,KodBot)

      logical lWat,lChem,lTemp,lVapor,lExtrap,lSat
      dimension hTemp(NumNP),hNew(NumNP),hOld(NumNP),ThOld(NumNP),
     !          ThNew(NumNP),vOld(NumNP),vNew(NumNP),ThVOld(NumNP),
     !          ThVNew(NumNP),vVOld(NumNP),vVNew(NumNP),ThOldIm(NumNP),
     !          ThNewIm(NumNP),TempO(NumNP),TempN(NumNP)

      lSat=.true.
      iBot=1
      if(KodBot.gt.0) iBot=2
      iTop=NumNP
      if(KodTop.gt.0) iTop=NumNP-1
      do 11 i=iBot,iTop
        if(lWat) then
          if(lExtrap.and.hNew(i).lt.0..and.hOld(i).lt.0.) then
            hTemp(i)=hNew(i)+(hNew(i)-hOld(i))*dt/dtOld
          else
            hTemp(i)=hNew(i)
          end if
          hOld(i) =hNew(i)
          hNew(i) =hTemp(i)
        end if
11    continue
      do 12 i=1,NumNP
        if(lWat) then
          ThOld(i)=ThNew(i)
          if(lTemp.or.lChem) vOld(i)=vNew(i)
          if(lVapor)         ThVOld(i)=ThVNew(i)
          if(lVapor)         vVOld(i) =vVNew(i)
          if(iDualPor.gt.0)  ThOldIm(i)=ThNewIm(i)
          if(hNew(i).lt.0.) lSat=.false.
        end if
        if(lTemp) TempO(i)=TempN(i)
12    continue
      if(lWat.and.lSat.and.(KodTop.eq.-1.or.KodTop.eq.4).
     !                                        and.rTop.ge.-ConSMax) then
        if(KodTop.eq.4) KodTop=-4
        hNew (NumNP)=-0.005*xConv      ! 0.5 cm
        hTemp(NumNP)=-0.005*xConv
        hOld (NumNP)=-0.005*xConv
      end if

      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||