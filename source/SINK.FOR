* Source file SINK.FOR |||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine SetSnk(N,NMat,MatNum,x,hRoot,vRoot,Sink,TPot,hNew,
     !                  lMoSink,lSolRed,lSolAdd,P0,POptm,P2H,P2L,P3,r2H,
     !                  r2L,aOsm,c50,P3c,Beta,lChem,NS,NSD,Conc,cRoot,
     !                  lMsSink,ThNew,ParD,dt,OmegaC,iModel,Con,lOmegaW,
     !                  OmegaW,rBot)

      logical lChem,lMoSink,lSolRed,lSolAdd,lMsSink,lHanks,lOmegaW
      dimension x(N),MatNum(N),hNew(N),POptm(NMat),Beta(N),Sink(N),
     !          Conc(NSD,N),cRoot(NS),aOsm(NS),ThNew(N),ParD(11,NMat),
     !          Con(N)

      Compen=1.
      nStep=1
      if(OmegaC.lt.1.) nStep=2
      Omega=0.

      vRoot=0.
      hRoot=0.
      ARoot=0.
      do 11 ii=1,NS
        cRoot(ii)=0.
11    continue

      lHanks=.false.
      if(lHanks) then
        hR=P3
        do 22 iter=1,5
          Sum1=0.
          Sum2=0.
          xConst=1.  ! Penalty for distance
          do 21 i=2,N
            if(Beta(i).gt.0.) then
              if(i.eq.N) then
                dxM=(x(i)-x(i-1))/2.
              else
                dxM=(x(i+1)-x(i-1))/2.
              end if
              if(hNew(i).gt.hR-xConst*x(i)) then
                Sum1=Sum1+Con(i)*Beta(i)*(hNew(i)+xConst*x(i))*dxM
                Sum2=Sum2+Con(i)*Beta(i)*dxM
              end if
            end if
21        continue
          hR=max((Sum1-TPot)/Sum2,P3)
22      continue
        do 23 i=1,N
          if(Beta(i).gt.0.) then
            if(i.eq.N) then
              dxM=(x(i)-x(i-1))/2.
            else
              dxM=(x(i+1)-x(i-1))/2.
            end if
            Sink(i)=max(-Con(i)*Beta(i)*(hR-xConst*x(i)-hNew(i)),0.)
c            Sink(i)=-Con(i)*Beta(i)*(hR-xConst*x(i)-hNew(i))
            vRoot=vRoot+Sink(i)*dxM
            hRoot=hRoot+hNew(i)*dxM
            ARoot=ARoot+dxM
          end if
23      continue
        if(ARoot.gt.0.001) hRoot=hRoot/ARoot
        return
      end if

      do 16 iStep=1,nStep
        do 13 i=2,N
          if(Beta(i).gt.0.) then
            if(i.eq.N) then
              dxM=(x(i)-x(i-1))/2.
            else
              dxM=(x(i+1)-x(i-1))/2.
            end if
            M=MatNum(i)
            hRed=hNew(i)
            SAlfa=1.
            if(lChem.and.lSolRed) then
              cRed=0.
              do 15 j=1,NS
                cRed=cRed+aOsm(j)*Conc(j,i)
15            continue
              if(lSolAdd) then
                hRed=hRed+cRed
              else
                SAlfa=FSAlfa(lMsSink,cRed,c50,P3c)
              end if
            end if
            Alfa=FAlfa(lMoSink,TPot,hRed,P0,POptm(M),P2H,P2L,P3,r2H,r2L)
            if(iStep.ne.nStep) then
              Omega=Omega+Alfa*SAlfa*Beta(i)*dxM
              goto 13
            else
              Compen=1.
              if(Omega.lt.OmegaC.and.Omega.gt.0.) Compen=OmegaC
              if(Omega.ge.OmegaC)                 Compen=Omega
            end if
            Sink(i)=Alfa*SAlfa*Beta(i)*TPot/Compen
            if(ThNew(i)-0.00025.lt.ParD(1,MatNum(i))) Sink(i)=0.
            if(lMoSink) PMin=P3
            if(.not.lMosink) PMin=10.*P0
            ThLimit=FQ(iModel,PMin,ParD(1,MatNum(i)))
c            Sink(i)=min(Sink(i),0.5*(ThNew(i)-ParD(1,MatNum(i)))/dt)
            Sink(i)=min(Sink(i),max(0.,0.5*(ThNew(i)-ThLimit)/dt))
            vRoot=vRoot+Sink(i)*dxM
            hRoot=hRoot+hNew(i)*dxM
            do 12 ii=1,NS
              if(lChem) cRoot(ii)=cRoot(ii)+Conc(ii,i)*dxM
12          continue
            ARoot=ARoot+dxM
          else
            Sink(i)=0.
          end if
          if(Beta(i).lt.0.) then ! Eddy Woehling's modification
            if(i.eq.N) then
              dxM=(x(i)-x(i-1))/2.
            else
              dxM=(x(i+1)-x(i-1))/2.
            end if
            Sink(i)=Beta(i)*rBot
            Sink(i)=max(Sink(i),0.5*(ThNew(i)-ParD(2,MatNum(i)))/dt)
          end if
13      continue
16    continue
      if(ARoot.gt.0.001) then
        hRoot=hRoot/ARoot
        do 14 ii=1,NS
          cRoot(ii)=cRoot(ii)/ARoot
14      continue
      end if
      if(lOmegaW.and.TPot.gt.0.) OmegaW=vRoot/TPot
      return
      end

************************************************************************

*     Subroutine calculating root solute uptake with and without compensation

      subroutine SetSSnk(jS,NS,N,t,x,Beta,Sink,SinkS,NSD,Conc,OmegaW,
     !                   cRootMax,lActRSU,OmegaS,SPot,rKM,cMin)

      double precision t
      logical lActRSU,lLast
      dimension x(N),Beta(N),Sink(N),SinkS(N),Conc(NSD,N)

*     Inputs:
*     SPot      - potential root solute uptake
*     OmegaS    - solute stress index
*     rKM       - Michaelis-Menten constant
*     lActRSU   - consider active root solute uptake
*     cRootMax  - maximum concentration for the passive solute uptake

*     From Water Flow
*     Sink(i)   - Root water uptake
*     OmegaW    - ratio of actual and potential transpiration

*     SPUptake  - passive root solute uptake (step 1)
*     SAUptakeP - potential active solute uptake (step 1)
*     SAUptakeA - uncompensated actual active solute uptake (step 2)
*     SAUptakeA - compensated actual active solute uptake (step 3)
*     SinkS(i)  - local active solute uptake

*     Initialization
      Compen=1.
      nStep=1
      if(lActRSU)                  nStep=2
      if(lActRSU.and.OmegaS.lt.1.) nStep=3
*     step 1: Passive uptake
*     step 2: Active uptake without compensation
*     step 3: Active uptake with compensation
      lLast=.false.                 ! Active uptake only for the last solute
      if(lLast.and.jS.lt.NS) nStep=1
      Omega=0.
      SPUptake=0.
      do 10 i=1,N
        SinkS(i)=0.
10    continue

      do 12 iStep=1,nStep
        SAUptakeA=0.
        do 11 i=1,N
          if(Beta(i).gt.0.) then
            if(i.eq.N) then
              dxM=(x(i)-x(i-1))/2.
            else if(i.eq.1) then
              dxM=(x(i)-x(i+1))/2.
            else
              dxM=(x(i+1)-x(i-1))/2.
            end if
            cc=amax1(Conc(jS,i)-cMin,0.)
            if(iStep.eq.1) then
              SinkS(i)=Sink(i)*amax1(amin1(Conc(jS,i),cRootMax),0.)
              SPUptake=SPUptake+SinkS(i)*dxM
*             This is needed only for the last node, but that node may not have beta
              SAUptakeP=amax1(SPot*OmegaW-SPUptake,0.)
            else if(iStep.eq.2) then
              AUptakeA=cc/(rKM+cc)*Beta(i)*SAUptakeP
              Omega=Omega+AUptakeA*dxM
              if(nStep.eq.2) SinkS(i)=SinkS(i)+AUptakeA
*             This is needed only for the last node, but that node may not have beta
              SAUptakeA =Omega
              SAUptakeAN=Omega
              if(SAUptakeP.gt.0.) Omega1=Omega/SAUptakeP
            else if(iStep.eq.3) then
*             This is needed only for the first node, but that node may not have beta
              if(Omega1.lt.OmegaS.and.Omega1.gt.0.) Compen=OmegaS
              if(Omega1.ge.OmegaS)                  Compen=Omega1
              if(Compen.gt.0.) AUptakeA=cc/(rKM+cc)*
     !                                          Beta(i)*SAUptakeP/Compen
              SinkS(i)=SinkS(i)+AUptakeA
              SAUptakeA=SAUptakeA+AUptakeA*dxM
            end if
          else
            SinkS(i)=0.
          end if
11      continue
        if(iStep.eq.nStep.and.jS.eq.NS) 
     !    write(78,100) t,SPUptake,SAUptakeP,SAUptakeA,SAUptakeAN ! the last is uncompensated
12    continue
      return

100   format(3x,e14.7,1x,4e12.4)
      end

************************************************************************

      real function FSAlfa(lMode,cRed,c50,P3c)

      logical lMode

      if(lMode) then
        FSAlfa=0.
        if(abs(c50).gt.0) FSAlfa=1./(1.+(cRed/c50)**P3c)
      else
        if(cRed.le.c50) then
          FSAlfa=1.
        else
          FSAlfa=max(0.,1.-(cRed-c50)*P3c*0.01)
        end if 
      end if
      return
      end

************************************************************************

      real function FAlfa(lMoSink,TPot,h,P0,P1,P2H,P2L,P3,r2H,r2L)

      logical lMoSink

      if(lMoSink) then
        if(TPot.lt.r2L) P2=P2L
        if(TPot.gt.r2H) P2=P2H
        if((TPot.ge.r2L).and.(TPot.le.r2H))
     !    P2=P2H+(r2H-TPot)/(r2H-r2L)*(P2L-P2H)
        FAlfa=0.0
        if((h.gt.P3).and.(h.lt.P2)) FAlfa=(h-P3)/(P2-P3)
        if((h.ge.P2).and.(h.le.P1)) FAlfa=1.0
        if((h.gt.P1).and.(h.lt.P0).and.P0-P1.gt.0.) FAlfa=(h-P0)/(P1-P0)
*       Uptake even at full saturation, when both P1 and P2 are equal to zero
        if(h.ge.P2.and.P1.eq.0..and.P0.eq.0.) FAlfa=1.0
      else
        FAlfa=1./(1.+(h/P0)**P3)
      end if
      return
      end

************************************************************************

      subroutine SetRG(NumNP,x,Beta,t,tRMin,tRHarv,xRMin,xRMax,RGR,
     !                 xRoot,lRoot,iRootIn,nGrowth,rGrowth,tRPeriod)

      dimension x(NumNP),Beta(NumNP),rGrowth(1000,5)
      double precision t
      logical lRoot

      if(lRoot.and.iRootIn.eq.1) then
        i=1000
        j=5
        call Table(nGrowth,rGrowth,i,j,t,rDummy,rDummy,rDummy,xRoot)
      end if

      if(lRoot.and.iRootIn.eq.2) then
        tRoot=amod(sngl(t),tRPeriod)
        if(tRoot.lt.tRMin.or.tRoot.gt.tRHarv) then
          do 11 i=1,NumNP
            Beta(i)=0.
11        continue
          return
        end if
        xR=xRMax
        if(xRMin.le.0.001) xRMin=0.001
        tt=tRoot-tRMin
        xR=(xRMax*xRMin)/(xRMin+(xRMax-xRMin)*exp(-RGR*tt))
      else
        xR=xRoot
      end if

      SBeta=0.
      do 12 i=2,NumNP-1
        if(x(i).lt.x(NumNP)-xR) then
          Beta(i)=0.
        else if(x(i).lt.x(NumNP)-0.2*xR) then
          Beta(i)=2.08333/xR*(1-(x(NumNP)-x(i))/xR)
        else
          Beta(i)=1.66667/xR
        end if
        if(i.ne.NumNP) then
          SBeta=SBeta+Beta(i)*(x(i+1)-x(i-1))/2.
        else
          SBeta=SBeta+Beta(i)*(x(i)-x(i-1))/2.
        end if
12    continue
      if(SBeta.lt.0.0001) then
        Beta(NumNP-1)=1./((x(NumNP)-x(NumNP-2))/2.)
      else
        do 13 i=2,NumNP-1
          Beta(i)=Beta(i)/SBeta
13      continue
      end if
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

*     Jasper Vrugt function

      subroutine RootInN(NumNP,Beta,z)

      dimension Beta(NumNP),z(NumNP)

*     read input
      read(44,*)
      read(44,*) Zm,Z0,Za

*     coordinate of the surface
      ZMax =z(NumNP)

*     calculate non-normalized uptake intensity
      r1=0.
      r2=0.
      do 11 i=1,NumNP
        if(abs(Zm).gt.1.e-5) then
          r1=(Zm-(ZMax-z(i)))/Zm
          rRootA=Za
          if((ZMax-z(i)).gt.Z0) rRootA=1.
          r2=rRootA/(Zm)*abs(Z0-(ZMax-z(i)))
          r2=exp(-r2)
        end if
        Beta(i)=amax1(r1*r2,0.)
11    continue

*     normalize uptake intensity
      SBeta=Beta(NumNP)*(z(NumNP)-z(NumNP-1))/2.
      do 12 i=2,NumNP-1
        SBeta=SBeta+Beta(i)*(z(i+1)-z(i-1))/2.
12    continue
      do 13 i=2,NumNP
        if(SBeta.gt.0.) then
          Beta(i)=Beta(i)/SBeta
        else
          Beta(i)=0.
        end if
13    continue

      return
      end

************************************************************************